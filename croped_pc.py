import os
import timeit
import numpy as np
import matplotlib.pyplot as plt

dirname = 'C:/Desktop/juicebox/dump/GEH/dencity/' # name of directory containing observed/expected matrices in dencity format
outname = 'C:/Users/admin/Documents/R/E/' # directory of corrected matrix output
picname = 'C:/Desktop/juicebox/scripts/pictures/' # directory of picture output

fname = 'geh.2.2.100.kr.oe' # name of matrix
resolution = 100000 # matrix resolution in bp
rge = (82000000,115000000) # genome coordinate in bp within matrix analysed
A = []
start_time = timeit.default_timer()

parse = fname.split('.')
elp = timeit.default_timer() - start_time
print 'start processing %s %.2f' % (fname, elp)
M = np.genfromtxt(dirname+fname,dtype=np.float32)
ln = len(M)
print ln

if rge[1] != 0: rge = rge[0]/resolution,rge[1]/resolution
else: rge = 0,ln

print rge
M = M[rge[0]:rge[1],rge[0]:rge[1]]
ln = len(M[0])
print ln
M = np.nan_to_num(M)
ma = set([])
for i in range(0,ln-80):
	if np.count_nonzero(M[i,i:(i+80)]) < 60: ma.add(i)
for i in range(80,ln):
	if np.count_nonzero(M[(i-80):i,i]) < 60: ma.add(i)
ma = list(ma)
ma.sort()
for i in range(1,len(ma)):
	if (ma[i] - ma[i-1]) < 4: ma.extend([j for j in range(ma[i-1]+1,ma[i])])
ma.sort()
print ma
elp = timeit.default_timer() - start_time
print '\treading oe %.2f' % elp
for i in range(ln):
	for j in range(i,ln):
		if i in ma or j in ma:
			M[i,j] = np.nan
			M[j,i] = np.nan
		else: pass
elp = timeit.default_timer() - start_time
print '\tresolution combined %.2f' % elp
M = M - np.nanmean(M,axis=0)
M = np.nan_to_num(M)
im = plt.imshow(M, cmap='coolwarm',vmin = -1.5, vmax = 1.5)
plt.colorbar(im)
plt.savefig(picname+fname[:-3]+'.cropped.norm.png',dpi=400)
plt.clf()
M = np.corrcoef(M)
M[ma,:] = np.nan
M[:,ma] = np.nan
for i in range(ln): 
	try:
		M[i,i] = np.nan
		M[i,i+1] = np.nan
		M[i+1,i] = np.nan
	except IndexError: pass
elp = timeit.default_timer() - start_time
print '\tbase transformations %.2f' % elp

im = plt.imshow(M, cmap='coolwarm',vmin = -.15, vmax = .15)
plt.colorbar(im)
plt.savefig(picname+fname[:-3]+'.cropped.prs.png',dpi=400)
plt.clf()
np.savetxt(outname+fname[:-3]+'.cropped.prs', M)
M = np.zeros([ln,ln])
elp = timeit.default_timer() - start_time
print '\t\t writing contracting matrix frame %.2f' % elp
elp = timeit.default_timer() - start_time
print 'end processing %s %.2f' % (fname, elp)