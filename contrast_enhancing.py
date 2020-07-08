import os
import timeit
import numpy as np
import matplotlib.pyplot as plt

dirname = 'C:/Desktop/juicebox/dump/GEH/dencity/' # name of directory containing observed/expected matrices in dencity format
outname = 'C:/Users/admin/Documents/R/E/' # directory of corrected matrix output
picname = 'C:/Desktop/juicebox/scripts/pictures/' # directory of picture output

fname = 'geh.2.2.100.kr.oe' # name of matrix
resolution = 100000 # matrix resolution in bp
frame = (0,1) # range of frame within contact smoothed (see error.py)
dist = (0,122) # distance in bins to change range of smoothing (see error.py)
# smoothing is equal resolution lowering
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
A = np.zeros([ln,ln])
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
#n = np.where((M==0))
#M[n[0],n[1]] = np.nan
# M[ma,:] = np.nan
# M[:,ma] = np.nan
# print ma
# im = plt.imshow( M, cmap='coolwarm',vmin = 0, vmax = 3)
# plt.colorbar(im)
# plt.savefig(fname+'.easy.png',dpi=400)
# plt.clf()
# plt.cla()
# C = np.corrcoef(M)
# im = plt.imshow(C, cmap='coolwarm',vmin = -1, vmax = 1)
# plt.colorbar(im)
# plt.savefig(fname+'.easy.prs.png',dpi=400)
# plt.clf()
elp = timeit.default_timer() - start_time
print '\treading oe %.2f' % elp
for d in range(len(dist)-1):
	if ln > d:
		f = frame[d]
		for k in range(dist[d],dist[d+1]):
			for i in range(ln):
				i0 = i-f
				i1 = i+f+1
				j0 = i0+k
				j1 = i1+k
				if i0 < 0: i0 = 0
					#i1 = 2*f + 1
				if i1 >= ln: i1 = ln
					#i0 = ln-2*f-1
				if j0 < 0: j0 = 0
					#j1 = 2*f + 1
				if j1 >= ln: j1 = ln
					#j0 = ln-2*f-1
				try:
					# if ((i1 - i0) * (j1-j0) >= f**2) and np.isnan(M[i,i+k]) != True:
					if i in ma or (i+k) in ma:
						A[i,(i+k)] = np.nan
						A[(i+k),i] = np.nan
					else:
						C = np.nanmedian(M[i0:i1,j0:j1])
						A[i,(i+k)] = C
						A[(i+k),i] = C
				except IndexError: pass
	else: pass
elp = timeit.default_timer() - start_time
print '\tresolution combined %.2f' % elp
A = np.nan_to_num(A)
im = plt.imshow(A, cmap='coolwarm',vmin = 0, vmax = 2.5)
plt.colorbar(im)
plt.savefig(fname+'.smooth.median.png',dpi=400)
plt.clf()
# np.savetxt(outname+fname+'.smouth.oe', A)
A = A - np.nanmean(A,axis=0)
im = plt.imshow(A, cmap='coolwarm',vmin = -2.5, vmax = 2.5)
plt.colorbar(im)
plt.savefig(fname+'.norm.mean.png',dpi=400)
plt.clf()
A = np.nan_to_num(A)
A = np.corrcoef(A)
A[ma,:] = np.nan
A[:,ma] = np.nan
im = plt.imshow(A, cmap='coolwarm',vmin = -1, vmax = 1)
plt.colorbar(im)
plt.savefig(fname+'.smouth.prs.png',dpi=400)
plt.clf()
# np.savetxt(outname+fname+'.smouth.prs', A)
S = np.zeros([ln,ln])
P = np.zeros([ln,ln])
for i in range(ln): 
	try:
		A[i,i] = np.nan
		A[i,i+1] = np.nan
		A[i+1,i] = np.nan
	except IndexError: pass
elp = timeit.default_timer() - start_time
print '\tbase transformations %.2f' % elp
for f in [1,3,5]: # range of frame within contacts contrasted
	elp = timeit.default_timer() - start_time
	print '\t\tstart contracting frame %i %.2f' % (f,elp)
	for i in range(ln):
		for j in range(i+1,ln):
			i0 = i-f
			i1 = i+f+1
			if i0 < 0: 
				i0 = 0
				i1 = 2*f + 1
			if i1 >= ln:
				i0 = ln-2*f-1
				i1 = ln
			j0 = j-f
			j1 = j+f+1
			if j0 < 0: 
				j0 = 0
				j1 = 2*f + 1
			if j1 >= ln:
				j0 = ln-2*f-1
				j1 = ln
			# try:
			R = np.nanpercentile(A[i0:i1,j0:j1],75),np.nanpercentile(A[i0:i1,j0:j1],25)
			C = (R[0] + R[1])/2.
			D = (R[0] - R[1])/2.
			P[i,j] = (A[i,j] - C)/D
			P[j,i] = P[i,j]
			R = np.nanmax(A[i0:i1,j0:j1]),np.nanmin(A[i0:i1,j0:j1])
			C = (R[0] + R[1])/2.
			D = (R[0] - R[1])/2.
			S[i,j] = (A[i,j] - C)/D
			S[j,i] = S[i,j]
	elp = timeit.default_timer() - start_time
	print '\t\tend contracting frame %i %.2f' % (f,elp)
	im = plt.imshow(P, cmap='coolwarm',vmin = -2, vmax = 2)
	plt.colorbar(im)
	plt.savefig(picname+fname+'.%i.smouth.prc.prs.png' % f,dpi=400)
	plt.clf()
	np.savetxt(outname+fname+'.%i.smouth.prc.prs' % f, P)
	P = np.zeros([ln,ln])
	im = plt.imshow(S, cmap='coolwarm',vmin = -1, vmax = 1)
	plt.colorbar(im)
	plt.savefig(picname+fname+'.%i.smouth.range.prs.png' % f,dpi=400)
	plt.clf()
	np.savetxt(outname+fname+'.%i.smouth.range.prs' % f, S)
	S = np.zeros([ln,ln])
	elp = timeit.default_timer() - start_time
	print '\t\t writing contracting matrix frame %i %.2f' % (f,elp)
del P
del S
elp = timeit.default_timer() - start_time
print 'end processing %s %.2f' % (fname, elp)
