import os
import timeit
import numpy as np
import matplotlib.pyplot as plt

# import oeData
# from oeData import *
# [2L,2R,3L,3R,X]
# [0,51,640],[0,60,700,1586],[0,57,716],[0,56,742],[0,36]
# [0,1,3],[0,1,3,9],[0,1,4],[0,1,4],[0,1]
dirname = 'C:/Arbeit/juicebox/dump/Anopheles/dencity/'
outname = 'C:/Arbeit/RD/E/'
picname = 'C:/Arbeit/juicebox/pictures/'
flist = os.listdir(dirname)# ['AatrE3_V3.3L.3L.25000.oe']#['AalbS2_V3.X.X.25000.oe','AalbS2_V3.2L.2L.25000.oe','AalbS2_V3.2R.2R.25000.oe','AalbS2_V3.3L.3L.25000.oe','AalbS2_V3.3R.3R.25000.oe']##os.listdir(dirname)#['Aatr.X.oe',]#['Aste.2L.prs']#['Acol.X.oe','Acol.2L.oe','Acol.2R.oe','Acol.3L.oe','Acol.3R.oe']#
# frag = [0,51,640,1586]
# frame = [0,1,3,9]
#flist = ['AcolNg_V4.2R.2R.25000.oe',]
for i in range(len(flist)-1,-1,-1):
	if flist[i].split('.')[1] in ['2','3']: pass
	else: del flist[i]
Hrange = {
	'AalbS2_V3':{'X':(0,11425000)},
	'AatrE3_V3':{'X':(0,17400000),'3L':(7300000,41225000)},
	'AcolNg_V4':{'X':(0,20300000),'2R':(0,58325000),'2L':(4500000,48125000),'3R':(0,50475000),'3L':(3250000,42075000)},
	'AmerR4_V4':{'X':(0,19700000),'2R':(0,57275000),'2L':(10300000,54175000),'3L':(4075000,39925000)},
	'Aste':{'X':(0,16150000)},
	'Aqmac_V3R':{'2L':(0,46650000), '3L':(34000000,48000000)}
}
resolution = 25000
A = []
start_time = timeit.default_timer()
for fname in flist:
	parse = fname.split('.')
	elp = timeit.default_timer() - start_time
	print 'start processing %s %.2f' % (fname, elp)
	M = np.genfromtxt(dirname+fname,dtype=np.float32)
	ln = len(M)
	print ln
	try: rge = Hrange[parse[0]][parse[1]][0]/resolution,Hrange[parse[0]][parse[1]][1]/resolution
	except KeyError: rge = 0,ln
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