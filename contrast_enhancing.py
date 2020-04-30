import os
import timeit
import numpy as np
import matplotlib.pyplot as plt

# import oeData
# from oeData import *
# [2L,2R,3L,3R,X]
# [0,51,640],[0,60,700,1586],[0,57,716],[0,56,742],[0,36]
# [0,1,3],[0,1,3,9],[0,1,4],[0,1,4],[0,1]
dirname = 'C:/Desktop/juicebox/dump/Anopheles/dencity/'
outname = 'C:/Users/admin/Documents/R/E/'
picname = 'C:/Desktop/juicebox/scripts/pictures/'
flist = os.listdir(dirname)# ['AatrE3_V3.3L.3L.25000.oe']#['AalbS2_V3.X.X.25000.oe','AalbS2_V3.2L.2L.25000.oe','AalbS2_V3.2R.2R.25000.oe','AalbS2_V3.3L.3L.25000.oe','AalbS2_V3.3R.3R.25000.oe']##os.listdir(dirname)#['Aatr.X.oe',]#['Aste.2L.prs']#['Acol.X.oe','Acol.2L.oe','Acol.2R.oe','Acol.3L.oe','Acol.3R.oe']#
# frag = [0,51,640,1586]
# frame = [0,1,3,9]
for i in range(len(flist)-1,-1,-1):
	if flist[i].split('.')[0] in ['AmerR4_V4']:pass#,'AcolNg_V4']: pass
	else: del flist[i]
Hfrag = {'AalbS2_V3':(0,46,510,1516),'AatrE3_V3':(0,79,916),'AcolNg_V4':(0,77,574),
	'AmerR4_V4':(0,45,525),'Aste':(0,72,814),'Aqmac_V3R':(0,50,520)}
Hframe = {'AalbS2_V3':(0,1,5,10),'AatrE3_V3':(0,1,4),'AcolNg_V4':(0,1,4),
	'AmerR4_V4':(0,1,3),'Aste':(0,1,5),'Aqmac_V3R':(0,1,3)}
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
	print ma
	frag = Hfrag[parse[0]] + (ln,)
	frame = Hframe[parse[0]]
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
	for fr in range(len(frag)-1):
		if ln > fr:
			f = frame[fr]
			for k in range(frag[fr],frag[fr+1]):
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
	# A = np.nan_to_num(A)
	# im = plt.imshow(A, cmap='coolwarm',vmin = 0, vmax = 2.5)
	# plt.colorbar(im)
	# plt.savefig(fname+'.smooth.median.png',dpi=400)
	# plt.clf()
	# np.savetxt(outname+fname+'.smouth.oe', A)
	A = A - np.nanmean(A,axis=0)
	# im = plt.imshow(A, cmap='coolwarm',vmin = -2.5, vmax = 2.5)
	# plt.colorbar(im)
	# plt.savefig(fname+'.norm.mean.png',dpi=400)
	# plt.clf()
	A = np.nan_to_num(A)
	A = np.corrcoef(A)
	A[ma,:] = np.nan
	A[:,ma] = np.nan
	# im = plt.imshow(A, cmap='coolwarm',vmin = -1, vmax = 1)
	# plt.colorbar(im)
	# plt.savefig(fname+'.smouth.prs.png',dpi=400)
	# plt.clf()
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
	for f in [5]:#,10,15]:
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
				# except IndexError: pass
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