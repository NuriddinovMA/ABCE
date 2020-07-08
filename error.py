import os
import timeit
import numpy as np
import matplotlib.pyplot as plt
#import oeData

dirname = 'C:/Desktop/juicebox/dump/GEH/none/'
flist = os.listdir(dirname)
file = open('error.txt','w')
start_time = timeit.default_timer()

for fname in flist:
	if fname.split('.')[-1] == 'none':
		print 'start processing %s' % fname
		print >> file, fname
		X,Y,H = [],[],[]
		f = 0
		M = np.genfromtxt(dirname+fname,dtype=np.float32)
		ln = len(M[0])
		for k in range(ln):
			for i in range(ln):
				try:
					C = np.nansum(M[(i-f):(i+f+1),(i+k-f):(i+k+f+1)])
					if C != 0: H.append( np.sqrt(1./C) )
				except IndexError: pass
			Y.append( np.nanmean(H) )
			if Y[-1] > 0.2: 
				f += 1
				print >> file, k,f
		X = list(range(len(Y)))
		elp = timeit.default_timer() - start_time
		print 'end processing %s %.2f' % (fname, elp)
file.close()
elp = timeit.default_timer() - start_time
print 'total end %.2f' % elp