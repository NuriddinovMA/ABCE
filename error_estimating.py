import os
import timeit
import numpy as np
import matplotlib.pyplot as plt

dirname = 'C:/Desktop/juicebox/dump/GEH/dencity/' #name of directory containing RAW CONTACT matrices in dencity format
fname = 'geh.2.2.100.none.oe' # name of matrix
outname = 'errors.txt' # name of output file

# this script analyzed dependence of contact variability from distance between contacted loci
# the output contains distance and frame with contact shall be united
# this parameters are used by contrast_enhanced.py script 
# frame as frame = (0,1,...)
# distance as dist = (0,122,...)

start_time = timeit.default_timer()

out = open(outname,'w')
print 'start processing %s' % fname
print >> out, fname
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
		print >> out, k,f
X = list(range(len(Y)))
elp = timeit.default_timer() - start_time

print 'end processing %s %.2f' % (fname, elp)
out.close()
elp = timeit.default_timer() - start_time
print 'total end %.2f' % elp