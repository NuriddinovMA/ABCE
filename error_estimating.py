import os
import argparse
import timeit
import numpy as np

parser = argparse.ArgumentParser(description='Calculate an area within the contacts are desirable united to avoid noise effects')
parser.add_argument('--input','-i', help='The name of input file. This file must contains a dense matrix with NONE-normalized contact counts. See juicebox dump')
parser.add_argument('--out','-o', help='The name of output file.')
parser.add_argument('--threshold','-t',type=int, default=25, help='The minimal mean contact number within some area for a given distance. If a contact number is less then the given threshold, the area are increased.')
args=parser.parse_args()


start_time = timeit.default_timer()

out = open(args.out,'w')
print 'start processing %s' % args.input
print >> out, args.input
print >> out, 'distance radius'
print >> out, '0 0'
X,Y,H = [],[],[]
f = 0
M = np.genfromtxt(args.input,dtype=np.float32)
ln = len(M[0])
for k in range(ln):
	if k % 10 == 0: print '.',
	for i in range(ln):
		try:
			C = np.nansum(M[(i-f):(i+f+1),(i+k-f):(i+k+f+1)])
			if C != 0: H.append( np.sqrt(1./C) )
		except IndexError: pass
	Y.append( np.nanmean(H) )
	if Y[-1] >= np.sqrt(1./args.threshold): 
		f += 1
		print >> out, k,f
X = list(range(len(Y)))
elp = timeit.default_timer() - start_time
print ''
print 'end processing %s %.2f' % (args.input, elp)
out.close()
elp = timeit.default_timer() - start_time
print 'total end %.2f' % elp