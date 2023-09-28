import os
import timeit
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Dump locus of interest from given oe-matrix and prepare it for compartment calculation')
parser.add_argument('--input','-i', help='The name of input file. This file must contains a dense observed/expected matrix. See juicebox dump')
parser.add_argument('--out','-o', default=False, help='The directory for a corrected matrix output. By default, the output matrix is placed with input matrix.')
parser.add_argument('--pic','-p', default=False, help='The directory for a picture output. By default, no picture output.')
parser.add_argument('--locus','-l', nargs=2, type=int, default=(0,0), help='Genomic coordinates of locus of interest in bp. For example "-l 1000000 10000000", if you wish analyze the locus between 1Mb and 10 Mb. By default, the script processes the full matrix.')
parser.add_argument('--resolution','-r',type=int, help='The matrix resolution in bp')
args=parser.parse_args()

fname = args.input.split('/')[-1]
dirname = args.input[:-len(fname)]
if args.out == False: args.out = dirname

A = []
start_time = timeit.default_timer()

parse = fname.split('.')
elp = timeit.default_timer() - start_time
print( 'start processing %s %.2f' % (fname, elp) )
M = np.genfromtxt(args.input,dtype=np.float32)
ln = len(M)
print( ln, args.locus, args.resolution )
if args.locus[1] != 0: args.locus = args.locus[0]//args.resolution,args.locus[1]//args.resolution
else: args.locus = 0,ln

print( args.locus )
M = M[args.locus[0]:args.locus[1],args.locus[0]:args.locus[1]]
ln = len(M[0])
print( ln )
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
elp = timeit.default_timer() - start_time
print( '\treading oe %.2f' % elp )
for i in range(ln):
	for j in range(i,ln):
		if i in ma or j in ma:
			M[i,j] = np.nan
			M[j,i] = np.nan
		else: pass
M = M - np.nanmean(M,axis=0)
M = np.nan_to_num(M)
if args.pic != False:
	im = plt.imshow(M, cmap='coolwarm',vmin = -1.5, vmax = 1.5)
	plt.colorbar(im)
	plt.savefig(args.pic+fname+'.cropped.norm.png',dpi=400)
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
print( '\tcorrelating and cleaning %.2f' % elp )

if args.pic != False:
	im = plt.imshow(M, cmap='coolwarm',vmin = -.15, vmax = .15)
	plt.colorbar(im)
	plt.savefig(args.pic+fname+'.cropped.prs.png',dpi=400)
	plt.clf()
np.savetxt(args.out+fname+'.cropped.prs', M)

M = np.zeros([ln,ln])
elp = timeit.default_timer() - start_time
print( '\t\t writing contracting matrix frame %.2f' % elp)
elp = timeit.default_timer() - start_time
print ('end processing %s %.2f' % (fname, elp))