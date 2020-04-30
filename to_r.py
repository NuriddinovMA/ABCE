import os
Rpath = "C:/Arbeit/R/R-3.6.3/bin/x64/Rscript.exe"
args_path = "AB.cat.txt"
f = open(args_path, 'r')
lines = f.readlines()
f.close()

for line in lines:
	args = line.strip()
	print line
	out = os.system( Rpath + " --vanilla eig_CRC.r --args " + args + " >> cor.tpm.framed.txt")
	print(out)
	# out = os.system( Rpath + " --vanilla eig_CR.r --args " + args + " >> cor.tpm.framed.txt")
	# print(out)
	# out = os.system( Rpath + " --vanilla eig_CR2.r --args " + args + " >> cor.tpm.framed.txt")
	# print(out)
	# out = os.system( Rpath + " --vanilla eig_FR.r --args " + args + " >> cor.tpm.framed.txt")
	# print(out)
	# out = os.system( Rpath + " --vanilla eig_CE.r --args " + args + " >> cor.tpm.framed.txt")
	# print(out)