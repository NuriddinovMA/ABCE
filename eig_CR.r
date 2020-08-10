args <- commandArgs(trailingOnly = T)
print(args)
print("start")
fname <- args[2]
pc <- as.integer(args[3])
res <- as.integer(args[4])
chr <- args[5]
st <- as.integer(as.integer(args[6])/res)+1
fi <- as.integer(as.integer(args[7])/res)
m.p <- c()

m.s <- scan(fname, sep=' ')
lnf <- as.integer( sqrt(length(m.s)) )
m.m <- matrix(m.s,nrow=lnf,byrow=T)
m.m[!is.finite(m.m)] <- 0

m.pca <- princomp(m.m)
m.pr <- predict(m.pca)[,pc]
m.p <- m.pr
ln <- length(m.p)

pos1 <- seq.int(0,ln-1,1)
pos1 <- (pos1 + st)*res
pos2 <- pos1 + res - 1
t.out <- format(data.frame(chr,pos1,pos2,m.p),trim=TRUE,scientific=FALSE)
fout <- paste0(fname, '.cropped.pc', pc, '.eig.bedGraph')

write.table(file=fout,t.out,dec=",",col.names=F,row.names=F,sep='\t',quote=FALSE)
