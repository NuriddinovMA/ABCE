args <- commandArgs(trailingOnly = T)
print(args)
print("start")
fname <- args[1]
corname <- args[2]
res <- as.integer(args[3])
pc <- as.integer(args[4])
chr <- args[5]
st <- as.integer(as.integer(args[6])/res)+1
fi <- as.integer(as.integer(args[7])/res)
m.p <- c()

m.s <- scan(fname, sep=' ')
lnf <- as.integer( sqrt(length(m.s)) )
m.m <- matrix(m.s,nrow=lnf,byrow=T)
m.m[!is.finite(m.m)] <- 0
gc.c <- read.table(corname)[,4]
gc.c[gc.c>20] <- 20

m.pca <- princomp(m.m[st:fi,st:fi])
m.pr <- predict(m.pca)[,pc]
m.h <- m.pr
m.h[m.h>10] <- 10
m.h[m.h<-10] <- -10
cr <- cor(m.h,gc.c[st:fi], method='pearson')
print(cr)

if (cr > 0) { m.p <- m.pr } else { m.p <- -1*m.pr }

ln <- length(m.p)
pos1 <- seq.int(0,ln-1,1)
pos1 <- (pos1 + st)*res
pos2 <- pos1 + res - 1
t.out <- format(data.frame(chr,pos1,pos2,m.p),trim=TRUE,scientific=FALSE)
fout <- paste0(fname, '.whole.pc', pc, '.eig.bedGraph')

write.table(file=fout,t.out,dec=",",col.names=F,row.names=F,sep='\t',quote=FALSE)
