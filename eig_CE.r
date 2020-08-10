args <- commandArgs(trailingOnly = T)
print(args)
a<-getwd()
print("111")
fname <- args[2]
corname <- args[3]
res <- as.integer(args[4])
fr <- as.integer(args[5])
chr <- args[6]
st <- as.integer(as.integer(args[7])/res)+1
fi <- as.integer(as.integer(args[7])/res)
m.p <- c()

m.s <- scan(fname, sep=' ')
lnf <- as.integer( sqrt(length(m.s)) )
m.m <- matrix(m.s,nrow=lnf,byrow=T)
m.m[!is.finite(m.m)] <- 0
gc.c <- read.table(corname)[,4]
gc.c[gc.c>20] <- 20

m.s <- m.m[st:fi,st:fi]
gc.s <- gc.c[st:fi]
gc.s[!is.finite(gc.s)] <- 0

ln <- (fi-st)
s <- seq.int(0,ln,fr)
step <- as.integer(ln/fr)
s <- cbind(s,s)
s[,2] <- s[,2]+fr-1
s[step,2] <- ln
print(s)

for(i in 1:step){
m.h <- m.s[s[i,1]:s[i,2],s[i,1]:s[i,2]]
m.pca <- princomp(m.h)
m.pr <- predict(m.pca)[,1]
m.h <- m.pr
m.h[m.h>10] <- 10
m.h[m.h<-10] <- -10

cr <- cor(m.h,gc.s[s[i,1]:s[i,2]], method='pearson')
print(cr)
if (cr > 0) { m.p <- c(m.p, m.pr ) }
else { m.p <- c(m.p, -1*m.pr) }
}

ln <- length(m.p)
pos1 <- seq.int(0,ln-1,1)
pos1 <- (pos1 + st)*res
pos2 <- pos1 + res - 1
t.out <- format(data.frame(chr,pos1,pos2,m.p),trim=TRUE,scientific=FALSE)
fout <- paste0(fname, '.eig.bedGraph')
write.table(file=fout,t.out,dec=",",col.names=F,row.names=F,sep='\t',quote=FALSE)
