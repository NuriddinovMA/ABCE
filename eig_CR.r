args <- commandArgs(trailingOnly = T)
print(args)
a<-getwd()
print("111")
sp <- args[2]
res <- as.integer(args[3])
fr <- as.integer(args[4])
chr <- args[5]
base <- as.integer(as.integer(args[6])/res)
n <- as.integer(length(args))
n <- seq.int(6,n,2)
m.p <- c()
N <- c(5)


for(k in N){
fn <- paste0('E/',sp)
m.s <- scan(fn, sep=' ')
lnf <- as.integer( sqrt(length(m.s)) )
m.m <- matrix(m.s,nrow=lnf,byrow=T)
m.m[!is.finite(m.m)] <- 0
print(lnf)
fn <- paste0('GC/chrmed/',sp,'.1.',chr,'.',res,'.rnaseq.bedGraph')#'.gene.dencity.bedGraph')#
#print(fn)
gc.c <- read.table(fn)[,4]
gc.c[gc.c>100] <- 100
for(i in n) {
print(args[i:(i+1)])
start <- as.integer(as.integer(args[i])/res)+1
end <- as.integer(as.integer(args[i+1])/res)
print(start)
print(end)
m.s <- m.m[(start-base):(end-base),(start-base):(end-base)]
l <- length(gc.c)
print(l)
gc.s <- gc.c[start:end]
gc.s[!is.finite(gc.s)] <- 0
l <- length(gc.s)
print(l)
m.pca <- princomp(m.s)
m.pr <- predict(m.pca)[,1]
m.h <- m.pr
l <- length(m.h)
print(l)
m.h[m.h>10] <- 10
m.h[m.h<-10] <- -10

cr <- cor(m.h,gc.s, method='pearson')
print(cr)
if (cr > 0) { m.p <- c(m.p, m.pr ) }
else { m.p <- c(m.p, -1*m.pr) }

}
ln <- length(m.p)
print(ln)
print(base)
pos1 <- seq.int(0,ln-1,1)
pos1 <- (pos1 + base)*res
pos2 <- pos1 + res - 1
t.out <- format(data.frame(chr,pos1,pos2,m.p),trim=TRUE,scientific=FALSE)
fout <- paste0('EIG/chrm/',spA,'.',chr,'.cropped.pc1.',res,'.eig.bedGraph')
#print(fn)
write.table(file=fout,t.out,dec=",",col.names=F,row.names=F,sep='\t',quote=FALSE)
}