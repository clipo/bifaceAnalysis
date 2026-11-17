options(width=100)
library(Momocs)
a <- seq(0,2*pi,length.out=101)[-101]
q <- cbind(cos(a),2*sin(a))
for(n in 1:3) for(k in 2:3) {
  cat('n=',n,'k=',k,'\n')
  print(tryCatch({
    o <- Out(lapply(seq_len(n),function(j) q %*% diag(c(1+j/10,1))))
    e <- efourier(o,nb.h=6,norm=TRUE)
    p <- PCA(e)
    cl <- hclust(dist(e$coe),method='ward.D2')
    z <- cutree(cl,k=k)
    list(scores=dim(p$x),PC2=p$x[,2],clusters=z)
  },error=function(e)conditionMessage(e)))
}
