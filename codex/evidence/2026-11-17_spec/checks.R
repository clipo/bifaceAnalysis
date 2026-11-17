options(width=140, digits=12)
library(imager)
cat('AUDIT environment\n',R.version.string,'\n')
for(p in c('Momocs','imager','shiny','bslib','ggplot2','cluster','testthat')) cat(p,as.character(packageVersion(p)),'\n')
cat('HEAD',system('git rev-parse HEAD',intern=TRUE),'\n')
cat('\nInstalled Momocs 1.4.1 source bodies\n')
for(f in c('efourier.default','efourier.Out','coo_align.default')) print(get(f,asNamespace('Momocs')))

# Independent implementation of the displayed equations. Closing segment comes first
# here solely to match Momocs' phase convention, which starts at the LAST vertex.
efa <- function(q,H=10) {
  prev <- q[c(nrow(q),seq_len(nrow(q)-1)),,drop=FALSE]
  d <- q-prev; dt <- sqrt(rowSums(d*d)); stopifnot(all(dt>0))
  t <- cumsum(dt); t0 <- c(0,head(t,-1)); T <- sum(dt)
  out <- t(vapply(seq_len(H),function(n) {
    dc <- cos(2*n*pi*t/T)-cos(2*n*pi*t0/T)
    ds <- sin(2*n*pi*t/T)-sin(2*n*pi*t0/T)
    T/(2*n*n*pi*pi)*c(sum(d[,1]/dt*dc),sum(d[,1]/dt*ds),sum(d[,2]/dt*dc),sum(d[,2]/dt*ds))
  },numeric(4)))
  colnames(out)<-c('a','b','c','d')
  list(coef=out, dc=colSums((q+prev)*dt)/(2*T),t=t/T)
}
ang <- (0:1999)*2*pi/2000
q <- cbind(3*cos(ang),sin(ang))
e <- efa(q,30); m <- Momocs::efourier(q,nb.h=30)
cat('\nELLIPSE coefficients and amplitudes\n');print(e$coef[1:5,])
print(c(first_x=sqrt(sum(e$coef[1,1:2]^2)),first_y=sqrt(sum(e$coef[1,3:4]^2)),
        largest_higher=max(abs(e$coef[-1,])),
        max_Momocs_error=max(abs(e$coef-cbind(m$an,m$bn,m$cn,m$dn)))))
cat('first five powers\n');print(rowSums(e$coef[1:5,]^2)/2)
# Arithmetic-centering is not arc-length centering.
tri <- rbind(c(0,3),c(-1,0),c(2,0)); tri<-sweep(tri,2,colMeans(tri))
cat('\nCentered triangle Fourier DC\n');print(efa(tri)$dc)
# Closed edge resampling with both endpoints repeats T and B.
edge1<-cbind(seq(0,1,length.out=100),rep(0,100))
edge2<-cbind(seq(1,0,length.out=100),sin(seq(0,pi,length.out=100)))
edge2[100,]<-edge1[1,]
cat('\nEndpoint-inclusive resampling rows/distinct\n');print(c(nrow(rbind(edge1,edge2)),nrow(unique(rbind(edge1,edge2)))))
# Close eigenvalues with an order exchange, but an unchanged two-dimensional plane.
X<-rbind(c(1.01,0),c(-1.01,0),c(0,1),c(0,-1))
p<-prcomp(X); p1<-prcomp(X[-1,])
cat('\nPCA eigenvalue/order counterexample\n');print(p$sdev^2);print(p1$sdev^2)
print(abs(diag(crossprod(p$rotation,p1$rotation))))
cat('Subspace principal-angle cosines\n');print(svd(crossprod(p$rotation,p1$rotation))$d)
cat('Three specimens leave-one-out centered ranks\n');print(c(qr(scale(X[1:3,],scale=FALSE))$rank,qr(scale(X[1:2,],scale=FALSE))$rank))
rv<-function(X,Y) {
 X<-scale(X,scale=FALSE);Y<-scale(Y,scale=FALSE)
 sum(crossprod(X,Y)^2)/sqrt(sum(crossprod(X)^2)*sum(crossprod(Y)^2))
}
set.seed(1711);X<-matrix(rnorm(4*40),4);Y<-matrix(rnorm(4*80),4)
perms<-replicate(999,rv(X,Y[sample(4),]))
cat('\nRV checks identical/independent observed/permutation mean/range\n')
print(c(identical=rv(X,X),observed=rv(X,Y),mean=mean(perms),quantile(perms,c(.025,.975))))

otsu<-function(z) {
 h<-hist(z,breaks=seq(0,1,length.out=257),plot=FALSE)$counts
 p<-h/sum(h); w<-cumsum(p); mu<-cumsum(p*(0:255)); mt<-sum(p*(0:255))
 score<-(mt*w-mu)^2/(w*(1-w));score[!is.finite(score)]<--Inf
 which.max(score)/256
}
components<-function(mask,high=FALSE) as.matrix(label(as.cimg(mask*1),high_connectivity=high))
polygon_info<-function(q) {
 r<-q[c(2:nrow(q),1),];cross<-q[,1]*r[,2]-r[,1]*q[,2];A<-sum(cross)/2
 cen<-colSums((q+r)*cross)/(6*A)
 xx<-sum((q[,1]^2+q[,1]*r[,1]+r[,1]^2)*cross)/(12*A)-cen[1]^2
 yy<-sum((q[,2]^2+q[,2]*r[,2]+r[,2]^2)*cross)/(12*A)-cen[2]^2
 xy<-sum((2*q[,1]*q[,2]+q[,1]*r[,2]+r[,1]*q[,2]+2*r[,1]*r[,2])*cross)/(24*A)-prod(cen)
 eg<-eigen(matrix(c(xx,xy,xy,yy),2));j<-which.max(abs(eg$vectors[1,]));u<-eg$vectors[,j];if(u[1]>0)u<--u
 v<-c(-u[2],u[1]);cents<-sweep(q,2,cen);a<-drop(cents%*%u);b<-drop(cents%*%v)
 bn<-b[c(2:length(b),1)];an<-a[c(2:length(a),1)];hit<-which((b<=0 & bn>0)|(bn<=0 & b>0))
 ints<-a[hit]-b[hit]*(an[hit]-a[hit])/(bn[hit]-b[hit])
 list(area=abs(A),ratio=eg$values[1]/eg$values[2],angle=atan2(u[2],u[1])*180/pi,
      extent_ratio=diff(range(a))/diff(range(b)),tip_ties=sum(abs(a-max(a))<1e-10),base_intersections=sum(ints<0),all_intersections=sort(ints))
}
cat('\nSEGMENTATION assumptions: grayscale Luma; first Otsu maximum; disk radius 2; 4-connectivity; default whole image.\n')
mask_disk<-outer(-2:2,-2:2,function(x,y) as.numeric(x*x+y*y<=4))
for(i in 1:4) {
 file<-sprintf('../bifaceAnalysis/www/artifact%d.jpg',i)
 img<-load.image(file);g<-as.matrix(grayscale(img)); th<-otsu(g)
 border<-unique(cbind(c(rep(1,ncol(g)),rep(nrow(g),ncol(g)),1:nrow(g),1:nrow(g)),
                       c(1:ncol(g),1:ncol(g),rep(1,nrow(g)),rep(ncol(g),nrow(g)))))
 above<-g>th; bgabove<-mean(above[border])>.5;raw<-if(bgabove)!above else above
 cat('IMAGE',i,'border light fraction',mean(above[border]),'light border pixels',sum(above[border]),'of',nrow(border),'\n')
 op<-as.matrix(mopening(as.cimg(raw*1),as.cimg(mask_disk),boundary_conditions=FALSE))>.5
 lab<-components(!op);bgids<-unique(lab[!op]);edgeids<-unique(lab[border]);holeids<-setdiff(bgids,edgeids)
 filled<-op | matrix(lab %in% holeids,nrow(op));labels<-components(filled)
 areas<-sort(table(labels[filled]),decreasing=TRUE);keep<-filled & labels==as.integer(names(areas)[1])
 cont<-contourLines(1:nrow(keep),1:ncol(keep),keep*1,levels=.5)
 closes<-vapply(cont,function(z) sqrt((z$x[1]-tail(z$x,1))^2+(z$y[1]-tail(z$y,1))^2)<1e-8,logical(1))
 cat('\nIMAGE',basename(file),'dimensions',nrow(g),ncol(g),'threshold',th,'background above?',bgabove,'\n')
 cat('top component areas',head(as.numeric(areas),4),'hole pixels filled',sum(filled)-sum(op),'\n')
 cat('retained area fraction',mean(keep),'border touches',any(keep[border]),'contours',length(cont),'closed',closes,'vertices',sapply(cont,function(z)length(z$x)),'\n')
 if(length(cont)==1) print(polygon_info(cbind(head(cont[[1]]$x,-1),-head(cont[[1]]$y,-1))))
 png(sprintf('evidence/2026-11-17_spec/data/segmentation%d.png',i),width=900,height=300)
 par(mfrow=c(1,3),mar=c(1,1,2,1));plot(img,main='Original');plot(as.cimg(raw),main='Threshold');plot(as.cimg(keep),main='Retained mask');dev.off()
 # Diagnostic only: manually force light foreground, holding every other step fixed.
 op<-as.matrix(mopening(as.cimg(above*1),as.cimg(mask_disk),boundary_conditions=FALSE))>.5
 lab<-components(!op);holeids<-setdiff(unique(lab[!op]),unique(lab[border]))
 filled<-op | matrix(lab %in% holeids,nrow(op));labels<-components(filled)
 areas<-sort(table(labels[filled]),decreasing=TRUE);keep<-filled & labels==as.integer(names(areas)[1])
 cont<-contourLines(1:nrow(keep),1:ncol(keep),keep*1,levels=.5)
 cat('Forced-light diagnostic: component areas',head(as.numeric(areas),4),'contours',length(cont),'border',any(keep[border]),'\n')
 if(length(cont)==1) print(polygon_info(cbind(head(cont[[1]]$x,-1),-head(cont[[1]]$y,-1))))
}
cat('\nDependency checks\n')
for(f in c('silhouette','hclust','prcomp')) {cat(f,': ');print(find(f))}
cat('silhouette supplied by cluster:',exists('silhouette',asNamespace('cluster')),'\n')
cat('\nComplete loaded-package versions\n');print(sessionInfo())
cat('Script finished\n')
