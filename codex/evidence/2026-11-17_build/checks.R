root <- normalizePath("..")
ev <- normalizePath("evidence/2026-11-17_build")
for (f in list.files(file.path(root,"bifaceAnalysis/R"),full.names=TRUE)) source(f)
ref <- readRDS(file.path(root,"tests/fixtures/reference.rds"))
cat("HEAD", system("git rev-parse HEAD",intern=TRUE),"\n")
cat("Fixture versions\n"); print(ref$versions)
attempt <- function(expr) tryCatch(expr,error=function(e) paste("ERROR",conditionMessage(e)))

# Independently project markers onto every segment of the plotted contour.
edge_distance <- function(p,o) {
  q <- o[c(2:nrow(o),1),]; v <- q-o
  t <- pmax(0,pmin(1,rowSums(sweep(o,2,p,"-") * -v)/rowSums(v*v)))
  min(sqrt(rowSums(sweep(o+t*v,2,p,"-")^2)))
}
files <- list.files(file.path(root,"bifaceAnalysis/www"),full.names=TRUE)
png(file.path(ev,"data/overlay.png"),1600,650)
par(mfrow=c(1,4),mar=c(0,0,2,0))
for(f in files) {
  im <- load_working_image(f); p <- process_image(im,direction="left")
  bb <- apply(p$outline_plot,2,range)
  crop <- c(bb[,1]+c(-15,15),bb[,2]+c(-15,15))
  pc <- process_image(im,crop=crop,direction="left")
  cat(basename(f),"dimensions",im$width,im$height,"tip",p$tip_plot,"base",p$base_plot,
      "distance",edge_distance(p$tip_plot,p$outline_plot),edge_distance(p$base_plot,p$outline_plot),
      "crop landmark max difference",max(abs(c(pc$tip_plot-p$tip_plot,pc$base_plot-p$base_plot))),"\n")
  plot(NA,xlim=c(0,im$width),ylim=c(im$height,0),asp=1,axes=FALSE,xlab="",ylab="",main=basename(f))
  rasterImage(im$raster,0,im$height,im$width,0); polygon(p$outline_plot,border="red")
  text(rbind(p$tip_plot,p$base_plot),c("T","B"),col="cyan",cex=1.5)
}
dev.off()

# Midline base at the second endpoint: a symmetric basal notch.
notch <- rbind(c(0,10),c(-3,0),c(-2,-3),c(0,-1),c(2,-3),c(3,0))
a <- standardize_outline(notch)
cat("Notch returned outline and base row\n"); print(a$outline); print(a$base_row)
cat("Notch consecutive zero segments",sum(rowSums(diff(a$outline)^2)==0),"\n")
man <- standardize_outline(notch,tip=c(0,10),base=c(0,-1))
cat("Manual vs auto point maximum error",max(abs(a$points-man$points)),"\n")
cat("Notch direct Fourier",attempt(efa_forward(a$outline,3)$x0),"\n")

# Independent quadrature of the linear polygon, with t=0 at vertex 1.
p <- ref$aligned[[1]]; q <- p[c(2:nrow(p),1),]; dt <- sqrt(rowSums((q-p)^2)); T <- sum(dt)
t0 <- c(0,head(cumsum(dt),-1)); nsub <- 1000
tt <- unlist(lapply(seq_along(dt),function(i) t0[i]+dt[i]*((1:nsub)-.5)/nsub))
xx <- do.call(rbind,lapply(seq_along(dt),function(i) outer(1-((1:nsub)-.5)/nsub,p[i,])+outer(((1:nsub)-.5)/nsub,q[i,])))
wt <- rep(dt/nsub,each=nsub)/T
quad <- c(sapply(1:20,function(h) 2*sum(wt*xx[,1]*cos(2*pi*h*tt/T))),
          sapply(1:20,function(h) 2*sum(wt*xx[,1]*sin(2*pi*h*tt/T))),
          sapply(1:20,function(h) 2*sum(wt*xx[,2]*cos(2*pi*h*tt/T))),
          sapply(1:20,function(h) 2*sum(wt*xx[,2]*sin(2*pi*h*tt/T))))
v <- efa_vector(efa_forward(p,20))
cat("Independent quadrature coefficient max error",max(abs(v-quad)),"\n")
suppressPackageStartupMessages(library(Momocs))
mom <- function(p) Momocs::efourier(Momocs::Out(list(p)),nb.h=20,norm=FALSE)$coe[1,]
cat("Momocs rolled error",max(abs(v-mom(p[c(2:nrow(p),1),]))),"unrolled error",max(abs(v-mom(p))),"\n")
g <- gpa(ref$points); rm <- apply(ref$gpagen,c(1,2),mean)
s <- svd(crossprod(g$mean,rm)); R <- s$u %*% diag(c(1,sign(det(s$u %*% t(s$v))))) %*% t(s$v)
cat("GPA geomorph maximum error",max(sapply(1:4,function(i) max(abs(g$aligned[[i]]%*%R-ref$gpagen[,,i])))),"\n")
x <- center_scale(p); y <- x; y[,1] <- -y[,1]; z <- rotate_onto(x,y)
cat("Reflection guard pair-distance max error",max(abs(as.matrix(dist(z))-as.matrix(dist(x)))),"\n")
if(exists(".Random.seed",.GlobalEnv)) rm(".Random.seed",envir=.GlobalEnv)
invisible(region_rv(ref$aligned,assign_regions(point_positions()),n_perm=9))
cat("RNG absent restored",!exists(".Random.seed",.GlobalEnv),"\n")
cat("RV scaled nonconstant block",rv_coefficient(matrix(1:20,4),matrix(1:20,4)),rv_coefficient(matrix(1:20,4)*1e-8,matrix(1:20,4)*1e-8),"\n")

# A numeric, finite, nonzero-length but collinear imported contour.
o <- cbind(x=c(0:4,4:0),y=0)
st <- list(new_specimen("line","line","",o,6))
f <- file.path(ev,"data/degenerate.csv"); write.csv(specimens_to_csv(st),f,row.names=FALSE)
loaded <- specimens_from_csv(f)
cat("Collinear CSV accepted count",length(loaded),"measurements\n"); print(attempt(specimen_measurements(loaded)))
for (case in c("huge","inf","zero_tip_base","empty_id")) {
  df <- specimens_to_csv(st)
  if(case=="huge") {df$x<-df$x*1e200; df$y<-c(1:5,5:1)*1e200}
  if(case=="inf") df$x[3]<-Inf
  if(case=="zero_tip_base") {df$x<-c(0,1,2,1,0,0,-1,-2,-1,0);df$y<-c(0,1,0,-1,-.5,0,1,0,-1,-.5)}
  if(case=="empty_id") df$id<-""
  write.csv(df,f,row.names=FALSE)
  ld <- attempt(specimens_from_csv(f))
  cat("CSV",case,"accepted",is.list(ld),"\n")
  if(is.list(ld)) {
    print(attempt(specimen_measurements(ld)))
    print(attempt({al<-gpa(rep(lapply(ld,function(s)s$points),4)); shape_pca(efa_matrix(al$aligned,6))$rank}))
    print(attempt(names(merge_specimens(list(),ld))))
  } else print(ld)
}

# Competing interior component of the opposite polarity.
gray <- matrix(.1,200,200); gray[20:90,20:90]<-.9
gray[110:190,110:190]<-.9; gray[112:188,112:188]<-.1
seg <- segment_image(gray,speck_side=0)
cat("Opposite polarity candidate warnings",paste(seg$warnings,collapse=" | "),"selected pixels",sum(seg$mask),"\n")

# Phase-free power of an actual sampled circle, not an exact synthetic vector.
theta <- 2*pi*(0:199)/200
pw <- harmonic_power(efa_matrix(list(cbind(cos(theta),sin(theta))),30))
cat("Sampled circle selected H",choose_harmonics(pw),"residual power",sum(pw$power[-1]),"\n")
disc<-function(cx,cy,r)outer(1:120,1:100,function(x,y)(x-cx)^2+(y-cy)^2<=r^2)
big<-disc(60,55,25); hole<-disc(60,55,6)
g<-matrix(.1,120,100);g[big]<-.9;g[hole]<-.1;g[disc(20,85,7)]<-.9;g[,1:8]<-.9
s<-segment_image(g)
cat("Disc note omitted cells",sum(big & !s$mask),"outside cells",sum(s$mask & !big),"unfilled hole cells",sum(hole & !s$mask),"\n")
