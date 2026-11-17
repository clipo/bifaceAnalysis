for(f in list.files("../bifaceAnalysis/R",full.names=TRUE)) source(f)
attempt<-function(e)tryCatch(e,error=function(e)conditionMessage(e))
poly <- function(L) rbind(c(-L,1),c(-L,-1),c(1,-1),c(1,-3),c(0,-4),c(2,-3),c(2,3),c(0,4),c(1,3),c(1,1))
L <- uniroot(function(L)polygon_moments(poly(L))$centroid[1],c(1,10),tol=1e-14)$root
o<-poly(L)
cat("L",format(L,digits=17),"centroid",format(polygon_moments(o)$centroid,digits=17),"\n")
for(k in 0:(nrow(o)-1)) {
  z<-o[c((k+1):nrow(o),seq_len(k)),]; a<-standardize_outline(z)
  cat("cyclic shift",k,"base",a$base_input,"\n")
}
for (eps in c(0,-1e-13,1e-13)) {
  a<-standardize_outline(poly(L+eps)); cat("base",a$base_input,"zero/nearzero segments",sum(sqrt(rowSums(diff(a$outline)^2))<1e-12),"\n")
  print(a$outline)
}
# Direct half-open crossing second-endpoint counterexample.
print(vertical_crossings(rbind(c(1,1),c(0,0),c(1,-1),c(2,0)),0))

# Crop coordinates are tested with a fixed threshold to isolate translations.
img<-load_working_image("../bifaceAnalysis/www/artifact1.jpg")
p<-process_image(img,direction="left")
bb<-apply(p$outline_plot,2,range)
q<-process_image(img,crop=c(max(0,bb[1,1]-10),min(img$width,bb[2,1]+10),max(0,bb[1,2]-10),min(img$height,bb[2,2]+10)),threshold=p$threshold,direction="left")
cat("Fixed-threshold crop marker max error",max(abs(c(q$tip_plot-p$tip_plot,q$base_plot-p$base_plot))),"\n")
cat("Stale crop beyond narrower image",attempt(process_image(load_working_image("../bifaceAnalysis/www/artifact4.jpg"),crop=c(160,240,70,170))),"\n")
# Tiny nonconstant coordinate blocks are well-defined under the stated RV formula.
x<-matrix(1:20,4); cat("RV rescaling",rv_coefficient(x,x),rv_coefficient(x*1e-8,x*1e-8),"\n")
