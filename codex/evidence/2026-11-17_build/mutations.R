root <- normalizePath("..")
ev <- normalizePath("evidence/2026-11-17_build")
testdir <- file.path(root,"tests/testthat")
library(testthat)
env <- new.env(parent=globalenv())
env$audit_test_path <- function(...) file.path(testdir,...)
eval(parse(text=gsub("testthat::test_path", "audit_test_path",readLines(file.path(testdir,"helper-app.R")),fixed=TRUE)),env)
for(f in list.files(file.path(root,"bifaceAnalysis/R"),full.names=TRUE)) sys.source(f,env)
run_mutation <- function(label,filename,mutate) {
  e <- new.env(parent=globalenv())
  for(n in ls(env,all.names=TRUE)) assign(n,get(n,env),e)
  # Rebind app functions so changed dependencies are used.
  for(n in ls(e)) if(is.function(e[[n]])) environment(e[[n]]) <- e
  mutate(e)
  for(n in ls(e,all.names=TRUE)) if(is.function(e[[n]])) environment(e[[n]]) <- e
  r <- test_file(file.path(testdir,filename),env=e,load_helpers=FALSE,reporter="summary")
  d <- as.data.frame(r)
  cat(label,"tests",nrow(d),"failures",sum(d$failed),"errors",sum(d$error),"\n")
  if(any(d$error)) print(d[d$error,c("test","result")])
}
run_mutation("stability plane always one","test-shape.R",function(e) {
  e$.original <- e$pca_stability
  e$pca_stability <- function(...) {z<-.original(...);if(!is.null(z))z$plane[]<-1;z}
})
run_mutation("PC2 never separable","test-shape.R",function(e) {
  e$.original <- e$pca_stability
  e$pca_stability <- function(...) {z<-.original(...);if(!is.null(z))z$separable[2]<-FALSE;z}
})
run_mutation("RV reference equals observation no permutations","test-shape.R",function(e) {
  e$.original <- e$region_rv
  e$region_rv <- function(...) {z<-.original(...);z$reference_mean<-z$reference_low<-z$reference_high<-z$rv;z}
})
run_mutation("ignore loaded metadata on identical contour","test-store.R",function(e) {
  e$.original <- e$merge_specimens
  e$merge_specimens <- function(store,loaded) {
    for(id in names(loaded)) if(!is.null(store[[id]]) && identical(store[[id]]$outline,loaded[[id]]$outline)) loaded[[id]] <- store[[id]]
    .original(store,loaded)
  }
})
run_mutation("manual landmarks reverse left and right edges","test-outline.R",function(e) {
  e$.original <- e$standardize_outline
  e$standardize_outline <- function(xy,direction="up",tip=NULL,base=NULL) {
    if(!is.null(tip)&&!is.null(base)) xy<-xy[nrow(xy):1,,drop=FALSE]
    .original(xy,direction,tip,base)
  }
})
run_mutation("disc returns smaller interior component","test-image.R",function(e) {
  e$.original<-e$segment_image
  e$segment_image<-function(gray,...) {
    z<-.original(gray,...)
    if(identical(dim(gray),c(120L,100L))) z$mask<-outer(1:120,1:100,function(x,y)(x-20)^2+(y-85)^2<=7^2)
    z
  }
})
run_mutation("no-op symmetrization","test-shape.R",function(e) e$symmetrize_shape<-function(m,pair)center_scale(m))
