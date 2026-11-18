# Run from repository root. All fixtures are synthetic.
for (f in c('shape', 'efa', 'collections', 'plots', 'interpret')) source(paste0('bifaceAnalysis/R/', f, '.R'))
options(digits = 16)
cat('HEAD:', system('git rev-parse HEAD', intern = TRUE), '\n', R.version.string, '\n')
for (p in c('testthat','shiny','bslib','ggplot2','maps','scales','imager','cluster')) cat(p, as.character(packageVersion(p)), '\n')
eq <- function(x,y,tol=1e-12) stopifnot(max(abs(x-y)) < tol)
X <- cbind(c(0,2,4,5,7,9,3,6,12), c(1,0,2,4,3,6,8,2,5))
g <- factor(rep(c('A','B','C'), each=3))
mu <- t(vapply(split(seq_len(nrow(X)),g),function(i) apply(X[i,,drop=FALSE],2,mean),numeric(2)))
centre <- apply(X,2,mean)
sst <- sum(vapply(seq_len(nrow(X)),function(i)sum((X[i,]-centre)^2),numeric(1)))
ssb <- sum(vapply(1:3,function(j)3*sum((mu[j,]-centre)^2),numeric(1)))
ssw <- sum(vapply(1:9,function(i)sum((X[i,]-mu[as.integer(g[i]),])^2),numeric(1)))
eq(sst,ssb+ssw); eq(share_between(X,g),ssb/sst)
cat('Independent partition: SST=',sst,' SSB=',ssb,' SSW=',ssw,' share=',ssb/sst,'\n')
oracle_share <- function(labels) {
  sum(vapply(levels(g),function(k) { z<-X[labels==k,,drop=FALSE]; nrow(z)*sum((colMeans(z)-centre)^2) },numeric(1)))/sst
}
set.seed(1); expected <- replicate(999,oracle_share(g[sample.int(9)]))
sr <- share_reference(X,g); eq(sr$reference$values,expected)
cat('Share random reference mean/low/high:', sr$reference$mean,sr$reference$low,sr$reference$high,' analytical expectation=',sr$expected,'\n')
pr <- pair_table(X,g)
for(i in seq_len(nrow(pr))) {
  A<-X[g==pr$collection_a[i],,drop=FALSE]; B<-X[g==pr$collection_b[i],,drop=FALSE]
  va<-sum(apply(A,2,var)); vb<-sum(apply(B,2,var))
  d2<-sum((apply(A,2,mean)-apply(B,2,mean))^2); Q<-d2-va/nrow(A)-vb/nrow(B)
  eq(pr$raw_distance[i],sqrt(d2)); eq(pr$corrected_square[i],Q); eq(pr$corrected_distance[i],sqrt(max(0,Q)))
  cat('Independent pair',pr$collection_a[i],pr$collection_b[i],': d2=',d2,' vA=',va,' vB=',vb,' Q=',Q,' dc=',sqrt(max(0,Q)),'\n')
}
# Independently enumerate combinations within the first pair (20 assignments).
pool<-X[1:6,]; all_d<-apply(combn(6,3),2,function(i)sqrt(sum((colMeans(pool[i,])-colMeans(pool[-i,]))^2)))
set.seed(1); rnd<-replicate(999,{ ii<-sample.int(6); sqrt(sum((colMeans(pool[ii[1:3],])-colMeans(pool[ii[4:6],]))^2)) })
eq(pr$reference_mean[1],mean(rnd)); eq(pr$reference_low[1],quantile(rnd,.025)); eq(pr$reference_high[1],quantile(rnd,.975))
stopifnot(all(vapply(rnd,function(z)any(abs(z-all_d)<1e-12),logical(1))))
cat('Pair reference independently reproduced; enumeration max=',max(all_d),'\n')

pos<-cbind(c(0,1,3,4,8,9,12,15),c(0,3,1,7,2,8,9,3))
# Great-circle oracle uses 3D unit vectors and atan2 cross/dot, not haversine.
unit<-cbind(cos(pos[,1]*pi/180)*cos(pos[,2]*pi/180),cos(pos[,1]*pi/180)*sin(pos[,2]*pi/180),sin(pos[,1]*pi/180))
geo<-outer(1:8,1:8,Vectorize(function(i,j) {
  a<-unit[i,]; b<-unit[j,]; cross<-c(a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
  6371.0088*atan2(sqrt(sum(cross^2)),sum(a*b))
}))
eq(geo,geo_distance_matrix(pos[,1],pos[,2],letters[1:8]),1e-9)
sh<-as.matrix(dist(cbind(c(0,2,1,5,3,8,6,4),c(2,0,4,1,3,7,5,9))))
for(n in c(5,7,8)) {
  a<-sh[1:n,1:n]; b<-geo[1:n,1:n]; p<-pos[1:n,]
  result<-shape_geo_correlation(a,b,p)
  eq(result$r,cor(a[lower.tri(a)],b[lower.tri(b)]))
  if(n==5) {
    perms<-as.matrix(expand.grid(rep(list(1:n),n)))
    perms<-perms[apply(perms,1,function(z)length(unique(z))==n),]
    oracle<-apply(perms,1,function(z)cor(a[lower.tri(a)],b[z,z][lower.tri(b)]))
    eq(sort(oracle),sort(result$reference$values))
  }
  if(n==8) {
    set.seed(1); oracle<-replicate(999,{z<-sample.int(n);cor(a[lower.tri(a)],b[z,z][lower.tri(b)])})
    eq(oracle,result$reference$values)
  }
  cat('Correlation G=',n,' r=',result$r,' assignments=',result$n_assignments,' ref mean/low/high=',result$reference$mean,result$reference$low,result$reference$high,'\n')
  cat(paste(text_collection_geo(result),collapse='\n'),'\n')
}

# Seed presence, absence, repeatability and early-return/error paths.
calls<-list(share=function()share_reference(X,g),pair=function()pair_table(X,g),
  geo5=function()shape_geo_correlation(sh[1:5,1:5],geo[1:5,1:5],pos[1:5,]),
  geo7=function()shape_geo_correlation(sh[1:7,1:7],geo[1:7,1:7],pos[1:7,]),
  geo8=function()shape_geo_correlation(sh,geo,pos),
  zero_share=function()share_reference(matrix(1,9,2),g),
  too_few=function()shape_geo_correlation(sh[1:2,1:2],geo[1:2,1:2],pos[1:2,]),
  coincident=function()shape_geo_correlation(sh,geo,matrix(0,8,2)),
  zero_shape=function()shape_geo_correlation(matrix(0,8,8),geo,pos),
  error=function()tryCatch(with_seed(1,stop('deliberate')),error=function(e)conditionMessage(e)))
for(k in names(calls)) {
  set.seed(864); old<-.Random.seed; a<-calls[[k]](); stopifnot(identical(old,.Random.seed)); b<-calls[[k]]();stopifnot(identical(a,b))
  rm('.Random.seed',envir=.GlobalEnv); calls[[k]](); stopifnot(!exists('.Random.seed',envir=.GlobalEnv))
  cat('RNG preserved present/absent; repeatable:',k,'\n')
}

fixtures<-list(coincident=rbind(c(42,-76),c(42,-76)), lon_only=rbind(c(42,-80),c(42,-70)),
 lat_only=rbind(c(40,-76),c(44,-76)), high=rbind(c(79,10),c(79.5,20)),
 clipping=rbind(c(-80,-179),c(80,-170)), outside=rbind(c(50,5),c(52,8)),
 straddle=rbind(c(48,-130),c(52,-120)))
for(k in names(fixtures)) {
 z<-fixtures[[k]]; mid<-colMeans(apply(z,2,range)); span<-pmax(1.4*apply(z,2,function(x)diff(range(x))),c(.5,.5/cos(mid[1]*pi/180)))
 bounds<-cbind(pmax(mid-span/2,c(-85,-180)),pmin(mid+span/2,c(85,180)))
 w<-map_window(z[,1],z[,2]);eq(w$ylim,bounds[1,]);eq(w$xlim,bounds[2,])
 bar<-scale_bar(w); target<-haversine_km(w$lat_mid,w$xlim[1],w$lat_mid,w$xlim[2])/4
 candidates<-as.vector(outer(c(1,2,5),10^(-8:8)))
 eq(bar$km,candidates[which.min(abs(candidates-target))]);eq(KM_PER_DEGREE*cos(bar$y*pi/180)*(bar$x1-bar$x0),bar$km,1e-8)
 loc<-data.frame(collection=c('A','B'),latitude=z[,1],longitude=z[,2])
 p<-plot_collection_map(loc,1:2,'n',w);b<-ggplot2::ggplot_build(p)
 eq(b$layout$panel_params[[1]]$x.range,w$xlim,1e-9);eq(b$layout$panel_params[[1]]$y.range,w$ylim,1e-9)
 cat('Map fixture',k,' bounds=',as.vector(t(bounds)),' backgrounds=',map_layers(w),' bar km=',bar$km,'\n')
}
delta_lon<-1000/(KM_PER_DEGREE*cos(pi/3))
cat('Build note (a): 1000 km along parallel 60 degrees; endpoint haversine=',haversine_km(60,0,60,delta_lon),' discrepancy percent=',100*(1-haversine_km(60,0,60,delta_lon)/1000),'\n')

# Labels are arbitrary strings. Two distinct pairs can get the same pasted key.
gg<-factor(rep(c('A','B and C','A and B','C'),each=3))
xx<-matrix(1:48,12,4); pp<-pair_table(xx,gg,n_perm=9)
bp<-ggplot2::ggplot_build(plot_pair_distances(pp))
cat('Pair label collision: pairs=',nrow(pp),' distinct plotted y rows=',length(unique(bp$data[[3]]$y)),' labels=',paste(paste(pp$collection_a,'and',pp$collection_b),collapse=' | '),'\n')

# In-memory mutation: zero pair-reference summaries pass every collection test.
# Evaluate test expressions only, retaining synthetic fixtures and required helpers.
test_env<-new.env(parent=.GlobalEnv)
library(testthat)
test_env$APP_DIR<-normalizePath('bifaceAnalysis')
test_env$DEMO_FILES<-list.files('bifaceAnalysis/www',pattern='(?i)\\.jpe?g$',full.names=TRUE)
names(test_env$DEMO_FILES)<-tools::file_path_sans_ext(basename(test_env$DEMO_FILES))
# App dependencies are needed for its existing label/outline test.
for(f in c('image','outline','store')) sys.source(paste0('bifaceAnalysis/R/',f,'.R'),envir=test_env)
test_env$KITE<-rbind(c(0,10),c(-3,2),c(-1.2,1.5),c(-1.2,-2),c(1.2,-2),c(1.2,1.5),c(3,2))
test_env$pair_table<-function(...) { p<-pair_table(...);p$reference_mean<-p$reference_low<-p$reference_high<-0;p }
cat('BEGIN MUTATION: pair references all zero\n')
testthat::with_reporter('summary',sys.source('tests/testthat/test-collections.R',envir=test_env))
cat('END MUTATION\nAll independent numeric checks completed.\n')
