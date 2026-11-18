library(shiny)
APP_DIR<-normalizePath('bifaceAnalysis')
for(f in list.files(file.path(APP_DIR,'R'),pattern='[.]R$',full.names=TRUE))source(f)
exprs<-parse('tests/testthat/test-server.R')
for(e in exprs)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])%in%c('defaults','run_app_test'))eval(e)
options(digits=17)
run_app_test({
 do.call(session$setInputs,defaults)
 for(rank in c(0,1,2)) {
  rv$specimens<-setNames(lapply(1:6,function(i) {
   t<-seq(0,2*pi,length.out=201)[-201];width<-if(rank==2)3+((i-1)%%3)/5 else if(i<=3||rank==0)3 else 4
   new_specimen(paste0('s',i),paste0('s',i),if(i<=3)'A' else 'B',cbind(-width*sin(t),10*cos(t)),101L,10)
  }),paste0('s',1:6))
  cat('PCA rank:',pca()$rank,'\n')
  cat('Observed share:',coll_stats()$share$share,'\n')
  for(g in c('A','B')) {
   x<-coll()$X[coll()$group==g,,drop=FALSE]
   cat('Group',g,'max absolute difference from its first row:',max(abs(sweep(x,2,x[1,]))),'no_variation:',no_variation(x),'row_variance:',row_variance(x),'\n')
  }
  cat('Full text:\n',paste(texts$coll(),collapse='\n'),'\n')
  pdf(tempfile(fileext='.pdf'))
  for(k in c('coll_means','coll_share','coll_pairs','coll_matrix','coll_within')) {
   withCallingHandlers(tryCatch({ggplot2::ggplotGrob(figures[[k]]());cat('Rendered',k,'\n')},error=function(e)cat('Message/error',k,conditionMessage(e),'\n')),
    warning=function(w){cat('DRAW WARNING',k,conditionMessage(w),'\n');invokeRestart('muffleWarning')})
  }
  dev.off()
 }
})
