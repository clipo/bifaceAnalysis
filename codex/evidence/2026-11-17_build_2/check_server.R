# Run from repository root. Synthetic outlines only, no real collections.
library(shiny)
APP_DIR <- normalizePath('bifaceAnalysis')
OUT <- normalizePath('codex/evidence/2026-11-17_build_2/data')
for (f in list.files(file.path(APP_DIR,'R'),pattern='[.]R$',full.names=TRUE)) source(f)
# Reuse only test infrastructure definitions, without running repository tests.
exprs <- parse('tests/testthat/test-server.R')
for (e in exprs) if (is.call(e) && identical(e[[1]],as.name('<-')) &&
  as.character(e[[2]]) %in% c('defaults','run_app_test','pdf_pages','try_figure','try_text')) eval(e)
is_validation <- function(e) inherits(e,'shiny.silent.error') || inherits(e,'validation')
eq <- function(a,b,tol=1e-11) stopifnot(isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE)))
make_s <- function(id,group,width=3,length=10,ppm=10) {
  theta<-seq(0,2*pi,length.out=201)[-201]
  # Valid CCW outline with tip in row 1 and base in row 101.
  o<-cbind(-width*sin(theta),length*cos(theta))
  new_specimen(id,id,group,o,101L,ppm)
}
put_csv<-function(loc) { f<-tempfile(fileext='.csv');write.csv(loc,f,row.names=FALSE);data.frame(name='locations.csv',datapath=f) }
copy_pdf<-function(f,name) {
  dest<-file.path(OUT,paste0(name,'.pdf'));file.copy(f,dest,overwrite=TRUE)
  system2('/opt/homebrew/bin/pdftotext',c('-layout',dest,file.path(OUT,paste0(name,'.txt'))))
  cat('PDF',name,'pages=',pdf_pages(dest),'\n')
  readLines(file.path(OUT,paste0(name,'.txt')),warn=FALSE)
}
run_app_test({
  do.call(session$setInputs,defaults)
  for(d in c('dl_coll_summary','dl_coll_dist')) {
    lines<-readLines(output[[d]]);stopifnot(length(lines)==1)
    cat('Empty CSV',d,':',lines,'\n')
  }
  for(rank in c(0,1)) {
    rv$specimens<-setNames(lapply(1:6,function(i)make_s(paste0('s',i),if(i<=3)'A' else 'B',if(i<=3 || rank==0)3 else 4)),paste0('s',1:6))
    session$setInputs(map_value='pc3_mean')
    stopifnot(pca()$rank==rank,map_value()=='length_width_mean')
    sm<-coll_summary();stopifnot(all(is.na(sm$pc2_mean)),all(is.na(sm$pc3_mean)))
    if(rank==0)stopifnot(all(is.na(sm$pc1_mean)))
    session$setInputs(locations_file=put_csv(data.frame(collection=c('A','B'),latitude=c(42,42),longitude=c(-76,-76))))
    for(k in names(figures)[startsWith(names(figures),'coll')]) {
      res<-try_figure(figures[[k]]);stopifnot(!startsWith(res,'ERROR'));cat('Rank',rank,k,res,'\n')
    }
    for(k in names(texts)[startsWith(names(texts),'coll')])stopifnot(!startsWith(try_text(texts[[k]]),'ERROR'))
    cat('Rank',rank,'within variances:',coll_stats()$within$variance,' map choices:',map_choices(),'\n')
    for(k in c('dl_coll_summary','dl_coll_dist','dl_loc_template'))stopifnot(file.exists(output[[k]]))
    cat('Rank',rank,'summary CSV:',paste(readLines(output$dl_coll_summary),collapse=' | '),'\n')
  }
  # Interleave eligible and excluded rows. Names exercise quoted CSV and display text.
  labels<-c('NA','007','A,B','Q"R','<A>&','P%','G7','G8','G9')
  sp<-list();k<-0
  for(j in 1:3)for(i in 1:9) {
    k<-k+1;sp[[paste0('n',k)]]<-make_s(paste0('n',k),labels[i],2+i/10+j/50,10+j/20,if(k==2)NA_real_ else 10)
    if(i==4){k<-k+1;sp[[paste0('n',k)]]<-make_s(paste0('n',k),'',7,13)}
  }
  sp[['small']]<-make_s('small','small',1,8)
  rv$specimens<-sp
  session$setInputs(locations_file=put_csv(data.frame(collection=rev(labels),latitude=30+(9:1),longitude=-80+(9:1)^1.1)))
  cd<-coll();eq(cd$X,coefs()[cd$rows,,drop=FALSE]);stopifnot(ncol(cd$X)==4*n_harm())
  sm<-read.csv(output$dl_coll_summary,colClasses='character',na.strings='',check.names=FALSE)
  m<-measurements(); scores<-pca()$scores
  for(i in 1:nrow(sm)) {
    rows<-which(groups()==sm$collection[i]);stopifnot(length(rows)==as.integer(sm$n[i]))
    for(j in 1:min(3,ncol(scores)))eq(as.numeric(sm[[paste0('pc',j,'_mean')]][i]),mean(scores[rows,j]))
    eq(as.numeric(sm$length_width_mean[i]),mean(m$length_width_ratio[rows]));eq(as.numeric(sm$max_width_position_mean[i]),mean(m$widest_position[rows]))
    eq(as.numeric(sm$within_variance[i]),sum(apply(coefs()[rows,,drop=FALSE],2,var)))
  }
  stopifnot(all(is.na(sm$length_mm_mean)),!('length_mm_mean'%in%map_choices()))
  stopifnot(setequal(sm$collection,labels),all(sm$located=='TRUE'))
  di<-read.csv(output$dl_coll_dist,colClasses='character',na.strings='',check.names=FALSE)
  stopifnot(nrow(di)==36,identical(names(di),c('collection_a','collection_b','n_a','n_b','raw_distance','reference_mean','reference_low','reference_high','corrected_square','corrected_distance','distance_km')))
  for(k in names(figures)[startsWith(names(figures),'coll')]) {
    res<-try_figure(figures[[k]]);stopifnot(res=='plot');cat('Nine special-name collections',k,res,'\n')
  }
  cat('Interleaved full-matrix/PC/measurement/download alignment verified; eligible rows=',paste(cd$rows,collapse=','),' H=',n_harm(),' rank=',pca()$rank,'\n')
  cat('Random geography sentence:',texts$coll_geo()[2],'\n')
  # No names match: all relevant outputs must remain available or validate cleanly.
  session$setInputs(locations_file=put_csv(data.frame(collection='unknown',latitude=1,longitude=2)))
  stopifnot(nrow(coll_loc()$located)==0,all(!coll_summary()$located),all(is.na(coll_pairs_km()$distance_km)))
  stopifnot(try_figure(figures$coll_map)=='message',try_figure(figures$coll_geo)=='message')
  cat('No matching names:',paste(texts$coll_map(),collapse=' | '),'\n')
  # Exact repeated means with different specimen order, all corrected distances zero.
  rv$specimens<-setNames(lapply(1:6,function(i)make_s(paste0('z',i),if(i<=3)'A' else 'B',3+((i-1)%%3)/5)),paste0('z',1:6))
  cat('Identical collection means corrected distance=',coll_stats()$pairs$corrected_distance,'\n')
  stopifnot(try_figure(figures$coll_pairs)=='plot',try_figure(figures$coll_matrix)=='plot')
  session$setInputs(locations_file=put_csv(data.frame(collection=c('A','B'),latitude=c(42,43),longitude=c(-76,-77))),map_in_report=TRUE)
  a<-copy_pdf(output$dl_report,'with-map')
  stopifnot(any(grepl('Map of collections',a,fixed=TRUE)))
  session$setInputs(map_in_report=FALSE)
  b<-copy_pdf(output$dl_report,'without-map')
  stopifnot(!any(grepl('Map of collections',b,fixed=TRUE)))
  cat('Map-off report has geographic-coordinate caution:',any(grepl('approximate coordinates',b,fixed=TRUE)),' sample-mean map caution:',any(grepl('colours show means',b,fixed=TRUE)),'\n')
  session$setInputs(clear_locations=1)
  c<-copy_pdf(output$dl_report,'no-locations')
  cat('No-location report has geographic-coordinate caution:',any(grepl('approximate coordinates',c,fixed=TRUE)),'\n')
  # Actual location rejection leaves previous table unchanged, tested by repository suite too.
  session$setInputs(locations_file=put_csv(data.frame(collection=c('A','B'),latitude=c(81,79),longitude=c(10,20))))
  stopifnot(try_figure(figures$coll_map)=='message')
  cat('Withheld map:',coll_window()$reason,'\n')
})
set.seed(20261117)
X<-matrix(rnorm(120*120),120);g<-factor(rep(letters[1:6],each=20))
cat('Performance 120 x 120 share + pair + within elapsed:',system.time({share_reference(X,g);pair_table(X,g);within_table(X,g)})[['elapsed']],'seconds\n')
cat('All additional server checks completed.\n')
