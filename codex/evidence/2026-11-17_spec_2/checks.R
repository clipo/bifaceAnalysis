# Independent specification checks. Run from codex/ with Rscript.
options(digits = 12)
cat("Audited HEAD:", system("git rev-parse HEAD", intern = TRUE), "\n")
cat(R.version.string, "\n")
for (p in c("ggplot2", "maps", "shiny", "testthat"))
  cat(p, if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "NOT INSTALLED", "\n")
source("../bifaceAnalysis/R/shape.R")
source("../bifaceAnalysis/R/efa.R")
source("../bifaceAnalysis/R/plots.R")

part <- function(x, g) {
  c0 <- colMeans(x)
  ss <- vapply(split(seq_len(nrow(x)), g), function(i) {
    m <- colMeans(x[i, , drop = FALSE])
    c(between = length(i) * sum((m - c0)^2),
      within = sum(sweep(x[i, , drop = FALSE], 2, m)^2))
  }, numeric(2))
  c(SST = sum(sweep(x, 2, c0)^2), SSB = sum(ss[1, ]), SSW = sum(ss[2, ]))
}
set.seed(117)
x <- matrix(rexp(15 * 12), 15, 12)
for (sizes in list(c(5, 5, 5), c(3, 5, 7))) {
  g <- rep(seq_along(sizes), sizes)
  pp <- part(x, g)
  draws <- replicate(50000L, { z <- part(x, sample(g)); z["SSB"] / z["SST"] })
  cat("\nPartition sizes", sizes, "\n")
  print(pp)
  cat("SSB+SSW-SST", pp["SSB"] + pp["SSW"] - pp["SST"], "\n")
  cat("Expected", 2/14, "simulation mean", mean(draws), "MCSE", sd(draws)/sqrt(length(draws)), "\n")
  fit <- lm(x ~ factor(g))
  cat("Pooled lm share", 1 - sum(residuals(fit)^2) / pp["SST"], "\n")
}
# Exhaustive labelled partitions into three groups of three.
z <- matrix(seq_len(36)^2, 9, 4)
shares <- numeric(0)
for (a in combn(9, 3, simplify = FALSE)) {
  for (b in combn(setdiff(seq_len(9), a), 3, simplify = FALSE)) {
    g <- rep(3L, 9); g[a] <- 1L; g[b] <- 2L
    p <- part(z, g); shares <- c(shares, p["SSB"] / p["SST"])
  }
}
cat("Exhaustive partition count", length(shares), "mean", mean(shares), "expected", 2/8, "\n")

variance <- function(x) sum(sweep(x, 2, colMeans(x))^2)/(nrow(x)-1)
set.seed(118)
for (delta in c(2, 0)) {
  sim <- replicate(4000L, {
    a <- matrix(rnorm(5*12), 5, 12)
    b <- matrix(rnorm(8*12), 8, 12); b[,1] <- b[,1]+delta
    raw2 <- sum((colMeans(a)-colMeans(b))^2)
    q <- raw2 - variance(a)/5 - variance(b)/8
    c(raw2 = raw2, q = q, clipped = max(0,q), dc = sqrt(max(0,q)))
  })
  cat("\nIndependent N(mu,I12), n=5,8, delta", delta, "draws", ncol(sim), "\n")
  print(data.frame(mean = rowMeans(sim), MCSE = apply(sim,1,sd)/sqrt(ncol(sim))))
  cat("Negative Q fraction", mean(sim["q",]<0), "raw2 expectation", delta^2+12*(1/5+1/8), "\n")
}
set.seed(119)
for (n in c(3,5,8,20)) {
  v <- replicate(20000, variance(matrix(rnorm(n*12), n,12)))
  cat("Variance n",n,"mean",mean(v),"MCSE",sd(v)/sqrt(length(v)),"expected",12,"\n")
}
hand <- matrix(c(0,1,2,7), ncol=1)
loo <- vapply(1:4, function(i) variance(hand[-i,,drop=FALSE]), numeric(1))
cat("LOO example rows 0,1,2,7; variance",variance(hand),"values",loo,"mean",mean(loo),"\n")

hav <- function(lat1,lon1,lat2,lon2) {
  r <- pi/180
  h <- sin((lat2-lat1)*r/2)^2 + cos(lat1*r)*cos(lat2*r)*sin((lon2-lon1)*r/2)^2
  2*6371.0088*asin(pmin(1,sqrt(h)))
}
cat("\nHaversine, radius 6371.0088 km\n")
print(c(one_degree=hav(0,0,0,1), equator_pole=hav(0,0,90,0), same=hav(35,-80,35,-80), meridian=hav(0,179.5,0,-179.5)))
cat("Pole expected error vs 10007.543:",hav(0,0,90,0)-10007.543,"\n")
perms <- function(v) if (length(v)==1) matrix(v,1) else do.call(rbind,lapply(v,function(i) cbind(i,perms(v[v!=i]))))
P <- perms(1:5)
lat <- c(30,31,35,40,43); lon <- c(-100,-90,-85,-80,-75)
D <- outer(1:5,1:5,Vectorize(function(i,j) hav(lat[i],lon[i],lat[j],lon[j])))
lt <- lower.tri(D)
ref <- apply(P,1,function(p) cor(D[lt],D[p,p][lt]))
cat("Matrix reference assignments",nrow(P),"identity count",sum(apply(P,1,function(p) identical(as.integer(p),1:5))),
    "observed",cor(D[lt],D[lt]),"mean",mean(ref),"maximum",max(ref),"\n")
cat("Reference quantiles",quantile(ref,c(.025,.975)),"\n")
pos <- c(1,1,2,2,2)
keys <- apply(P,1,function(p) paste(pos[p],collapse=","))
cat("Coincident positions: 5 labelled rows, multiplicities 2 and 3; assignments",nrow(P),
    "distinct coordinate arrangements",length(unique(keys)),"weights",unique(as.integer(table(keys))),"\n")

cat("\nMinimum map window (minimum final latitude span 0.5 degrees)\n")
for (lat in c(0,45,89.9,90))
  cat("middle latitude",lat,"latitude limits",lat+c(-.25,.25),"longitude span",.5/cos(lat*pi/180),"\n")
cat("Quarter-width scale at latitudes 45 and 60 for same longitude span:",cos(60*pi/180)/cos(45*pi/180),"ratio\n")
# A full polygon retained and clipped by the coordinates builds in ggplot2.
# This is a synthetic geometry check, not a check of absent maps databases.
poly <- data.frame(long=c(-130,-60,-60,-130,-130),lat=c(20,20,55,55,20),group=1)
pg <- ggplot2::ggplot(poly,ggplot2::aes(long,lat,group=group)) + ggplot2::geom_polygon() +
  ggplot2::coord_quickmap(xlim=c(-125.2,-124.8),ylim=c(39.8,40.2),expand=FALSE)
gb <- ggplot2::ggplot_build(pg)
cat("Synthetic quickmap polygon: retained vertices",nrow(gb$data[[1]]),"aspect",pg$coordinates$aspect(gb$layout$panel_params[[1]]),"\n")
cat("category_scale nine levels unique colours",length(unique(category_scale(letters[1:9])$palette(9))),"\n")
cat("category_scale fill aesthetic:",category_scale(letters[1:9],"fill")$aesthetics,"\n")
cat("Preserved groups vs trimmed collections:\n")
labels <- c("A"," A "," ","", "NA")
print(data.frame(stored=labels, collection=trimws(labels), nonempty_stored=nzchar(labels), nonempty_trimmed=nzchar(trimws(labels))))

# Synthetic shapes demonstrate an automatic-H change with an unlabelled row.
theta <- 2*pi*(0:199)/200
shapes <- lapply(seq(.30,.36,length.out=6), function(w) cbind(w*sin(theta),cos(theta)))
names(shapes) <- paste0("s",1:6)
added <- cbind(.35*sin(theta),cos(theta))*(1+.20*cos(14*theta))
Hof <- function(s) choose_harmonics(harmonic_power(efa_matrix(gpa(s)$aligned,30)))
cat("Synthetic eligible set automatic H",Hof(shapes),"with unlabelled lobed shape",Hof(c(shapes,list(extra=added))),"\n")
cat("Constant coefficient matrix rank",shape_pca(matrix(rep(1:8,each=6),6,8))$rank,
    "rank-one matrix rank",shape_pca(outer(1:6,1:8))$rank,"\n")
cat("\nNo app implementation or project test suite executed. No maps data rendered.\n")
