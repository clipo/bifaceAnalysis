# Checks run by the author's assistant against Codex findings F1 to F12 of the
# 2026-11-17 spec audit. Run from the repository root:
#   Rscript docs/audits/2026-11-17_spec_checks.R > docs/audits/2026-11-17_spec_checks.log 2>&1
suppressMessages({library(imager); library(Momocs)})
border <- function(m) { nx <- nrow(m); ny <- ncol(m); c(m[1, ], m[nx, ], m[-c(1, nx), 1], m[-c(1, nx), ny]) }
otsu <- function(v, bins = 256L) {
  br <- seq(0, 1, length.out = bins + 1L); ct <- tabulate(findInterval(v, br, rightmost.closed = TRUE), bins)
  p <- ct / sum(ct); mids <- (br[-1] + br[-(bins + 1L)]) / 2; w0 <- cumsum(p); m0 <- cumsum(p * mids); mt <- m0[bins]
  bt <- ifelse(w0 > 0 & w0 < 1, (mt * w0 - m0)^2 / (w0 * (1 - w0)), NA); br[median(which(bt == max(bt, na.rm = TRUE))) + 1L]
}
cat("== F1: spec 6.1 step 4 (border-majority polarity) on the four demo images\n")
for (f in list.files("bifaceAnalysis/www", full.names = TRUE)) {
  g <- as.matrix(load.image(f)[, , 1, 1]); above <- g > otsu(g)
  light_border <- mean(border(above))
  fg <- if (light_border > 0.5) !above else above       # majority border class is background
  cat(sprintf("%s: fraction of border pixels above threshold %.3f -> artifact taken as %s; foreground touches border: %s; foreground share of image %.2f\n",
              basename(f), light_border, if (light_border > 0.5) "dark" else "light", any(border(fg)), mean(fg)))
}
cat("\n== F2: 3-by-1 ellipse, 2000 vertices equally spaced in angle, Momocs efourier norm = FALSE\n")
th <- 2 * pi * (0:1999) / 2000
e <- efourier(Out(list(cbind(3 * cos(th), sin(th)))), nb.h = 5, norm = FALSE)$coe[1, ]
print(round(e[c("A1", "D1", "A3", "D3")], 6))
cat("circle radius 2:\n"); print(round(efourier(Out(list(cbind(2 * cos(th), 2 * sin(th)))), nb.h = 3, norm = FALSE)$coe[1, ], 6))
cat("\n== F3: arc-length mean of a vertex-centred triangle is not zero\n")
tri <- rbind(c(0, 0), c(4, 0), c(0, 3)); tri <- sweep(tri, 2, colMeans(tri))
nx <- tri[c(2, 3, 1), ]; dt <- sqrt(rowSums((nx - tri)^2))
print(colSums(dt * (tri + nx) / 2) / sum(dt))
cat("\n== F5: 100 points per edge including both ends of each edge\n")
sq <- rbind(c(0, 1), c(-1, 0), c(0, -1), c(1, 0)); left <- sq[1:3, ]; right <- sq[c(3, 4, 1), ]
rs <- function(p, fr) { s <- c(0, cumsum(sqrt(rowSums(diff(p)^2)))); cbind(approx(s, p[, 1], fr * max(s))$y, approx(s, p[, 2], fr * max(s))$y) }
both <- rbind(rs(left, seq(0, 1, length.out = 100)), rs(right, seq(0, 1, length.out = 100)))
cat("rows:", nrow(both), " distinct:", nrow(unique(round(both, 10))), "\n")
cat("\n== F6: order-swap fixture\n")
x <- rbind(c(1.01, 0), c(-1.01, 0), c(0, 1), c(0, -1))
V <- prcomp(x)$rotation; W <- prcomp(x[-1, ])$rotation
cat("same-rank |cosine| PC1, PC2:", round(abs(colSums(V * W)), 6), " principal-angle cosines of the planes:", round(svd(crossprod(V, W))$d, 6), "\n")
cat("components available with 3 specimens:", sum(prcomp(matrix(rnorm(30), 3))$sdev > 1e-12), "\n")
cat("\n== F10: silhouette conventions in cluster::silhouette\n")
d <- dist(c(0, 0, 1, 5)); print(unclass(cluster::silhouette(c(1, 1, 1, 2), d))[, 3]); print(unclass(cluster::silhouette(c(1, 1, 2, 2), dist(c(0, 0, 3, 3))))[, 3])
cat("\n== F12: current pipeline (Out -> align -> efourier -> PCA -> hclust -> cutree) at small n\n")
data(bot); run <- function(n, k) tryCatch({
  o <- Out(bot$coo[seq_len(n)]) %>% coo_center() %>% coo_scale() %>% coo_align()
  ef <- suppressMessages(efourier(o, nb.h = 10, norm = TRUE)); p <- PCA(ef); cl <- cutree(hclust(dist(ef$coe), method = "ward.D2"), k = k)
  sprintf("n=%d k=%d: ran; PCs %d", n, k, ncol(p$x)) }, error = function(e) sprintf("n=%d k=%d: ERROR %s", n, k, conditionMessage(e)))
for (a in list(c(1, 3), c(2, 3), c(2, 2), c(3, 3))) cat(run(a[1], a[2]), "\n")
