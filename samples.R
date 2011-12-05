source('patchCluster.R')
reps <- 2000
nombres <- c('var', 'skew', 'kurt', 'udist', 'slope', 'frac')
tb2   <- matrix(0, reps, length(nombres))
tb2   <- as.data.frame(tb2)
names(tb2) <- nombres
for (i in 1:reps) {
  a <- patchCluster(20, puntos=rnorm, doplot=F)
  b <- a$dists - punif(a$dists, min=a$dmin, max=a$dmax)
  b <- sqrt(mean(b ^ 2))
  tb2[i,] <- c(var(a$dists),
              skewness(a$dists),
              kurtosis(a$dists),
              b,
              a$slope,
              a$frac)
}
