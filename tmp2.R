plot.patchClusterAnalysis <- function(x) {
  require(vegan)
  
  for (i in 1:length(names(x)))
    assign(names(x)[i], x[[i]])
  for (i in 1:length(names(ca)))
    assign(names(ca)[i], ca[[i]])
  dis    <- dist(pts)
  sptr   <- spantree(dis)
  dists  <- get.edge.attribute(mst, 'weight')
  npatch <- nrow(pts)
  
  op <- par(mfrow=c(2,2), mar=c(5, 4, 4, 1) + 0.1)
  plot(pts, cex.lab=1.4, cex=0.9, pch=19, col='#5D5D5D',
       asp=1, main='Paisaje')
  points(pts, cex=0.9)
  plot(sptr, cmdscale(dis), type='b', xlab='', ylab='',
       cex=1, cex.lab=1.4, bty='n', xaxt='n', yaxt='n', 
       main='Árbol Recubridor Mínimo')
#   par(mar=c(5, 4, 4, 1) + 0.1)
#   den <- density(dists)
#   hist(dists, plot=FALSE) -> histo
#   m <- max(max(den$y), max(histo$intensities))
#   plot(histo, freq=F, main='Histograma de distancias del MST',
#        xlab='Dist. en el MST', ylim=c(0, m))
#   lines(den, lwd=2)
  par(mar=c(5, 5, 4, 2) + 0.1)
  plot(ncl ~ movDist, type='n', xlab='D',
       main='No. de componentes & o(M)', bty='l', ylab='')
  points(ncl ~ movDist, type='s', lwd=5, col='#848484')
  segments(dmin, 1, dmin, npatch, lwd=2, lty=3)
  segments(dmax, 1, dmax, npatch, lwd=2, lty=3)
  lines(epatch ~ movDist, type='s', lwd=5)
  text(dmax + .0005, npatch / 2, expression(d[max]), pos=4)
  text(dmin - .0005, npatch / 2, expression(d[min]), pos=2)
  curve(stepFun(x, a=fstep$minimum, min=min(epatch), max=max(epatch)),
        n=2001, add=FALSE, col='#C2C2C2', lwd=2, main='Hill & Escalón ajustadas',
        ylab='', xlab='D')
  if (class(fhill) %in% c("nls", "list")) {
    ep <- epatch / max(epatch)
    curve(hill(x, a=cf(fhill)[1], k=cf(fhill)[2], ymin=min(ep)) * 
      max(epatch), n=2001, add=TRUE, col=1, lwd=2, lty=2)
  }
  
  par(mfcol=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}
patchCluster(60, rnorm)