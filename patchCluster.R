## Clustering y E(parches)
# ejemplos:
# patchCluster(60, runif)
# patchCluster(60, rnorm)
# patchCluster(60, rnorm, sd=3)
read.csv('ptos.csv') -> pts
# patchCluster(puntos=ej.pts)

patchCluster <- function (npatch=10, puntos=runif,
                          doplot=TRUE, polar=FALSE, ...) {
  require(vegan)
  require(igraph)
  if (is.function(puntos)) {
    if (polar) {
      # Los puntos se crean con coordenadas polares y luego se convierte a un
      # sistema cartesiano.
      r <- puntos(npatch, ...)
      theta <- runif(npatch, 0, 2 * pi)
      x <- r * cos(theta)
      y <- r * sin(theta)
    } else {
      x <- puntos(npatch, ...)
      y <- puntos(npatch, ...)
    }
  } else {
     x <- puntos[,1]
     y <- puntos[,2]
     npatch <- length(x)
  }
  pts    <- data.frame(x=x, y=y)
  rownames(pts) <- 0:(npatch - 1)
  dis    <- dist(pts)
  dismat <- as.matrix(dist(pts))
  grafo  <- graph.adjacency(dismat, mode='undirected',
                            weighted=TRUE, diag=FALSE)
  mst    <- minimum.spanning.tree(grafo)
  sptr   <- spantree(dis)
  dists  <- get.edge.attribute(mst, 'weight')
  if (doplot) {
    op <- par(mfrow=c(2,2), mar=c(5, 4, 4, 1) + 0.1)
    plot(pts, cex.lab=1.4, cex=2, pch=19, col='#5D5D5D',
         asp=1, main='Paisaje')
    points(pts, cex=2)
  #   if(npatch < 200)
  #     mstlayout <- layout.fruchterman.reingold
  #   else
  #     mstlayout <- NULL
  #   par(mar=c(1, 1, 3, 1) + 0.1)
  #   plot(mst, main='MST', layout=mstlayout,
  #     vertex.color='#5D5D5D', vertex.frame.color='1', edge.width=3,
  #     edge.color=1, vertex.size=12, vertex.label.family='sans',
  #     vertex.label.color=0, vertex.label.cex=.7)
    plot(sptr, cmdscale(dis), type='b', xlab='', ylab='',
         cex=1, cex.lab=1.4, bty='n', xaxt='n', yaxt='n', main='MST')
    par(mar=c(5, 4, 4, 1) + 0.1)
  }
  # text(.28, -.15, expression(d[max]))
  # text(.115, .028, expression(d[min]))
  
  d_min <- min(dists)
  d_max <- max(dists)
  dife  <- d_max - d_min
  d_i <- d_min - .2 * dife
  d_f <- d_max + .2 * dife
#   d_  <- seq(d_i, d_f, , nd_)
  nd_ <- length(dists) + 2
  ncl <- numeric(nd_)
  cls <- vector('list', nd_)

  lineaE <- as.data.frame(matrix(0, nd_, 2))
  lineaE[1,] <- c(d_i, 1)
  lineaE[nd_,] <- c(d_f, npatch)
  d_ <- c(d_i, sort(dists), d_f)
  for(i in 1:nd_) {
    salen <- which(dists > d_[i])
    roto  <- delete.edges(mst, salen - 1)
    cl    <- clusters(roto)
    ncl[i]   <- cl$no
    cls[[i]] <- cl$csize
  }
  expectPatches <- sapply(cls, expval)
  if (doplot) {
    par(mar=c(5, 4, 4, 2) + 0.1)
    den <- density(dists)
    hist(dists, plot=FALSE) -> histo
    m <- max(max(den$y), max(histo$intensities))
    plot(histo, freq=F, main='Histograma de distancias del MST',
         xlab='Dist. en el MST', ylim=c(0, m))
    lines(den, lwd=2)
    par(mar=c(5, 5, 4, 2) + 0.1)
    plot(ncl ~ d_, type='n', xlab='Dist. de movimiento',
         main=expression(paste(N^{o}, ' de componentes & ',
         E(parches))), bty='l', ylab='')
    points(ncl ~ d_, type='s', lwd=5, col='#848484')
    segments(d_min, 1, d_min, npatch, lwd=2, lty=3)
    segments(d_max, 1, d_max, npatch, lwd=2, lty=3)
    lines(expectPatches ~ d_, type='s', lwd=5)
    text(d_max + .0005, npatch / 2, expression(d[max]), pos=4)
    text(d_min - .0005, npatch / 2, expression(d[min]), pos=2)
    
    par(mfcol=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
  }
  frac <- (d_max - d_min) / max(dis)
  slp  <- 1 / (d_max - d_min)
#   cat(' - (d_max - d_min) / max.dist =', round(frac, 4), '\n')
#   cat(' - Relative slope             =', round(slp, 4), '\n')
  invisible(list(genera=puntos, puntos=pts, dists=dists, frac=frac,
                 mst=mst, ncl=ncl, cls=cls, slope=slp, dmin=d_min, dmax=d_max,
                 epatch=expectPatches))
}


makegraph <- function(dis, pts, d_) {
   if(missing(dis)) {
		if(missing(pts))
			stop("Both dis & pts are empty!!")
		dis <- dist(pts)
	}
	m <- as.matrix(dis)
	m <- ifelse(m <= d_, 1, 0)
	0 -> diag(m)
	g <- graph.adjacency(m, mode='undirected', diag=FALSE)
	return(g)
}
expval <- function (x) {
# Extracts expected value from a vector of ocurrences
  table(x) -> tab
  sort(unique(x)) -> values
  values * as.vector(tab) -> suma
  suma / sum(suma) -> pr
  sum(values * pr) -> expvalue
  return(expvalue)
}
