## Clustering y E(parches)
# ejemplos:
# patchCluster(60, runif)
# patchCluster(60, rnorm)
# patchCluster(60, rnorm, sd=3)
# read.csv('ptos.csv') -> pts
# patchCluster(puntos=ej.pts)

mstClusterAnalysis <- function(mst, verbose=TRUE) {
# Toma el mst y hace el análisis del número de clusters y el valor esperado 
# de número de parches a los que se accede en base a un gradiente de 
# distancias de movimiento.
  require(igraph)
  dists <- get.edge.attribute(mst, 'weight')
  d_min <- min(dists)
  d_max <- max(dists)
  dife  <- d_max - d_min
  d_i   <- d_min - .2 * dife
  d_f   <- d_max + .2 * dife
#   movDist  <- seq(d_i, d_f, , nd_)
  nd_ <- length(dists) + 2
  ncl <- numeric(nd_)
  cls <- vector('list', nd_)

#   lineaE       <- as.data.frame(matrix(0, nd_, 2))
#   lineaE[1,]   <- c(d_i, 1)
#   lineaE[nd_,] <- c(d_f, npatch)
  movDist      <- c(d_i, sort(dists), d_f)
  for(i in 1:nd_) {
    salen <- which(dists > movDist[i])
    roto  <- delete.edges(mst, salen - 1)
    cl    <- clusters(roto)
    ncl[i]   <- cl$no
    cls[[i]] <- cl$csize
  }
  epatch <- sapply(cls, expval)
  # epatch: Expected number of patches
  ajuste <- fit2step(epatch, movDist)
  ajuste$objective <- ajuste$objective / length(epatch)
#   epSlps <- diff(epatch) / diff(movDist)
#   frac   <- (d_max - d_min) / max(dis)
#   
#   slp    <- 1 / (d_max - d_min)
  if (verbose) {
    cat(' - Square Sum from a step function =', round(ajuste$objective, 4), '\n')
#     cat(' - Fraction of all max. dist. (d_max - d_min) / max.dist =', round(frac, 4), '\n')
#     cat(' - Relative slope: 1 / (d_max - d_min)                   =', round(slp, 4), '\n')
#     cat(' - Slopes summary:\n')
#     print(summary(epSlps))
  }
  
  out <- list(cls=cls, dmin=d_min, dmax=d_max, epatch=epatch, ncl=ncl,
#               epSlps=epSlps, frac=frac,
              movDist=movDist, fit=ajuste)
#               slope=slp)
  invisible(out)  
}

patchCluster <- function (npatch=10, puntos=runif,
                          doplot=TRUE, polar=FALSE,
                          verbose=TRUE, ...) {
  require(vegan)
  require(igraph)

  pts   <- makePointsLocation(npatch, puntos, polar, ...)
  grafo <- makeGraph(pts)
  mst   <- minimum.spanning.tree(grafo)
  ca    <- mstClusterAnalysis(mst, verbose)
  # ca: Cluster Analysis

  out <- list(ca=ca, genera=puntos, mst=mst,
              pts=pts)
  class(out) <- 'patchClusterAnalysis'
  if (doplot) {
    plot(out)
  }
  invisible(out)
}

makeGraph <- function(pts) {
  dis    <- dist(pts)
  dismat <- as.matrix(dist(pts))
  grafo  <- graph.adjacency(dismat, mode='undirected',
                            weighted=TRUE, diag=FALSE)
  return(grafo)
}

makeGraphOnDists <- function(dis, pts, movDist) {
# Makes a unweighted & undirected graph based on the locations of nodes and the 
# maximum distance between each posible pair.
	if(missing(dis)) {
		if(missing(pts))
			stop("Both dis & pts are empty!!")
		dis <- dist(pts)
	}
	m <- as.matrix(dis)
	m <- ifelse(m <= movDist, 1, 0)
	0 -> diag(m)
	g <- graph.adjacency(m, mode='undirected', diag=FALSE)
	return(g)
}

makeIsoGrid <- function(n_, dist_, noiseFactor=0.1) {
  xcor <- 0:(n_ - 1) * dist_ + dist_ / 4
  altura <- sin(pi / 3) * dist_
  xycoords <- matrix(0, n_ ** 2, 2)
  for (i in 1:n_) {
    ycor <- (i - 1) * altura 
    xycoords[((i - 1) * n_ + 1):(i * n_), ] <-
      cbind(xcor + dist_ * .25 * (-1) ^ i, ycor)
  }
  noise <- dist_ * noiseFactor
  ruido <- matrix(rnorm(length(xycoords), 0, noise), ncol=2)
  xycoords <- xycoords + ruido
  return(xycoords)
}

makePointsLocation <- function(npatch, puntos, polar, ...) {
  if (is.function(puntos)) {
    if (!identical(makeIsoGrid, puntos)) {
      if (polar) {
        # Los puntos se crean con coordenadas polares y luego se convierte a un
        # sistema cartesiano, siguiendo estos pasos:
        #  1. Usar el generador de números aleatorios para obtener la coordenada
        #     radial ("r") + runif para obtener la coordenada angular ("theta")
        #     (entre 2 y 2*pi).
        r <- puntos(npatch, ...)
        theta <- runif(npatch, 0, 2 * pi)
        #  2. Convertir r y theta a valores de x e y utilizando las fórmulas:
        x <- r * cos(theta)
        y <- r * sin(theta)
      } else {
        x <- puntos(npatch, ...)
        y <- puntos(npatch, ...)
      }
    } else {
      np <- round(sqrt(npatch))
      npatch <- np ** 2 
      pts <- makeIsoGrid(np, ...)
    }
  } else {
     x <- puntos[,1]
     y <- puntos[,2]
     npatch <- length(x)
  }
  if (!exists('pts'))
    pts <- data.frame(x=x, y=y)
  rownames(pts) <- 0:(npatch - 1)
  
  return(pts)
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

fit2step <- function(epatch, movDist) {
  ep <- epatch / max(epatch)
  e1 <- min(ep)
  e2 <- max(ep)
  md <- movDist
  if (length(unique(md)) > 1) {
    op <- optimize(stepSS, interval=range(md), x=md, y=ep, min=e1, max=e2)
  } else {
    op <- list(minimum=unique(md), objective=0)
  }
  return(op)
}

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
    den <- density(dists)
    hist(dists, plot=FALSE) -> histo
    m <- max(max(den$y), max(histo$intensities))
    plot(histo, freq=F, main='Histograma de distancias del MST',
         xlab='Dist. en el MST', ylim=c(0, m))
    lines(den, lwd=2)
    par(mar=c(5, 5, 4, 2) + 0.1)
    plot(ncl ~ movDist, type='n', xlab='Dist. de movimiento',
         main=expression(paste(N^{o}, ' de componentes & ',
                         E(parches))), bty='l', ylab='')
    points(ncl ~ movDist, type='s', lwd=5, col='#848484')
    segments(dmin, 1, dmin, npatch, lwd=2, lty=3)
    segments(dmax, 1, dmax, npatch, lwd=2, lty=3)
    curve(stepFun(x, a=fit$minimum, min=min(epatch), max=max(epatch)),
          n=2001, add=TRUE, col='#C2C2C2', lwd=2)          
    lines(epatch ~ movDist, type='s', lwd=5)
    text(dmax + .0005, npatch / 2, expression(d[max]), pos=4)
    text(dmin - .0005, npatch / 2, expression(d[min]), pos=2)
    
    par(mfcol=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}

sampler <- function(reps=100, npatch=20, ptsFun=runif, ...) {
  require(moments)
  require(igraph)
  nombres <- c('var', 'skew', 'kurt', 'udist', 'stepFit')
  tabla   <- matrix(0, reps, length(nombres))
  tabla   <- as.data.frame(tabla)
  names(tabla) <- nombres
  for (i in 1:reps) {
    pc <- patchCluster(npatch, puntos=ptsFun, doplot=FALSE, verbose=FALSE, ...)
    dists <- get.edge.attribute(pc$mst, 'weight')
    distFromUniform <- dists - punif(dists, min=pc$ca$dmin, max=pc$ca$dmax)
    distFromUniform <- sqrt(mean(distFromUniform ^ 2))
    tabla[i,] <- c(var(dists),
                   skewness(dists),
                   kurtosis(dists),
                   distFromUniform,
                   pc$ca$fit$objective)
  }
  print(summary(tabla))
  print(cor(tabla))
  pairs(tabla)
  invisible(tabla)
}

stepFun <- function(x, a=0, min=0, max=1) {
# Step Function
  y <- x
  y[x <= a] <- min
  y[x >  a] <- max
  return(y)
}

stepSS <- function(a, x=movDist, y=epatch, ...) {
# Step Function Sum of Squares
  out <- sum(sqrt((stepFun(x, a, ...) - y) ** 2))
  return(out)
}
