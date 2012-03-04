## Clustering y E(parches)
# ejemplos:
# patchCluster(60, runif)
# patchCluster(60, rnorm)
# patchCluster(60, rnorm, sd=3)
# read.csv('ptos.csv') -> pts
# patchCluster(puntos=ej.pts)

mstClusterAnalysis <- function(mst, verbose=TRUE, ...) {
# Toma el mst y hace el análisis del número de clusters y el valor esperado 
# de número de parches a los que se accede en base a un gradiente de 
# distancias de movimiento.
  require(igraph)
  dists <- get.edge.attribute(mst, 'weight')
  dists <- dists / max(dists)
  d_min <- min(dists)
  d_max <- max(dists)
  dife  <- d_max - d_min
  d_i   <- d_min - .2 * dife
  d_f   <- d_max + .2 * dife
  nd_ <- length(dists) + 2
  ncl <- numeric(nd_)
  cls <- vector('list', nd_)
  movDist      <- c(d_i, sort(dists), d_f)
  if (movDist[1] < 0)
    movDist[1] <- 1e-6
  for(i in 1:nd_) {
    salen <- which(dists > movDist[i])
    roto  <- delete.edges(mst, salen - 1)
    cl    <- clusters(roto)
    ncl[i]   <- cl$no
    cls[[i]] <- cl$csize
  }
  epatch <- sapply(cls, expval)
  # epatch: Expected number of patches

  # Fitting models:
  fstep  <- fit2step(epatch, movDist) # Step function
  fhill  <- fit2hill(epatch, movDist, ...) # Hill function
  
  # Correction:
  if (class(fhill) %in% c("nls", "list")) {
    hill.h <- cf(fhill)[1]
    hill.k <- cf(fhill)[2]
  } else {
    hillSS <- hill.h <- hill.k <- NA
  }
  
  if (verbose) {
    cat(' - R Squared for Step Function   =', round(fstep$rsq, 4), '\n')
    cat(' - Step location parameter       =', round(fstep$minimum, 4), '\n')
    cat(' - R Squared for Hill Function   =', round(fhill$rsq, 4), '\n')
    cat(' - Hill Function coefficient     =', round(hill.h, 4), '\n')
    cat(' - Hill Function location param  =', round(hill.k, 4), '\n')
  }
  
  out <- list(cls=cls, dmin=d_min, dmax=d_max, epatch=epatch, fhill=fhill,
              fstep=fstep, movDist=movDist, ncl=ncl)
  invisible(out)  
}

patchCluster <- function (npatch=10, puntos=runif, doplot=TRUE, polar=FALSE,
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

makeIsoGrid <- function(n_, dist_=1, noiseFactor=0.1) {
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

cf <- function(x) {
  if (class(x) %in% c('lm', 'nls', 'glm'))
    return(coef(x))
  if (is.list(x) && 'par' %in% names(x))
	return(x$par)

}

expval <- function (x)
# Weighted mean with weights vector being the same as the input vector
  weighted.mean(x, x)

fit2hill <- function(epatch, movDist, ...) {
  ep <- epatch / max(epatch)
#   ep <- epatch
  e1 <- min(ep)
  e2 <- max(ep)
  md <- movDist
  if (length(unique(md)) > 1) {
    out <- optim(c(1, .5), hillSumSq, ymin=e1, ymax=e2, x=md, y=ep)
  } else { 
    out <- list(par=c(1e6, unique(md)), value=0, rsq=1)
  }
  out$rsq <- rsquared(ep, out$value)
  return(out)
}

fit2step <- function(epatch, movDist) {
  ep <- epatch / max(epatch)
  e1 <- min(ep)
  e2 <- max(ep)
  md <- movDist
  if (length(unique(md)) > 1) {
    out <- optimize(stepSumSq, interval=range(md), x=md, y=ep, min=e1, max=e2)
  } else {
    out <- list(minimum=unique(md), objective=0, rsq=1)
  }
  browser()
  out$rsq <- rsquared(ep, out$objective)
  return(out)
}

hill <- function(x, h, k, ymin=0, ymax=1) {
  y <- x ^ h / (x ^ h + k ^ h)
  y <- ymax * (y + ymin * (1 - y))
  # ymax * y + ymax * ymin - ymax * ymin * y
  return(y)
}

hillderiv <- function(x, h, k, ymin, ymax=1) {
  denom <- (x ^ h + k ^ h)
  y <- (h * x ^ (h - 1)) / denom - (h * x ^ (2 * h - 1)) / denom ^ 2
  ymax * y - ymax * ymin * y
}

hillSumSq <- function(pm=c(h=1, k=.5), x, y, ymin, ymax) {
# Hill Function Sum of Squares
  out <- sum((hill(x, h=pm[1], k=pm[2], ymin=ymin, ymax=ymax) - y) ** 2)
  return(out)
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
    curve(stepFun(x, a=fstep$minimum, min=min(epatch), max=max(epatch)),
          n=2001, add=TRUE, col='#C2C2C2', lwd=2)
    if (class(fhill) %in% c("nls", "list")) {
      ep <- epatch / max(epatch)
      curve(hill(x, h=cf(fhill)[1], k=cf(fhill)[2], ymin=min(ep)) * 
            max(epatch), n=2001, add=TRUE, col=1, lwd=2, lty=2)
    }
    lines(epatch ~ movDist, type='s', lwd=5)
    text(dmax + .0005, npatch / 2, expression(d[max]), pos=4)
    text(dmin - .0005, npatch / 2, expression(d[min]), pos=2)
    
    par(mfcol=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}

sampler <- function(reps=100, npatch=20, ptsFun=runif, intMax=30, 
                    verbose=TRUE, ...) {
  require(moments)
  require(igraph)
  nombres <- c('var', 'skew', 'kurt', 'udist', 'stepRsq', 'hillRsq',
               'hillExp', 'hillLoc', 'hillSlope')
  tabla   <- matrix(nrow=reps, ncol=length(nombres))
  tabla   <- as.data.frame(tabla)
  names(tabla) <- nombres
  i <- 1
  int <- 1
  while (i < reps + 1) {
    cat('Rep:', i, '  intento:', int, '\n')
    pc <- patchCluster(npatch, puntos=ptsFun, doplot=FALSE, verbose=FALSE, ...)
    stepRsq <- pc$ca$fstep$rsq
    if (class(pc$ca$fhill) %in% c('nls', 'list')) {
      dists <- get.edge.attribute(pc$mst, 'weight')
      distFromUniform <- dists - punif(dists, min=pc$ca$dmin, max=pc$ca$dmax)
      distFromUniform <- sqrt(mean(distFromUniform ^ 2))
      hillRsq <- pc$ca$fhill$rsq
      hillExp <- cf(pc$ca$fhill)[1]
      hillLoc <- cf(pc$ca$fhill)[2]
      hillSlope <- hillExp / (4 * hillLoc)
      
      tabla[i,] <- c(var(dists),
                     skewness(dists),
                     kurtosis(dists),
                     distFromUniform,
                     stepRsq,
                     hillRsq,
                     hillExp,
                     hillLoc,
		     hillSlope)
      i <- i + 1
      int <- 1
    } else {
      int <- int + 1
      if (int > intMax) {
        i <- i + 1
        int <- 1
      }
    }
  }
  if (verbose) {
    print(summary(tabla))
    print(cor(tabla))
    pairs(tabla)
  }
  invisible(tabla)
}

stepFun <- function(x, a=0, min=0, max=1) {
# Step Function
  y <- x
  y[x <= a] <- min
  y[x >  a] <- max
  return(y)
}

stepSumSq <- function(a, x, y, ...) {
# Step Function Sum of Squares
  out <- sum((stepFun(x, a, ...) - y) ** 2)
  return(out)
}
# 
rsquared <- function(ep, SSerr) {
# Works only inside fit2hill and fit2step
  SStot <- sum((ep - mean(ep)) ** 2)
  return(1 - SSerr / SStot)
}

#      ans$r.squared <- mss/(mss + rss)
#         ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
#             df.int)/rdf)
