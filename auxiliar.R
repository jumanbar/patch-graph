mklands <- function(dim_=2, dist_=1, lmax_=2, n_=3, rdist_=3, type='fractal') {
# ejemplo:
#   plot(mklands())
  require(splancs, quietly=TRUE)
  require(ellipse, quietly=TRUE)
  n_ <- round(n_)
  npatch <- (n_ ^ dim_) ^ lmax_

  areas  <- vector('list', lmax_)
  belong <- as.data.frame(matrix(1:npatch, nrow=npatch, ncol=lmax_))
  names(belong) <- names(areas) <- paste('l', 0:(lmax_ - 1), sep='')

  if (lmax_ <= 1) {
    areas <- vector('list', 1)
    if (npatch == 1) {
      coordsAll <- matrix(0, 1, 2)
      pos       <- 0
    }
  }

  if (type == 'fractal' && npatch > 1) {
    dists_ <- NULL
    seg <- dists_
    for (l_ in 1:lmax_) {
      for (i in 1:(n_ - 1)) {
        dists_ <- c(dists_, c(dist_ * rdist_ ^ (l_ - 1), seg))
      }
      seg <- dists_
    }
    pos <- c(0, cumsum(dists_))
  }

  if (type == 'regular' && npatch > 1) {
    n_ <- (n_ ^ dim_) ^ lmax_
    pos <- 0:n_ * dist_
  }

  if (type == 'randUnif' && npatch > 1) {
    x  <- sort(runif (npatch, 0, (n_ - 1) * dist_))
    y  <- runif (npatch, 0, (n_ - 1) * dist_)
    coordsAll <- cbind(x=x, y=y)
    pos <- NULL
  }

  if (type != 'randUnif' && npatch > 1) {
    npos      <- length(pos)
    coordsAll <- matrix(nrow=npatch, ncol=2)

    p <- 0:npos * npos + 1
    for (i in 1:npos) {
      coordsAll[p[i]:(p[i + 1] - 1),] <- cbind(pos[i], pos)
    }

    if (type == 'fractal' && lmax_ > 1) {

      for (l_ in 1:(lmax_ - 1)) {
        p   <- 1
        nl  <- n_ ^ l_
        b   <- (dist_ * rdist_ ^ l_) * .15
        ini <- 0:(npos / nl - 1) * nl + 1
        level <- vector('list', (n_ ^ dim_) ^ (lmax_ - l_))

        for (j in ini) {
          x <- c(pos[j], pos[j + nl - 1])
          for (i in ini) {
            y <- c(pos[i], pos[i + nl - 1])
            linea <- rbind(c(x[1] - b, y[1] - b),
                      c(x[1] - b, y[2] + b),
                      c(x[2] + b, y[2] + b),
                      c(x[2] + b, y[1] - b),
                      c(x[1] - b, y[1] - b))
            level[[p]] <- linea
            id <- inpip(coordsAll, linea, bound=TRUE)
            belong[id, l_ + 1] <- p
            p <- p + 1
          }
        }
        areas[[l_ + 1]] <- level
      }
    }
  }

  patchAreas <- vector('list', nrow(coordsAll))
  for (i in 1:nrow(coordsAll)) {
    patchAreas[[i]] <- ellipse(0, centre=c(coordsAll[i,1], coordsAll[i,2]), t=dist_ * .5, npoints=15)
  }
  areas[[1]] <- patchAreas
  
  coordsAll <- as.data.frame(coordsAll)
  names(coordsAll) <- c('x', 'y')

  frm <- formals()
  parms <- lapply(names(frm), get, envir=sys.parent(0))
  names(parms) <- names(frm)
  out <- list(areas=areas, belong=belong, coordsAll=coordsAll,
              parms=parms, pos=pos)
  class(out) <- 'lands'
  return(out)

}
# Ejemplo x <- mklands() (paisaje fractal):
# x
#   areas
#     l0: areas al nivel 0, circulares
#     l1: areas al nivel 1, cuadradas
#     ...
#     ln: areas al nivel n, cuadradas
#   belong
#     l0: pertenencia de parches respecto al nivel 0
#     l1: pertenencia de parches respecto al nivel 1
#     ...
#     ln: pertenencia de parches respecto al nivel n
#   coordsAll: coordenadas x y de todos los parches nivel 0
#   parms
#     dim_:  dimensión del paisaje (2)
#     dist_: distancia entre parches nivel 0
#     lmax_: número de niveles del paisaje
#     n_: número de parches 0 en un parche nivel 1
#     rdist_: razón de distancias entre niveles consecutivos
#     type: tipo de paisaje ('fractal' por defecto).
#   pos: posiciones de los parches en la proyección sobre cualquiera de los ejes
#        de coordenadas.


n_ <- 10001
out <- numeric(n_)
a <- seq(1.00091, 1.00093, len=n_)
for (i in 1:n_)
  out[i] <- stepSS(a[i], md, ep, min=1/196, max=1) / 196
plot(out ~ a, type='l')







