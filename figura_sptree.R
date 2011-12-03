## MST Y E(parches)

require(vegan)
require(igraph)
curr <- getwd()
ptos <- read.csv('ptos.csv')
rownames(ptos) <- paste('p', 1:12, sep='')
dis  <- dist(ptos)
sptr <- spantree(dis)
par(mfcol=c(1,2), mar=c(1, 1, 4, 1) + 0.1)
#~ par(mar=c(5, 4, 4, 2) + 0.1)
plot(sptr, cmdscale(dis), type='t', xlab='', ylab='', cex=1, cex.lab=1.4,
  bty='n', xaxt='n', yaxt='n')
text(.28, -.15, expression(d[max]))
text(.115, .028, expression(d[min]))

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

nd_ <- 250
d_min <- min(dis)
d_max <- max(sptr$dist)
d_i <- d_min - .09
d_f <- d_max + .09
d_  <- seq(d_i, d_f, , nd_)

ncl <- numeric(nd_)
cls <- vector('list', nd_)

for(i in 1:nd_) {
	g   <- makegraph(dis, ptos, d_[i])
	cl  <- clusters(g)
	ncl[i]   <- cl$no
	cls[[i]] <- cl$csize
}
par(mar=c(5, 4.5, 4, 2) + 0.1)
plot(ncl ~ d_, type='n', xlab='d', ylab='', main=expression(paste(N^{o},
  ' de componentes & ', E(parches))), bty='l')
points(ncl ~ d_, type='s', lwd=5, col='#848484')
# abline(v=c(d_min, d_max), col='#5C5C5C', lwd=2, lty=3)
segments(d_min, 1, d_min, 12, lwd=2, lty=3)
segments(d_max, 1, d_max, 12, lwd=2, lty=3)
# for(j in 1:nd_) points(cbind(d_[j], cls[[j]]), pch=20, cex=.8)
expectPatches <- sapply(cls, expval)
lines(expectPatches ~ d_, type='s', lwd=5)
text(d_max + .0005, 6, expression(d[max]), pos=4)
text(d_min - .0005, 6, expression(d[min]), pos=2)

par(mfcol=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

## IMAGEN DEL PAISAJE (PARCHES)

## dist mínima: entre los parches 6 y 10
## dist máxima: entre los parches 7 y 8
# 
# par(mfcol=c(1,1), mar=c(5, 4.5, 4, 2) + 0.1)
# plot(ptos, cex.lab=1.4, cex=4, pch=19, col='#7F7F7F', asp=1)             
# points(ptos, cex=4)
# with(ptos, segments(x[6], y[6], x[10], y[10]))
# with(ptos, segments(x[7], y[7], x[8], y[8]))
# 
# l1 <- ptos[6,] + (ptos[10,] - ptos[6,]) / 2
# l2 <- ptos[7,] + (ptos[8,] - ptos[7,]) / 2
# text(l1, labels=expression(d[min]), pos=3)
# text(l2, expression(d[max]), pos=2)

