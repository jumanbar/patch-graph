require(manipulate)

hill <- function(x, a, k) {
  y <- x ^ a / (x ^ a + k ^ a)
}

manipulate(
  curve(hill(x, a=a,  k=1), n=1e4+1, from=0, to=2, ylab='y', lwd=2, ylim=0:1)
  , a=slider(0, 160, 30))
