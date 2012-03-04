sigmas <- function(tope, r, tf=1) {
  cr <- c(.12, .12, rep(0, 6))
  text(1:8 - 0.48 - cr, tope, labels=expression(sigma ^ 2 == ''), pos=4,
       cex=tf, family='serif')
  text(1:8 + 0.02 - cr, tope - r, labels=nf, pos=4, cex=tf, family='serif')
}
nf <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8)
ths <- thillSlope
the <- thillExp
tsr <- tstepRsq
names(ths) <- c(paste('I.', 1:8, sep=''), 'U', 'N', 'lnN')
names(the) <- names(ths)
names(tsr) <- names(the)
png('~/Dropbox/jmb-maestria/tesis/resultados-cm/figuras/boxplot-hill-exp-step-rsq.png',
    width=600, height=1500)
op <- par(mfrow=c(3, 1), oma=c(5, 0, 0, 0), mar=c(3, 3.5, 6, 2), cex.axis=2, 
          cex.main=3)
tf <- 3
boxplot(ths[-12], log='y', main='Pendiente máxima', xlab='',
        boxwex=0.6)
mtext(' a', line=1.5, adj=0, font=2, cex=3)
# sigmas(max(ths[-12]), 4)
boxplot(the[-12], log='y', main='Coeficiente de Hill (h)', xlab='', boxwex=0.6)
# sigmas(0.07, 0, tf)
text(1:8, 2, labels=nf, cex=tf, family='serif')
mtext(' b', line=1.5, adj=0, font=2, cex=3)
boxplot(tsr[-12], xlab='', boxwex=0.6, main='Coeficiente de determinación')
mtext(' c', line=1.5, adj=0, font=2, cex=3)
# sigmas(-0.65, 0)
# mtext(expression(paste(R^2, ' respecto a función escalón')), line=1, font=2, cex=1)
mtext('Tratamiento', side=1, line=3, cex=par('cex.main') * 1.2, outer=TRUE)
par(op)
dev.off()
