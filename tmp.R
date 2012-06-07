# path <- '~/Dropbox/jmb-maestria/tesis/resultados-cm/R/tablas'
path <- '~/projects/patch-graph/'
load(file.path(path, '/thillSlope.RData'))
load(file.path(path, '/thillLoc.RData'))
load(file.path(path, '/thillExp.RData'))
load(file.path(path, '/tstepRsq.RData'))
nf <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8)
ths <- thillSlope
thl <- thillLoc
tsr <- tstepRsq
the <- thillExp
names(ths) <- c(paste('I.', 1:8, sep=''), 'U', 'N', 'lnN')
names(thl) <- names(ths)
names(tsr) <- names(thl)
names(the) <- names(ths)
png('~/Dropbox/jmb-maestria/tesis/resultados-cm/figuras/boxplot-slope-exp.png',
    width=600, height=800)
op <- par(mfrow=c(2, 1), oma=c(4, 0, 0, 0), mar=c(3, 3.5, 3, 2), cex.axis=1.5, 
          cex.main=2)
boxplot(ths, log='y', main='Pendiente máxima', xlab='',
        boxwex=0.6)
mtext(' a', line=1.5, adj=0, font=2, cex=2)
# sigmas(max(ths[-12]), 4)
boxplot(the, log='y', main='h', xlab='y', boxwex=0.6)
# sigmas(0.07, 0, tf)
text(1:8, 1.2, labels=nf, cex=1.6, family='serif')
mtext(' b', line=1.5, adj=0, font=2, cex=2)
# boxplot(tsr[-12], xlab='', boxwex=0.6, main='Coeficiente de determinación')
# mtext(' c', line=1.5, adj=0, font=2, cex=3)
# sigmas(-0.65, 0)
# mtext(expression(paste(R^2, ' respecto a función escalón')), line=1, font=2, cex=1)
mtext('Tratamiento', side=1, line=2, cex=par('cex.main') * 1.2, outer=TRUE)
par(op)
dev.off()
