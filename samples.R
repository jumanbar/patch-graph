nreps  <- 2000
npatch <- 200

source('patchCluster.R')

nf <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8)
isoSamp <- vector('list', length(nf))

names(isoSamp) <- paste('nf_', nf, sep='')

for (i in 1:length(nf)) {
  isoSamp[[i]] <- sampler(reps=nreps, npatch=npatch, ptsFun=makeIsoGrid,
                          intMax=100, dist_=1, noiseFactor=nf[i])
}
save(isoSamp, file='isoSamp.RData')

samples <- vector('list', 4)
names(samples) <- list('runifSamp', 'rnormSamp', 'rlnormSamp', 'rparetoSamp')
runifSamp   <- sampler(reps=nreps, ptsFun=runif,  intMax=50, npatch=npatch)
rnormSamp   <- sampler(reps=nreps, ptsFun=rnorm,  intMax=50, npatch=npatch)
rlnormSamp  <- sampler(reps=nreps, ptsFun=rlnorm, intMax=50, npatch=npatch)
rparetoSamp <- sampler(reps=nreps, ptsFun=rpareto, intMax=50, polar=TRUE,
                       location=0.01, shape=0.2, npatch=npatch)
save(isoSamp, file='isoSamp.RData')                   

samples$runifSamp   <- runifSamp
samples$rnormSamp   <- rnormSamp
samples$rlnormSamp  <- rlnormSamp
samples$rparetoSamp <- rparetoSamp
save(samples,  file='samples.RData')

### HACER TABLAS
makeTabla <- function(variable) {
# Antes de definir esta función en la sesión hay que tener los objetos
# isoSamp y samples en el espacio de trabajo.
  tabla <- data.frame(iso0.1=isoSamp$nf_0.1[variable],
                      iso0.25=isoSamp$nf_0.25[variable],
                      iso0.5=isoSamp$nf_0.5[variable],
                      iso1=isoSamp$nf_1[variable],
                      iso2=isoSamp$nf_2[variable],
                      iso4=isoSamp$nf_4[variable],
                      iso6=isoSamp$nf_6[variable],
                      iso8=isoSamp$nf_8[variable],
                      runif=samples$runifSamp[variable],
                      rnorm=samples$rnormSamp[variable],
                      rlnorm=samples$rlnormSamp[variable],
                      rpareto=samples$rparetoSamp[variable])
  names(tabla) <- c('iso0.1', 'iso0.25', 'iso0.5', 'iso1', 'iso2', 'iso4', 'iso6', 'iso8', 'runif', 'rnorm', 'rlnorm', 'rpareto')
  return(tabla)
}

vnames <- c("var", "skew", "kurt", "udist", "stepRsq", "hillRsq", "hillExp",
            "hillLoc", "hillSlope")

for (i in 1:length(vnames)) {
  name <- paste('t', vnames[i], sep='')
  assign(name, makeTabla(vnames[i]))
  save(list=name, file=paste(name, 'RData', sep='.'))
}

### ###
# sigmas <- function(tope, r, tf=1) {
#   cr <- c(.12, .12, rep(0, 6))
#   text(1:8 - 0.48 - cr, tope, labels=expression(sigma ^ 2 == ''), pos=4,
#        cex=tf, family='serif')
#   text(1:8 + 0.02 - cr, tope - r, labels=nf, pos=4, cex=tf, family='serif')
# }

# ths <- thillSlope
# thl <- thillLoc
# tsr <- tstepRsq
# names(ths) <- c(paste('Iso', 1:8), 'U', 'N', 'lnN')
# names(thl) <- names(ths)
# names(tsr) <- names(thl)
# png('boxplot-hill-exp-step-rsq.png', width=900, height=2100)
# op <- par(mfrow=c(3, 1), oma=c(5, 0, 0, 0), mar=c(2, 2.5, 4, 2), cex.axis=1.2, 
#           cex.main=1.5)
# tf <- 1
# boxplot(ths[-12], log='y', main='Pendiente máxima', xlab='',
#         pch=4, boxwex=0.4)
## sigmas(max(ths[-12]), 4)
# boxplot(thl[-12], log='y', main='K', xlab='', pch=4, boxwex=0.4)
# sigmas(0.07, 0, tf)
# boxplot(tsr[-12], xlab='', pch=4, boxwex=0.4, main='R2 respecto a la función escalón')
## sigmas(-0.65, 0)
## mtext(expression(paste(R^2, ' respecto a función escalón')), line=1, font=2, cex=1)
# mtext('Tratamiento', side=1, line=2, cex=1.5, outer=TRUE)
# par(op)
# dev.off()








