nreps  <- 2000
npatch <- 200

source('patchCluster.R')

require(VGAM)

nf <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8)
isoSamp <- vector('list', length(nf))

names(isoSamp) <- paste('nf_', nf, sep='')

for (i in 1:length(nf)) {
  isoSamp[[i]] <- sampler(reps=nreps, npatch=npatch, ptsFun=makeIsoGrid, intMax=100,
	                      dist_=1, noiseFactor=nf[i])
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

png('boxplot-hill-exp-step-rsq.png', width=700*.65, height=1650*.65)
par(mfrow=c(3, 1))
boxplot(thillSlope[-12], log='y', main='Pendiente de Hill | x=K', xlab='Tratamiento')
boxplot(thillExp[-12], log='y', main='Coeficiente de Hill', xlab='Tratamiento')
boxplot(tstepRsq[-12], main='R^2 respecto a función escalón', xlab='Tratamiento')
dev.off()







