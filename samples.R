source('patchCluster.R')

require(VGAM)

system('mv isoSamp.RData isoSamp_0.1.RData')

nf <- c(0.25, 0.5, 1, 2, 4)
isoSamp <- vector('list', length(nf))

names(isoSamp) <- paste('nf_', nf, sep='')

for (i in 1:length(nf)) {
  isoSamp[[i]] <- sampler(reps=1e4, ptsFun=makeIsoGrid, intMax=50,
	                      dist_=1, noiseFactor=0.05)
}
save(isoSamp, file='isoSamp.RData')
rm(isoSamp)

samples <- vector('list', 4)
names(samples) <- list('runifSamp', 'rnormSamp', 'rlnormSamp', 'rparetoSamp')
runifSamp   <- sampler(reps=1e4, ptsFun=runif,  intMax=50)
rnormSamp   <- sampler(reps=1e4, ptsFun=rnorm,  intMax=50)
rlnormSamp  <- sampler(reps=1e4, ptsFun=rlnorm, intMax=50)
rparetoSamp <- sampler(reps=1e4, ptsFun=rparetoSamp, intMax=50, polar=TRUE,
                       location=0.01, shape=0.2)

save(samples, file='samples.RData')
