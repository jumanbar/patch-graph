source('patchCluster.R')

require(VGAM)

nf <- c(0.25, 0.5, 1, 2, 4, 6, 8)
isoSamp <- vector('list', length(nf))

names(isoSamp) <- paste('nf_', nf, sep='')

for (i in 1:length(nf)) {
  isoSamp[[i]] <- sampler(reps=1e4, npatch=200, ptsFun=makeIsoGrid, intMax=50,
	                      dist_=1, noiseFactor=0.05)
}
save(isoSamp, file='isoSamp.RData')
rm(isoSamp)

samples <- vector('list', 4)
names(samples) <- list('runifSamp', 'rnormSamp', 'rlnormSamp', 'rparetoSamp')
runifSamp   <- sampler(reps=1e4, ptsFun=runif,  intMax=50)
rnormSamp   <- sampler(reps=1e4, ptsFun=rnorm,  intMax=50)
rlnormSamp  <- sampler(reps=1e4, ptsFun=rlnorm, intMax=50)
rparetoSamp <- sampler(reps=1e4, ptsFun=rpareto, intMax=50, polar=TRUE,
                       location=0.01, shape=0.2)
# 
# 
# ### Nuevo
load('isoSamp.RData')
load('samples.RData')
# 
isoSamp$nf_8 <- sampler(reps=1e4, ptsFun=makeIsoGrid, intMax=50,
	                    dist_=1, noiseFactor=8)
save(isoSamp, file='isoSamp.RData')                   
rparetoSamp  <- sampler(reps=1e4, ptsFun=rpareto, intMax=50, polar=TRUE,
                       location=0.01, shape=0.2)
samples$runifSamp   <- runifSamp
samples$rnormSamp   <- rnormSamp
samples$rlnormSamp  <- rlnormSamp
samples$rparetoSamp <- rparetoSamp

thillExp <- data.frame(iso0.25=isoSamp$nf_0.25$hillExp,
                         iso0.5=isoSamp$nf_0.5$hillExp,
				         iso1=isoSamp$nf_1$hillExp,
				         iso2=isoSamp$nf_2$hillExp,
                         iso4=isoSamp$nf_4$hillExp,
                         iso8=isoSamp$nf_8$hillExp,
                         runif=samples$runifSamp$hillExp,
                         rnorm=samples$rnormSamp$hillExp,
                         rlnorm=samples$rlnormSamp$hillExp,
                         rpareto=samples$rparetoSamp$hillExp)

thillSS <- data.frame(iso0.25=isoSamp$nf_0.25$hillSS,
                 iso0.5=isoSamp$nf_0.5$hillSS,
                 iso1=isoSamp$nf_1$hillSS,
                 iso2=isoSamp$nf_2$hillSS,
                 iso4=isoSamp$nf_4$hillSS,
                 iso8=isoSamp$nf_8$hillSS,
                 runif=samples$runifSamp$hillSS,
                 rnorm=samples$rnormSamp$hillSS,
                 rlnorm=samples$rlnormSamp$hillSS,
                 rpareto=samples$rparetoSamp$hillSS)

save(thillExp, file='tillParam.RData')
save(thillSS, file='tillSS.RData')
save(samples, file='samples.RData')
rm(samples, isoSamp)

boxplot(thillExp, ylim=c(0, 140))
boxplot(thillExp, log='y')
boxplot(thillSS, ylim=c(0, 0.0115))
boxplot(thillSS, log='y')

pruebita <- sampler(reps=1e4, ptsFun=makeIsoGrid, intMax=50,
	                dist_=1, noiseFactor=8)







