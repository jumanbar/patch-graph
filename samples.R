source('patchCluster.R')

require(VGAM)

rm(uniSamp)
rm(rnoSamp)
rm(rnlSamp)
rm(parSamp)
rm(isoSamp)
system('mv isoSamp.RData isoSamp_0.1.RData')

isoSamp_0.05 <- sampler(reps=1e4, ptsFun=makeIsoGrid, dist_=1, noiseFactor=0.05)
save(isoSamp_0.05, file='isoSamp_0.05.RData')
rm(isoSamp_0.05)

isoSamp_0.2 <- sampler(reps=1e4, ptsFun=makeIsoGrid, dist_=1, noiseFactor=0.2)
save(isoSamp_0.2, file='isoSamp_0.2.RData')
rm(isoSamp_0.2)

isoSamp_0.3 <- sampler(reps=1e4, ptsFun=makeIsoGrid, dist_=1, noiseFactor=0.3)
save(isoSamp_0.3, file='isoSamp_0.3.RData')
rm(isoSamp_0.3)

isoSamp_0.4 <- sampler(reps=1e4, ptsFun=makeIsoGrid, dist_=1, noiseFactor=0.4)
save(isoSamp_0.4, file='isoSamp_0.4.RData')
rm(isoSamp_0.4)

isoSamp_0.5 <- sampler(reps=1e4, ptsFun=makeIsoGrid, dist_=1, noiseFactor=0.5)
save(isoSamp_0.5, file='isoSamp_0.5.RData')
rm(isoSamp_0.5)

isoSamp_0.8 <- sampler(reps=1e4, ptsFun=makeIsoGrid, dist_=1, noiseFactor=0.8)
save(isoSamp_0.8, file='isoSamp_0.8.RData')
rm(isoSamp_0.8)

isoSamp_1 <- sampler(reps=1e4, ptsFun=makeIsoGrid, dist_=1, noiseFactor=1)
save(isoSamp_1, file='isoSamp_1.RData')
rm(isoSamp_1)

isoSamp_2 <- sampler(reps=1e4, ptsFun=makeIsoGrid, dist_=1, noiseFactor=2)
save(isoSamp_2, file='isoSamp_2.RData')
rm(isoSamp_2)

isoSamp_3 <- sampler(reps=1e4, ptsFun=makeIsoGrid, dist_=1, noiseFactor=3)
save(isoSamp_3, file='isoSamp_3.RData')
rm(isoSamp_3)

# uniSamp <- sampler(reps=1e4, ptsFun=runif)
# save(uniSamp, file='uniSamp.RData')
# 
# rnoSamp <- sampler(reps=1e4, ptsFun=rnorm)
# save(rnoSamp, file='rnoSamp.RData')
# 
# rnlSamp <- sampler(reps=1e4, ptsFun=rlnorm, polar=TRUE)
# save(rnlSamp, file='rnlSamp.RData')
# 
# parSamp <- sampler(reps=1e4, ptsFun=rpareto, polar=TRUE, location=0.01,
#                    shape=0.2)
# save(parSamp, file='parSamp.RData')

