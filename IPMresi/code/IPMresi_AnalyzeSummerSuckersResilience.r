###########################################################################
### Exploring number of bins and resilience metrics
### Christina Hernandez
### Sept 2025
###########################################################################

require(popdemo)
require(here)
require(lattice)

require(doParallel)
require(doSNOW)
cl <- makeCluster(10)
registerDoSNOW(cl)

source(here('IPMresi', 'code', 'MatrixImage.R'))
source(here('IPMresi', 'code', 'IPMresi_SummerSuckersIPMfunctions.R'))

## Calculate resilience metrics across a range of number of mesh points: -------
min.size <- 0.00
max.size <- 600.00
# model parameters for building the kernel are stored in ss_m_par

# test set:
n_bins<- seq(2,82, by=5)
# initialize output:
output<- data.frame(Nbins=n_bins,
                    AmplificationFirst=vector(length=length(n_bins)),
                    AmplificationMax=vector(length=length(n_bins)),
                    AttenuationFirst=vector(length=length(n_bins)),
                    AttenuationMax=vector(length=length(n_bins)),
                    RecoveryTime=vector(length=length(n_bins)))

for (i in 1:length(n_bins)){
  cat(i, '\n')
  # generate the kernel at this number of bins:
  thisK<- make_K_suckers(n_bins[i], ss_m_par, min.size, max.size)

  # Calculate the resilience metrics from popdemo:
  output$AmplificationFirst[i]<- popdemo::reac(thisK, bound='upper')
  output$AmplificationMax[i]<- popdemo::maxamp(thisK)
  output$AttenuationFirst[i]<- popdemo::reac(thisK, bound='lower')
  output$AttenuationMax[i]<- popdemo::maxatt(thisK)
  output$RecoveryTime[i]<- max(popdemo::convt(thisK))
}

png(filename=here('IPMresi', 'results', 'SummerSuckers_NbinsResilience_linegraph.png'),
    height=10, width=7, units='in', res = 300)
par(mfrow=c(3,2))

plot(AmplificationFirst~Nbins, data=output, type='l', log='y',
     xlab='Number of bins', ylab='Amplification (1st timestep)')
plot(AmplificationMax~Nbins, data=output, type='l', log='y',
     xlab='Number of bins', ylab='Amplification (max)')
plot(AttenuationFirst~Nbins, data=output, type='l', log='y',
     xlab='Number of bins', ylab='Attenuation (1st timestep)')
plot(AttenuationMax~Nbins, data=output, type='l', log='y',
     xlab='Number of bins', ylab='Attenuation (max)')
plot(RecoveryTime~Nbins, data=output, type='l',
     xlab='Number of bins', ylab='Recovery Time')
dev.off()


## What would it mean to do a Gaussian?
initialvec<- rep(NA, 50)
par(mfrow=c(2,2))
plot(1:50, dnorm(1:50, mean=4, sd=1), type='l')
plot(1:50, dnorm(1:50, mean=4, sd=2), type='l')
plot(1:50, dnorm(1:50, mean=4, sd=0.5), type='l')
plot(1:50, dnorm(1:50, mean=10, sd=2), type='l')

reac_with_gaussian_vectors<- function(Amat, nbins, sd){

  vectors<- sapply(1:nbins, function(x){dnorm(1:nbins, mean=x, sd=sd)})
  # Each row is an initial vector centered at one of the bins, with the given sd

  reacz<- sapply(1:nbins, function(x){popdemo::reac(Amat, vector = vectors[x,])})

  return(data.frame(upper=max(reacz), lower=min(reacz)))
}

# Heat map plots:

# test set:
n_bins<- seq(2,300, by=5)
sd_vals<- seq(0.25, 2, by=0.25)
testgrid<- expand.grid(n_bins, sd_vals)
names(testgrid)<- c("n_bins", "sd_vals")

pb <- txtProgressBar(max = length(testgrid$n_bins), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

results<- foreach(i=1:length(testgrid$n_bins), .combine = rbind,
                  .options.snow = opts) %dopar% {

  # generate the kernel at this number of bins:
  thisK<- make_K_suckers(testgrid$n_bins[i], ss_m_par, min.size, max.size)

  thisoutput<- reac_with_gaussian_vectors(thisK, testgrid$n_bins[i], testgrid$sd_vals[i])
}

temp<- as.data.frame(results)
testgrid$upper<- unlist(temp$upper)
testgrid$lower<- unlist(temp$lower)

png(here('IPMresi', 'results', 'SummerSuckersExample_AmplificationFirst_levelplot.png'),
    height=5, width=5, units='in', res=300)

lattice::levelplot(upper~n_bins+sd_vals, data=testgrid, xlab='N bins', ylab='SD of Gaussian',
                   main='Amplification (1st time step)')

dev.off()

# For attenuation, values greater than 1 mean that the algorithm doesn't work at all - let's make those NA
testgrid$lower[testgrid$lower>1]<- NA

png(here('IPMresi', 'results', 'SummerSuckersExample_AttenuationFirst_levelplot.png'),
    height=5, width=5, units='in', res=300)

lattice::levelplot(log10(lower)~n_bins+sd_vals, data=testgrid,
                   xlab='N bins', ylab='SD of Gaussian', main='log10[Attenuation (1st time step)]')
dev.off()

## Extend the plots for the first-time-step metrics: ---------------------------

min.size <- 0.00
max.size <- 600.00
# model parameters for building the kernel are stored in ss_m_par

# test set:
n_bins<- seq(2,500, by=5)
# initialize output:
output<- data.frame(Nbins=n_bins,
                    AmplificationFirst=vector(length=length(n_bins)),
                    AttenuationFirst=vector(length=length(n_bins)))

for (i in 1:length(n_bins)){
  cat(i, '\n')
  # generate the kernel at this number of bins:
  thisK<- make_K_suckers(n_bins[i], ss_m_par, min.size, max.size)

  # Calculate the resilience metrics from popdemo:
  output$AmplificationFirst[i]<- popdemo::reac(thisK, bound='upper')
  output$AttenuationFirst[i]<- popdemo::reac(thisK, bound='lower')
}

png(filename=here('IPMresi', 'results', 'SummerSuckers_NbinsExtended_ResilienceFirstTimeStep.png'),
    height=5, width=10, units='in', res = 300)
par(mfrow=c(1,2))

plot(AmplificationFirst~Nbins, data=output, type='l', log='y',
     xlab='Number of bins', ylab='Amplification (1st timestep)')
plot(AttenuationFirst~Nbins, data=output, type='l', log='y',
     xlab='Number of bins', ylab='Attenuation (1st timestep)')
dev.off()





