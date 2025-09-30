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

source(here('IPMresi', 'code', 'IPMresi_SoaySheepIPMfunctions.R'))

## Calculate resilience metrics across a range of number of mesh points: -------
min.size <- 1.60
max.size <- 3.70
# model parameters for building the kernel are stored in m.par.true

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
  thisK<- mk_K(n_bins[i], m.par.true, min.size, max.size)

  # Calculate the resilience metrics from popdemo:
  output$AmplificationFirst[i]<- popdemo::reac(thisK$K, bound='upper')
  output$AmplificationMax[i]<- popdemo::maxamp(thisK$K)
  output$AttenuationFirst[i]<- popdemo::reac(thisK$K, bound='lower')
  output$AttenuationMax[i]<- popdemo::maxatt(thisK$K)
  output$RecoveryTime[i]<- max(popdemo::convt(thisK$K))
}

png(filename=here('IPMresi', 'results', 'SoaySheep_NbinsResilience2.png'),
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

min.size <- 1.60
max.size <- 3.70
# model parameters for building the kernel are stored in m.par.true

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
  cat(j, '\n')

  # generate the kernel at this number of bins:
  thisK<- mk_K(testgrid$n_bins[i], m.par.true, min.size, max.size)

  thisoutput<- reac_with_gaussian_vectors(thisK$K, testgrid$n_bins[i], testgrid$sd_vals[j])
}

temp<- as.data.frame(results)
upper<- unlist(temp$upper)
lower<- unlist(temp$lower)

png(here('IPMresi', 'results', 'SoaySheepExample_levelplots.png'),
    height=5, width=10, units='in', res=300)
par(mfrow=c(1,2))

subfigz<- list()
subfigz[[1]]<- lattice::levelplot(upper~output$n_bins+output$sd_vals, xlab='N bins', ylab='SD of Gaussian',
                   main='Amplification (1st time step)')
subfigz[[2]]<- lattice::levelplot(lower~output$n_bins+output$sd_vals, xlab='N bins', ylab='SD of Gaussian',
                   main='Attenuation (1st time step)')

# Plot them:
lay <- c(1,2)
gridExtra::grid.arrange(grobs = subfigz, layout_matrix = lay)

dev.off()

