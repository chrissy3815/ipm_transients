###########################################################################
### Comparing resilience metrics from MPMs and IPMs built for the same population
### Christina Hernandez
### Sept 2025
###########################################################################

library(popdemo)
library(here)

# Initialize resilience metrics dataframe:
resilience_metrics<- data.frame(species=c("Opuntia", "Opuntia", "Lupine", "Lupine"),
                                MPMorIPM=c("MPM", "IPM", "MPM", "IPM"),
                                Amp_first=vector("numeric", 4),
                                Amp_max=vector("numeric", 4),
                                Att_first=vector("numeric", 4),
                                Att_max=vector("numeric", 4),
                                Recov_t=vector("numeric", 4),
                                Lambda=vector("numeric", 4))

###########################################################################
## Opuntia
###########################################################################
opuntia_MPM<- readRDS(here("IPMresi", "data", "opuntia_MPM.rds"))
opuntia_IPM<- readRDS(here("IPMresi", "data", "Opuntia_IPM.rds"))

I<- which(resilience_metrics$species=="Opuntia" & resilience_metrics$MPMorIPM=="MPM")
resilience_metrics$Amp_first[I]<- popdemo::reac(opuntia_MPM, bound='upper')
resilience_metrics$Att_first[I]<- popdemo::reac(opuntia_MPM, bound='lower')
resilience_metrics$Amp_max[I]<- popdemo::maxamp(opuntia_MPM)
resilience_metrics$Att_max[I]<- popdemo::maxatt(opuntia_MPM)
resilience_metrics$Recov_t[I]<- max(popdemo::convt(opuntia_MPM))
resilience_metrics$Lambda[I]<- max(Re(eigen(opuntia_MPM, only.values = TRUE)$values))

I<- which(resilience_metrics$species=="Opuntia" & resilience_metrics$MPMorIPM=="IPM")
resilience_metrics$Amp_first[I]<- popdemo::reac(opuntia_IPM, bound='upper')
resilience_metrics$Att_first[I]<- popdemo::reac(opuntia_IPM, bound='lower')
resilience_metrics$Amp_max[I]<- popdemo::maxamp(opuntia_IPM)
resilience_metrics$Att_max[I]<- popdemo::maxatt(opuntia_IPM)
resilience_metrics$Recov_t[I]<- max(popdemo::convt(opuntia_IPM))
resilience_metrics$Lambda[I]<- max(Re(eigen(opuntia_IPM, only.values = TRUE)$values))

###########################################################################
## Lupine
###########################################################################
lupine_MPM<- readRDS(here("IPMresi", "data", "lupine_MPM.rds"))
lupine_IPM<- readRDS(here("IPMresi", "data", "lupine_IPM.rds"))

I<- which(resilience_metrics$species=="Lupine" & resilience_metrics$MPMorIPM=="MPM")
resilience_metrics$Amp_first[I]<- popdemo::reac(lupine_MPM, bound='upper')
resilience_metrics$Att_first[I]<- popdemo::reac(lupine_MPM, bound='lower')
resilience_metrics$Amp_max[I]<- popdemo::maxamp(lupine_MPM)
resilience_metrics$Att_max[I]<- popdemo::maxatt(lupine_MPM)
resilience_metrics$Recov_t[I]<- max(popdemo::convt(lupine_MPM))
resilience_metrics$Lambda[I]<- max(Re(eigen(lupine_MPM, only.values = TRUE)$values))

I<- which(resilience_metrics$species=="Lupine" & resilience_metrics$MPMorIPM=="IPM")
resilience_metrics$Amp_first[I]<- popdemo::reac(lupine_IPM, bound='upper')
resilience_metrics$Att_first[I]<- popdemo::reac(lupine_IPM, bound='lower')
resilience_metrics$Amp_max[I]<- popdemo::maxamp(lupine_IPM)
resilience_metrics$Att_max[I]<- popdemo::maxatt(lupine_IPM)
resilience_metrics$Recov_t[I]<- max(popdemo::convt(lupine_IPM))
resilience_metrics$Lambda[I]<- max(Re(eigen(lupine_IPM, only.values = TRUE)$values))

###########################################################################
## Scatter plot with 1:1 line
###########################################################################

colz<- c('pink', 'darkred', 'lightblue', 'darkblue', 'black', 'grey')

png(filename=here("IPMresi","results","CompareMPMvIPM_OpuntiaLupine.png"),
    height=5, width=10, units='in', res=300)
par(mfrow=c(1,2))

I<- which(resilience_metrics$species=="Opuntia" & resilience_metrics$MPMorIPM=="MPM")
exes<- as.vector(resilience_metrics[I,3:ncol(resilience_metrics)])
I<- which(resilience_metrics$species=="Opuntia" & resilience_metrics$MPMorIPM=="IPM")
whys<- as.vector(resilience_metrics[I,3:ncol(resilience_metrics)])
plot(exes, whys, pch=19, col=colz,
     xlab='MPM metrics', ylab='IPM metrics', main="Opuntia")
abline(0,1)

I<- which(resilience_metrics$species=="Lupine" & resilience_metrics$MPMorIPM=="MPM")
exes<- as.vector(resilience_metrics[I,3:ncol(resilience_metrics)])
I<- which(resilience_metrics$species=="Lupine" & resilience_metrics$MPMorIPM=="IPM")
whys<- as.vector(resilience_metrics[I,3:ncol(resilience_metrics)])
plot(exes, whys, pch=19, col=colz,
     xlab='MPM metrics', ylab='IPM metrics', main="Lupine")
abline(0,1)

dev.off()


###########################################################################
## bar plot of ratios
###########################################################################

colz<- c('pink', 'darkred', 'lightblue', 'darkblue', 'black', 'grey')
barlabels<- names(resilience_metrics)[3:ncol(resilience_metrics)]

png(filename=here("IPMresi","results","CompareMPMvIPM_OpuntiaLupine_Barplot.png"),
    height=5, width=10, units='in', res=300)
par(mfrow=c(1,2))

I<- which(resilience_metrics$species=="Opuntia" & resilience_metrics$MPMorIPM=="MPM")
denom<- as.numeric(resilience_metrics[I,3:ncol(resilience_metrics)])
I<- which(resilience_metrics$species=="Opuntia" & resilience_metrics$MPMorIPM=="IPM")
num<- as.numeric(resilience_metrics[I,3:ncol(resilience_metrics)])
toplot<- num/denom
# flip the attenuation items so that more extreme attenuation looks bigger:
toplot[c(3,4)]<- 1/toplot[c(3,4)]

barplot(toplot, col=colz, names.arg=barlabels,
     xlab='Resilience metrics', ylab='Ratio IPM vs MPM', main="Opuntia")

I<- which(resilience_metrics$species=="Lupine" & resilience_metrics$MPMorIPM=="MPM")
denom<- as.numeric(resilience_metrics[I,3:ncol(resilience_metrics)])
I<- which(resilience_metrics$species=="Lupine" & resilience_metrics$MPMorIPM=="IPM")
num<- as.numeric(resilience_metrics[I,3:ncol(resilience_metrics)])
toplot<- num/denom
# flip the attenuation items so that more extreme attenuation looks bigger:
toplot[c(3,4)]<- 1/toplot[c(3,4)]

barplot(toplot, col=colz, names.arg=barlabels,
        xlab='Resilience metrics', ylab='Ratio IPM vs MPM', main="Lupine")

dev.off()

## log-scaled version:
colz<- c('pink', 'darkred', 'lightblue', 'darkblue', 'black', 'grey')
barlabels<- names(resilience_metrics)[3:ncol(resilience_metrics)]

png(filename=here("IPMresi","results","CompareMPMvIPM_OpuntiaLupine_BarplotLOG10.png"),
    height=5, width=10, units='in', res=300)
par(mfrow=c(1,2))

I<- which(resilience_metrics$species=="Opuntia" & resilience_metrics$MPMorIPM=="MPM")
denom<- as.numeric(resilience_metrics[I,3:ncol(resilience_metrics)])
I<- which(resilience_metrics$species=="Opuntia" & resilience_metrics$MPMorIPM=="IPM")
num<- as.numeric(resilience_metrics[I,3:ncol(resilience_metrics)])
toplot<- num/denom
# flip the attenuation items so that more extreme attenuation looks bigger:
toplot[c(3,4)]<- 1/toplot[c(3,4)]

barplot(log10(toplot), col=colz, names.arg=barlabels,
        xlab='Resilience metrics', ylab='log10[Ratio IPM vs MPM]', main="Opuntia")

I<- which(resilience_metrics$species=="Lupine" & resilience_metrics$MPMorIPM=="MPM")
denom<- as.numeric(resilience_metrics[I,3:ncol(resilience_metrics)])
I<- which(resilience_metrics$species=="Lupine" & resilience_metrics$MPMorIPM=="IPM")
num<- as.numeric(resilience_metrics[I,3:ncol(resilience_metrics)])
toplot<- num/denom
# flip the attenuation items so that more extreme attenuation looks bigger:
toplot[c(3,4)]<- 1/toplot[c(3,4)]

barplot(log10(toplot), col=colz, names.arg=barlabels,
        xlab='Resilience metrics', ylab='log10[Ratio IPM vs MPM]', main="Lupine")

dev.off()

