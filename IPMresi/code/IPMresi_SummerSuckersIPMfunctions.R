###########################################################################
### Modeling life history of white suckers and summer suckers: summer sucker version
### Christina Hernandez
### July 2023
### Reusing IPM for other analysis - September 2025
###########################################################################
## Change log
# survival model changed to 4 param with shallow slope
# egg viability changed from 0.1 to 0.0035 -> this led to bimodal size distribution

# Load necessary packages:
library(here)

###########################################################################
### Mostly following Pierce et al. 2023, model building:
###########################################################################

# IPM parameters (from summer suckers life history modeling, take fitted values
# and a few literature values)
ss_m_par <- list(
  ## Growth parameters
  grow_rate = 0.1384936, # growth rate
  Linf  = 350, # ss_growth_params$Linf, # maximum length in mm
  grow_sd   = 25,  # growth sd
  ## Survival parameters a
  surv_min =  0.003,
  surv_mid = 0.62,
  surv_max = 0.65,
  surv_alpha = 112,
  surv_alpha2 = 220,
  surv_beta = -25,
  surv_beta2 = -25,
  ## Size of age-1 individuals:
  recruit_mean = 112, # mean size of age-1 individuals
  recruit_sd = 25, # same as grow_sd
  ## PLACEHOLDER:
  egg_viable = 0.022,
  ## Estimated from fecundity data
  egg_logslope = 3.4,
  egg_logintercept = -11.35,
  ## Spawning Probability
  pb_max = 0.6524654, # maximum probability of spawning
  pb_k = 0.1834677, # rate of increase of spawning probability with size
  pb_midsize = 163.5221, # size at which 50% of individuals spawn
  ## YOY survival probability:
  s0= 0.1 #
)

##########################
## Section 2: Model Set-up
##########################

## Growth function
# given you are size z now returns the pdf of size z1 next time
# computed from von Bertanaffy equation z(t) = L_inf(1-e^K(t-t0))
# to find z(t+1) = L_inf*(1-e^(-K)) + e^(-K)*z(t)

ss_g_z1z <- function(z1, z, ss_m_par) {
  mu <- ss_m_par$Linf * (1 - exp(- ss_m_par$grow_rate)) +
    exp(-ss_m_par$grow_rate) * z           # mean size next year
  sig <- ss_m_par$grow_sd                       # sd about mean
  p_den_grow <- dnorm(z1, mean = mu, sd = sig)    # pdf that you are size z1
  # given you were size z
  return(p_den_grow)
}

# Adult Survival function a, 4-parameter logistic
# ss_s_z <- function(z, ss_m_par) {
#   ss_m_par$surv_min + (ss_m_par$surv_max - ss_m_par$surv_min) /
#     (1 + exp(ss_m_par$surv_beta * (log(z) - log(ss_m_par$surv_alpha))))
# }


## Adult Survival function b, 4-parameter double logistic
ss_s_z <- function(z, ss_m_par) {
  ss_m_par$surv_min + ((ss_m_par$surv_mid - ss_m_par$surv_min) /
                         (1 + exp(ss_m_par$surv_beta * (log(z) - log(ss_m_par$surv_alpha))))) +
    ((ss_m_par$surv_max - ss_m_par$surv_mid) /
       (1 + exp(ss_m_par$surv_beta2 * (log(z) - log(ss_m_par$surv_alpha2)))))
}

## Reproduction, log-linear
ss_eggs_z <- function(z, ss_m_par) { # Eggs produced (note: data are in thousands)
  eggz<- exp(ss_m_par$egg_logslope*log(z) + ss_m_par$egg_logintercept)
  return(eggz)
}

## Probability of spawning, 3-parameter logistic
ss_pb_z<- function(z, ss_m_par){
  ss_m_par$pb_max / (1 + exp(-ss_m_par$pb_k * (z - ss_m_par$pb_midsize)))
}

## Recruit size pdf
ss_c_1z1 <- function(z1, ss_m_par) {
  mu <- ss_m_par$recruit_mean
  sig <- ss_m_par$recruit_sd
  p_den_recruit <- dnorm(z1, mean = mu, sd = sig)
  return(p_den_recruit)
}

#####################################################
## Section 3 - Build IPM kernels F and P
#####################################################

## Fecundity Kernel
ss_f_z1z <- function(z1, z, ss_m_par) {
  age1_dist <- ss_pb_z(z, ss_m_par) * ss_eggs_z(z, ss_m_par) *
    ss_m_par$egg_viable * ss_m_par$s0
  #returns fecundity kernel (as a matrix). Recruits= F.dot(n*delta_z)
  return(outer(ss_c_1z1(z1, ss_m_par), age1_dist))
}

## Growth and Survival Kernel
ss_p_z1z <- function(z1, z, ss_m_par) {
  N<- length(z1)
  delta_z<- z1[2]-z1[1]
  g_matrix <- matrix(0, N, N)
  for (x in 1:N) {
    g_matrix[, x] <- ss_g_z1z(z, rep(z[x], times = N), ss_m_par)
    g_matrix[, x] <- g_matrix[, x] / (sum(g_matrix[, x]) * delta_z)
  }
  return(g_matrix %*% diag(ss_s_z(z, ss_m_par)))
}

make_K_suckers<- function(m, m.par, L, U){
  h <- (U-L)/m # integration bin width
  ss_meshpts <-  L + (1:m)*h - h/2

  ss_Pmat<- h*ss_p_z1z(ss_meshpts, ss_meshpts, m.par)
  ss_Fmat<- h*(ss_f_z1z(ss_meshpts, ss_meshpts, m.par))
  ss_Kmat<- ss_Pmat+ss_Fmat
}

## Build the deterministic kernels ##
m <- 300 # number of meshpoints: bins for the integration
L <- 0.00   # lower size limit in mm
U <- 600.00    # upper size limit in mm - must be larger than Linf

K_suckers<- make_K_suckers(m, ss_m_par, L, U)


