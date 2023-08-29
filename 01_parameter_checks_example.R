# Basic checks of the parameter values in the VL model

# Author: Luc Coffeng
# Date created: 4 Aug 2023 (updated by Ananthu James on 10 Aug 2023) 

# Prep session ----
rm(list = ls())
library(rstudioapi)  # requires RStudio 

base_dir <- dirname(getActiveDocumentContext()$path) 

setwd(file.path(base_dir, "00_Source"))
source("param_values.r")


# Function to analytically calculate the average detection delay (dwell time)
# for a population of undetected symptomatic cases subject to some background
# mortality rate, excess mortality due to VL, the detection rate, and the shape
# of Erlang distribution for time until death. The function returns a vector of
# values, which if summed up represent the average undetected dwell time in
# years. The individual vector elements represent the relative size of the
# population being detection in the different Erlang stages.
calc_undetect_dwell <- function(mu, muVL, rho, shape = 3) {
  
  # Calculate dwell time per compartment
  dwell_compartment <- 1 / (mu + shape * muVL + rho)
  
  # Calculate proportions making it to specific exits
  prop <- matrix(NA,
                 nrow = shape,
                 ncol = 3,
                 dimnames = list(paste0("Erlang_", 1:shape),
                                 c("Detected",
                                   "Progressed",
                                   "Death due to other causes")))
  par_vector <- c(rho, shape * muVL, mu)
  
  prop[1, ] <- par_vector / sum(par_vector)
  if (shape > 1) {
    for(i in 2:shape) {
      prop[i, ] <- par_vector / sum(par_vector) * prop[(i - 1), 2]  
    }
  }
  
  # sum(prop) - sum(prop[1:(shape - 1), "Progressed"])  # should be 1.0
  
  prop[, "Detected"] * 1:shape * dwell_compartment / sum(prop[, "Detected"])
  
}


# Function to numerically calibrate detection rate, given some desired undetected
# dwell time (inverse function of calc_detect_dwell()).
# N.B. this process is sensitive to initial values!
calibrate_rho <- function(target_dur,  # duration in days!
                          rho_init = param_values["rhoI1"] * 2,
                          mu = param_values["mu"],
                          muVL = param_values["muVL"],
                          shape = 3) {
  
  temp <- optim(par = log(rho_init),
                fn = function(x, target) {
                  (sum(calc_undetect_dwell(mu = mu,
                                           muVL = muVL,
                                           rho = exp(x),
                                           shape = shape)) * 365 - target)^2
                }, method = "BFGS", target = target_dur)
  
  # if(temp$value > 1) stop("Non-convergence; provide better inits")
  # temp$par <- exp(temp$par)
  # temp
  if(temp$value > 1) return(NA)
  exp(temp$par)
  
}


# Function to calculate average duration of asymptomatic infection
calc_dur_asymp <- function(fh, fs, rhoL1, rhoL2, rhoH) {
  1/rhoL1 + 1/rhoL2 + (fh * (1 - fs) * 1/rhoH) / (1 - fh * fs)
}


# Function to calculate duration of L2, given duration of asymptomatic infection
calc_rhoL2_inv <- function(dur_asymp_target, fh, fs, rhoL1, rhoH) {
  dur_asymp_target - 1 / rhoL1 - fh*(1-fs) / ((1 - fh * fs) * rhoH)
}


# Example calibration of detection rate ----
# Note: the detection rates rhoI1 and rhoI2 were chosen in two steps. First,
# parameters were calibrated to a target detection delay, using the function
# 'calibrate_rho()', of which you will find an example below. Second, the
# calibrated rates were re-expressed as the closest rational number (i.e., a
# fraction of integers).

## Example of calibrating detection rate rho to some target detection delay
## among detected cases, given mortality and excess mortality among not yet
## detected cases
rho_proposal <- with(as.list(param_values),
                       calibrate_rho(target_dur = 90,  # duration in days!
                                     rho_init = 1,  # may need to be adjusted
                                     mu = mu,
                                     muVL = muVL,
                                     shape = 3))
dur_days <- round(365 / rho_proposal, 0)

sum(with(as.list(param_values),
         calc_undetect_dwell(mu = mu,
                             muVL = muVL,
                             rho = 365 / dur_days,
                             shape = 3))) * 365  # should be ~90 days


## Check that default parameter values make sense ----
# Average duration of asymptomatic infection
with(as.list(param_values),
     calc_dur_asymp(fh = fh,
                    fs = fs,
                    rhoL1 = rhoL1,
                    rhoL2 = rhoL2,
                    rhoH = rhoH)) * 365  # should be about 200 days

# Average duration of L2, given duration of asymptomatic infection
with(as.list(param_values),
     calc_rhoL2_inv(dur_asymp_target = 200 / 365,
                    fh = fh,
                    fs = fs,
                    rhoL1 = rhoL1,
                    rhoH = rhoH)) * 365


### END OF CODE ### ----
