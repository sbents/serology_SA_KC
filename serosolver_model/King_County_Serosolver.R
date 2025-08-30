##################################################################
# Project: Multiplex serology reveals age-specific immunodynamics 
# of respiratory pathogens in the wake of the COVID-19 pandemic 
# Author: Samantha Bents, sjbents@stanford.edu

##################################################################
# install packages 
install.packages("ggplot2")
install.packages("plyr") # Note plyr must be loaded before dplyr. 
install.packages("dplyr")
install.packages("tidyr")
install.packages("data.table")
install.packages("doParallel")
install.packages("coda")
install.packages("tidyverse")
install.packages("ggpubr")
install.package("doRNG")

# Install serosolver. 
devtools::install_github("seroanalytics/serosolver@e53b06335f6f8c4d4dd7ee390d04bfe411c662b9")

# load packages 
library(ggplot2)
library(plyr) 
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)
library(tidyverse)
library(serosolver)
library(ggpubr)
library(doRNG)

##################################################################
#################### Pediatric Samples   #########################
##################################################################

# Set MCMC parameters. 
##################################################################
set.seed(1)
mcmc_pars_use <- c(adaptive_iterations=10000, iterations=50000,proposal_ratio=1,thin_inf_hist=200)
prior_version = 2 # Indicates individual probability of infection is not independent of other individuals in the population. 

#################################################################
# Set working directory to where you have saved your antibody data. 
setwd("/Users/sambents/Desktop/serology")

################################################################
# Read in simulated data and start code here:
antibody_data <- read.csv("simulated_hcov_hku1_pediatric_serology.csv", row.names = 1)
head(antibody_data)
# Set up model. 
################################################################
time_per_year <- 1 # Samples per year 

## Read in the parameter table to be modified
par_tab <- read.csv(system.file("extdata", "par_tab_base.csv", package="serosolver"))

## Read in antigenic map and translate to quarterly time. Note this analysis does not 
## use an antigenic map, but the input is still required. 
antigenic_map <- read.csv("antigenic_map_quarters.csv") %>%
  mutate(inf_times = inf_times/4) %>%
  filter(inf_times <= max(antibody_data$sample_time))

## If estimating infection histories back to birth: 
possible_exposure_times <-  seq(min(antibody_data$birth), 2022,by=1)
print(possible_exposure_times)

## If estimating only recent times: 
#possible_exposure_times <- print(unique(antibody_data$sample_time))

# Filter antigenic map to align with exposure times.
antigenic_map <- antigenic_map %>% filter(inf_times %in% possible_exposure_times)

## Fix cross-reactivity parameters as this is a single-antigen analysis
par_tab[par_tab$names %in% c("cr_short","cr_long"),c("fixed")] <- 1
par_tab[par_tab$names %in% c("cr_short","cr_long"),c("values")] <- c(0.03,0.1)

## Fix antigenic seniority parameter as not used
par_tab[par_tab$names == "antigenic_seniority","fixed"] <- 1
par_tab[par_tab$names == "antigenic_seniority","values"] <- 0

## Fitting short-term boosting parameters
par_tab[par_tab$names == "boost_short","fixed"] <- 0
par_tab[par_tab$names == "boost_short","values"] <- 2

## Fitting short-term waning parameters
par_tab[par_tab$names == "wane_short","fixed"] <- 0
par_tab[par_tab$names == "wane_short","values"] <- 0.25/time_per_year

## NOT fitting boost long 
par_tab[par_tab$names %in% c("boost_long"),"fixed"] <- 1
par_tab[par_tab$names %in% c("boost_long"),"values"] <- 0

##  NOT fitting wane long 
par_tab[par_tab$names %in% c("wane_long"),"fixed"] <- 1
par_tab[par_tab$names %in% c("wane_long"),"values"] <- 0

## Infection model prior space parameters 
par_tab[par_tab$names %in% c("infection_model_prior_shape1","infection_model_prior_shape2"),"values"] <- c(.37, 3.5)

# Set upper bounds for boost and wane short 
par_tab[par_tab$names %in% c("boost_short"),"upper_bound"] <- max(antibody_data$measurement)
par_tab[par_tab$names %in% c("wane_short"),"upper_bound"] <-1

n_indiv = length(unique(antibody_data$individual))

## This places a standard normal prior on each of the coefficients. 
prior_func <- function(par_tab){
  coef_pars <- which(par_tab$names %like% "coef")
  par_names <- par_tab$names
  f <- function(pars){
    prior_p <- sum(dnorm(pars[coef_pars],0,1,log=TRUE))
    
    ## You can add your own priors on the model parameters here, for example, this places a log-normal prior with mean 2 on the boost_short parameter:
    # prior_p_boost <- dlnorm(pars[which(par_names == "boost_short")], log(2.7), 0.25, log=TRUE)
    # prior_p <- prior_p + prior_p_boost
    
    p1 <- dlnorm(pars["boost_short"],log(2), 0.5,log=TRUE)
    p2 <- dbeta(pars["wane_short"],1, 1,log=TRUE)
    p3 <- dbeta(pars["antigenic_seniority"],1, 1,log=TRUE)
    p4 <- dbeta(pars["cr_long"],1, 1,log=TRUE)
    p5 <- dbeta(pars["cr_short"],1, 1,log=TRUE)
    p6 <- dnorm(pars["obs_sd"],0, 100,log=TRUE)
    prior_p  + p1 + p2 + p3 + p4 + p5 + p6
  }
}


plot_all_outputs <- function(res, chain_wd, save_name, save_wd){
  
  ## Posterior estimates
  par_traces <- res$all_diagnostics$p_thetas[[1]]
  par_densities <- res$all_diagnostics$p_thetas[[2]]
  rho_traces <- res$all_diagnostics$p_thetas[[5]]
  rho_densities <- res$all_diagnostics$p_thetas[[6]]
  coef_traces <- res$all_diagnostics$p_thetas[[7]]
  coef_densities <- res$all_diagnostics$p_thetas[[8]]
  
  ######################################################################
  ## ANALYSES
  ######################################################################
  ## Read in MCMC chains
  chains <- load_mcmc_chains(chain_wd,par_tab,burnin=mcmc_pars_use["adaptive_iterations"])
  
  ## Compare parameter estimates to true values
  ## True values used for simulation shown with dashed lines, shaded regions show posterior densities
  p_par_ests <- chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no)) %>%
    filter(name %in% par_tab[par_tab$fixed == 0, "names"]) %>%
    rename(names=name,est=value) %>%
    left_join(par_tab %>% select(names,values)) %>%
    ggplot() + geom_density(aes(x=est,fill="Posterior"),alpha=0.5) + 
    scale_color_viridis_d(name="") +
    scale_fill_viridis_d(name="") +
    scale_linetype_manual(name="",values=c("True value"="dashed")) +
    facet_wrap(~names,scales="free") +
    xlab("Value") +
    ylab("Density") +
    theme_minimal()
  
  ## Compare model-estimated antibody levels to observations
  p_ab_predictions <- plot_antibody_predictions(chains$theta_chain,chains$inf_chain,settings=res$settings)
  p_ab_predictions1 <- p_ab_predictions[[4]]
  p_ab_predictions2 <- p_ab_predictions[[5]]
  
  
  ## Plot individual fits to data
  ## Orange shows posterior probability of infection in a given year
  ## Dots show data
  ## Purple shows model-predicted antibody levels
  ## Vertical dashed lines show true infection times
  p_fits <- plot_model_fits(chains$theta_chain,chains$inf_chain,individuals = 1:9,
                            known_infection_history= NULL, ## Set this to NULL for real data
                            settings=res$settings,orientation="longitudinal",expand_to_all_times = TRUE) + facet_wrap(~individual,ncol=3)
  
  ## Plot estimated attack rates vs. truth
  p_ar <- plot_attack_rates_pointrange(chains$inf_chain,settings = res$settings,
                                       # true_ar=all_simulated_data$attack_rates, ## Set this to NULL for real data
                                       true_ar= NULL,
                                       by_group=TRUE, 
                                       plot_den = FALSE,
                                       prior_pars = list(prior_version=1,infection_model_prior_shape1=3.5,infection_model_prior_shape2=.37))
  
  pdf(paste0(save_wd,"/",save_name,".pdf"),width=12,height=9)
  print(par_traces)
  print(par_densities)
  print(rho_traces)
  print(rho_densities)
  print(coef_traces)
  print(coef_densities)
  print(p_par_ests)
  print(p_ab_predictions1)
  print(p_ab_predictions2)
  print(p_fits)
  print(p_ar)
  dev.off()
}

# Set working directory where you want to save results. Note that you must create 
# an entirely separate folder (i.e. "results") to properly same model outputs 
# and figures. 

setwd("/Users/sambents/Desktop/NIH/serology/results")
# Create folders for outputs. 
dir.create("hku1_simulated_output") # Folder to save chain outputs
dir.create("plot_wd") # Folder to save model plot diagnositcs

######################################################################
## RUN MODEL
######################################################################
# Note that a browser will interrupt to the run, allowing for potential debugging.
# Entering 'c'  and then pressing "enter" will allow the model to keep running to completion. 
res <- serosolver(par_tab, 
                  antibody_data, 
                  #  demographics=demographics, # not using demographics 
                  antigenic_map= NULL,       # not using antigenic map 
                  fixed_inf_hists= NULL, ## Set fixed infection states
                  prior_func=prior_func,
                  possible_exposure_times= possible_exposure_times,
                  measurement_bias = NULL,
                  filename="hku1_simulated_output/phirst_sim_recovery", # change file name based on antigen of interest
                  n_chains=4, ## Run 4 chains
                  parallel=TRUE, ## Run in parallel
                  mcmc_pars=mcmc_pars_use, ## Some MCMC control settings
                  verbose=TRUE,
                  data_type=1,
                  start_level="none")

# Print results and model diagnostics. 
print(res$all_diagnostics)

################################################################
## SAVE MODEL DIAGNOSTICS: 
plot_all_outputs(res, "hku1_simulated_output", "hku1_simulated_output",  "plot_wd")

