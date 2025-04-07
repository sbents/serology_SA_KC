# Load packages
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)
library(tidyverse)
library(serosolver)

# Set MCMC parameters. 
##########################################################################
set.seed(1)
mcmc_pars_use <- c(adaptive_iterations=100000, iterations=300000,proposal_ratio=1,thin_inf_hist=200)
prior_version = 2 # Indicates individual probability of infection is not independent of other individuals in the population. 
##########################################################################
setwd("/Users/sambents/Desktop/NIH/serology/unprocessed")


# Load unprocessed data in. 
##########################################################################
hai = read.csv("HAI.CSV") # serology data
meta = read.csv("meta.csv") # meta data 

sero_samp = left_join(hai, meta, by = c("ind_id", "year")) %>%
  mutate(titre = log2(flub_victoria_1/5)) %>%
  mutate(collection_date = as.Date(collectiondate, "%m/%d/%Y") , collection_year = year(collection_date), collection_month = month(collection_date)) %>%
  mutate(age =  collection_year - dob_year) %>%
  mutate(half = ifelse(collection_month < 7, 1, 2)) %>%
  mutate(quarter = ifelse(collection_month < 4, 1,
                          ifelse(collection_month >= 4 & collection_month < 7, 2,
                                 ifelse(collection_month >= 7 & collection_month < 10,3, 4)))) %>%
  dplyr::select(ind_id, year, age, titre, dob_year, collection_year, half, quarter) %>%
# filter(age > 11) %>% # adults 
# filter(age >= 5 & age < 11) %>% # middle-aged children
  filter(age < 5) %>% 
  drop_na()

ind_id= print(unique(sero_samp$ind_id))
n_ind_id = seq(1, length(ind_id), 1)
ids = data.frame(ind_id, n_ind_id)
sero_samp = left_join(sero_samp, ids, by = "ind_id")

antibody_data = sero_samp %>% 
  mutate(sample_time = as.numeric(collection_year )) %>%
  mutate(virus = 2016, biomarker_group = 1, 
         obs_type = 1, individual = n_ind_id, biomarker_id = 2016, repeat_number = 1, measurement = titre) %>%
  drop_na(measurement) %>%
  mutate(birth = dob_year) %>% # This model does not use information on DOB so this column can be arbitrarily filled in. 
  distinct(sample_time, individual, .keep_all = TRUE) %>% # one sample per time period 
  dplyr::select(individual, sample_time, birth, repeat_number, biomarker_id, biomarker_group, measurement) %>%
  as.data.frame()
head(antibody_data)

# Set up model. 
##########################################################################
time_per_year <- 1
min(antibody_data$birth)

## Read in the parameter table to be modified
par_tab <- read.csv(system.file("extdata", "par_tab_base.csv", package="serosolver"))
antigenic_map <- read.csv("antigenic_map_quarters.csv") %>%
  mutate(inf_times = inf_times/4) %>%
  filter(inf_times <= max(antibody_data$sample_time))

min(sero_samp$dob_year)
## If estimating back to birth: 
possible_exposure_times <- c(seq(min(sero_samp$dob_year), 2015,by= 12), seq(2016, 2018,by=1))

possible_exposure_times <- c(seq(min(sero_samp$dob_year), 2003,by= 12), seq(2004, 2018,by=1))
possible_exposure_times <-  seq(2012, 2018,by=1)

possible_exposure_times <- c(seq(min(sero_samp$dob_year), 2005,by= 12), seq(2006, 2018,by=1))
print(possible_exposure_times)

## If estimating only recent times: 
#possible_exposure_times <- print(unique(antibody_data$sample_time))
antigenic_map <- antigenic_map %>% filter(inf_times %in% possible_exposure_times)

## Fix cross-reactivity parameters as this is a single-antigen analysis
par_tab[par_tab$names %in% c("cr_short","cr_long"),c("fixed")] <- 1
par_tab[par_tab$names %in% c("cr_short","cr_long"),c("values")] <- c(0.03,0.1)


## Fix antigenic seniority parameter as not used
par_tab[par_tab$names == "antigenic_seniority","fixed"] <- 1
par_tab[par_tab$names == "antigenic_seniority","values"] <- 0

## Using both long- and short-term antibody kinetics
par_tab[par_tab$names == "boost_short","fixed"] <- 0
par_tab[par_tab$names == "boost_short","values"] <- 2

par_tab[par_tab$names == "wane_short","fixed"] <- 0
par_tab[par_tab$names == "wane_short","values"] <- 0.25/time_per_year

#  we are NOT fitting boost long but not wane long 
par_tab[par_tab$names %in% c("boost_long"),"fixed"] <- 1
par_tab[par_tab$names %in% c("boost_long"),"values"] <- 0

par_tab[par_tab$names %in% c("wane_long"),"fixed"] <- 1
par_tab[par_tab$names %in% c("wane_long"),"values"] <- 0

## Infection model prior -- remember you can check what the attack rate prior from this looks like with
par_tab[par_tab$names %in% c("infection_model_prior_shape1","infection_model_prior_shape2"),"values"] <- c(.37, 3.5)
n_indiv = length(unique(antibody_data$individual))

## IMPORTANT -- prior function for the model coefficients. This places a standard normal prior on each of the coefficients. 
## This is quite strong but shrinks their effects to 0
prior_func <- function(par_tab){
  coef_pars <- which(par_tab$names %like% "coef")
  par_names <- par_tab$names
  f <- function(pars){
    prior_p <- sum(dnorm(pars[coef_pars],0,1,log=TRUE))
    
    ## You can add your own priors on the model parameters here, for example, this places a log-normal prior with mean 2 on the boost_short parameter:
    # prior_p_boost <- dlnorm(pars[which(par_names == "boost_short")], log(2.7), 0.25, log=TRUE)
    # prior_p <- prior_p + prior_p_boost
    # p1 <- dlnorm(pars["boost_long"],log(2), 0.5,log=TRUE)
    p2 <- dlnorm(pars["boost_short"],log(2), 0.5,log=TRUE)
    p3 <- dbeta(pars["wane_short"],1, 1,log=TRUE)
    p4 <- dbeta(pars["antigenic_seniority"],1, 1,log=TRUE)
    p5 <- dbeta(pars["cr_long"],1, 1,log=TRUE)
    p6 <- dbeta(pars["cr_short"],1, 1,log=TRUE)
    p7 <- dnorm(pars["obs_sd"],0, 100,log=TRUE)
    prior_p  + p2 + p3 + p4 + p5 + p6 + p7
  }
}


plot_all_outputs <- function(res, chain_wd, save_name, save_wd){
  ## Have a look at the posterior estimates
  ## there are now extensive MCMC diagnostics in res$all_diagnostics
  ## this includes all parameter estimates, ESS, Rhat, trace plots, posterior densities etc
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
    # geom_vline(aes(xintercept=values,linetype="True value"),col="black") + ## Comment out for real data
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
  
  ## Plot model-estimated antibody kinetics for a single infection, stratified by group
  p_model <- plot_estimated_antibody_model(chains$theta_chain, 
                                           solve_times=possible_exposure_times,
                                           settings=res$settings,
                                           nsamp = 1000,by_group = TRUE)
  
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
  print(p_model)
  dev.off()
}

setwd("/Users/sambents/Desktop/NIH/serology/final_phirst")
getwd()

dir.create("chains_ah1_under5")
dir.create("chains_ah1_mid")
dir.create("chains_ah1_adult")

dir.create("chains_ah3_young")
dir.create("chains_ah3_mid")
dir.create("chains_ah3_adult")

dir.create("chains_byam_young")
dir.create("chains_byam_mid")
dir.create("chains_byam_adult")

dir.create("chains_bvic_adult")
dir.create("chains_bvic_mid")
dir.create("chains_bvic_young")

dir.create("plot_wd")

######################################################################
## RUN 1: no antigenic map, no rho, estimating back to birth, AH3, < 3 yo is group 1 
######################################################################

res = NULL
res <- serosolver(par_tab, 
                  antibody_data, 
                 #  demographics=demographics,
                  antigenic_map= NULL,
                  fixed_inf_hists= NULL, ## Set fixed infection states
                  prior_func=prior_func,
                  #  possible_exposure_times= antigenic_map$inf_times,
                  possible_exposure_times= possible_exposure_times,
                  measurement_bias = NULL,
                  filename="chains_bvic_young/phirst_sim_recovery", 
                  n_chains=3, ## Run 3 chains
                  parallel=TRUE, ## Run in parallel
                  mcmc_pars=mcmc_pars_use, ## Some MCMC control settings
                  verbose=TRUE,
                  data_type=1,
                  start_level="none")
print(res$all_diagnostics)

plot_all_outputs(res, "chains_bvic_young","chains_bvic_young","plot_wd")


