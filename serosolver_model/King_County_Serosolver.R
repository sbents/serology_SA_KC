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
mcmc_pars_use <- c(adaptive_iterations=100000, iterations=500000,proposal_ratio=1,thin_inf_hist=200)
prior_version = 2 # Indicates individual probability of infection is not independent of other individuals in the population. 
##########################################################################
setwd("/Users/sambents/Desktop/NIH/serology/unprocessed")

# Load unprocessed data in. 
##########################################################################
cov = read.csv("serosamples_CoV_22Mar2024.csv") 
resp = read.csv("serosamples_Resp_22Mar2024.csv")
print(unique(adult$Assay))
meta = read.csv("demos_children.csv")  %>%
  mutate(Sample = sample_id) # meta data 

sero_samp = left_join(cov, meta, by = c("Sample")) %>%
  filter(Assay == "CoV-HKU1 Spike") %>%
  mutate(titre = log10(titer_mean)) %>%
  mutate(age = age_raw) %>%
  mutate(collection_date = as.Date(blood_date) , collection_year = year(blood_date), collection_month = month(blood_date)) %>%
  dplyr::select(Sample, age, titre, collection_year) %>%
  # filter(age > 11) %>% # specify age group
 #  filter(age >= 5 & age < 11) %>%
   filter(age < 5) %>% 
  filter(age > .5) %>% 
  filter(titre > 0) %>%
  drop_na()
head(sero_samp)

Sample = print(unique(sero_samp$Sample))

n_ind_id = seq(1, length(Sample), 1)
ids = data.frame(Sample, n_ind_id)
sero_samp = left_join(sero_samp, ids, by = "Sample")

antibody_data = sero_samp %>% 
  mutate(sample_time = as.numeric(collection_year )) %>%
  mutate(virus = 2020, biomarker_group = 1, 
         obs_type = 1, individual = n_ind_id, biomarker_id = 2020, repeat_number = 1, measurement = titre) %>%
  drop_na(measurement) %>%
  mutate(birth = collection_year - age) %>% # This model does not use information on DOB so this column can be arbitrarily filled in. 
  distinct(sample_time, individual, .keep_all = TRUE) %>% # one sample per time period 
  dplyr::select(individual, sample_time, birth, repeat_number, biomarker_id, biomarker_group, measurement) %>%
  mutate(birth = round(birth, digits = 0)) %>%
  group_by(individual) %>%
  mutate(birth = min(birth)) %>%
  ungroup() %>%
  mutate(measurement = round(measurement, digits = 0)) %>% 
  mutate(measurement = measurement - min(measurement)) %>%
  as.data.frame() %>%
  drop_na()
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

## If estimating back to birth: 
possible_exposure_times <-  seq(min(antibody_data$birth), 2022,by=1)
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

par_tab[par_tab$names %in% c("boost_short"),"upper_bound"] <- max(antibody_data$measurement)
par_tab[par_tab$names %in% c("wane_short"),"upper_bound"] <-1

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

setwd("/Users/sambents/Desktop/NIH/serology/final_kc_log10")
getwd()

dir.create("hku1_young")
dir.create("oc43_young")
dir.create("229e_young")
dir.create("nl63_young")
dir.create("sarsn_young")
dir.create("sarsrbd_young")
dir.create("sarsspike_young")
dir.create("rsv_young")
dir.create("fluah3_young")
dir.create("fluah1_young")
dir.create("flubvic_young")
dir.create("flubyam_young")

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
                  filename="hku1_young/phirst_sim_recovery", 
                  n_chains=4, ## Run 4 chains
                  parallel=TRUE, ## Run in parallel
                  mcmc_pars=mcmc_pars_use, ## Some MCMC control settings
                  verbose=TRUE,
                  data_type=1,
                  start_level="none")
print(res$all_diagnostics)

plot_all_outputs(res, "hku1_young","hku1_young","plot_wd")




##################################################################
##################################################################
#### Adults ######################################################

# Set MCMC parameters. 
##########################################################################
set.seed(1)
mcmc_pars_use <- c(adaptive_iterations=200000, iterations=700000,proposal_ratio=1,thin_inf_hist=200)
prior_version = 2 # Indicates individual probability of infection is not independent of other individuals in the population. 

##########################################################################
setwd("/Users/sambents/Desktop/NIH/serology/unprocessed")


# Load unprocessed data in. 
##########################################################################
adult = read.csv("serology_samples_all.csv") %>%
  filter(age > 18) %>%
  mutate(Sample = sample_ID)
print(unique(adult$Assay))

sero_samp = adult %>%
  filter(Assay == "SARS-CoV-2 RBD") %>%
#  mutate(titre = log2(as.numeric(mean_titer)/5)) %>%
  mutate(titre = log(mean_titer)) %>%
  dplyr::select(Sample, age, titre, year) %>%
  filter(	Sample != 1302773) %>%
  filter(	Sample != 988161) %>%
  drop_na()
head(sero_samp)

count = sero_samp %>%
  count(Sample) %>%
  filter(n == 2)

Sample = print(unique(count$Sample))
Sample1 = print(unique(count$Sample))

n_ind_id = seq(1, length(Sample), 1)
ids = data.frame(Sample, n_ind_id)


sero_samp = left_join(sero_samp %>% filter(Sample %in% Sample1), ids, by = "Sample")

antibody_data = sero_samp %>% 
  mutate(sample_time = year) %>%
  mutate(virus = 2021, biomarker_group = 1, 
         obs_type = 1, individual = n_ind_id, biomarker_id = 2021, repeat_number = 1, measurement = titre) %>%
  drop_na(measurement) %>%
  mutate(birth = year - age) %>% # This model does not use information on DOB so this column can be arbitrarily filled in. 
  distinct(sample_time, individual, .keep_all = TRUE) %>% # one sample per time period 
  dplyr::select(individual, sample_time, birth, repeat_number, biomarker_id, biomarker_group, measurement) %>%
  mutate(birth = round(birth, digits = 0)) %>%
  group_by(individual) %>%
  mutate(birth = min(birth)) %>%
  ungroup() %>%
  arrange(individual) %>%
  mutate(individual = rep(1:length(Sample1), each = 2)) %>%
  mutate(measurement = round(measurement, digits = 0)) %>% 
  mutate(measurement = measurement - min(measurement)) %>%
   as.data.frame() 
head(antibody_data)

# Set up model. 
##########################################################################
time_per_year <- 1

print(min(antibody_data$birth))

possible_exposure_times <- c(seq(min(antibody_data$birth), 2016, by = 44), seq(2017, 2022,by=1)) 
possible_exposure_times <- c(2021, 2022) 
print(possible_exposure_times)



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
par_tab[par_tab$names == "wane_short","values"] <- 0.03

par_tab[par_tab$names == "wane_short","values"] <- 0.25


#  we are NOT fitting boost long but not wane long 
par_tab[par_tab$names %in% c("boost_long"),"fixed"] <- 1
par_tab[par_tab$names %in% c("boost_long"),"values"] <- 0
par_tab[par_tab$names %in% c("wane_long"),"fixed"] <- 1
par_tab[par_tab$names %in% c("wane_long"),"values"] <- 0

## Infection model prior -- remember you can check what the attack rate prior from this looks like with
par_tab[par_tab$names %in% c("infection_model_prior_shape1","infection_model_prior_shape2"),"values"] <- c(.37, 3.5)

par_tab[par_tab$names %in% c("boost_short"),"upper_bound"] <- max(antibody_data$measurement)
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

setwd("/Users/sambents/Desktop/NIH/serology/final_kc_log10")
getwd()

dir.create("sarsspike_adult")
dir.create("sarsrbd_adult")

dir.create("hku1_adult")
dir.create("oc43_adult")
dir.create("nl63_adult")
dir.create("229e_adult")
dir.create("flubyam_adult")
dir.create("fluah1_adult")
dir.create("fluah1_adults")
dir.create("flubvic_adult")
dir.create("fluah3_adult")
dir.create("rsv_adult")
dir.create("sarsn_adult")
dir.create("sarsnc_adult")


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
                  filename="sarsrbd_adult/phirst_sim_recovery", 
                  n_chains=4, ## Run 3 chains
                  parallel=TRUE, ## Run in parallel
                  mcmc_pars=mcmc_pars_use, ## Some MCMC control settings
                  verbose=TRUE,
                  data_type=1,
                  start_level="none")
print(res$all_diagnostics)

plot_all_outputs(res, "sarsrbd_adult","sarsrbd_adult","plot_wd")




##########################
#####
# with rho
#######################
res = NULL

par_tab_rhos <- data.frame(names="rho",values= 6,fixed=0,
                           lower_bound= 0, upper_bound= 12, lower_start= 2,
                           upper_start=7,biomarker_group = 1,
                           stratification= NA,
                           par_type=3) # lower bound is first quadrant
par_tab_rhos$values <- rnorm(1, 1)
measurement_bias <- data.frame(biomarker_id = unique(antibody_data$biomarker_id),
                               obs_type = rep(1:1, each=length(unique(antibody_data$biomarker_id))),
                               rho_index=1:(length(unique(antibody_data$biomarker_id))*1))
par_tab_with_rho <- bind_rows(par_tab, par_tab_rhos)


res <- serosolver(par_tab_with_rho, 
                  antibody_data, 
                  #  demographics=demographics,
                  antigenic_map= NULL,
                  fixed_inf_hists= NULL, ## Set fixed infection states
                  prior_func=prior_func,
                  #  possible_exposure_times= antigenic_map$inf_times,
                  possible_exposure_times= possible_exposure_times,
                  measurement_bias = measurement_bias,
                  filename="oc43_adult/phirst_sim_recovery", 
                  n_chains=4, ## Run 3 chains
                  parallel=TRUE, ## Run in parallel
                  mcmc_pars=mcmc_pars_use, ## Some MCMC control settings
                  verbose=TRUE,
                  data_type=1,
                  start_level="none")
print(res$all_diagnostics)

plot_all_outputs(res, "oc43_adult","oc43_adult","plot_wd")
