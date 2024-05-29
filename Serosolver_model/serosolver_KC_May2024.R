# May 29, 2024
library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(readr)
library(viridis)
library(ggpubr)
library(tidyverse)
#devtools::install_github("seroanalytics/serosolver", ref= "multiple_obs_types") 
library(serosolver)

setwd("C:/Users/bentssj/OneDrive - National Institutes of Health/Year_2024/serology/seattle")

# Run pediatric samples 
#___________________________________________________________________________
# Read in pediatric sero samples and demographic data.
cov = read.csv("serosamples_CoV_22Mar2024.csv") %>%
  drop_na(titer_mean)
resp = read.csv("serosamples_Resp_22Mar2024.csv")  %>%
  drop_na(titer_mean)
dem = read.csv("demos_children.csv") %>%
  mutate(Sample = ï..sample_id)


sero_samp = left_join(cov, dem, by = "Sample") %>%
  mutate(titre = log(titer_mean)) %>%
  filter(titre > 0 ) %>%
  mutate(sample_ID = Sample) %>%
  mutate(blood_date = as.Date(blood_date), year = year(blood_date)) %>%
  filter(Assay == "SARS-CoV-2 N") %>%
  filter(age_raw < 5) %>%
  filter(age_raw > .25) %>%
  mutate(individual = seq(1, length(unique(sero_samp$sample_ID)), 1) )

setwd("C:/Users/bentssj/OneDrive - National Institutes of Health/Year_2024/serology/seattle")
serosolver_wd = setwd("C:/Users/bentssj/OneDrive - National Institutes of Health/Year_2024/serology/seattle")


# Prepare data format for serosolver ___________________________________________________________
# singular observation per person.

titre_dat = sero_samp %>% 
  mutate(samples = as.numeric(year), virus = 2020, run = 1, group = 1, 
         obs_type = 1) %>%
  mutate(DOB = 2015) %>% # This model does not use information on DOB so this column can be arbitrarily filled in. 
  dplyr::select(individual, samples, virus, obs_type, titre, run, group, DOB)

ggplot(data = titre_dat) + geom_violin(aes(x=samples,y=titre,fill=obs_type,group=samples)) + facet_wrap(~obs_type)+
  theme_bw()
#__________________________________________________________________________________________________

# run serosolver model 
run_name <- "covN_under5_noMI_seattle_may29"
main_wd <- "C:/Users/bentssj/OneDrive - National Institutes of Health/Year_2024/serology/seattle"
chain_wd <- paste0(main_wd,"/chains/",run_name)
save_wd <- paste0(main_wd,"/figures/")

if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)

buckets = 1  ## Time resolution, meaning 1 year
prior_version <- 2 ## Which prior on the infection histories to use. 
n_chains <- 5 ## Number of MCMC chains to run

rerun <- TRUE
setwd(main_wd)
print(paste0("In directory: ", main_wd))
print(paste0("Saving to: ", save_wd))
set.seed(1)


## Run multiple chains in parallel
cl <- makeCluster(n_chains)
registerDoParallel(cl)

## MCMC settings, not super important but can tweak number of iterations
mcmc_pars <- c("save_block"=100,"thin"=10,"thin_hist"=50,
               "iterations"=700000,
               "adaptive_period"=50000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.05,
               "year_swap_propn"=0.8,"swap_propn"=0.5,
               "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=3, "hist_opt"=0,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)



n_obs_types <- 1 ## 1 observation type
obs_type_dists <- c(1) # # number of samples per individual 
sampled_viruses <- c(2020)*buckets

# define null antigenic map
antigenic_map <- data.frame(inf_times= seq(2020,  2022, 1),x_coord=1,y_coord=1)
strain_isolation_times <- antigenic_map$inf_times
n_times <- length(strain_isolation_times)

## Set up parameter table
par_tab <- read.csv("par_tab_base_new.csv",stringsAsFactors=FALSE)
par_tab <- par_tab[par_tab$names != "phi",]
par_tab[par_tab$names %in% c("alpha","beta"),c("values")] <- c(.37,3.5) # prior on the per capita attack rate
par_tab$fixed <- 1
par_tab[par_tab$names %in% c("obs_sd","mu_short","wane"),"fixed"] <- 0 # parameters we want to fit 
par_tab[par_tab$names %in% c("MIN_TITRE"),"values"] <- 0
par_tab[par_tab$names %in% c("MAX_TITRE"),"values"] <-16
par_tab[par_tab$names %in% c("mu"),"values"] <- 0 # not fitting long-term boost 


## Create some random measurement offsets -- these add a systematic shift to each observation
## These are unknown parameters to be estimated, assuming the offsets are drawn from ~ norm(0, 1)
par_tab_rhos <- data.frame(names="rho",values= 6,fixed=0,steps=0.1,
                           lower_bound= 2, upper_bound= 15.4, lower_start= 5,
                           upper_start=12,type=3) # lower bound is first quadrant

par_tab_rhos$values <- rnorm(length(sampled_viruses), 1)
measurement_indices <- data.frame(virus = sampled_viruses, 
                                  obs_type = rep(1:n_obs_types, each=length(sampled_viruses)),
                                  rho_index=1:(length(sampled_viruses)*n_obs_types))

par_tab <- bind_rows(par_tab, par_tab_rhos)


## Extend parameter table for each additional observation type
par_tab$obs_type <- 1
par_tab_tmp <- par_tab
if(n_obs_types > 1){
  for(i in 2:n_obs_types){
    par_tab_tmp2 <- par_tab_tmp
    antigenic_map_tmp2 <- antigenic_map_tmp
    antigenic_map_tmp2$obs_type <- i 
    par_tab_tmp2$obs_type <- i
    par_tab <- bind_rows(par_tab, par_tab_tmp2 %>% filter(!(names %in% c("alpha","beta"))))
    antigenic_map <- bind_rows(antigenic_map, antigenic_map_tmp2)
  }
}

## Randomize model parameters
par_tab <- generate_start_tab(par_tab)
head(par_tab)

# option to put prior on other variables 
#prior_func <- function(par_tab){
#  f <- function(pars){
#    pr <- dlnorm(pars["wane"],log(.32),0.34,log=TRUE) #1.8
#    pr
#  }
#}


# create posterior function 
f <- create_posterior_func(par_tab,titre_dat,antigenic_map=antigenic_map,
                           strain_isolation_times = strain_isolation_times,
                           version=prior_version,solve_likelihood=TRUE,n_alive=NULL,
                           measurement_indices_by_time=measurement_indices,
                           data_type=obs_type_dists)

## Time runs and use dopar to run multiple chains in parallel
t1 <- Sys.time()
filenames <- paste0(chain_wd, "/",run_name, "_",1:n_chains)



if(rerun){
  res <- foreach(x = filenames, .packages = c('data.table','plyr',"dplyr","tidyverse","serosolver")) %dopar% {
    #devtools::load_all(serosolver_wd)
    index <- 1
    lik <- -Inf
    inf_hist_correct <- 1
    while((!is.finite(lik) || inf_hist_correct > 0) & index < 100){
      start_tab <- generate_start_tab(par_tab)
      start_inf <- setup_infection_histories_total(titre_dat,strain_isolation_times,2,3) #2,3 
      
      inf_hist_correct <- sum(check_inf_hist(titre_dat, strain_isolation_times, start_inf))
      y <- f(start_tab$values, start_inf)
      lik <- sum(y[[1]])
      index <- index + 1
    }
    
    write.csv(start_tab, paste0(x, "_start_tab.csv"))
    write.csv(start_inf, paste0(x, "_start_inf_hist.csv"))
    write.csv(antigenic_map, paste0(x, "_antigenic_map.csv"))
    write.csv(titre_dat, paste0(x, "_titre_dat.csv"))
    
    res <- run_MCMC(start_tab, 
                    titre_dat, 
                    antigenic_map, 
                    start_inf_hist=start_inf,
                    filename=x,
                    CREATE_POSTERIOR_FUNC=create_posterior_func, 
                    CREATE_PRIOR_FUNC = NULL, # was null, CAN BE PRIOR FUNC
                    version=prior_version,
                    mcmc_pars=mcmc_pars,
                    measurement_indices= measurement_indices, ## measurement_indices, ## NULL
                    measurement_random_effects = FALSE, ## TRUE, ## FALSE
                    solve_likelihood=TRUE,
                    data_type=obs_type_dists)
  }
}


run_time_fast <- Sys.time() - t1

run_time_fast


## load in chains when done running. 
chains <- load_mcmc_chains(chain_wd,convert_mcmc=TRUE,burnin = mcmc_pars["adaptive_period"],unfixed = TRUE)

list_chains = chains$theta_list_chains
list_chains1 <- lapply(list_chains, function(x) x[,c("mu_short", "wane", "total_infections", "rho")])

#check diagnositcs 
print(gelman.diag(as.mcmc.list(list_chains1)))
effectiveSize(as.mcmc.list(list_chains1))

# view parameter estimates 
print(summary(as.mcmc.list(list_chains1)))






#############################################################################
# Run adult samples
serology_samples = read.csv("serology_samples_all.csv") %>%
  filter(Assay== "SARS-CoV-2 N")%>%
  filter(logtiter > 0)  %>%
  filter(age >17.999) %>% 
  group_by(sample_ID) %>% mutate(count = count(sample_ID)) %>%
  filter(count$freq == 2)

sample_ID = unique(serology_samples$sample_ID)
individual = seq(1, length(sample_ID), 1)
in_dat = data.frame(sample_ID, individual)

titre_dat = left_join(serology_samples, in_dat, by = "sample_ID") %>%
  mutate(samples = as.numeric(year), titre = logtiter, virus = 2020, run =1, 
         group =1, DOB =2015, obs_type =1) %>%
  dplyr::select(individual, samples, virus, obs_type, titre, run, group, DOB)

titre_dat = data.frame(titre_dat[c(2:9)])
head(titre_dat)

ggplot(data = titre_dat) + geom_violin(aes(x=samples,y=titre,fill=obs_type,group=samples)) + facet_wrap(~obs_type)+
  theme_bw()

#_________________________________

# run serosolver model 
run_name <- "cov4_adult_seattle_may29"
main_wd <- "C:/Users/bentssj/OneDrive - National Institutes of Health/South Africa"
chain_wd <- paste0(main_wd,"/chains/",run_name)
save_wd <- paste0(main_wd,"/figures/")

if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)

buckets = 1  ## Time resolution
prior_version <- 2 ## Which prior on the infection histories to use? Prior version 2 is generally preferred
n_chains <- 5 ## Number of MCMC chains to run

rerun <- TRUE
setwd(main_wd)
print(paste0("In directory: ", main_wd))
print(paste0("Saving to: ", save_wd))
set.seed(1)

## Run multiple chains in parallel
cl <- makeCluster(n_chains)
registerDoParallel(cl)

## MCMC settings, not super important but can tweak number of iterations
mcmc_pars <- c("save_block"=100,"thin"=10,"thin_hist"=50,
               "iterations"=600000,
               "adaptive_period"=50000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.05,
               "year_swap_propn"=0.8,"swap_propn"=0.5,
               "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=3, "hist_opt"=0,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)



n_obs_types <- 1 ## 1 observation type
obs_type_dists <- c(2) # # number of samples per individual 
sampled_viruses <- c(2020)*buckets

# set up null antigenic map
antigenic_map <- data.frame(inf_times= seq(2020,  2022, 1),x_coord=1,y_coord=1)
strain_isolation_times <- antigenic_map$inf_times
n_times <- length(strain_isolation_times)


## Set up parameter table
par_tab <- read.csv("par_tab_base_new.csv",stringsAsFactors=FALSE)
par_tab <- par_tab[par_tab$names != "phi",]
par_tab[par_tab$names %in% c("alpha","beta"),c("values")] <- c(.37,3.5) # prior on the per capita attack rate
par_tab$fixed <- 1
par_tab[par_tab$names %in% c("obs_sd","mu_short","wane"),"fixed"] <- 0 # parameters we want to fit 
par_tab[par_tab$names %in% c("MIN_TITRE"),"values"] <- 0
par_tab[par_tab$names %in% c("MAX_TITRE"),"values"] <-16
par_tab[par_tab$names %in% c("mu"),"values"] <- 0 # not fitting long-term boost 


## Create some random measurement offsets -- these add a systematic shift to each observation
## These are unknown parameters to be estimated, assuming the offsets are drawn from ~ norm(0, 1)
par_tab_rhos <- data.frame(names="rho",values= 8,fixed=0,steps=0.1,
                           lower_bound= 2, upper_bound= 15.406, lower_start= 4,
                           upper_start=12,type=3) # lower bound is first quadrant

par_tab_rhos$values <- rnorm(length(sampled_viruses), 1)
measurement_indices <- data.frame(virus = sampled_viruses, 
                                  obs_type = rep(1:n_obs_types, each=length(sampled_viruses)),
                                  rho_index=1:(length(sampled_viruses)*n_obs_types))

par_tab <- bind_rows(par_tab, par_tab_rhos)

## Extend parameter table for each aditional observation type
par_tab$obs_type <- 1
par_tab_tmp <- par_tab
if(n_obs_types > 1){
  for(i in 2:n_obs_types){
    par_tab_tmp2 <- par_tab_tmp
    antigenic_map_tmp2 <- antigenic_map_tmp
    antigenic_map_tmp2$obs_type <- i 
    par_tab_tmp2$obs_type <- i
    par_tab <- bind_rows(par_tab, par_tab_tmp2 %>% filter(!(names %in% c("alpha","beta"))))
    antigenic_map <- bind_rows(antigenic_map, antigenic_map_tmp2)
  }
}

## Randomize model parameters
par_tab <- generate_start_tab(par_tab)
head(par_tab)


# create posterior function 
f <- create_posterior_func(par_tab,titre_dat,antigenic_map=antigenic_map,
                           strain_isolation_times = strain_isolation_times,
                           version=prior_version,solve_likelihood=TRUE,n_alive=NULL,
                           measurement_indices_by_time=measurement_indices,
                           data_type=obs_type_dists)

## Time runs and use dopar to run multiple chains in parallel
t1 <- Sys.time()
filenames <- paste0(chain_wd, "/",run_name, "_",1:n_chains)



if(rerun){
  res <- foreach(x = filenames, .packages = c('data.table','plyr',"dplyr","tidyverse","serosolver")) %dopar% {
    #devtools::load_all(serosolver_wd)
    index <- 1
    lik <- -Inf
    inf_hist_correct <- 1
    while((!is.finite(lik) || inf_hist_correct > 0) & index < 100){
      start_tab <- generate_start_tab(par_tab)
      start_inf <- setup_infection_histories_total(titre_dat,strain_isolation_times,2,3) #2,3 
      
      inf_hist_correct <- sum(check_inf_hist(titre_dat, strain_isolation_times, start_inf))
      y <- f(start_tab$values, start_inf)
      lik <- sum(y[[1]])
      index <- index + 1
    }
    
    write.csv(start_tab, paste0(x, "_start_tab.csv"))
    write.csv(start_inf, paste0(x, "_start_inf_hist.csv"))
    write.csv(antigenic_map, paste0(x, "_antigenic_map.csv"))
    write.csv(titre_dat, paste0(x, "_titre_dat.csv"))
    
    res <- run_MCMC(start_tab, 
                    titre_dat, 
                    antigenic_map, 
                    start_inf_hist=start_inf,
                    filename=x,
                    CREATE_POSTERIOR_FUNC=create_posterior_func, 
                    CREATE_PRIOR_FUNC = NULL, # was null, CAN BE PRIOR FUNC
                    version=prior_version,
                    mcmc_pars=mcmc_pars,
                    measurement_indices= measurement_indices, ## measurement_indices, ## NULL
                    measurement_random_effects = FALSE, ## TRUE, ## FALSE
                    solve_likelihood=TRUE,
                    data_type=obs_type_dists)
  }
}


run_time_fast <- Sys.time() - t1
run_time_fast


## load in chains 
chains <- load_mcmc_chains(chain_wd,convert_mcmc=TRUE,burnin = mcmc_pars["adaptive_period"],unfixed = TRUE)
list_chains = chains$theta_list_chains
list_chains1 <- lapply(list_chains, function(x) x[,c("mu_short", "wane", "total_infections", "rho")])

#check diagnositcs 
print(gelman.diag(as.mcmc.list(list_chains1)))
effectiveSize(as.mcmc.list(list_chains1))
print(summary(as.mcmc.list(list_chains1)))




