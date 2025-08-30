library(dplyr)
library(lubridate)
library(MMWRweek)
library(ggplot2)
library(odin)
library(readr)
library(tidyr)
library(coda)
library(pracma)
library(stringr)
library(zoo)
library(MetBrewer)
library(cowplot)
require("pals")
require("lubridate")
require("fields")

# Set working directory 
setwd("/Users/sambents/Desktop/NIH/serology/unprocessed")

# Load in odin file to initialize model 
x <- odin::odin("toy_odin.R") # save odin as function

init_conds_from_file <- 0  # Choose whether to read in some existing ICs- we don't read any ICs in
save_init_conds <- 1 # Choose whether to save final model state as ICs for next time
max_t <- 2000 # max time point (in months)

# Adapted mixing matrix from Mossong et al expanded for 0-5 yr converted to monthly
mixing <- as.matrix(read.csv("mixing_75.csv", header = TRUE), header = TRUE)*365/12 


# Model function
##################################################################
run_imm_model_fit <- function(b0 = .8,
                              b1 =   .5,
                              phi =  0.98456,
                              prop_detected_1 = 0.001381, # proportion detected, hospitalization rate for age group 1
                              prop_detected_2 = 0.000398, # hos rate 2 
                              prop_detected_3 = 0.000541, # hos rate 3
                              prop_detected_4 = 0.001369,# hos rate 4
                              prop_detected_5 = 0.003164,# hos rate 5
                              imm_detected_1 = .0417, # proportion detected, hospitalization rate for age group 1
                              imm_detected_2 = .0417, # hos rate 2 
                              imm_detected_3 = .0417, # hos rate 3
                              imm_detected_4 = .0417,# hos rate 4
                              imm_detected_5 = .0417,
                              mixing = mixing,
                              max_t = 2100,
                              init_conds_from_file = 0, 
                              timestart = 0, 
                              timeend = 0, 
                              timestart2 = 0,
                              timeend2 = 0,
                              timestart3 = 0,
                              timeend3 = 0,
                              timestart4 = 0,
                              timeend4 = 0,
                              betared = 0, 
                              betared2 = 0,
                              betared3 = 0,
                              betared4 = 0,
                              maternalimmunity = FALSE,
                              popInput = 737000,# population size 
                              age_threshhold = 5){ # was 10 in the Hogan model, should be adpated if dealing with older ages 
  

  omega <- .5   # reduced infectiousness from older age groups
  delta <- 15.21 # latency rate 
  gamma <- 13.41  # infectious rate 
  dt <- .033  
  T0 <- 0
  
  
  # maternal protection parameters - represent scaled susceptibility to infection in first few months of life
  sigma_1 <- 0.08
  sigma_2 <- 0.45
  sigma_3 <- 0.45
  sigma_4 <- 1
  sigma_5 <- 1
  sigma_6 <- 1
  if(maternalimmunity==FALSE){
    sigma_1 <- 1
    sigma_2 <- 1
    sigma_3 <- 1
    sigma_4 <- 1
    sigma_5 <- 1
    sigma_6 <- 1
  }
  
  ##########################################################
  # ageing-related parameters
  seattle_pop <- popInput 
  age_vect_years <- c(seq(0,5,1/12), seq(10,75,5)) # age structure defined here
  nAges <- length(age_vect_years) #number of age cohorts
  final_age <- 80
  age_vect_months <- age_vect_years*12
  size_cohorts_months <- c(diff(age_vect_months), final_age*12 - age_vect_months[length(age_vect_months)])
  trans_rate <- 1/size_cohorts_months
  rel_sizes <- size_cohorts_months/sum(size_cohorts_months)
  
  ##########################################################
  #  reduced infectiousness
  omega_vect <- as.vector(rep(1, nAges))
  omega_vect[!(age_vect_years < age_threshhold)] <- omega
  
  
  #Pproportion hospitalizations detected, this needs to be changed if age structure changes 
  prop_detected_vect <- as.vector(c(rep(prop_detected_1, 60),#0-4 yo
                                    rep(prop_detected_2, 3), #5-10 yo
                                    rep(prop_detected_3, 6), #10-20 yo
                                    rep(prop_detected_4, 3), #20-64 yo
                                    rep(prop_detected_5, (nAges-72))))  #65+
  
  # immunity age structure, this needs to be changed if age structure changes 
  imm_detected_vect  <- as.vector(c(rep(imm_detected_1, 60),#0-4 yo
                                    rep(imm_detected_2, 3), #5-10 yo
                                    rep(imm_detected_3, 6), #10-20 yo
                                    rep(imm_detected_4, 3),  #20-64 yo
                                    rep(imm_detected_5, (nAges-72)))) #65+

  # Reduced susceptibility in youngest cohorts to reflect natural maternally-derived immunity
  sigma_vect <- matrix(1, 1, nAges)
  sigma_vect[1] <- sigma_1
  sigma_vect[2] <- sigma_2
  sigma_vect[3] <- sigma_3
  sigma_vect[4] <- sigma_4
  sigma_vect[5] <- sigma_5
  sigma_vect[6] <- sigma_6
  sigma_vect <- as.vector(c(sigma_1, sigma_2, sigma_3, sigma_4, sigma_5, sigma_6, rep(1, (nAges-6))))
  
  ##################################################################
  # initial conditions: choose whether to reset here, or read in from csv file
  if (init_conds_from_file == 1) {
    ic <- readRDS("init_conds.rds")
    S0 <- ic$S0
    E0 <- ic$E0
    I0 <- ic$I0
    R0 <- ic$R0
    Incidence0 <- rel_sizes * 0
    DetIncidence0 <- rel_sizes * 0
  } else {
    I0 <- rel_sizes * seattle_pop * 0.01
    S0 <- rel_sizes * seattle_pop * 0.99
    E0 <- rel_sizes * seattle_pop * 0
    R0 <- rel_sizes * seattle_pop * 0
    DetIncidence0 <- rel_sizes * seattle_pop * 0
    Incidence0 <- R0
  }
  
  ##########################################################
  
  pars <- list(
    b0 = b0,
    b1 = b1,
    phi = phi,
    delta = delta,
    gamma = gamma,
    prop_detected_vect = prop_detected_vect,
    imm_detected_vect = imm_detected_vect,
    sigma_vect = sigma_vect,
    omega_vect = omega_vect,
    mixing = mixing,
    timestart = timestart,
    timeend = timeend, 
    timestart2 = timestart2,
    timeend2 = timeend2, 
    timestart3 = timestart3,
    timeend3 = timeend3, 
    timestart4 = timestart4,
    timeend4 = timeend4, 
    betared = betared,
    betared2 = betared2,
    betared3 = betared3,
    betared4 = betared4,
    S0 = S0,
    E0 = E0,
    I0 = I0,
    R0 = R0,
    Incidence0 = Incidence0,
    DetIncidence0 = DetIncidence0
  )
  
  ##########################################################
  # run model with cohort ageing
  
  mod <- x(user = pars)
  
  pop_out <- NULL
  while (T0 <= max_t){
    
    # solve the odes first
    t <- seq(from = T0, to = T0 + 1, by = dt)
    m <- mod$run(t)
    pop <- mod$transform_variables(m)
    
    if (T0 == 0){
      pop_out <- pop
    } else {
      pop_out$time <- c(pop_out$time, pop$time[5])
      pop_out$S <- rbind(pop_out$S, pop$S[5,])
      pop_out$E <- rbind(pop_out$E, pop$E[5,])
      pop_out$I <- rbind(pop_out$I, pop$I[5,])
      pop_out$R <- rbind(pop_out$R, pop$R[5,])
      pop_out$Incidence <- rbind(pop_out$Incidence, pop$Incidence[5,])
      pop_out$DetIncidence <- rbind(pop_out$DetIncidence, pop$DetIncidence[5,])
    }
    
    # cohort ageing
    
    # extract the final state from pop
    S <- as.vector(t(data.table::last(pop$S)))
    E <- as.vector(t(data.table::last(pop$E)))
    I <- as.vector(t(data.table::last(pop$I)))
    R <- as.vector(t(data.table::last(pop$R)))
    Incidence <- as.vector(t(data.table::last(pop$Incidence)))
    DetIncidence <- as.vector(t(data.table::last(pop$DetIncidence)))
    
    # initialise the new initial condition vectors
    I0 <- rel_sizes*0
    S0 <- rel_sizes*0
    E0 <- rel_sizes*0
    R0 <- rel_sizes*0
    Incidence <- rel_sizes*0
    DetIncidence <- rel_sizes*0
    
    # then fill them in
    S0[1] <- seattle_pop * rel_sizes[1]
    
    for(i in c(2:nAges)){
      S0[i] = S[(i - 1)] * trans_rate[(i-1)] + S[i] - S[i] * trans_rate[i]
      E0[i] = E[(i - 1)] * trans_rate[(i-1)] + E[i] - E[i] * trans_rate[i]
      I0[i] = I[(i - 1)] * trans_rate[(i-1)] + I[i] - I[i] * trans_rate[i]
      R0[i] = R[(i - 1)] * trans_rate[(i-1)] + R[i] - R[i] * trans_rate[i]
      Incidence0[i] = Incidence[(i - 1)] * trans_rate[(i-1)] + Incidence[i] - Incidence[i] * trans_rate[i]
      DetIncidence0[i] = DetIncidence[(i - 1)] * trans_rate[(i-1)] + DetIncidence[i] - DetIncidence[i] * trans_rate[i]
    }
    
    pars <- list(
      b0 = b0,
      b1 = b1,
      phi = phi,
      delta = delta,
      gamma = gamma,
      prop_detected_vect = prop_detected_vect,
      imm_detected_vect = imm_detected_vect,
      sigma_vect = sigma_vect,
      omega_vect = omega_vect,
      mixing = mixing,
      timestart = timestart,
      timeend = timeend, 
      timestart2 = timestart2,
      timeend2 = timeend2, 
      timestart3 = timestart3,
      timeend3 = timeend3, 
      timestart4 = timestart4,
      timeend4 = timeend4, 
      betared = betared,
      betared2 = betared2,
      betared3 = betared3,
      betared4 = betared4,
      S0 = S0,
      E0 = E0,
      I0 = I0,
      R0 = R0,
      Incidence0 = Incidence0,
      DetIncidence0 = DetIncidence0
    )
    # code in NPIs
    if(T0 >= timestart & T0 <= timeend){
      pars <- list(
        b0 = b0*betared,
        b1 = b1,
        phi = phi,
        delta = delta,
        gamma = gamma,
        prop_detected_vect = prop_detected_vect,
        imm_detected_vect = imm_detected_vect,
        sigma_vect = sigma_vect,
        omega_vect = omega_vect,
        mixing = mixing,
        timestart = timestart,
        timeend = timeend, 
        timestart2 = timestart2,
        timeend2 = timeend2, 
        timestart3 = timestart3,
        timeend3 = timeend3, 
        timestart4 = timestart4,
        timeend4 = timeend4, 
        betared = betared,
        betared2 = betared2,
        betared3 = betared3,
        betared4 = betared4,
        S0 = S0,
        E0 = E0,
        I0 = I0,
        R0 = R0,
        Incidence0 = Incidence0,
        DetIncidence0 = DetIncidence0
      )
    }
    if(T0 >= timestart2 & T0 <= timeend2) {
      pars <- list(
        b0 = b0*betared2,
        b1 = b1,
        phi = phi,
        delta = delta,
        gamma = gamma,
        prop_detected_vect = prop_detected_vect,
        imm_detected_vect = imm_detected_vect,
        sigma_vect = sigma_vect,
        omega_vect = omega_vect,
        mixing = mixing,
        timestart = timestart,
        timeend = timeend, 
        timestart2 = timestart2,
        timeend2 = timeend2, 
        timestart3 = timestart3,
        timeend3 = timeend3, 
        timestart4 = timestart4,
        timeend4 = timeend4, 
        betared = betared,
        betared2 = betared2,
        betared3 = betared3,
        betared4 = betared4,
        S0 = S0,
        E0 = E0,
        I0 = I0,
        R0 = R0,
        Incidence0 = Incidence0,
        DetIncidence0 = DetIncidence0
      )
    }
    if(T0 >= timestart3 & T0 <= timeend3) {
      pars <- list(
        b0 = b0*betared3,
        b1 = b1,
        phi = phi,
        delta = delta,
        gamma = gamma,
        prop_detected_vect = prop_detected_vect,
        imm_detected_vect = imm_detected_vect,
        sigma_vect = sigma_vect,
        omega_vect = omega_vect,
        mixing = mixing,
        timestart = timestart,
        timeend = timeend, 
        timestart2 = timestart2,
        timeend2 = timeend2, 
        timestart3 = timestart3,
        timeend3 = timeend3, 
        timestart4 = timestart4,
        timeend4 = timeend4, 
        betared = betared,
        betared2 = betared2,
        betared3 = betared3,
        betared4 = betared4,
        S0 = S0,
        E0 = E0,
        I0 = I0,
        R0 = R0,
        Incidence0 = Incidence0,
        DetIncidence0 = DetIncidence0
      )
    }
    if(T0 >= timestart4 & T0 <= timeend4) {
      pars <- list(
        b0 = b0*betared4,
        b1 = b1,
        phi = phi,
        delta = delta,
        gamma = gamma,
        prop_detected_vect = prop_detected_vect,
        imm_detected_vect = imm_detected_vect,
        sigma_vect = sigma_vect,
        omega_vect = omega_vect,
        mixing = mixing,
        timestart = timestart,
        timeend = timeend, 
        timestart2 = timestart2,
        timeend2 = timeend2, 
        timestart3 = timestart3,
        timeend3 = timeend3, 
        timestart4 = timestart4,
        timeend4 = timeend4, 
        betared = betared,
        betared2 = betared2,
        betared3 = betared3,
        betared4 = betared4,
        S0 = S0,
        E0 = E0,
        I0 = I0,
        R0 = R0,
        Incidence0 = Incidence0,
        DetIncidence0 = DetIncidence0
      )
    }
    mod <- x(user = pars) # still in while loop
    T0 <- T0 + 1
    
  }
  pop_out <- pop_out[2:8]
  return(pop_out)
}



####################################################################
# Run mode output 
# Fitted uniform immunity scenario
output = run_imm_model_fit(b0 = .53, 
                           b1 = .280, 
                           phi = 4.59, 
                           prop_detected_1 = .124, # hospitalization rate 1 for age group 1
                           prop_detected_2 = .089, # hr2 
                           prop_detected_3 = .062, # hr3 
                           prop_detected_4 = .062, # hr4
                           prop_detected_5 =  .104, # hr4 
                           imm_detected_1 = .0208,  # immunity for age group 1
                           imm_detected_2 = .0208,  
                           imm_detected_3 = .0208,  
                           imm_detected_4 = .0208, 
                           imm_detected_5 = .0208,
                           max_t = max_t,
                           mixing = mixing,
                           timestart = 1885, # NPI 1 period start time, note time is given in months since 1880 
                           timeend = 1886,   #NPI 1 period end time
                           timestart2 = 1888, #NPI period 2 start time
                           timeend2 = 1900,  #NPI period 2 end time    
                           timestart3 = 1901,               
                           timeend3 = 1913,                 
                           timestart4 = 1965,              
                           timeend4 = 1966,                 
                           betared = 1,  # percent beta reduction for NPI period 1, leave at 1 for no NPI         
                           betared2 = .57,                 
                           betared3 = .73,             
                           betared4 = 1,                  
                           init_conds_from_file = 0, maternalimmunity = TRUE, popInput = 2280540)

all_inc <- output$DetIncidence
all_inc <- all_inc[1880:2000,]
all_inc_sum <- as.data.frame(all_inc)%>%
  summarize(yr0_5 = V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9+ V10 + V11 + V12+
              V13 + V14 + V15 + V16 + V17 + V18+ V19 + V20 + V21+ V22 + V23 + V24 +
              V25 + V26 + V27 + V28 + V29 + V30 + V31 + V32 + V33 + V34 + V35 + V36 + V37 + V38 + V39 + V40 + 
              V41 + V42 + V43 + V44 + V45 + V46 + V47 + V48 + V49 + V50 + 
              V51 + V52 + V53 + V54 + V55 + V56 + V57 + V58 + V59 + V60 ,
            yr5_10 = V61 ,
            yr10_20 = V62 + V63,  
            yr20_65 =  V64 + V65 + V66 + V67 + V68 + V69 + V70 + V71 + V72, 
            yr_65 = V73+ V74+ V75) %>%
  mutate(total = (yr0_5 +  yr5_10 + yr10_20 + yr20_65 +  yr_65)) %>% 
  mutate(t = seq(1880, 2000, 1)) %>% 
  mutate(time = seq(2017, 2027, 1/12 )) %>%
  mutate(pct = total/2280540*100)


# Fitted tiered immunity scenario 
output = run_imm_model_fit(b0 = .51,  
                           b1 = .32,  
                           phi = 4.61, 
                           prop_detected_1 = .130,   # hospitalization rate 1 for age group 1
                           prop_detected_2 = .195,   # hr2 
                           prop_detected_3 = .080,   # hr3 
                           prop_detected_4 = .081,   # hr4
                           prop_detected_5 = .08,    # hr5
                           imm_detected_1 = .0208*2.05,  # immunity for age group 1
                           imm_detected_2 = .0208,  
                           imm_detected_3 = .0208,  
                           imm_detected_4 = .0208,  
                           imm_detected_5 = .0208,  
                           max_t = max_t,
                           mixing = mixing,
                           timestart = 1885, # NPI 1 period start time, note time is given in months since 1880 
                           timeend = 1886, #NPI 1 period end time
                           timestart2 = 1888,#NPI period 2 start time
                           timeend2 = 1900, #NPI period 2 end time    
                           timestart3 = 1901,               
                           timeend3 = 1913,                 
                           timestart4 = 1965,              
                           timeend4 = 1966,                 
                           betared = 1,  # percent beta reduction for NPI period 1, leave at 1 for no NPI         
                           betared2 =  .50,  #.75 #.79                   
                           betared3 = .605,     # .83             
                           betared4 = 1,                  
                           init_conds_from_file = 0, maternalimmunity = TRUE, popInput = 2280540)

all_inc <- output$DetIncidence
all_inc <- all_inc[1880:2000,]

all_inc_sum <- as.data.frame(all_inc)%>%
  summarize(yr0_5 = V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9+ V10 + V11 + V12+
              V13 + V14 + V15 + V16 + V17 + V18+ V19 + V20 + V21+ V22 + V23 + V24 +
              V25 + V26 + V27 + V28 + V29 + V30 + V31 + V32 + V33 + V34 + V35 + V36 + V37 + V38 + V39 + V40 + 
              V41 + V42 + V43 + V44 + V45 + V46 + V47 + V48 + V49 + V50 + 
              V51 + V52 + V53 + V54 + V55 + V56 + V57 + V58 + V59 + V60 ,
            yr5_10 = V61 ,
            yr10_20 = V62 + V63,  
            yr20_65 =  V64 + V65 + V66 + V67 + V68 + V69 + V70 + V71 + V72, 
            yr_65 = V73+ V74+ V75) %>%
  mutate(total = (yr0_5 +  yr5_10 + yr10_20 + yr20_65 +  yr_65)) %>% 
  mutate(t = seq(1880, 2000, 1)) %>% 
  mutate(time = seq(2017, 2027, 1/12 )) %>%
  mutate(pct = total/2280540*100)


#####################################################
# Pull in observed data 
########################
kingco = read.csv("kingco_flu_weekly.csv") # King County weekly flu data 

kc_monthly = kingco %>%
  filter(demo_age =="All") %>%
  mutate(first_day_of_week = ymd(paste(year_mmwr, "-01-01")) + weeks(week_mmwr - 1)) %>%
  mutate(month = month(first_day_of_week)) %>%
  mutate(hosp_combined = hosp_inf + ed_inf) %>%
  mutate(year = year(first_day_of_week)) %>%
  filter(first_day_of_week < "2021-01-01") %>%
  mutate(hosp = hosp_combined*1) %>%
  group_by(year, month)%>%
  dplyr::mutate(monthly_hosp = sum(hosp)) %>%
  mutate(time = month/12 + year) %>%
  distinct(time, monthly_hosp)

kc_monthly_2021 = kingco %>%
  filter(demo_age =="All")%>%
  mutate(first_day_of_week = ymd(paste(year_mmwr, "-01-01")) + weeks(week_mmwr - 1)) %>%
  mutate(month = month(first_day_of_week)) %>%
  mutate(hosp_combined = hosp_inf + ed_inf) %>%
  mutate(month = month(first_day_of_week)) %>%
  mutate(year = year(first_day_of_week)) %>%
  filter(first_day_of_week> "2021-01-01") %>%
  mutate(hosp = hosp_combined) %>%
  group_by(year, month)%>%
  dplyr::mutate(monthly_hosp = sum(hosp)) %>%
  mutate(time = month/12 + year) %>%
  distinct(time, monthly_hosp)

full_time_series = rbind(kc_monthly, kc_monthly_2021)
head(full_time_series)

# Control for number of reporting hospitals 
multiplier = c(rep(21/9, 12), rep(21/10, 3), rep(21/15, 4), rep(21/16, 2),
               rep(21/17, 2), rep(21/18, 2), rep(21/19, 2), rep(21/20, 40),
               rep(1, 8)) # removed one week to end was rep(1, 9))
time = print(full_time_series$time)
dates = data.frame(time, multiplier) 

full_time_series_fin = left_join(full_time_series, dates, by = "time") %>%
  mutate(monthly_hosp_fin = monthly_hosp*multiplier)%>%
  ungroup() %>%
  mutate(t = seq(1880, 1954,1)) # removed one was 1955

# Age structured hospitalizations 
kc_monthly_age = kingco %>%
  filter(demo_age !="Unknown")%>%
  filter(demo_age !="All")%>%
  mutate(first_day_of_week = ymd(paste(year_mmwr, "-01-01")) + weeks(week_mmwr - 1)) %>%
  mutate(month = month(first_day_of_week)) %>%
  mutate(hosp_combined = hosp_inf + ed_inf) %>%
  mutate(month = month(first_day_of_week)) %>%
  mutate(year = year(first_day_of_week)) %>%
  filter(first_day_of_week < "2021-01-01") %>%
  mutate(hosp = hosp_combined*1) %>%
  mutate(age_band = ifelse(demo_age == "<1 year", "0-4",
                           ifelse(demo_age == "1 years", "0-4",
                                  ifelse(demo_age == "2-4 years", "0-4",
                                         ifelse(demo_age == "5-9 years","5-19",
                                                ifelse(demo_age == "10-19 years", "5-19",
                                                       ifelse(demo_age == "20-39 years", "20-64",
                                                              ifelse(demo_age == "40-64 years", "20-64", "65+"
                                                              )))))))) %>%
  group_by(year, month, age_band)%>%
  dplyr::mutate(monthly_hosp = sum(hosp)) %>%
  mutate(time = month/12 + year) %>%
  distinct(time, monthly_hosp, age_band)


kc_monthly_2021_age = kingco %>%
  filter(demo_age !="Unknown")%>%
  filter(demo_age !="All")%>%
  mutate(first_day_of_week = ymd(paste(year_mmwr, "-01-01")) + weeks(week_mmwr - 1)) %>%
  mutate(month = month(first_day_of_week)) %>%
  mutate(hosp_combined = hosp_inf + ed_inf) %>%
  mutate(year = year(first_day_of_week)) %>%
  filter(first_day_of_week > "2021-01-01") %>%
  mutate(age_band = ifelse(demo_age == "<1 year", "0-4",
                           ifelse(demo_age == "1 years", "0-4",
                                  ifelse(demo_age == "2-4 years", "0-4",
                                         ifelse(demo_age == "5-9 years","5-19",
                                                ifelse(demo_age == "10-19 years", "5-19",
                                                       ifelse(demo_age == "20-39 years", "20-64",
                                                              ifelse(demo_age == "40-64 years", "20-64", "65+"
                                                              )))))))) %>%
  mutate(hosp = hosp_combined) %>%
  group_by(year, month, age_band )%>%
  dplyr::mutate(monthly_hosp = sum(hosp)) %>%
  mutate(time = month/12 + year) %>%
  distinct(time, monthly_hosp, age_band)


full_time_series_age = rbind(kc_monthly_age, kc_monthly_2021_age)
head(full_time_series_age)

time = print(unique(full_time_series_age$time))
multiplier = c(rep(21/9, 12), rep(21/10, 3), rep(21/15, 4), rep(21/16, 2),
               rep(21/17, 2), rep(21/18, 2), rep(21/19, 2), rep(21/20, 40),
               rep(1, 8)) # was 9

mult = data.frame(time, multiplier)
full_time_series_age_mult = left_join(full_time_series_age, mult, by = "time")


rates = c(127000, 354700, 1344600, 252800 ) # Age group sizes in King County, Seattle 
age_band = c("0-4", "5-19", "20-64", "65+")
rayband = data.frame(rates, age_band)

age_plot = left_join(full_time_series_age_mult, rayband, by  = c("age_band")) %>%
  mutate(hosp_rate_per_1000 = ((monthly_hosp*multiplier)/rates)*1000) %>%
  mutate(hosp_rate_per_1000 = round(hosp_rate_per_1000, digits = 2))
head(age_plot)

# Save WDOH hospitalization rates 
wdoh = age_plot %>% dplyr::select(year, month, age_band, hosp_rate_per_1000)
write.csv(wdoh, "RHINO_hospitalization_rates.csv")

# Plot age structure hospitalizations 
re_pal = met.brewer(name = "Homer1", n = 4)
a = ggplot(data = age_plot %>% filter(time > 2017.6), aes(x = time, y = ((monthly_hosp*multiplier)/rates)*1000, col  = factor(age_band, levels = c("0-4", "5-19", "20-64", "65+")))) +
  geom_line(lwd = 1.3, aes(y = rollmean(((monthly_hosp*multiplier)/rates)*1000, 2, na.pad = TRUE))) +
  ylab("Monthly healthcare encounters per 1000")+
  theme_light()+
  xlab("Date") +
  scale_color_manual(values = re_pal) +
  labs(color = "Age group") +
  theme(legend.position="bottom") +
  theme(axis.title.y = element_text(size = 8)) +
  theme(axis.title.x = element_blank() ) +
  theme(axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 9, color = "black"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))

a

# Pre-pandemic age structure 
before = full_time_series_age %>%
  filter(time < 2020.5) %>%
  group_by(age_band) %>%
  dplyr::mutate(age_sum = sum(monthly_hosp)) %>%
  distinct(age_sum, age_band) %>%
  ungroup() %>%
  mutate(tot = sum(age_sum)) %>%
  mutate(prop = age_sum/tot) %>%
  mutate(period = "Pre-pandemic")

# Post-pandemic age structure
after = full_time_series_age %>%
  filter(time > 2022) %>%
  group_by(age_band) %>%
  dplyr::mutate(age_sum = sum(monthly_hosp)) %>%
  distinct(age_sum, age_band) %>%
  ungroup() %>%
  mutate(tot = sum(age_sum)) %>%
  mutate(prop = age_sum/tot) %>%
  mutate(period = "Post-pandemic")

age = rbind(before, after)

# Plot age structure 
b = ggplot(data = age, aes(x = factor(age_band, levels = c("0-4", "5-19", "20-64", "65+")), y = prop, fill = factor(period, levels = c("Pre-pandemic", "Post-pandemic"))))+
  geom_col(position = "dodge")+
  scale_fill_manual(values = re_pal) +
  labs(fill = "Time period") +
  ylab("Proportion of healthcare encounters") +
  xlab('Age group') + 
  theme_light()+
  theme(legend.position="bottom") +
  theme(axis.title.y = element_text(size = 8)) +
  theme(axis.text = element_text(size = 7, color = "black"),
      axis.title = element_text(size = 7, color = "black"),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
      plot.title.position = "plot",
      plot.subtitle = element_text(colour = "black", size = 12.5),
      legend.position = "bottom",
      legend.key.width = unit(0.5, "cm"),
      legend.text = element_text(size = 8, color = "black"),
      legend.title = element_text(size =9, color = "black"),
      strip.text = element_text(colour = "black", size = 15, hjust = 0),
      strip.background = element_rect(colour="white", fill="white"),
      panel.border = element_rect(colour = "black", fill=NA))
b

real = plot_grid(a,b, rel_widths = c(1/2, 1/4), labels = c("a", "b"))
real



###################################################################
# Plot model results 
plot_un = all_inc_sum # after running uniform scenario
plot_ti = all_inc_sum # after running tiered scenario

c = ggplot(data = plot_un) +
  annotate("rect", xmin = c(2020.33, 2021.33), xmax = c(2021.33, 2022.33), 
           ymin = 0, ymax = Inf, alpha = .5, fill = c("gray65", "gray80"))+
  geom_line(aes(x = time, y = rollmean(total, k = 2, na.pad = TRUE)), lwd = 1.5, col = "skyblue",  alpha = .9)+
  geom_area(aes(x = time, rollmean(total, k = 2, na.pad = TRUE), fill = "lightgray"), alpha = .5) +
  geom_line(data = plot_ti, aes(x = time, y = rollmean(total, k = 2, na.pad = TRUE)), lwd = 1.5, col = "darkred",  alpha = .9)+
  geom_area(data = plot_ti, aes(x = time,rollmean(total, k = 2, na.pad = TRUE), fill = "darkred"), alpha = .5) +
  scale_fill_manual(values = c(  "darkred", "skyblue"), labels = c("Tiered immunity","Uniform immunity"))+
  xlim(c(2017.75, 2023.5))+
  geom_line(data = full_time_series_fin, aes(x  = time, y = monthly_hosp_fin, col = "black"), lwd = 1.5)+
  scale_color_manual(values = c("black"), labels = c("Observed")) +
  labs(color=NULL, fill = NULL) +
  theme_bw() +
  theme(axis.title.x = element_blank() ) +
  ylab("Monthly healthcare encounters") +
  xlab("Time") +
  theme(legend.position="bottom") +
  theme(axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size =9, color = "black"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  theme(axis.title.y = element_text(size = 7)) 
c

# Predicted age structure - uniform immunity 
pre_uniform = plot_un %>%
  filter(time < 2020.5) %>%
  dplyr::summarize(across(yr0_5:yr_65, ~sum(.x, na.rm = TRUE))) %>%
  mutate(total = yr0_5 + yr5_10 + yr10_20 + yr20_65 + yr_65) %>%
  mutate(prop0_5 = yr0_5/total , prop5_20 = (yr5_10 +yr10_20) /total,
         prop20_65 = yr20_65/total, prop65 = yr_65 /total ) %>%
  mutate(time_period = "Pre-pandemic")

post_uniform = plot_un %>%
  filter(time > 2020.5) %>%
  filter(time < 2023.5) %>%
  dplyr::summarize(across(yr0_5:yr_65, ~sum(.x, na.rm = TRUE))) %>%
  mutate(total = yr0_5 + yr5_10 + yr10_20 + yr20_65 + yr_65) %>%
  mutate(prop0_5 = yr0_5/total , prop5_20 = (yr5_10 +yr10_20) /total,
         prop20_65 = yr20_65/total, prop65 = yr_65 /total ) %>%
  mutate(time_period = "Post-pandemic")

uniform_age = rbind(pre_uniform, post_uniform)
head(uniform_age)

ag = c("0-4", "5-19", "20-64", "65+" ,"0-4", "5-19", "20-64", "65+")
pb = c(.157, .251, .413, .178, .155, .231, .415, .198)
tp = c(rep("Pre-pandemic", 4), rep("Post-pandemic", 4))
uniform = data.frame(ag, pb, tp) %>% mutate(Immunity = "Uniform immunity")

model_uni = uniform %>% mutate(tp = factor(tp, levels = c("Pre-pandemic", "Post-pandemic"))) %>%
  mutate(ag = factor(ag, levels = c("0-4", "5-19", "20-64", "65+"))) %>%
  ggplot(aes(x =ag, y = pb, fill = tp)) +
  geom_col(position = "dodge") +
  # facet_wrap(vars(Immunity)) +
  theme_bw() +
  xlab("Age group (years)") + ylab("Proportion of hospitalizations") +
  scale_fill_manual(values = re_pal) +
  labs(fill = "Time period") + 
  theme(legend.position="bottom") +
  theme(text = element_text(size=10))  +
  theme(strip.background = element_rect(fill="gray94")) 

model_uni

# Predicted age structure - tiered immunity 
pre_tiered = plot_ti %>%
  filter(time < 2020.5) %>%
  dplyr::summarize(across(yr0_5:yr_65, ~sum(.x, na.rm = TRUE))) %>%
  mutate(total = yr0_5 + yr5_10 + yr10_20 + yr20_65 + yr_65) %>%
  mutate(prop0_5 = yr0_5/total , prop5_20 = (yr5_10 +yr10_20) /total,
         prop20_65 = yr20_65/total, prop65 = yr_65 /total ) %>%
  mutate(time_period = "Pre-pandemic")

post_tiered = plot_ti %>%
  filter(time > 2020.5) %>%
  filter(time < 2023.5) %>%
  dplyr::summarize(across(yr0_5:yr_65, ~sum(.x, na.rm = TRUE))) %>%
  mutate(total = yr0_5 + yr5_10 + yr10_20 + yr20_65 + yr_65) %>%
  mutate(prop0_5 = yr0_5/total , prop5_20 = (yr5_10 +yr10_20) /total,
         prop20_65 = yr20_65/total, prop65 = yr_65 /total ) %>%
  mutate(time_period = "Post-pandemic")

tiered_age = rbind(pre_tiered, post_tiered)
head(tiered_age)  

ag = c("0-4", "5-19", "20-64", "65+" ,"0-4", "5-19", "20-64", "65+")
pb = c(.158, .25, .42, .17, .157, .36, .37, .11)
tp = c(rep("Pre-pandemic", 4), rep("Post-pandemic", 4))
tiered = data.frame(ag, pb, tp) %>% mutate(Immunity = "Tiered immunity")

model_t = tiered %>% mutate(tp = factor(tp, levels = c("Pre-pandemic", "Post-pandemic"))) %>%
  mutate(ag = factor(ag, levels = c("0-4", "5-19", "20-64", "65+"))) %>%
  ggplot(aes(x =ag, y = pb, fill = tp)) +
  geom_col(position = "dodge") +
  theme_bw() +
  xlab("Age group (years)") + ylab("Proportion of hospitalizations") +
  scale_fill_manual(values = re_pal) +
  labs(fill = "Time period") + 
  theme(legend.position="bottom") +
  theme(text = element_text(size=10))  +
  theme(strip.background = element_rect(fill="gray94"))  +
  theme(legend.position="none") 
model_t

# Combine tiered and uniform 
model_ageplot = rbind(uniform, tiered) %>% mutate(tp = factor(tp, levels = c("Pre-pandemic", "Post-pandemic"))) %>%
  mutate(ag = factor(ag, levels = c("0-4", "5-19", "20-64", "65+"))) %>%
  ggplot(aes(x =ag, y = pb, fill = tp)) +
  geom_col(position = "dodge") +
  facet_wrap(vars(Immunity), nrow = 2) +
  theme_bw() +
  xlab("Age group") + ylab("Proportion of healthcare encounters") +
  scale_fill_manual(values = re_pal) +
  labs(fill = "Time period") + 
  theme(axis.title.y = element_text(size = 8)) +
  theme(legend.position="none") +
  theme(text = element_text(size=7))  +
  theme(strip.background = element_rect(fill="gray94")) +
  theme(axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 9, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 9),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size =9, color = "black"),
        strip.text = element_text(colour = "black", size = 9, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))
model_ageplot

# Plot altogether 
figure3 = plot_grid(a, b, c, model_ageplot, rel_widths = c(2/3, 1/3, 2/3, 1/3),
          rel_heights = c(1, 1, 3, 3), labels = c("a", "b", "c", "d"))
figure3

ggsave(
  "Figure3.pdf",
  plot = figure3,
  device = "pdf",   
  width = 250,
  height = 178,
  units = "mm"
)



