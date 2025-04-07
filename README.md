Multiplex serology reveals age-specific immunodynamics of respiratory pathogens in the wake of the COVID-19 pandemic 

Overview 
The rebound of endemic respiratory viruses following the COVID-19 pandemic was marked by atypical transmission dynamics, with children experiencing increased disease burden and out-of-season epidemics as restrictions relaxed. Here we used serology from a newly developed quantitative multiplex assay to assess the post-pandemic immunity debt, a drop in immunity due to a lack of endemic virus circulation during COVID-19. We assessed age-specific antibody dynamics in Seattle, Washington, US, across a broad range of respiratory viruses, including influenza, respiratory syncytial virus, seasonal coronaviruses, and SARS-CoV-2. We found that respiratory virus immunodynamics differed between individuals <5 years of age and older individuals, with young children experiencing both larger boosts and quicker waning of antibodies across pathogens. We confirmed that these patterns are upheld in a non-pandemic setting by analyzing influenza serological data from South Africa. We incorporated our serological insights into an influenza transmission model calibrated to epidemiological data from Seattle and show that consideration of age-specific immunodynamics may be important to anticipate the effects of pandemic perturbations. 

Software requirements 
All analyses were conducted in R Studio Version 2024.12.1+563 (2024.12.1+563). 

Installation Guide 
R Studio can be installed at this link: https://posit.co/download/rstudio-desktop/. Typical installation is under ten minutes using an updated operating system.

Demo 
The analysis is separated into two phases, 1) influenza tranmission model and 2) serosolver model for individual-level serology data.

1) influenza_model: All code needed to generate age-structured SEIRS influenza transmission model and produce relevant figures. Expected outputs and run time: Model generated for tiered and uniform immunity structure time series and age structure compared to observed hospitalization data, 25 min run time. 
   
  - seattle_flu_model_outputs.R: Provides model code for age-structured compartmental Susceptible-Exposed-Infected-Recovered influenza       transmission model.

  - toy_odin.R: Code required to integrate ordinary differential equations in seattle_flu_model_outputs.R. Must be present in R       
    environment to run analysis in seattle_flu_model_outputs.R.
  
  - mixing_75.csv: Characterizes mixing patterns between 75 age classes, in 60 1-month age classes from ages 0-5, and in 5-year age   
    bands from ages 5-80. Mixing matrix was adapted from Mossing et al. 2008 (10.1371/journal.pmed.0050074). Must be present in R 
    environment to run analysis in seattle_flu_model_outputs.R.

2) serosolver_model: All code needed to produce estimates of boosting and waning rates from South Africa and King County serological data. Expected outputs and run time: Figures related to serology and relevant statistical analyses.

  - King_County_Serosolver.R: Code to produce boosting and waning estimates for one sample antigen using the Serosolver model on King        County serology data. Example code provided for cross-sectional pediatric serology and adult paired serology separately. Simulated   
    data provided to run the model in King_County_simulated_dat.csv.
    
  - South_Africa_Serosolver.R --> Set up for running serosolver model on one antigen using PHIRST South Africa data 
    serology_manuscript_analysis_plots.R --> Code to produce all serology figures and run statistical analyses 
