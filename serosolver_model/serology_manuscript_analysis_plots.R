# Make clean figures for serosolver outputs 
library(ggplot2)
library(cowplot)
library(dplyr)
library(MetBrewer)
library(lubridate)
library(tidyr)
library(cowplot)  

# Set working directory 
setwd("/Users/sambents/Desktop/NIH/serology/unprocessed")

# Figure 1A
#_________________________________________________________________________

#Load in endemic respiratory monthly incidence data
endemic_incidence = read.csv("monthly_incidence_2018_2023.csv") %>%
  mutate(pct_pos = num_pos/num_tested, time = as.numeric(collection_month)/12 + as.numeric(collection_year) )%>%
  filter(virus != "Flu A") %>% # non-subtyped Flu A 
  mutate(virus = replace(virus, virus == "Corona 229E", "HCoV 229E")) %>%
  mutate(virus = replace(virus, virus == "Corona HKU1", "HCoV HKU1")) %>%
  mutate(virus = replace(virus, virus == "Corona NL63", "HCoV NL63")) %>%
  mutate(virus = replace(virus, virus == "Corona OC43", "HCoV OC43")) %>%
  mutate(virus = replace(virus, virus == "Flu A H1 2009", "Flu A/H1")) %>%
  mutate(virus = replace(virus, virus == "Flu A H3", "Flu A/H3")) %>%
  dplyr::select(time, virus, pct_pos)

# Number tested total 
end =  read.csv("monthly_incidence_2018_2023.csv")  %>%
  distinct(collection_month, collection_year, num_tested)
sum(end$num_tested)
print(unique(endemic_incidence$virus))

# Add COVID-19 time series 
covid_incidence = read.csv("cov_kingcount.csv") %>%
  mutate(date = as.Date(date, "%m/%d/%Y"),month = month(date), year = year(date))%>%
  mutate(time = month/12 + year)%>%
  group_by(time) %>% 
  dplyr::summarize(across(pcr_test_count:pcr_test_pos_count, ~sum(.x, na.rm = TRUE))) %>%
  mutate(pct_pos = pcr_test_pos_count/pcr_test_count, virus = "CoV N") %>%
  dplyr::select(time, virus, pct_pos)

# Combine covid and endemic pathogens
all_path = rbind(endemic_incidence, covid_incidence) %>%
  mutate(Virus = factor(virus, levels = c("Flu A/H3", "Flu A/H1", "Flu B", "RSV", 
                                          "HCoV 229E", "HCoV NL63", "HCoV OC43", "HCoV HKU1", "CoV N")))

re_pal <- met.brewer(name="Homer1", n=9)
a = ggplot(data = all_path, aes(x = time, y = pct_pos, col = Virus))+
  annotate("rect", xmin = c(2020.58, 2021.5, 2022.33), xmax = c(2020.92, 2021.95, 2022.67), 
           ymin = 0, ymax = Inf, alpha = .12, fill = c("darkred", "darkred", "darkred"))+
  annotate("rect", xmin = c(2021.83, 2022.25), xmax = c(2021.95, 2022.75), 
           ymin = 0, ymax = Inf, alpha = .3, fill = c("skyblue", "skyblue"))+
  geom_line(lwd = 1.2) +
 # scale_color_manual(values = c( "salmon1","darkred","skyblue1", "royalblue3", "plum1", "plum4","lightgreen", "darkgreen", "black" )) +
  scale_color_manual(values = re_pal) +
  theme_bw() +
  ylab("PCR positive (%)") +
  xlab("Year") +
  theme(legend.position = "bottom")+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))

a

#_________________________________________________________________________
# Fig 1B 

### Children MSD serology samples
# Covid
cov = read.csv("serosamples_CoV_22Mar2024.csv")
print(length(unique(cov$Sample)))
cov_df = data.frame(print(unique(cov$Sample)))

# Respiratory 
resp = read.csv("serosamples_Resp_22Mar2024.csv") 
print(length(unique(resp$Sample)))

# Demographic 
dem = read.csv("demos_children.csv") %>%
  mutate(Sample = sample_id)
print(length(unique(dem$Sample)))

# Combine pediatric serology 
ped_serology = rbind(cov, resp) %>%
  mutate(log_titer = log(titer_mean))
print(length(unique(ped_serology$Sample)))
ped_dat = left_join(ped_serology, dem, by = "Sample") %>%
  mutate(blood_date = as.Date(blood_date), year = year(blood_date)) %>%
  dplyr::select(Assay, log_titer, year) %>%
  mutate(Population = "Children <11 yo") 

# Adult samples, use paired sera
adult_serologyc = read.csv("adult_serology_cov.csv") 
adult_serologyr = read.csv("adult_serology_res.csv")
adult_serology = rbind(adult_serologyc, adult_serologyr) %>%
  mutate(date = as.Date(Collection.Date), year = year(date), 
         log_titer = log(titer_mean)) %>%
  dplyr::select(Assay, log_titer, year) %>%
  mutate(Population ="Adult >18 yo")

# Join pediatric and adult serology
all_age_serology = rbind(ped_dat, adult_serology)
Assay = print(unique(all_age_serology$Assay))
Virus = c("HCoV HKU1", "HCoV OC43", "HCoV 229E", "HCoV NL63", NA, NA, "CoV N",
          "CoV RBD", "CoV S", "Flu B/Yam", "Flu A/H1", "Flu A/H3", NA, "Flu B/Vic", "RSV", NA)
recode_assay = data.frame(Assay, Virus)

all_age_dat = left_join(all_age_serology, recode_assay, by = "Assay") %>%
  drop_na(Virus, log_titer, year) 
head(all_age_dat)

b = ggplot(data = all_age_dat, aes(x = as.character(year), y = log_titer, fill = Population )) +
  geom_violin( ) +
  facet_wrap(vars(Virus)) +
  theme_bw() +
  scale_fill_manual(values = c("skyblue", "darkred")) +
  theme(legend.position = "bottom") +
  xlab("Year") + 
  ylab("Log titer") +
  theme(strip.background =element_rect(fill="lightgray"))+
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))
b

# ____________________________________________________________________________
# Fig 1c-d
# Compare mean titer values throughout time 
ped_dat_age = left_join(ped_serology, dem, by = "Sample") %>%
  mutate(blood_date = as.Date(blood_date), year = year(blood_date), age = age_raw) %>%
  dplyr::select(age, year, Assay, log_titer)

adult_serology_age = rbind(adult_serologyc, adult_serologyr) %>%
  mutate(date = as.Date(Collection.Date), year = year(date), 
         log_titer = log(titer_mean), age = Age) %>%
  dplyr::select(age, year, Assay, log_titer)

# Join 
age_year_pop = rbind(ped_dat_age, adult_serology_age ) %>%
  mutate(age_band = ifelse(age < 1, "<1",
                           ifelse(age >= 1 & age < 3, "1-2",
                                  ifelse(age >= 3 & age < 5, "3-4", 
                                         ifelse(age >= 5 & age < 11, "5-10", 
                                                ifelse(age >= 20 & age < 49, "18-49",
                                                       ifelse(age >= 50 & age < 65, "50-64", "65+"))))))) %>%
  filter(log_titer > 0) %>%
  group_by(year, age_band, Assay) %>%
  dplyr::mutate(mean = exp(mean(log(log_titer)))) %>%
  distinct(year, age_band, Assay, mean)


Assay = print(unique(age_year_pop$Assay))
Virus = c("HCoV HKU1", "HCoV OC43", "HCoV 229E", "HCoV NL63", NA, NA, "CoV N",
          "CoV RBD", "CoV S", "Flu B/Yam", "Flu A/H1", "Flu A/H3", NA, "Flu B/Vic", "RSV")
recode_assay = data.frame(Assay, Virus)
hm_dat = left_join(age_year_pop, recode_assay, by = "Assay") %>%
  drop_na(Virus)
head(hm_dat)

yr_2020 = hm_dat  %>%
  filter(year == 2020) %>%
  mutate(mean_2020 = mean)
yr_2021 = hm_dat %>%
  filter(year == 2021)%>%
  mutate(mean_2021 = mean)
hm1 = left_join(yr_2020, yr_2021, by = c("age_band", "Virus")) %>%
  mutate(change = (mean_2021 - mean_2020)/mean_2020)

white <- "#FFFFFF"
skyblue1 <- "#87CEEB"

# Extract RGB components of each color
rgb_white <- col2rgb(white)
rgb_skyblue1 <- col2rgb(skyblue1)
# Blend the two colors: 50% white, 50% skyblue1
blended_rgb <- (0.25 * rgb_white + 0.75 * rgb_skyblue1) / 255
# Create the blended color
blended_color <- rgb(blended_rgb[1], blended_rgb[2], blended_rgb[3])
# Show the blended color
blended_color

# 2021 compared to 2020
c = ggplot(data = hm1, aes(x = age_band, y = Virus)) +
  geom_tile(aes(fill = change)) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "#A5DAF0") +
  theme_bw() +
  xlab("Age band (years)") +
  ylab("") +
  theme(legend.position = "bottom") + 
  labs(fill = "Percent change 2020-2021 (%)")+ 
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 7, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))

c

# 2022 compared to 2021, add adults since there is data available 
yr_2022 = hm_dat %>%
  filter(year == 2022)%>%
  mutate(mean_2022 = mean)
hm2 = left_join(yr_2021, yr_2022, by = c("age_band", "Virus")) %>%
  mutate(change = (mean_2022 - mean_2021)/mean_2021)
hm2$age_band = factor(hm2$age_band, levels = c("<1", "1-2", "3-4", "5-10",
                                               "18-49", "50-64", "65+" ))

d = ggplot(data = hm2, aes(x = age_band, y = Virus)) +
  geom_tile(aes(fill = change)) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "skyblue1") +
  theme_bw() +
  xlab("Age band (years)") +
  ylab("") +
  theme(legend.position = "bottom") + 
  labs(fill = "Percent change 2021-2022 (%)")  + 
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 7, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))
d

plot_grid(a, b, c, d, nrow =2, labels = c("a", 'b', 'c', 'd'))
# width 1000. height 800


###########################################################################################
## Figure 2, all outputs from serology models  

# King County RSV/HCoV boosting
mu_lower = c(3.68, 3.91, 3.85, 3.94,     4.29, 5.15, 4.76, 4.60,    2.40, 1.12, 1.14, 1.14,     2.50, 3.84, 8.95, 9.34)
mu_short = c( 4.06, 4.36, 4.09, 4.73,      4.67, 5.37, 4.97, 4.86,   2.55,  1.23, 1.25, 1.62,   2.66, 4.02, 9.35, 9.99)
mu_upper = c(4.46, 4.86, 4.38, 5.33,      5.09, 5.59, 5.21, 5.11,     2.80, 1.73, 1.42, 1.55,   2.72, 4.21, 9.84, 10.77)
wane_lower = c(.01, .09, .00, .09,        .01, .03, .00, .05,        .01, .00, .00, .00,       .01, .01, .00, .00)
wane_upper =  c(.13, .19, .08, .26,      .12, .13, .10, .19,          .02, .04, .09,  .02,     .01, .16, .07 ,.29)
wane_short = c(.06, .14, .03, .15,       .07, .08, .05, .11,        .01, .01, .03,  .01,       .01, .08, .02, .08 )
sub = c ("HCoV HKU1", "HCoV OC43", "HCoV NL63", "HCoV 229E", 'RSV', 'CoV N', 'CoV S', 'CoV RBD',
         "HCoV HKU1", "HCoV OC43", "HCoV NL63", "HCoV 229E", 'RSV', 'CoV N', 'CoV S', 'CoV RBD')
age = c(rep("<5 yo", 8), rep( "> 18 yo", 8))

seattle_hrsv = data.frame(mu_lower, mu_short, mu_upper, 
                          wane_lower, wane_short, wane_upper, sub, age)
seattle_hrsv$sub = factor(seattle_hrsv$sub, levels = c ("HCoV HKU1", "HCoV OC43", "HCoV NL63", "HCoV 229E", 'RSV', 'CoV N', 'CoV S', 'CoV RBD'))

re_pal <- met.brewer(name="Homer1", n=8)
b_kc1 = ggplot(data = seattle_hrsv, aes(x = factor(age, level = c("<5 yo",  "> 18 yo")), y = mu_short, fill = sub)) + 
  geom_col(position = "dodge")+
  scale_fill_manual(
    values = re_pal,            # Custom color palette
    name = "Subtype"            # Custom legend title for 'fill'
  ) +
  geom_errorbar(
    aes(ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "gray60", position = position_dodge(.9))+ 
  theme_light()+
  ylab("Log titer boost")+
  theme(axis.text=element_text(size=12)) +
  xlab("Age group") +
  ylim(c(0, 11.0)) +
  geom_errorbar(
    aes(ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "gray60", position = position_dodge(.9)) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 7, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA)) 
b_kc1

# King County RSV/HCoV waning 
w_kc1 = ggplot(data = seattle_hrsv, aes(x = factor(age, level = c("<5 yo",  "> 18 yo")), y = wane_short, fill =sub)) + 
  geom_errorbar(
       aes(ymin = wane_lower, 
         ymax = wane_upper), lwd = .8, width = .2, position = position_dodge(.9))+ 
  theme_light()+
  geom_point(position = position_dodge(.9), cex = 6,
             shape = 23)+
  ylab("Percent boost lost in 1 year")+
  xlab("Age group")+
  theme(axis.text=element_text(size=12)) +
  scale_fill_manual(
    values = re_pal,            # Custom color palette
    name = "Subtype"            # Custom legend title for 'fill'
  ) +
  theme(legend.position = "bottom") +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        #  plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 7, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))  +
  ylim(c(0, .35)) 
w_kc1

# King County influenza 
mu_lower = c(4.13, 3.72, 3.42, 3.84,1.26, 1.57, 1.62, 1.53) 
mu_short = c(4.47, 4.11, 3.97, 4.47,  1.36, 1.90, 1.66, 1.88)
mu_upper = c(4.84, 4.57, 4.61, 5.72,2.50, 2.08, 1.69, 2.05)
wane_lower = c(.18, .13, .08, .05,.02, .05, .00, .00 )
wane_short = c(.25, .21, .16, .16, .03, .07, .01, .01 )
wane_upper = c(.32, .28, .24, .27, .06, .08, .01, .01)

sub = c("Flu B/Yam", "Flu A/H1", "Flu B/Vic",
        "Flu A/H3", "Flu B/Yam", "Flu A/H1", "Flu B/Vic",
        "Flu A/H3")
age = c(rep("<5 yo", 4), rep( "> 18 yo", 4))

seattle_flu = data.frame(mu_lower, mu_short, mu_upper, 
                         wane_lower, wane_short, wane_upper, sub, age)
seattle_flu$sub = factor(seattle_flu$sub, levels = c ("Flu A/H3", "Flu A/H1", "Flu B/Vic", "Flu B/Yam"))

# King County influenza boosting
b_kc2 = ggplot(data = seattle_flu, aes(x = factor(age, level = c("<5 yo",  "> 18 yo")), y = mu_short, fill = sub)) + 
  geom_col(position = "dodge")+
  scale_fill_manual(
    values = re_pal,            # Custom color palette
    name = "Subtype"            # Custom legend title for 'fill'
  ) +
  geom_errorbar(
    aes(ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "gray60", position = position_dodge(.9))+ 
  theme_light()+
  ylab("Log titer boost")+
  theme(axis.text=element_text(size=12)) +
  xlab("Age group") +
  # scale_fill_manual(name = "Subtype", values = c("lightsteelblue1", "lightsteelblue2", "lightsteelblue3", "lightsteelblue4", "gray30", "black")) +
  ylim(c(0, 6.5)) +
  geom_errorbar(
    aes(ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "gray60", position = position_dodge(.9)) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        #  plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 7, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        #  legend.margin=margin(1,1.5,0.5,0.5,unit = "line"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA)) 
b_kc2


# King County influenza waning
w_kc2 = ggplot(data = seattle_flu, aes(x = factor(age, level = c("<5 yo",  "> 18 yo")), y = wane_short, fill =sub)) + 
  geom_errorbar(
    aes(ymin = wane_lower, 
        ymax = wane_upper), lwd = .8, width = .2, position = position_dodge(.9))+ 
  theme_light()+
  geom_point(position = position_dodge(.9), cex = 6,
             shape = 23)+
  ylab("Percent boost lost in 1 year")+
  xlab("Age group")+
  theme(axis.text=element_text(size=12)) +
  #geom_col(position = "dodge")+
  scale_fill_manual(
    values = re_pal,            # Custom color palette
    name = "Subtype"            # Custom legend title for 'fill'
  )  +
  theme(legend.position = "bottom") +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        #  plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 7, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        #  legend.margin=margin(1,1.5,0.5,0.5,unit = "line"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))  +
  ylim(c(0, .35)) 
  # scale_fill_manual(name = "Subtype", values = c("lightsalmon", "salmon1", "tomato1", "tomato3",  "tomato4", "black")) +
 # ggtitle("King County")
w_kc2

# South Africa boosting and waning 
mu_lower = c(5.09, 5.31, 3.12, 4.17 ,     2.67, 2.67, 1.85, 2.75,   3.53, 3.48,2.56, 2.95) 
mu_short = c( 5.50, 5.86, 3.87, 4.66,    2.90, 2.88, 2.14, 3.17,    3.73, 3.62, 2.77, 3.13)
mu_upper = c(5.87, 6.40 , 4.71, 5.14,    3.13, 3.09, 2.45, 3.83,    3.93, 3.80, 2.99, 3.31)
wane_lower = c(.08, 0 , .03, .07,    .06, 0, .01, .07,             .04, .04, .05, .05 )
wane_short = c(.13, .09 , .13, .17,  .08, .04, .04, .09,            .05, .05, .05, .06)
wane_upper = c(.16, .15 , .21, .25,  .09, .07, .07, .11,             .05, .06, .06, .07)
sub = c("Flu A/H3", "Flu A/H1", "Flu B/Vic", "Flu B/Yam", "Flu A/H3", "Flu A/H1", "Flu B/Vic", "Flu B/Yam", "Flu A/H3", "Flu A/H1", "Flu B/Vic", "Flu B/Yam")
age = c(rep("<5 yo", 4),rep("5-10 yo", 4), rep( "> 18 yo", 4))

sa_flu = data.frame(mu_lower, mu_short, mu_upper, 
                    wane_lower, wane_short, wane_upper, sub, age)
sa_flu$sub = factor(sa_flu$sub, levels = c ("Flu A/H3", "Flu A/H1", "Flu B/Vic", "Flu B/Yam"))

# South Africa influenza boosting
b_sa = ggplot(data = sa_flu, aes(x = factor(age, level = c("<5 yo","5-10 yo",  "> 18 yo")), y = mu_short, fill = sub)) + 
  geom_col(position = "dodge")+
  scale_fill_manual(
    values = re_pal,            # Custom color palette
    name = "Subtype"            # Custom legend title for 'fill'
  ) +
  geom_errorbar(
    aes(ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "gray60", position = position_dodge(.9))+ 
  theme_light()+
  ylab("Log titer boost")+
  theme(axis.text=element_text(size=12)) +
  xlab("Age group") +
  ylim(c(0, 6.5)) +
  geom_errorbar(
    aes(ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "gray60", position = position_dodge(.9)) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 7, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA)) 
b_sa

# South Africa influenza waning
w_sa= ggplot(data = sa_flu, aes(x = factor(age, level = c("<5 yo", "5-10 yo", "> 18 yo")), y = wane_short, fill =sub)) + 
  geom_errorbar(
    aes(ymin = wane_lower, 
        ymax = wane_upper), lwd = .8, width = .2, position = position_dodge(.9))+ 
  theme_light()+
  ylim(c(0, .35)) +
  geom_point(position = position_dodge(.9), cex = 6,
             shape = 23)+
  ylab("Percent boost lost in 1 year")+
  xlab("Age group")+
  theme(axis.text=element_text(size=12)) +
  scale_fill_manual(
    values = re_pal,            # Custom color palette
    name = "Subtype"            # Custom legend title for 'fill'
  )  +
  theme(legend.position = "bottom")+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 7, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))  
w_sa

# Full Figure 2 
fig2 = plot_grid(b_kc1, b_kc2, b_sa, w_kc1, w_kc2, w_sa, nrow = 2, rel_widths = c(1, .75, 1, 1, .75, 1), labels = c("a", "b", "c", "d", "e", "f"))
fig2

# 1500,750
jpeg("Figure_2_kinetics.jpeg", units = "in", width=14, height=7, res=300)
ggarrange(fig2)
dev.off()



######################################
# End of figures 
#####################################

# Table S1: Counts and Kolmogorov-Smirnov tests by pathogen and age group comparing average 
# antibody concentration levels in the 2020 baseline year to 2021 and 2022.


ped_dat_age = left_join(ped_serology, dem, by = "Sample") %>%
  mutate(blood_date = as.Date(blood_date), year = year(blood_date), age = age_raw, sample_id = Sample) %>%
  dplyr::select(age, year, Assay, log_titer, sample_id)
adult_serology_age = rbind(adult_serologyc, adult_serologyr) %>%
  mutate(date = as.Date(Collection.Date), year = year(date), 
         log_titer = log(titer_mean), age = Age, sample_id  = Donor.ID) %>%
  dplyr::select(age, year, Assay, log_titer, sample_id)

# Join 
age_year_pop = rbind(ped_dat_age, adult_serology_age ) %>%
  mutate(age_band = ifelse(age < 1, "<1",
                           ifelse(age >= 1 & age < 3, "1-2",
                                  ifelse(age >= 3 & age < 5, "3-4", 
                                         ifelse(age >= 5 & age < 11, "5-10", 
                                                ifelse(age >= 20 & age < 49, "18-49",
                                                       ifelse(age >= 50 & age < 65, "50-64", "65+"))))))) %>%
  filter(log_titer > 0) %>%
  group_by(sample_id, Assay) %>%
  sample_n(1, replace = FALSE)  %>%
  drop_na(age_band)
age_group = print(unique(age_year_pop$age_band))

table = age_year_pop %>%
  group_by(year, age_band) %>%
  summarise(num_samp = n_distinct(sample_id))

# 
### Compare time points 
baseline  = age_year_pop %>%
  filter(year == 2020) %>%
  filter(Assay != "Flu-H7-Shanghai" ) %>%
  filter(Assay != "MERS Spike" ) %>%
  filter(Assay != "SARS-CoV-1 Spike" ) 
compare  = age_year_pop %>%
  filter(year == 2022) %>%
  filter(Assay != "Flu-H7-Shanghai" ) %>%
  filter(Assay != "MERS Spike" ) %>%
  filter(Assay != "SARS-CoV-1 Spike" ) 

# Run 
assay = print(unique(baseline$Assay)) 
age_band = print(unique(baseline$age_band))

for(i in age_band){
  for(j in assay) { 
    
    baseline_i = baseline %>%  filter(age_band == i) %>%
      filter(Assay == j) 
    compare_i = compare %>% filter(age_band == i) %>%
      filter(Assay == j) 
    
    ks_run = (ks.test(baseline_i$log_titer, compare_i$log_titer))
    ks_stat = ks_run$p.value
    
    combine = data.frame(cbind(i, j, ks_stat))
    print(combine)
    
  } 
  
} 

fin_ks  = do.call(rbind, combine)
head(fin_ks)

##########################
# Vaccination rates in King County 
dem = read.csv("demos_children.csv") %>%
  mutate(Sample = ï..sample_id) %>%
  filter(flu_vx == TRUE)
head(dem)
391/999 # 39% vaccinated 

#load in meta data to get DOB
meta  = read.csv("meta.csv") 
head(meta)    
length(unique(meta$ind_id))
table(meta$receivefluvacc)


##########################
# Number of samples for South Africa by age group 
hai = read.csv("hai.csv") 
print(length(unique(hai$ind_id)))
meta = read.csv("meta.csv")  
head(meta)

full_sa = left_join(hai, meta, by = "ind_id") %>%
  mutate(date_format = as.Date(collectiondate, format = "%m/%d/%Y")) %>%
  mutate(year = year(date_format)) %>%
  dplyr::select(ind_id, dob_year, year ) %>%
  distinct(ind_id, dob_year, year) %>%
  count(ind_id) %>% filter(n == 1)
no = print(full_sa$ind_id)

full_sa = left_join(hai, meta, by = "ind_id") %>%
  mutate(date_format = as.Date(collectiondate, format = "%m/%d/%Y")) %>%
  mutate(year = year(date_format)) %>%
  dplyr::select(ind_id, dob_year, year ) %>%
  distinct(ind_id, dob_year, year) %>%
  filter(!ind_id %in% no) %>%
  mutate(age = year - dob_year) %>%
  group_by(ind_id) %>%
  mutate(max_age = max(age)) %>%
  mutate(adult = ifelse(max_age >= 18, 3,
                        ifelse(max_age < 5, 1, 2) )) %>%
  filter(year == 2018) %>%
  distinct(ind_id, adult) %>%
  drop_na()





