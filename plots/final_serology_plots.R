## Make final clean figures for serology paper 
library(ggplot2)
library(cowplot)
library(dplyr)


setwd("C:/Users/bentssj/OneDrive - National Institutes of Health/Year_2024/serology/seattle")


# Figure 1A
### load in endemic respiratory monthly incidence data
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

sum(endemic_incidence$num_tested)

### adds covid king county 
covid_incidence = read.csv("cov_kingcount.csv") %>%
  mutate(date = as.Date(date, "%m/%d/%Y"),month = month(date), year = year(date))%>%
  mutate(time = month/12 + year)%>%
  group_by(time) %>% 
  dplyr::summarize(across(pcr_test_count:pcr_test_pos_count, ~sum(.x, na.rm = TRUE))) %>%
  mutate(pct_pos = pcr_test_pos_count/pcr_test_count, virus = "CoV N") %>%
  dplyr::select(time, virus, pct_pos)

### all pathogens incidence curve
all_path = rbind(endemic_incidence, covid_incidence) %>%
  mutate(Virus = factor(virus, levels = c("Flu A/H3", "Flu A/H1", "Flu B", "RSV", 
                                          "HCoV 229E", "HCoV NL63", "HCoV OC43", "HCoV HKU1", "CoV N")))


a = ggplot(data = all_path, aes(x = time, y = pct_pos, col = Virus))+
  annotate("rect", xmin = c(2020.58, 2021.5, 2022.33), xmax = c(2020.92, 2021.95, 2022.67), 
           ymin = 0, ymax = Inf, alpha = .15, fill = c("orangered", "orangered", "orangered"))+
  annotate("rect", xmin = c(2021.83, 2022.25), xmax = c(2021.95, 2022.75), 
           ymin = 0, ymax = Inf, alpha = .5, fill = c("lightblue", "lightblue"))+
  geom_line(lwd = 1.2) +
  scale_color_manual(values = c( "salmon1","tomato","darkred", "royalblue2", "lightblue","darkslategray3", "darkslategray4","darkslategray", "black" )) +
  theme_bw() +
  ylab("PCR positive (%)") +
  theme(legend.position = "bottom")

#_________________________________________________________________________
# Fig 1B 
### children MSD serology samples
cov = read.csv("serosamples_CoV_22Mar2024.csv")
resp = read.csv("serosamples_Resp_22Mar2024.csv") 
dem = read.csv("demos_children.csv") %>%
  mutate(Sample = ï..sample_id)

ped_dat = left_join(ped_serology, dem, by = "Sample") %>%
  mutate(blood_date = as.Date(blood_date), year = year(blood_date)) %>%
  dplyr::select(Assay, log_titer, year) %>%
  mutate(Population = "Pediatric")


### adult samples 
# use paired sera to get it 
adult_serologyc = read.csv("C:/Users/bentssj/OneDrive - National Institutes of Health/Scenario Hub/Seattle Flu Study/adult_serology_cov.csv") 
adult_serologyr = read.csv("C:/Users/bentssj/OneDrive - National Institutes of Health/Scenario Hub/Seattle Flu Study/adult_serology_res.csv")

adult_serology = rbind(adult_serologyc, adult_serologyr) %>%
  mutate(date = as.Date(Collection.Date), year = year(date), 
         log_titer = log(titer_mean)) %>%
  dplyr::select(Assay, log_titer, year) %>%
  mutate(Population ="Adult")

# Join children and adults 
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
  scale_fill_manual(values = c("lightblue2", "salmon1")) +
  theme(legend.position = "bottom") +
  xlab("Year") + 
  ylab("Log titer") +
  theme(strip.background =element_rect(fill="lightgray"))


# ____________________________________________________________________________
# Fig 1c-d
library(dplyr)

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

c = ggplot(data = hm1, aes(x = age_band, y = Virus)) +
  geom_tile(aes(fill = change)) +
  scale_fill_gradient2(low = "tomato", mid = "white", high = "royalblue1") +
  theme_bw() +
  xlab("Age band (years)") +
  ylab("") +
  theme(legend.position = "bottom") + 
  labs(fill = "Percent change (%)")


yr_2022 = hm_dat %>%
  filter(year == 2022)%>%
  mutate(mean_2022 = mean)

hm2 = left_join(yr_2021, yr_2022, by = c("age_band", "Virus")) %>%
  mutate(change = (mean_2022 - mean_2021)/mean_2021)
hm2$age_band = factor(hm2$age_band, levels = c("<1", "1-2", "3-4", "5-10",
                                              "18-49", "50-64", "65+" ))

d = ggplot(data = hm2, aes(x = age_band, y = Virus)) +
  geom_tile(aes(fill = change)) +
  scale_fill_gradient2(low = "tomato", mid = "white", high = "royalblue3") +
  theme_bw() +
  xlab("Age band (years)") +
  ylab("") +
  theme(legend.position = "bottom") + 
  labs(fill = "Percent change (%)")


plot_grid(a, b, c, d, nrow =2, labels = c("a", 'b', 'c', 'd'))
  
  


## Figure 2 
### kING County rsv/hcov
# plot
mu_lower = c(3.72, 4.31, .37, 3.79, 4.38, 5.41, 1.29, 1.15, 1.13, 1.26, 3.49, 1.18) 
mu_short = c(4.36, 4.61, 2.70, 4.21, 4.77, 5.79, 1.45, 1.27, 1.27, 1.49, 3.85, 1.42)
mu_upper = c(5.00, 4.92, 4.23, 4.62, 5.17, 6.15, 1.61, 1.39, 1.41, 1.73, 4.20, 1.70)
wane_lower = c(.72, .87, 0, .02, .01, .90, 0, 0, 0, 0, .13, 0)
wane_short = c(.92, .96, .37, .11, .15, .96, .02, .03, .01, .02, .20, .03)
wane_upper = c(.99, .99, .98, .29, .36, .99, .06, .08, .05, .06, .26, .08)
sub = c("HCoV 229E", "HCoV NL63", "HCoV OC43", "HCoV HKU1", "CoV N", "RSV","HCoV 229E", "HCoV NL63", "HCoV OC43", "HCoV HKU1", "CoV N", "RSV" )
age = c(rep("<5 yo", 6), rep( "> 18 yo", 6))

seattle_hrsv = data.frame(mu_lower, mu_short, mu_upper, 
                         wane_lower, wane_short, wane_upper, sub, age)
seattle_hrsv$sub = factor(seattle_hrsv$sub, levels = c ("HCoV 229E", "HCoV NL63", "HCoV OC43", "HCoV HKU1", "CoV N", "RSV"))



b_kc1 = ggplot(data = seattle_hrsv, aes(x = factor(age, level = c("<5 yo",  "> 18 yo")), y = mu_short, fill = sub)) + 
  geom_col(position = "dodge")+
  geom_errorbar(
    aes(ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "gray60", position = position_dodge(.9))+ 
  theme_light()+
  ylab("Short-term log titer boost")+
  theme(axis.text=element_text(size=12)) +
  xlab("Age group") +
  scale_fill_manual(name = "Subtype", values = c("lightsteelblue1", "lightsteelblue2", "lightsteelblue3", "lightsteelblue4", "gray30", "black")) +
  ylim(c(0, 6.5)) +
  geom_errorbar(
    aes(ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "gray60", position = position_dodge(.9))
b_kc1

w_kc1 = ggplot(data = seattle_hrsv, aes(x = factor(age, level = c("<5 yo",  "> 18 yo")), y = wane_short, fill =sub)) + 
  geom_errorbar(
    aes(ymin = wane_lower, 
        ymax = wane_upper), col = c("lightsalmon", "salmon1", "tomato1", "tomato3",  "tomato4", "black", "lightsalmon", "salmon1", "tomato1", "tomato3",  "tomato4", "black"), lwd = .8, width = .2, position = position_dodge(.9))+ 
  theme_light()+
  geom_point(position = position_dodge(.9), cex = 6,
             shape = 23)+
  ylab("Percent boost lost in one year")+
  xlab("Age group")+
  theme(axis.text=element_text(size=12)) +
  scale_fill_manual(name = "Subtype", values = c("lightsalmon", "salmon1", "tomato1", "tomato3",  "tomato4", "black"))
w_kc1


## seattle flu 
mu_lower = c(1.88, .42, .12, .10, 1.66, 1.71, 1.34, 1.47 ) 
mu_short = c(3.39, 2.48, 2.11, 1.80, 1.84, 1.93, 1.52, 1.65 )
mu_upper = c(4.92, 4.73, 4.16, 4.03, 2.03, 2.22, 1.75, 1.84)
wane_lower = c(.39, .12, .05, .06, .05, .05, .03, .02)
wane_short = c(.82, .63, .52, .59, .12, .12, .09, .08 )
wane_upper = c(.99, .99, .99, .99, .17, .18, .14, .14 )
sub = c("Flu A/H3", "Flu A/H1", "Flu B/Vic", "Flu B/Yam", "Flu A/H3", "Flu A/H1", "Flu B/Vic", "Flu B/Yam")
age = c(rep("<5 yo", 4), rep( "> 18 yo", 4))

seattle_flu = data.frame(mu_lower, mu_short, mu_upper, 
                          wane_lower, wane_short, wane_upper, sub, age)
seattle_flu$sub = factor(seattle_flu$sub, levels = c ("Flu A/H3", "Flu A/H1", "Flu B/Vic", "Flu B/Yam"))


b_kc2 = ggplot(data = seattle_flu, aes(x = factor(age, level = c("<5 yo",  "> 18 yo")), y = mu_short, fill = sub)) +
  geom_col(position = "dodge")+
  geom_errorbar(
    aes(ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "gray60", position = position_dodge(.9))+ 
  theme_light()+
  ylab("Short-term log titer boost")+
  theme(axis.text=element_text(size=12)) +
  xlab("Age group") +
  scale_fill_manual(name= "Subtype", values = c("lightsteelblue1", "lightsteelblue2", "lightsteelblue3", "lightsteelblue4")) +
  guides(fill="none") +
  ylim(c(0, 6.5))
b_kc2

w_kc2 = ggplot(data = seattle_flu, aes(x = factor(age, level = c("<5 yo",  "> 18 yo")), y = wane_short, fill =sub)) + 
  geom_errorbar(
    aes(ymin = wane_lower, 
        ymax = wane_upper), col = c( "salmon1", "tomato1", "tomato3",  "tomato4",  "salmon1", "tomato1", "tomato3",  "tomato4"), lwd = .8, width = .2, position = position_dodge(.9))+ 
  theme_light()+
  geom_point(position = position_dodge(.9), cex = 6,
             shape = 23)+
  ylab("Percent boost lost in one year")+
  xlab("Age group")+
  theme(axis.text=element_text(size=12)) +
  scale_fill_manual(name = "Subtype", values = c("salmon1", "tomato1", "tomato3",  "tomato4")) +
  guides(fill="none")
w_kc2

######### South Africa 

mu_lower = c(1.30, 1.34, 1.03, 1.14, .77, .88, .86, .77, .89, .93, .77, .83) 
mu_short = c(1.43, 1.47, 1.21, 1.26, .86, .97, .93, .86, .94, .98, .82, .88 )
mu_upper = c(1.56, 1.60, 1.45, 1.37, .95, 1.06, 1, .95, 1.00, 1.02, .86, .93)
wane_lower = c(.13, .07, .11, .28, .04, 0, .08, .09, .04,.10, .11, .10)
wane_short = c(.20, .16, .23, .34, .12, .01, .13, .15, .08, .14, .15, .13)
wane_upper = c(.26, .24, .37, .39, .19, .05, .18, .20, .13, .17, .20, .17 )
sub = c("Flu A/H3", "Flu A/H1", "Flu B/Vic", "Flu B/Yam", "Flu A/H3", "Flu A/H1", "Flu B/Vic", "Flu B/Yam", "Flu A/H3", "Flu A/H1", "Flu B/Vic", "Flu B/Yam")
age = c(rep("<5 yo", 4),rep("5-10 yo", 4), rep( "> 18 yo", 4))


sa_flu = data.frame(mu_lower, mu_short, mu_upper, 
                         wane_lower, wane_short, wane_upper, sub, age)
sa_flu$sub = factor(sa_flu$sub, levels = c ("Flu A/H3", "Flu A/H1", "Flu B/Vic", "Flu B/Yam"))


b_sa= ggplot(data = sa_flu, aes(x = factor(age, level = c("<5 yo", "5-10 yo", "> 18 yo")), y = mu_short, fill = sub)) + 
  geom_col(position = "dodge")+
  geom_errorbar(
    aes(ymin = mu_lower, 
        ymax = mu_upper), lwd = .8, width = .2, col = "gray60", position = position_dodge(.9))+ 
  theme_light()+
  ylab("Short-term log titer boost")+
  theme(axis.text=element_text(size=12)) +
  xlab("Age group") +
  scale_fill_manual(name = "Subtype", values = c("lightsteelblue1", "lightsteelblue2", "lightsteelblue3", "lightsteelblue4"))

w_sa = ggplot(data = sa_flu, aes(x = factor(age, level = c("<5 yo", "5-10 yo", "> 18 yo")), y = wane_short, fill =sub)) + 
  geom_errorbar(
    aes(ymin = wane_lower, 
        ymax = wane_upper), col = c("salmon1", "tomato1", "tomato3",  "tomato4", "salmon1", "tomato1", "tomato3",  "tomato4", "salmon1", "tomato1", "tomato3",  "tomato4"), lwd = .8, width = .2, position = position_dodge(.9))+ 
  theme_light()+
  geom_point(position = position_dodge(.9), cex = 6,
             shape = 23)+
  ylab("Percent boost lost in one year")+
  xlab("Age group")+
  theme(axis.text=element_text(size=12)) +
  scale_fill_manual(name = "Subtype", values = c("salmon1", "tomato1", "tomato3",  "tomato4"))
w_sa


##### Figure

library(cowplot)

fig2 = plot_grid(b_kc1, b_kc2, b_sa, w_kc1, w_kc2, w_sa, nrow = 2, rel_widths = c(1, .75, 1, 1, .75, 1), labels = c("a", "b", "c", "d", "e", "f"))
fig2


jpeg("Figure_2_kinetics.jpeg", units = "in", width=14, height=7, res=300)
ggarrange(fig2)

dev.off()

