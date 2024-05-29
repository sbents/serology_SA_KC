## Make final clean figures for serology paper 
library(ggplot2)
library(cowplot)
library(dplyr)



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
  scale_fill_manual(name = "Subtype", values = c("lightsteelblue1", "lightsteelblue2", "lightsteelblue3", "lightsteelblue4", "gray40", "gray30"))
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
  scale_fill_manual(name = "Subtype", values = c("lightsteelblue1", "lightsteelblue2", "lightsteelblue3", "lightsteelblue4")) +
  guides(fill="none")
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
plot_grid(b_kc1, b_kc2, b_sa, w_kc1, w_kc2, w_sa, nrow = 2, rel_widths = c(1, .75, 1, 1, .75, 1), labels = c("A", "C", "E", "B", "D", "F"))
