library(MuMIn)
library(ggplot2)
library(lme4)
library(tidyverse)
library(GGally)
library(xtable)
library(gamm4)
library(bbmle)
library(lubridate)
library(mgcv)
library(nlme)
library(ggeffects)
library(visreg)
library(scales)
library(zoo)
library(weathercan)
library(sp)
library(mapview)
library(mgcViz)
library(DHARMa)

# ggplot custom  ----------------------------------------------------------
theme_set(theme_bw())
scale_colour_discrete <- function(...,palette="Set1") {
  scale_colour_brewer(...,palette=palette)
}
scale_colour_orig <- ggplot2::scale_colour_discrete
scale_fill_discrete <- function(...,palette="Set1") {
  scale_fill_brewer(...,palette=palette)
}
# clean up 
rm(list = ls())

# what are the predictors of parturition date and mismatch ----------------------------
rm(list = ls())

trends_df <- read.csv2("data/mine/trends_data_noID.csv")

m1.g <- gamm4(birthdate ~ s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m2.g <- gamm4(birthdate ~ s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m3.g <- gamm4(birthdate ~ s(evi_tm1, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m4.g <- gamm4(birthdate ~ s(spring_prec, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m5.g <- gamm4(birthdate ~ s(spring_temp, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)

AICtab(m1.g$mer, m2.g$mer, m3.g$mer,m4.g$mer, m5.g$mer)
# dAIC df
# m2.g$mer  0.0 6 
# m3.g$mer  8.2 6 
# m1.g$mer  9.0 6 
# m4.g$mer 11.5 6 
# m5.g$mer 12.3 6 

plot(m2.g$gam) # this is linear
summary(m2.g$gam) # this is linear  1      1 18.91 2.01e-05 ***

m1.l <- lmer(birthdate ~ fall_temp + (1|mom_id) + (1|year), data = dt4)
m2.l <- lmer(birthdate ~ fall_prec + (1|mom_id) + (1|year), data = dt4)
m3.l <- lmer(birthdate ~ evi_tm1 + (1|mom_id) + (1|year), data = dt4)
# m4.g <- lmer(birthdate ~ spring_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
# m5.g <- lmer(birthdate ~ pring_temp, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)

AICtab(m1.l, m2.l, m3.l)
# dAIC df
# m2.l  0.0 5 
# m1.l 10.5 5 
# m3.l 16.3 5 

t=AICtab(m1.l, m2.l, m3.l)

t=kable(t) %>%
  kable_styling(font_size = 10) %>%
  kable_styling("bordered") %>%
  #save_kable(file = "Graphs/TableS1_Appendix2.html", self_contained = T)



summary(m2.l)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# mom_id   (Intercept)  75.07    8.664  
# year     (Intercept)  10.45    3.232  
# Residual             121.74   11.034  
# Number of obs: 261, groups:  mom_id, 75; year, 17
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  147.104      2.424  60.698
# fall_prec     15.045      3.460   4.348
# 
# Correlation of Fixed Effects:
#   (Intr)
# fall_prec -0.766

round(confint(m2.l),3)
# 2.5 %  97.5 %
#   .sig01        4.445  12.230
# .sig02        0.594   5.516
# .sigma        9.857  12.547
# (Intercept) 142.155 151.942
# fall_prec     8.220  22.041

summary(m1.l)
# fall_temp     -4.528      2.355  -1.922

round(confint(m1.l),3)
# 2.5 %  97.5 %
#   .sig01        3.423  11.710
# .sig02        2.923   8.005
# .sigma        9.961  12.794
# (Intercept) 154.413 195.018
# fall_temp    -9.223   0.107

MuMIn::r.squaredGLMM(m2.l) #    R2m      R2c [1,] 0.09478338 0.468289



# mismatch 
m1 <- gamm4(mismatch ~ s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m2 <- gamm4(mismatch ~ s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m3 <- gamm4(mismatch ~ s(evi_tm1, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m4 <- gamm4(mismatch ~ s(spring_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m5 <- gamm4(mismatch ~ s(spring_temp, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
dev.off()

AICtab(m1$mer, m2$mer, m3$mer,m4$mer, m5$mer)
# dAIC df
# m5$mer  0.0 6 
# m2$mer  7.8 6 
# m3$mer  8.0 6 
# m4$mer 10.1 6 
# m1$mer 10.3 6



# mismatch # REVISED FOLLOWING REV 3 ROUND OF MAJOR REVISIONS
# dt4$res.temp <- resid(lm(dt4$spring_temp ~ dt4$evi_up))

# replace spring temp by residual spring temp? nope
m1 <- gamm4(mismatch ~ s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m2 <- gamm4(mismatch ~ s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m3 <- gamm4(mismatch ~ s(evi_tm1, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m4 <- gamm4(mismatch ~ s(spring_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m5 <- gamm4(mismatch ~ s(spring_temp, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m6 <- gamm4(mismatch ~ s(spring_temp, k=10) + s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m7 <- gamm4(mismatch ~ s(spring_prec, k=10) + s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m8 <- gamm4(mismatch ~ s(spring_prec, k=10) + s(evi_tm1, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m9 <- gamm4(mismatch ~ s(spring_temp, k=10) + s(evi_tm1, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
#m10 <- gamm4(mismatch ~ s(spring_temp, k=10) + s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m10 <- gamm4(mismatch ~ s(spring_temp, k=10) + s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m11 <- gamm4(mismatch ~ s(spring_prec, k=10) + s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m12 <- gamm4(mismatch ~ s(evi_tm1, k=10) + s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m13 <- gamm4(mismatch ~ s(evi_tm1, k=10) + s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
#m14 <- gamm4(mismatch ~ s(spring_temp, k=10) + fall_prec, random = ~ (1|mom_id) + (1|year), data = dt4)


gam.check(m10$gam)
plot(m10$gam, pages=1)
#s(spring_temp) 9.00 1.97    0.63  <2e-16 ***

gam.check(m10_l$gam)
# end check
dev.off()

# summary(m10$mer)
# Linear mixed model fit by REML ['lmerMod']
# 
# REML criterion at convergence: 2081
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.4664 -0.4756 -0.1603  0.1947  5.4191 
# 
# Random effects:
#   Groups   Name           Variance Std.Dev.
# mom_id   (Intercept)     54.03    7.350  
# year     (Intercept)     24.58    4.958  
# Xr.0     s(fall_prec)     0.00    0.000  
# Xr       s(spring_temp)  33.17    5.759  
# Residual                129.17   11.365  
# Number of obs: 261, groups:  mom_id, 75; year, 17; Xr.0, 8; Xr, 8
# 
# Fixed effects:
#   Estimate Std. Error t value
# X(Intercept)          8.618      1.739   4.955
# Xs(spring_temp)Fx1    9.184      4.144   2.216
# Xs(fall_prec)Fx1      5.083      1.456   3.492
# 
# Correlation of Fixed Effects:
#   X(Int) Xs(s_)F1
# Xs(sprn_)F1  0.042         
# Xs(fll_p)F1 -0.077  0.007  


summary(m10$gam)

# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   mismatch ~ s(spring_temp, k = 10) + s(fall_prec, k = 10)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    8.618      1.739   4.955 1.31e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
# s(spring_temp) 2.447  2.447 15.99 2.18e-07 ***
#   s(fall_prec)   1.000  1.000 12.20 0.000564 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.406   
# lmer.REML =   2081  Scale est. = 129.17    n = 261
# 
confint(m10$mer)
# 2.5 %    97.5 %
#   .sig01              0.1316327 11.020699
# .sig02              2.2974823  7.840790
# .sig03              0.0000000 11.119779
# .sig04              0.0000000 15.767131
# .sigma             10.0920091 13.216284
# X(Intercept)        5.2331932 12.180876
# Xs(spring_temp)Fx1  1.2625998 18.599179
# Xs(fall_prec)Fx1   -0.7668800  8.030466

m10_l <- gamm4(mismatch ~ s(spring_temp, k=10) + fall_prec, random = ~ (1|mom_id) + (1|year), data = dt4)
summary(m10_l$gam)


# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   mismatch ~ s(spring_temp, k = 10) + fall_prec
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.09218    3.10448   0.030 0.976335    
# fall_prec   16.44865    4.71009   3.492 0.000564 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
# s(spring_temp) 2.447  2.447 15.99 2.18e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.406   
# lmer.REML = 2078.6  Scale est. = 129.17    n = 261

summary(m10_l$mer)
# Linear mixed model fit by REML ['lmerMod']
# 
# REML criterion at convergence: 2078.6
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.4664 -0.4756 -0.1603  0.1947  5.4191 
# 
# Random effects:
#   Groups   Name           Variance Std.Dev.
# mom_id   (Intercept)     54.02    7.350  
# year     (Intercept)     24.58    4.958  
# Xr       s(spring_temp)  33.17    5.759  
# Residual                129.17   11.365  
# Number of obs: 261, groups:  mom_id, 75; year, 17; Xr, 8
# 
# Fixed effects:
#   Estimate Std. Error t value
# X(Intercept)        0.09218    3.10448   0.030
# Xfall_prec         16.44865    4.71009   3.492
# Xs(spring_temp)Fx1  9.18410    4.14414   2.216
# 
# Correlation of Fixed Effects:
#   X(Int) Xfll_p
# Xfall_prec  -0.829       
# Xs(sprn_)F1  0.018  0.007

round(confint(m10_l$mer), 2)
# 2.5 % 97.5 %
#   .sig01              0.20  11.02
# .sig02              2.30   7.87
# .sig03              0.00  15.77
# .sigma             10.09  13.22
# X(Intercept)       -5.92   6.41
# Xfall_prec          7.03  25.44
# Xs(spring_temp)Fx1  1.26  18.60
round(confint(m10_l$gam), 2)





AICtab(m1$mer, m2$mer, m3$mer,m4$mer, m5$mer, m6$mer, m7$mer, m8$mer, m9$mer, m10$mer, m11$mer, m12$mer, m13$mer)
# dAIC df
# dAIC df
# m10$mer  0.0 8 
# m9$mer   6.6 8 
# m5$mer   7.7 6 
# m6$mer   8.3 8 
# m13$mer 14.5 8 
# m7$mer  15.5 8 
# m2$mer  15.6 6 
# m8$mer  15.6 8 
# m3$mer  15.7 6 
# m12$mer 15.8 8 
# m11$mer 17.7 8 
# m4$mer  17.9 6 
# m1$mer  18.1 6


library(kableExtra)

#t=AICtab(m1$mer, m2$mer, m3$mer, m4$mer, m5$mer)
t=AICtab(m1$mer, m2$mer, m3$mer,m4$mer, m5$mer, m6$mer, m7$mer, m8$mer, m9$mer, m10$mer, m11$mer, m12$mer, m13$mer)
t=kable(t, digits =2) %>%
  kable_styling(font_size = 10) %>%
  kable_styling("bordered") %>%
  save_kable(file = "Graphs/revisions2/Revised_TableS2_Appendix2.html", self_contained = T)

# summary(m5$gam) # s(spring_temp) 1.972  1.972 8.43 0.000215 ***
# # check 
# plot(m5$gam, pages=1)
# gam.check(m5$gam)
# #s(spring_temp) 9.00 1.97    0.63  <2e-16 ***
# # end check 

# check 
pdf("Graphs/revisions2/test.pdf", width = 8, height = 5)
par(mfrow=c(1, 2)) 
#plot.gam(m10$gam, select=1, xlabs = c("Autumn precipitation (mm)", "XX") )
# xlab=""
plot.gam(m10$gam, select=1, all.terms=T, shade=T, xlab="Spring temperature (°C)", ylab="Smoothed term for spring temperature", cex=0.8) 
# xlabs=""
plot.gam(m10$gam, select=2, all.terms=T, shade=T, xlab="Autumn precipitation (mm)", ylab="Smoothed term for autumn precipitation",cex=0.8)
dev.off()




# showing raw plot in supp
# par(mfrow=c(1,2))
# plot(m2.g$gam, xlab = "Autumn precipitation (mm)", ylab = "Smoothed term for autumn precipitation")
# plot(l$m5$gam, xlab = "Spring temperature (°C)", ylab = "Smoothed term for spring temperature")
# 
# 

# get estimates and variance - need mer only 
m10 <- lmer(mismatch ~ spring_temp + fall_prec + (1|mom_id) + (1|year), data = dt4)
summary(m10)
r.squaredGLMM(m10) 
# R2m       R2c
#0.2741147 0.5869986
m10a <- lmer(mismatch ~ fall_prec + (1|mom_id) + (1|year), data = dt4) # rm : 0.05998344
r.squaredGLMM(m10a) 
#0.2741147 -0.05998344 : 0.2141313

m10b <- lmer(mismatch ~ spring_temp + (1|mom_id) + (1|year), data = dt4)
summary(m10b)
r.squaredGLMM(m10b) # 0.1847434
#0.2741147 - 0.1847434=0.0893713

# venn diagram for response to reviewer 
library(partR2)
m10 <- lmer(mismatch ~ spring_temp + fall_prec + (1|mom_id) + (1|year), data = dt4)
res <- partR2(m10, partvars = c("spring_temp", "fall_prec"), nboot=100)

print(res)
# R2 (marginal) and 95% CI for the full model: 
#   R2     CI_lower CI_upper nboot ndf
# 0.2741 0.1253   0.4167   100   3  
# 
# Part (semi-partial) R2:
#   Predictor(s)          R2     CI_lower CI_upper nboot ndf
# Model                 0.2741 0.1253   0.4167   100   3  
# spring_temp           0.2134 0.0505   0.3635   100   2  
# fall_prec             0.0919 0.0000   0.2581   100   2  
# spring_temp+fall_prec 0.2741 0.1253   0.4167   100   1  

summary(res, round_to = 2)

# Structure coefficients range from −1 to 1 with their absolute value expressing the correlation relative to a perfect correlation 
# if a single predictor explains as much as the total fixed part of the model.
res$R2
res$R2[1,2]

p1 <- forestplot(res, type = "IR2")
# inclusive is calculated from structure coefficients - generates total variance explained by a predictor. 
# calculated as the squared structure coefficient, i.e., its contribution to the linear predictor independent of other predictors 
# times the proportion of variance explained by the linear predictor (which is the ‘total’ marginal R2 of the model)

# Inclusive R2 (SC^2 * R2):
#   Predictor   IR2  CI_lower CI_upper
# spring_temp 0.20 0.07     0.34    
# fall_prec   0.08 0.01     0.16   

# Inclusive R2 as we define it here, complements part R2 by giving additional insights.
# While part R2 quantifies the variance uniquely explained by a predictor (or set of predictors), 
# inclusive R2 quantifies the total proportion of variance explained in the model, both uniquely and jointly with other predictors. 

library("VennDiagram")

# move to new plotting page
grid.newpage()

# create pairwise Venn diagram
draw.pairwise.venn(area1=21, area2=9,cross.area=3,
                   category=c("Spring T°","Autumn P"),fill=c("Red","Yellow"))






library(kableExtra)

#t=AICtab(m1$mer, m2$mer, m3$mer, m4$mer, m5$mer)
t=AICtab(m1$mer, m2$mer, m3$mer,m4$mer, m5$mer, m6$mer, m7$mer, m8$mer, m9$mer, m10$mer, m11$mer, m12$mer, m13$mer)
t=kable(t) %>%
  kable_styling(font_size = 10) %>%
  kable_styling("bordered") %>%
  save_kable(file = "Graphs/revisions2/Revised_TableS2_Appendix2.html", self_contained = T)

summary(m5$gam) # s(spring_temp) 1.972  1.972 8.43 0.000215 ***
# check 
plot(m5$gam, pages=1)
gam.check(m5$gam)
#s(spring_temp) 9.00 1.97    0.63  <2e-16 ***
# end check 


# showing raw plot in supp
par(mfrow=c(1,2))
plot(m2.g$gam, xlab = "Autumn precipitation (mm)", ylab = "Smoothed term for autumn precipitation")
plot(l$m5$gam, xlab = "Spring temperature (°C)", ylab = "Smoothed term for spring temperature")

# Generating some new data for which you'd like predictions:

# 
# newdat <- data.frame(spring_temp=seq(min(dt4$spring_temp), max(dt4$spring_temp), length = 100),
#                      mother_id="A44",
#                      year="2010")
# table(dt4$mother_id)
# table(dt4$year)
# 
# predictions = predict(m5$gam, newdata=newdat, se.fit = TRUE) # previously m1 but too circular... 
# 
# # Consolidating new data and predictions
# newdat = cbind(newdat, predictions)
# 
# # If you want CIs 
# newdat <- within(newdat, {
#   lower = fit-1.96*se.fit
#   upper = fit+1.96*se.fit
# })
# 
# # Plot, for example, the predicted outcomes as a function of x1...
# egplot <- ggplot(newdat, aes(x=spring_temp, y=fit)) + 
#   geom_point(data=dt4, aes(x = spring_temp, y = mismatch), colour = "grey") +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)+ 
#   ylab("Mismatch (number of days)")+ xlab("Spring temperature (°C)") + 
#   theme_cowplot()
# egplot


# REVISIONS residuals of spring temp given certain level of green-up

# no difference if take model with linear term for fall_prec 


m10 <- gamm4(mismatch ~ s(spring_temp, k=10) + s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)

newdat <- data.frame(spring_temp=seq(min(dt4$spring_temp), max(dt4$spring_temp), length = 100),
                     fall_prec=mean(dt4$fall_prec, na.rm = T), 
                     mother_id="A44",
                     year="2010")
table(dt4$mother_id)
table(dt4$year)

predictions = predict(m10$gam, newdata=newdat, se.fit = TRUE) # previously m1 but too circular... 

# Consolidating new data and predictions
newdat = cbind(newdat, predictions)

# If you want CIs 
newdat <- within(newdat, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

# Plot, for example, the predicted outcomes as a function of x1...
p1 <- ggplot(newdat, aes(x=spring_temp, y=fit)) + 
  geom_point(data=dt4, aes(x = spring_temp, y = mismatch), colour = "grey") +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)+ 
  ylab("Mismatch (number of days)")+ xlab("Spring temperature (°C)")+
theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p1


newdat <- data.frame(fall_prec=seq(min(dt4$fall_prec), max(dt4$fall_prec), length = 100),
                     spring_temp =mean(dt4$spring_temp, na.rm = T), 
                     mother_id="A44",
                     year="2010")
table(dt4$mother_id)
table(dt4$year)

predictions = predict(m10$gam, newdata=newdat, se.fit = TRUE) # previously m1 but too circular... 

# Consolidating new data and predictions
newdat = cbind(newdat, predictions)

# If you want CIs 
newdat <- within(newdat, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

# Plot, for example, the predicted outcomes as a function of x1...
p2 <- ggplot(newdat, aes(x=fall_prec, y=fit)) + 
  geom_point(data=dt4, aes(x = fall_prec, y = mismatch), colour = "grey") +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)+ 
  ylab("Mismatch (number of days)")+ xlab("Autumn precipitation (mm)")+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
p2



# panel 
library(cowplot)
FIG2_mismatch= cowplot::plot_grid(p1,p2,
                                  ncol=2,
                                  align = "vh",
                                  labels = c("a)", "b)"),
                                  rel_widths = c(1,1))
ggsave("Graphs/revisions2/FIG2_mismatch_revised.png", width = 150, height = 100, units="mm")

# ggsave("FIGS1_detrendedCorr.png")
# ggsave("FIGS1_detrendedCorr.tiff")
# #
cowplot::save_plot("Graphs/revisions2/FIG2_mismatch_revised_2.png", FIG2_mismatch,
                   ncol = 2, #
                   nrow = 1, #
                   # each individual subplot should have an aspect ratio of 1.3
                   base_aspect_ratio = 1)


# # for SCEE - background only 
# p_spring<- ggplot(newdat, aes(y=fit, x = spring_temp)) +
#  # geom_ribbon(aes(ymin = min, ymax = max),  fill="white", alpha=0.2) +
#   #geom_line(aes(y = predi), size=1.2, color = "White") +
#   geom_point(data=dt4, aes(x = spring_temp, y = mismatch), colour = "transparent") +
#   ylab("Mismatch (number of days)")+ xlab("Spring temperature (°C)") + 
#   #scale_y_continuous(limits = c(0, 120),breaks=seq(0, 120,40)) +
#   #scale_x_continuous(limits = c(-2, 4),breaks=seq(-2, 4,1)) +
#   # theme_pander(12) +
#   # theme(legend.title=element_blank(), legend.position = c(0.3, 0.4)) +
#   theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), legend.position = "none", 
#         axis.text=element_text(size=12, color="white"),
#         axis.title=element_text(size=14, color="white"),
#         axis.line = element_line(color="white"),
#         axis.ticks = element_line(color = "white"),
#         rect = element_rect(fill = "white"),
#         panel.grid.major = element_line(color="white"), 
#         panel.grid.minor = element_line(color="transparent"),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA)) 
# 
# ggsave(p_spring, filename = "spring_temp_mis_bg.png",  bg = "transparent", width = 150, height = 100, units="mm")
# 
# # for SCEE 
# p_spring<- ggplot(newdat, aes(y=fit, x = spring_temp)) +
#   geom_point(data=dt4, aes(x = spring_temp, y = mismatch), alpha = .5, colour = "white") +
#   geom_line(size=1.2, color = "White") +
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill="white", alpha=0.2)+ 
#   ylab("Mismatch (number of days)")+ xlab("Spring temperature (°C)") + 
#   # scale_y_continuous(limits = c(10, 40),breaks=seq(10, 40,10)) +
#   #scale_x_continuous(limits = c(-2, 4),breaks=seq(-2, 4,1)) +
#   # theme_pander(12) +
#   # theme(legend.title=element_blank(), legend.position = c(0.3, 0.4)) +
#   theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), legend.position = "none", 
#         axis.text=element_text(size=12, color="white"),
#         axis.title=element_text(size=14, color="white"),
#         axis.line = element_line(color="white"),
#         axis.ticks = element_line(color = "white"),
#         rect = element_rect(fill = "white"),
#         panel.grid.major = element_line(color="white"), 
#         panel.grid.minor = element_line(color="transparent"),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA)) 
# 
# ggsave(p_spring, filename = "spring_temp_mis.png",  bg = "transparent", width = 150, height = 100, units="mm")
# 
# # fall -  background only 
# p_fall<- ggplot(newdat, aes(x=fall_prec, y=fit)) +
#  # geom_line() +
#   #geom_ribbon(aes(ymin = lower, ymax = upper), fill="white", alpha=0.2)+ 
#   geom_point(data=dt4, aes(x = fall_prec, y = mismatch), alpha = 0.2, colour = "transparent") +
#   ylab("Mismatch (number of days)")+ xlab("Fall precipitation (mm)") +
# #scale_y_continuous(limits = c(0, 120),breaks=seq(0, 120,40)) +
#   #scale_x_continuous(limits = c(-2, 4),breaks=seq(-2, 4,1)) +
#   # theme_pander(12) +
#   # theme(legend.title=element_blank(), legend.position = c(0.3, 0.4)) +
#   theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), legend.position = "none", 
#         axis.text=element_text(size=12, color="white"),
#         axis.title=element_text(size=14, color="white"),
#         axis.line = element_line(color="white"),
#         axis.ticks = element_line(color = "white"),
#         rect = element_rect(fill = "white"),
#         panel.grid.major = element_line(color="white"), 
#         panel.grid.minor = element_line(color="transparent"),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA)) 
# 
# ggsave(p_fall, filename = "fall_prec_mis_bg.png",  bg = "transparent", width = 150, height = 100, units="mm")
# 
# p_fall<- ggplot(newdat, aes(x=fall_prec, y=fit)) +
#   geom_line(size=1.2, color = "White") +
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill="white", alpha=0.2)+ 
#   geom_point(data=dt4, aes(x = fall_prec, y = mismatch), alpha = .5, colour = "white") +
#   ylab("Mismatch (number of days)")+ xlab("Fall precipitation (mm)") +
#   #scale_y_continuous(limits = c(0, 120),breaks=seq(0, 120,40)) +
#   #scale_x_continuous(limits = c(-2, 4),breaks=seq(-2, 4,1)) +
#   # theme_pander(12) +
#   # theme(legend.title=element_blank(), legend.position = c(0.3, 0.4)) +
#   theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), legend.position = "none", 
#         axis.text=element_text(size=12, color="white"),
#         axis.title=element_text(size=14, color="white"),
#         axis.line = element_line(color="white"),
#         axis.ticks = element_line(color = "white"),
#         rect = element_rect(fill = "white"),
#         panel.grid.major = element_line(color="white"), 
#         panel.grid.minor = element_line(color="transparent"),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA)) 
# 
# ggsave(p_fall, filename = "fall_prec_mis.png",  bg = "transparent", width = 150, height = 100, units="mm")



# newdat <- data.frame(spring_temp=seq(min(dt4$spring_temp), max(dt4$spring_temp), length = 100), 
#                      mother_id="A44", 
#                      year="2010")
# table(dt4$mother_id)
# table(dt4$year)

# Getting predicted outcomes for new data
# # These include the splines but ignore other REs
# predictions = predict(m5$gam, newdata=newdat, se.fit = TRUE) # previously m1 but too circular... 
# 
# # Consolidating new data and predictions
# newdat = cbind(newdat, predictions)
# 
# # If you want CIs 
# newdat <- within(newdat, {
#   lower = fit-1.96*se.fit
#   upper = fit+1.96*se.fit
# })
# 
# # Plot, for example, the predicted outcomes as a function of x1...
# egplot <- ggplot(newdat, aes(x=spring_temp, y=fit)) + 
#   geom_point(data=dt4, aes(x = spring_temp, y = mismatch), colour = "grey") +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)+ 
#   ylab("Mismatch (number of days)")+ xlab("Spring temperature (°C)") + 
#   theme_cowplot()
# egplot
# ggsave("Graphs/FIG2_mismatchSpringTemp.png", width = 140, height = 140, units = "mm")


# export fig2 -------------------------------------------------------------

# FIG2_mismatch= cowplot::plot_grid(yr,egplot,
#                                   ncol=2, 
#                                   align = "vh",
#                                   labels = c("a)", "b)"),
#                                   rel_widths = c(1,1))
# 
# # ggsave("FIGS1_detrendedCorr.png")
# # ggsave("FIGS1_detrendedCorr.tiff")
# # # 
# cowplot::save_plot("Graphs/FIG2_mismatch.png", FIG2_mismatch,
#                    ncol = 2, #
#                    nrow = 1, #
#                    # each individual subplot should have an aspect ratio of 1.3
#                    base_aspect_ratio = 1)
# 


# # verification : remove extreme high values  or late birthdate 
tmp <- dt4 %>%filter(birthdate<180)

tmp <- na.omit(tmp)
m1 <- gamm4(mismatch ~ s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = tmp)
m2 <- gamm4(mismatch ~ s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = tmp)
m3 <- gamm4(mismatch ~ s(evi_tm1, k=10), random = ~ (1|mom_id) + (1|year), data = tmp)
m4 <- gamm4(mismatch ~ s(spring_prec, k=10), random = ~ (1|mom_id) + (1|year), data = tmp)
m5 <- gamm4(mismatch ~ s(spring_temp, k=10), random = ~ (1|mom_id) + (1|year), data = tmp)

AICctab(m1$mer, m2$mer, m3$mer, m4$mer, m5$mer)
# dAICc df
# m5$mer  0.0  6 
# m2$mer  8.0  6 
# m3$mer  8.7  6 
# m4$mer 10.4  6 
# m1$mer 10.6  6 


par(mfrow=c(2,2))
    
plot(m5$gam, page =1)
summary(m5$gam) # this is closer to a linear relationship 1.232  1.232 12.03 0.000213 ***
summary(m5$mer)

# end verification 




# correlations on detrended variables -------------------------------------


# simple correlations are based on unique values per year, and for the case of birthdate this is median birthdates 
# the fall weather conditions were 'windowed' based on plasticity analyses, and reduced dataset just taken 

mod.sel <- list()
mod.sel$mod1 <- lm(bd.r~fP.r, data = detrended.fall.pheno)
mod.sel$mod2 <- lm(bd.r~fT.r, data = detrended.fall.pheno) # 
mod.sel$mod3 <- lm(bd.r~fEVI.r, data = detrended.fall.pheno)

AICtab(mod.sel) # very few points 
# dAIC df
# mod1  0.0 3 
# mod2  5.8 3 
# mod3 10.8 3 
summary(mod.sel$mod1) # MAR fP.r         1.172e+01  3.119e+00   3.756  0.00191 **
summary(mod.sel$mod2) # Y fT.r        -5.207e+00  2.185e+00  -2.383   0.0308 *
summary(mod.sel$mod3) # N fEVI.r       -1.066e-01  1.708e-01  -0.624    0.542


confint(mod.sel$mod1)# fP.r    5.067662 18.365147
confint(mod.sel$mod2) #fT.r          -9.863935 -0.5498044
confint(mod.sel$mod3)#fEVI.r    -0.4705747 0.2574405


# correlate raw variables  -----------------------------------------------------

colnames(detrended.fall.pheno)


mod.raw <- list()
mod.raw$mod1 <- lm(medBD~fall_prec, data = detrended.fall.pheno)
mod.raw$mod2 <- lm(medBD~fall_temp, data = detrended.fall.pheno) # 
mod.raw$mod3 <- lm(medBD~evi_tm1, data = detrended.fall.pheno)

AICtab(mod.raw) # very few points 
# dAIC df
# mod1 0.0  3 
# mod2 1.1  3 
# mod3 9.5  3
summary(mod.raw$mod1) # Y
summary(mod.raw$mod2) # Y
summary(mod.raw$mod3) # N


confint(mod.raw$mod1)
#fall_prec     4.860142  20.23597

confint(mod.raw$mod2)
#fall_temp    -9.871136  -2.027796

confint(mod.raw$mod3)
#evi_tm1     -0.5380959   0.2751294


# prediction figures raw variables  ---------------------------------------


# mod1
range(detrended.fall.pheno$fall_prec)
xnew <- data.frame(fall_prec=seq(0, 1.5, 0.5))
fitmod1 <- predict(mod.raw$mod1, xnew, type="response") 
new = cbind(xnew, fitmod1)

# now compare with raw values and put in appendix 
p.raw= ggplot(data = detrended.fall.pheno, aes(x = fall_prec,y=medBD))+
  geom_point(size =3, color = "grey") + 
  geom_line(data=new, aes(x=fall_prec, y=fitmod1),size = 1, linetype = "solid") + 
  labs(x = "Precipitation (mm)", y = "Median parturition (Julian day)")+
  cowplot::theme_cowplot()+ theme(plot.margin = unit(c(0.5,0.5,1,1), "cm"))

# 2
range(detrended.fall.pheno$fall_temp)
xnew <- data.frame(fall_temp=seq(2.5, 5.5, 0.5))
fitmod2 <- predict(mod.raw$mod2, xnew, type="response") 
new = cbind(xnew, fitmod2)

# now compare with raw values and put in appendix 
t.raw = ggplot(data = detrended.fall.pheno, aes(x = fall_temp,y=medBD))+
  geom_point(size =3, color = "grey") + 
  geom_line(data=new, aes(x=fall_temp, y=fitmod2),size = 1) + 
  labs(x = "Temperature (°C)", y = "Median parturition (Julian day)")+
  cowplot::theme_cowplot()+ theme(plot.margin = unit(c(0.5,0.5,1,1), "cm"))

# evi
range(detrended.fall.pheno$evi_tm1)
xnew <- data.frame(evi_tm1=seq(260, 295, 5))
fitmod3 <- predict(mod.raw$mod3, xnew, type="response") 
new = cbind(xnew, fitmod3)
summary(mod.raw$mod3) # marginal, -0.1

# now compare with raw values and put in appendix 
evi.raw=ggplot(data = detrended.fall.pheno, aes(x = evi_tm1,y=medBD))+
  geom_point(size =3, color = "grey") + 
  geom_line(data=new, aes(x=evi_tm1, y=fitmod3),size = 1, linetype = "dashed") + 
  labs(x = "Green-down (Julian day)", y = "Median parturition (Julian day)")+
  cowplot::theme_cowplot()+ theme(plot.margin = unit(c(0.5,0.5,1,1), "cm"))

# prediction figures detrended variables -------------------------------------------------------------


# mod1
range(detrended.fall.pheno$fP.r)
xnew <- data.frame(fP.r=seq(-0.5, 1, 0.5))
fitmod1 <- predict(mod.sel$mod1, xnew, type="response") 
new = cbind(xnew, fitmod1)

# now compare with raw values and put in appendix 
p.res= ggplot(data = detrended.fall.pheno, aes(x = fP.r,y=bd.r))+
  geom_point(size =3, color = "grey") + 
  geom_line(data=new, aes(x=fP.r, y=fitmod1),size = 1, linetype = "solid") + 
  labs(x = "Precipitation residuals", y = "Median parturition residuals")+
  cowplot::theme_cowplot()+ theme(plot.margin = unit(c(0.5,0.5,1,1), "cm"))

# 2
range(detrended.fall.pheno$fT.r)
xnew <- data.frame(fT.r=seq(-1, 1.5, 0.5))
fitmod2 <- predict(mod.sel$mod2, xnew, type="response") 
new = cbind(xnew, fitmod2)
summary(mod.sel$mod2)

# now compare with raw values and put in appendix 
t.res = ggplot(data = detrended.fall.pheno, aes(x = fT.r,y=bd.r))+
  geom_point(size =3, color = "grey") + 
  geom_line(data=new, aes(x=fT.r, y=fitmod2),size = 1) + 
  labs(x = "Temperature residuals", y = "Median parturition residuals")+
cowplot::theme_cowplot()+ theme(plot.margin = unit(c(0.5,0.5,1,1), "cm"))

# evi
range(detrended.fall.pheno$fEVI.r)
xnew <- data.frame(fEVI.r=seq(-15, 15, 1))
fitmod3 <- predict(mod.sel$mod3, xnew, type="response") 
new = cbind(xnew, fitmod3)
summary(mod.sel$mod3) # marginal, -0.1

# now compare with raw values and put in appendix 
evi.r=ggplot(data = detrended.fall.pheno, aes(x = fEVI.r,y=bd.r))+
  geom_point(size =3, color = "grey") + 
  geom_line(data=new, aes(x=fEVI.r, y=fitmod3),size = 1, linetype = "dashed") + 
  labs(x = "Green-down residuals", y = "Median parturition residuals")+
cowplot::theme_cowplot()+ theme(plot.margin = unit(c(0.5,0.5,1,1), "cm"))



# export fig1 -------------------------------------------------------------

FIG1_detrendedCorr= cowplot::plot_grid(p.raw, p.res,t.raw, t.res, evi.raw, evi.r,ncol=2, align = "vh", 
                   labels = c("a)", "", "b)", "", "c)", ""),
                   rel_widths = c(1,1))

# cowplot::save_plot("Graphs/FIG1_detrendedCorr20210421.png", FIG1_detrendedCorr,
#           ncol = 2, #
#           nrow = 3, #
#           # each individual subplot should have an aspect ratio of 1.3
#           base_aspect_ratio = 1)


rm(p.res, t.res,evi.r, ndvi.r,snow.r, lai.r, gpp.r,  pss.r,fpar.r,
   prec.raw,temp.raw, evi.raw, ndvi.raw,snow.raw, lai.raw, gpp.raw,pss.raw,fpar.raw)
rm(bdGAM, eviGAM, eviLM, ndviGAM, ndviLM, new, pGAM, pLM, snowGAM, snowLM, tLM, xnew)





# year effect 
# myear <- gamm4(mismatch ~ s(as.numeric(as.character(year)), k = 10), random = ~ (1|mom_id), data = dt4) #statistically significant, in the sense that the fitted smoother is different from a null model.
# summary(myear$gam) # R-sq.(adj) =  0.355   
# 
# # k'  edf k-index p-value    
# # s(as.numeric(as.character(year))) 8.211  8.211 16.31  <2e-16 ***
# 
# summary(myear$mer)$coef # useless since Non param
# 
# # Estimate Std. Error   t value
# # X(Intercept)                            7.192521  0.9435312  7.622982
# # Xs(as.numeric(as.character(year)))Fx1 -58.146268 16.1220294 -3.606635
# 
# # checkup model 
# plot(myear$gam, page =1)
# gam.check(myear$gam)
# # end checkup 
# 
# 
# # predict on year 
# xnew <- data.frame(year=seq(2000, 2017, 1))
# fit <- predict(myear$gam, xnew, type="response", se.fit=T, re.form=NA)
# new = cbind(xnew, fit)
# new <- within(new, {
#   lower = fit-1.96*se.fit
#   upper = fit+1.96*se.fit
# })
# 
# yr <- ggplot(new, aes(x=year, y=fit)) + 
#   geom_point(data=dt4, aes(x = year, y = mismatch), colour = "grey") +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) + 
#   ylab("Mismatch (number of days)") + xlab("Year") 
# yr
# 
