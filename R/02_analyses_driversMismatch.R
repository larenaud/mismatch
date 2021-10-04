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
library(kableExtra)

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

# load clean data
trends_df <- read.csv2("data/mine/trends_data_noID.csv")

# model sel. What are the predictors of parturition date and mismatch ----------------------------
m1.g <- gamm4(birthdate ~ s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m2.g <- gamm4(birthdate ~ s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m3.g <- gamm4(birthdate ~ s(evi_tm1, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
#m4.g <- gamm4(birthdate ~ s(spring_prec, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
#m5.g <- gamm4(birthdate ~ s(spring_temp, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)

AICtab(m1.g$mer, m2.g$mer, m3.g$mer)
# dAIC df
# m2.g$mer 0.0  6 
# m3.g$mer 8.2  6 
# m1.g$mer 9.0  6 

t=AICtab(m1.g$mer, m2.g$mer, m3.g$mer)

# t=kable(t, digits = 3) %>%
#   kable_styling(font_size = 10) %>% 
#   kableExtra::kable_styling("bordered") %>% 
#   save_kable(file = "output/table/APPENDIX2_Table_S1_modSelBD.html", self_contained = T)

plot(m2.g$gam) # this is linear
summary(m2.g$gam) # this is linear  1      1 18.91 2.01e-05 ***

# for supp
plot(m1.g$gam) 
summary(m1.g$gam) 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(fall_temp) 1.274  1.274 3.568  0.0762 .
# -


# linear estimates 
m1.l <- lmer(birthdate ~ fall_temp + (1|mom_id) + (1|year), data = trends_df)
m2.l <- lmer(birthdate ~ fall_prec + (1|mom_id) + (1|year), data = trends_df)
m3.l <- lmer(birthdate ~ evi_tm1 + (1|mom_id) + (1|year), data = trends_df)
# m4.g <- lmer(birthdate ~ spring_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
# m5.g <- lmer(birthdate ~ pring_temp, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)

AICtab(m1.l, m2.l, m3.l)
# dAIC df
# m2.l  0.0 5 
# m1.l 10.5 5 
# m3.l 16.3 5 

# effect of fall prec
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

MuMIn::r.squaredGLMM(m2.l) #    R2m      R2c [1,] 0.09478338 0.468289


# effect of fall temp
summary(m1.l)
# fall_temp     -4.528      2.355  -1.922

round(confint(m1.l),3)
# 2.5 %  97.5 %
#   .sig01        3.423  11.710
# .sig02        2.923   8.005
# .sigma        9.961  12.794
# (Intercept) 154.413 195.018
# fall_temp    -9.223   0.107


# mismatch # REVISED 
m1 <- gamm4(mismatch ~ s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m2 <- gamm4(mismatch ~ s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m3 <- gamm4(mismatch ~ s(evi_tm1, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m4 <- gamm4(mismatch ~ s(spring_prec, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m5 <- gamm4(mismatch ~ s(spring_temp, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m6 <- gamm4(mismatch ~ s(spring_temp, k=10) + s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m7 <- gamm4(mismatch ~ s(spring_prec, k=10) + s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m8 <- gamm4(mismatch ~ s(spring_prec, k=10) + s(evi_tm1, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m9 <- gamm4(mismatch ~ s(spring_temp, k=10) + s(evi_tm1, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m10 <- gamm4(mismatch ~ s(spring_temp, k=10) + s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m11 <- gamm4(mismatch ~ s(spring_prec, k=10) + s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m12 <- gamm4(mismatch ~ s(evi_tm1, k=10) + s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)
m13 <- gamm4(mismatch ~ s(evi_tm1, k=10) + s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = trends_df)

AICtab(m1$mer, m2$mer, m3$mer,m4$mer, m5$mer, m6$mer, m7$mer, m8$mer, m9$mer, m10$mer, m11$mer, m12$mer, m13$mer)
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

# make some checkup
gam.check(m10$gam)
plot(m10$gam, pages=1)
#s(spring_temp) 9.00 1.97    0.63  <2e-16 ***
gam.check(m10_l$gam)
# end check

summary(m10$mer)
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


summary(m10$gam) # R-sq.(adj) =  0.406   
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
# s(spring_temp) 2.447  2.447 15.99 2.18e-07 ***
#   s(fall_prec)   1.000  1.000 12.20 0.000564 ***
#   ---


# to get parametric coefficient of linear term fall prec:
m10_l <- gamm4(mismatch ~ s(spring_temp, k=10) + fall_prec, random = ~ (1|mom_id) + (1|year), data = trends_df)
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

round(confint(m10_l$mer), 3)
# 2.5 % 97.5 %
#   .sig01              0.205 11.021
# .sig02              2.299  7.870
# .sig03              0.000 15.768
# .sigma             10.091 13.222
# X(Intercept)       -5.917  6.414
# Xfall_prec          7.025 25.436
# Xs(spring_temp)Fx1  1.263 18.599
round(confint(m10_l$gam), 3)

#t=AICtab(m1$mer, m2$mer, m3$mer, m4$mer, m5$mer)
t2=AICtab(m1$mer, m2$mer, m3$mer,m4$mer, m5$mer, m6$mer, m7$mer, m8$mer, m9$mer, m10$mer, m11$mer, m12$mer, m13$mer)
# t2=kable(t2, digits =3) %>%
#   kable_styling(font_size = 10) %>%
#   kable_styling("bordered") %>%
#   save_kable(file = "output/table/APPENDIX2_Table_S2_modSelMismatch.html", self_contained = T)

# summary(m5$gam) # s(spring_temp) 1.972  1.972 8.43 0.000215 ***
# # check 
# plot(m5$gam, pages=1)
# gam.check(m5$gam)
# #s(spring_temp) 9.00 1.97    0.63  <2e-16 ***
# # end check 

# save models for figure S2 A2
# save(list = ls(), file = "cache/modSel_mismatch.RData")

# variance partitioning with R2part---------------------------------------------------
# get estimates and variance - need mer only - MANUALLY
m10 <- lmer(mismatch ~ spring_temp + fall_prec + (1|mom_id) + (1|year), data = trends_df)
summary(m10)
r.squaredGLMM(m10) 
# R2m       R2c
#0.2741147 0.5869986
m10a <- lmer(mismatch ~ fall_prec + (1|mom_id) + (1|year), data = trends_df) # rm : 0.05998344
r.squaredGLMM(m10a) 
#0.2741147 -0.05998344 : 0.2141313

m10b <- lmer(mismatch ~ spring_temp + (1|mom_id) + (1|year), data = trends_df)
summary(m10b)
r.squaredGLMM(m10b) # 0.1847434
#0.2741147 - 0.1847434=0.0893713

# venn diagram - LIBRARY R2PART
library(partR2)
m10 <- lmer(mismatch ~ spring_temp + fall_prec + (1|mom_id) + (1|year), data = trends_df)
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

summary(res, round_to = 3)

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

