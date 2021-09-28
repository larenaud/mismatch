# script modified 2020-01-21 by L. Renaud to reproduce analyses in: 
# Renaud et al. GCB , section 'Temporal trends in autumn and spring phenology' + supplementary results in Appendix 1

library (tidyverse)
library(lme4)
library (ggplot2)
library(xtable)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(pander)
library(cowplot)
library(forcats)
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
library(ggpubr)
library(mgcv)
library(gamm4)


getwd()
rm(list = ls ())


# ggplot custom  ----------------------------------------------------------
theme_set(theme_bw())
scale_colour_discrete <- function(...,palette="Set1") {
    scale_colour_brewer(...,palette=palette)
}
scale_colour_orig <- ggplot2::scale_colour_discrete
scale_fill_discrete <- function(...,palette="Set1") {
    scale_fill_brewer(...,palette=palette)
}


# load df as .csv ------------------------------------------------------------------
dat.trend <- read.csv2("data/mine/trends_data.csv")

# create dat.trend.yr for analyses at the annual scale --------------------
dat.trend.yr <- dat.trend[, c("year", "spring_temp", "fall_temp", "evi_up",  "evi_tm1","birthdate")] %>%
    group_by(year) %>%
    summarise_all(mean)


# fit models of temporal trends  ----------------------------------------------
# spring 
mod.ts.l <- lm(spring_temp~year,data=unique(dat.trend[, c("year", "spring_temp")]) )
mod.ts.gam <- gam(spring_temp~s(year), gamma = 1.4, data=unique(dat.trend[, c("year", "spring_temp")]) )
summary(mod.ts.gam)
# parametric coefficients : (Intercept)   4.3396     0.3197   13.57 7.87e-10 ***
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(year)   1      1 2.308   0.149



# without penalty 
mod.ts.gam <- gam(spring_temp~s(year), data=unique(dat.trend[, c("year", "spring_temp")]) ) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   4.3396     0.2136   20.32 9.82e-09 ***
#   edf Ref.df     F p-value  
# s(year) 7.146  8.185 3.226  0.0502 .
summary(mod.ts.gam)
summary(mod.ts.l) # year: 0.09913    0.06525   1.519    0.149
round(confint(mod.ts.l),3) # year: -0.040  0.238




# autumn
mod.tf.l <- lm(fall_temp~year,data=unique(dat.trend[, c("year", "fall_temp")]) )
mod.tf.gam <- gam(fall_temp~s(year),gamma = 1.4, data=unique(dat.trend[, c("year", "fall_temp")]) ) # Woods 
summary(mod.tf.gam)
# Parametric coefficients:  #4.3023     0.1397    30.8 5.62e-15 ***       
#   s(year)   1      1 6.682  0.0207 *
summary(mod.tf.l) # 0.07369    0.02851   2.585   0.0207 *
round(confint(mod.tf.l),3)  #0.013   0.134


# supp - removed penalty gamma 1.4
mod.ts.gam <- gam(spring_temp~s(year), data=unique(dat.trend[, c("year", "spring_temp")]))
mod.tf.gam <- gam(fall_temp~s(year),data=unique(dat.trend[, c("year", "fall_temp")])) 
summary(mod.ts.gam)
# param est: 4.3396, SE, 0.2136, p: 9.82e-09 ***, s(year) 7.146  8.185 3.226  0.0502 .
summary(mod.tf.gam)
# param: 4.3023, se: 0.1397,p: 5.62e-15, *** edf 1      1 6.682  0.0207 *




# green-up 
mod.gu.l <- lm(evi_up~year,data=unique(dat.trend[, c("year", "evi_up")]) )
mod.gu.gam <- gam(evi_up~s(year),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_up")]))  
#mod.gu.gam2 <- gam(evi_up~s(spring_temp),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_up", "spring_temp")])) 
#mod.gu.lm2 <- lm(evi_up~spring_temp,data=unique(dat.trend[, c("year", "evi_up", "spring_temp")]))  

summary(mod.gu.l) # year          -1.0049     0.4107  -2.447   0.0272 *
round(confint(mod.gu.l), 3) # year         -1.880   -0.130
summary(mod.gu.gam) #s(year)   1      1 5.987  0.0272 *
#summary(mod.gu.gam2) ##s(spring_temp) 1.56  1.934 18.26 0.000225 ***
#summary(mod.gu.lm2) 
#spring_temp   -5.588      1.058   -5.28 9.26e-05 ***
round(confint(mod.gu.lm2),3) # spring_temp  -7.844  -3.332




# useful later 
mod.gu.gam <- gam(evi_up~s(year),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_up")]) ) # Woods 
mod.gu.gam2 <- gam(evi_up~s(spring_temp),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_up", "spring_temp")]) ) # Woods 
mod.gu.lm2 <- lm(evi_up~spring_temp,data=unique(dat.trend[, c("year", "evi_up", "spring_temp")]) ) # Woods 

dat.trend.yr$pred.gu <- predict(mod.gu.gam2, newdata = dat.trend.yr) # this is to make figure 1 version 2




# revisions - added brown down
mod.gd.l <- lm(evi_tm1~year,data=unique(dat.trend[, c("year", "evi_tm1")]) )
mod.gd.gam <- gam(evi_tm1~s(year),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_tm1")]) ) 
mod.gd.gam2 <- gam(evi_tm1~s(year),data=unique(dat.trend[, c("year", "evi_tm1")]) ) 
#mod.gd.gam2 <- gam(evi_tm1~s(spring_temp),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_tm1", "spring_temp")]) ) # Woods 
#mod.gd.lm2 <- lm(evi_tm1~spring_temp,data=unique(dat.trend[, c("year", "evi_tm1", "spring_temp")]) ) # Woods 

summary(mod.gd.l) # year      year          0.1225     0.4228   0.290    0.776
round(confint(mod.gd.l), 3) # year         -0.779    1.024

summary(mod.gd.gam) #  s(year)   1      1 0.084   0.776
summary(mod.gd.gam2) # s(year)   1      1 0.084   0.776



# parturition date 
mod.bd.gam<- gamm4(birthdate ~ s(year), random=~(1|mom_id), REML=T, data=dat.trend)
mod.bd.lmer<- lmer(birthdate ~ year + (1|mom_id), data=dat.trend)

summary(mod.bd.gam$gam) # s(year) 5.898  5.898 7.466 5.33e-07 ***
summary(mod.bd.lmer) # year          -0.5932     0.2145  -2.765
round(confint(mod.bd.lmer2),3) #fall_temp    -9.223   0.107

# quick check up
plot(mod.bd.gam2$gam)
plot(mod.bd.gam$gam)

# useful later for supp, Appendix 2 Figure S3
mod.bd.gam2 <- gamm4(birthdate~s(fall_temp), random=~(1|mom_id) + (1|year), REML=T, data=dat.trend)
mod.bd.lmer2 <- lmer(birthdate~fall_temp +(1|mom_id) + (1|year), REML=T, data=dat.trend)

pred=predict(mod.bd.gam2$gam, newdata =dat.trend.yr, se.fit=T)
dat.trend.yr=dat.trend.yr %>% mutate(pred.bd=pred$fit, ci.h=pred$se.fit*1.96+pred$fit, ci.l=-pred$se.fit*1.96+pred$fit) 


# mismatch
mod.mis.gam2 <- gamm4(mismatch~s(year), random=~(1|mom_id), REML=T, data=dat.trend)
plot(mod.mis.gam2$gam)
summary(mod.mis.gam2$gam)

mod.mis.l <- lmer(mismatch~year +I(year^2) + (1|mom_id), REML=T, data=dat.trend)# convergence issue
summary(mod.mis.l)

mod.mis.l2 <- lmer(mismatch~year+ (1|mom_id), REML=T, data=dat.trend)
summary(mod.mis.l2)


# predicted changes in conception and green-up date - correlations --------
cor.test(dat.trend.yr$pred.gu, dat.trend.yr$pred.bd)
# Pearson's product-moment correlation
# 
# data:  dat.trend.yr$pred.gu and dat.trend.yr$pred.bd
# t = -0.30516, df = 15, p-value = 0.7644
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5388492  0.4178730
# sample estimates:
#       cor 
# -0.078548 

summary(lm(pred.bd~pred.gu,data=dat.trend.yr)) # pred.gu      -0.03104    0.10172  -0.305    0.764  
confint(lm(pred.bd~pred.gu,data=dat.trend.yr)) # pred.gu      -0.247861   0.1857774




