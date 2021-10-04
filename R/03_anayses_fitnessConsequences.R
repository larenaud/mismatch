
# pkgs_CRAN <- c("lme4","MCMCglmm","blme",
#                "pbkrtest","coda","aods3","bbmle","ggplot2",
#                "reshape2","plyr","numDeriv","Hmisc",
#                "plotMCMC","gridExtra","R2admb",
#                "broom.mixed","dotwhisker")
# install.packages(pkgs_CRAN)
# rr <- "http://www.math.mcmaster.ca/bolker/R"
# 
# install.packages("glmmADMB",type="source",repos=rr)
# install.packages("devtools")
library("devtools")
library("lme4")
library("tidyverse")
library("glmmADMB")      ## (not on CRAN)
library("glmmTMB")
library("MCMCglmm")
library("blme")
library("bbmle")     ## AICtab
library("MASS")          ## for glmmPQL (base R)
## auxiliary
library("ggplot2")
## to squash facets together ...
zmargin <- theme(panel.spacing=grid::unit(0,"lines"))
library("gridExtra")     ## for grid.arrange()
library("broom.mixed")
## n.b. as of 25 Sep 2018, need bbolker github version of dotwhisker ...
library("dotwhisker")
library("coda")      ## MCMC diagnostics
library("aods3")     ## overdispersion diagnostics
library("plotMCMC") ## pretty plots from MCMC fits

library("pbkrtest")  ## parametric bootstrap
library("Hmisc")
## for general-purpose reshaping and data manipulation:
library("reshape2")

## for illustrating effects of observation-level variance in binary data:
library("numDeriv")
library(kableExtra)

sessionInfo()
packageVersion("bbmle") # ‘1.0.23.1’

# clean up 
rm(list = ls())

# load data
fitness.data <- read.csv2("data/mine/fitness_data_noID.csv")
# load("Cache/20210426fitnessData_centeredDetrended.RData")


# retidy in case 
colnames(fitness.data)
fitness.data$year <- as.factor(fitness.data$year)
fitness.data$mom_id <- as.factor(fitness.data$mom_id)
fitness.data$lamb_id <- as.factor(fitness.data$lamb_id)
fitness.data$prs <- as.factor(fitness.data$prs)
fitness.data$newRS <- as.factor(fitness.data$newRS)

# check out distribution of response variable 
fitness.data$mom_id<-droplevels(fitness.data$mom_id) # important because large amount of 0s
fitness.data$lamb_id<-droplevels(fitness.data$lamb_id) # important because large amount of 0s

# remove residuals
data.frame(colnames(fitness.data))

fitness.data$weanMass_z <- scale(fitness.data$weanMass)
fitness.data$birthdate_z <- scale(fitness.data$birthdate)
fitness.data$fem_z <- scale(fitness.data$fem)
fitness.data$fem_tm1_z <- scale(fitness.data$fem_tm1)
fitness.data$fallMass_tm1_z <- scale(fitness.data$fallMass_tm1)
fitness.data$evi_pm_z <- scale(fitness.data$evi_pm)
fitness.data$evi_im_z <- scale(fitness.data$evi_im)
fitness.data$bd_ind_z <- scale(fitness.data$bd_ind)
fitness.data$bd_pop_z <- scale(fitness.data$bd_pop)

fitness.data <- fitness.data %>% 
  mutate_at(vars(contains("_z")), as.numeric)

# visual check 
hist(fitness.data$evi_im_z)
hist(fitness.data$evi_pm_z)
hist(fitness.data$fallMass_tm1_z)
hist(fitness.data$birthdate_z)
hist(fitness.data$bd_ind_z)
hist(fitness.data$bd_pop_z)

t =table(fitness.data$mom_id, fitness.data$year) # this is for year 2001 to 2017
t=rowSums(t)
hist(t)
mean(t) # 3.5125
range(t) # 1 12
sd(t) #  2.916534
table(fitness.data$year, fitness.data$neonatal) # this is for year 2001 to 2017

#  perfect separation in 2 years
#        0  1
# 2002  0  7
# 2007  0 16
# end check

# sample sizes and data/graph exploration ---------------------------------------------------------
fitness.data$mom_id %>%
  droplevels() %>%
  table %>%
  mean()
# 3.5125

t = fitness.data %>% 
  group_by(mom_id) %>% 
  droplevels () %>%
  summarise(lamb = length(unique(lamb_id)))

mean(t$lamb) # 3.2875


# confirm the design 
with(fitness.data, table(year, mom_id)) # a few ID have 0 everywhere - 

# neonatal model fitting - test random structure -----------------------------------------------------------
str(fitness.data)
colnames(fitness.data)
cat(colnames(fitness.data),sep="\n")

neonatal.data <- fitness.data[, c("year","mom_id","lamb_id", "neonatal","newRS", "prs","sex", "birthdate_z","birthdate", "fem_z","fem_tm1_z","fallMass_tm1_z",
                                  "evi_pm_z","evi_im_z", "bd_ind_z", "bd_pop_z", "evi_pm","evi_im", "bd_ind", "bd_pop")] 
neonatal.data <- neonatal.data[!is.na(neonatal.data$birthdate_z),] #  281
neonatal.data <- neonatal.data[!is.na(neonatal.data$prs),] # makes it drop to 250 
neonatal.data <- neonatal.data[!is.na(neonatal.data$fem_z),] #
neonatal.data <- neonatal.data[!is.na(neonatal.data$evi_pm_z),] #
neonatal.data <- neonatal.data[!is.na(neonatal.data$evi_im_z),] #
neonatal.data <- neonatal.data[!is.na(neonatal.data$bd_ind_z),] #
neonatal.data <- neonatal.data[!is.na(neonatal.data$bd_pop_z),] # n = 250 observation 
neonatal.data <- neonatal.data[!is.na(neonatal.data$evi_pm),] #
neonatal.data <- neonatal.data[!is.na(neonatal.data$evi_im),] #
neonatal.data <- neonatal.data[!is.na(neonatal.data$bd_ind),] #
neonatal.data <- neonatal.data[!is.na(neonatal.data$bd_pop),] # n = 250 observation 
neonatal.data$mom_id <- droplevels(neonatal.data$mom_id)

t =table(neonatal.data$mom_id, neonatal.data$year) # this is for year 2001 to 2017
t=rowSums(t)
hist(t)
mean(t) # 3.623188
range(t) # 1 12
sd(t) # 2.859975

# sample sizes 
neonatal.data$mom_id %>%
  droplevels() %>%
  table %>%
  mean()
# 3.623188
neonatal.data$mom_id %>%
  droplevels() %>%
  nlevels() # 69 females 

# number of lambs 
t = neonatal.data %>% 
  group_by(mom_id) %>% 
  droplevels () %>%
  summarise(lamb = length(unique(lamb_id)))
mean(t$lamb) # 3.376812
range(t$lamb) #  1 11

# test random effect structure 
cmod_blme_Lyrmom_id <- bglmer(neonatal~fallMass_tm1_z+(1|year) + (1|mom_id),data=neonatal.data, 
                        family=binomial,
                        fixef.prior = normal(sd = c(10, 10)))

cmod_blme_Lmom_id <- bglmer(neonatal~fallMass_tm1_z+ (1|mom_id),data=neonatal.data,
                        family=binomial,
                        fixef.prior = normal(sd = c(10, 10)))
cmod_blme_Lyr <- bglmer(neonatal~fallMass_tm1_z+ (1|year),data=neonatal.data, # convergence issues
                        family=binomial,
                        fixef.prior = normal(sd = c(10, 10)))

anova(cmod_blme_Lmom_id,cmod_blme_Lyrmom_id)
#                         npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
# cmod_blme_Lmom_id      3 218.33 228.90 -106.17   212.33                    
# cmod_blme_Lyrmom_id    4 219.72 233.81 -105.86   211.72 0.612  1      0.434

anova(cmod_blme_Lyr,cmod_blme_Lyrmom_id)
#                       par    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
# cmod_blme_Lyr          3 216.35 226.91 -105.17   210.35                    
# cmod_blme_Lyrmom_id    4 219.72 233.81 -105.86   211.72     0  1          1

summary(cmod_blme_Lyrmom_id)
AICtab(cmod_blme_Lmom_id, cmod_blme_Lyr,cmod_blme_Lyrmom_id) # yr only is best 
summary(cmod_blme_Lyr)

# neonatal model selection by AIC ----------------------------------
# change sd value of 10 to less if take detrended values 
# weak prior is one with a standard deviation that is large relative
# to the expected scale of the effect (e.g. we might say sigma=10 is
#                                      large, but it won't be if the units of measurement are very small
# so that a typical value of the mean is 100,000 ...

# # first select best variables to include as control 
# base.neo<-list()
# base.neo$bm1 <- glmer(neonatal ~1 + (1|yr), data = neonatal.data, family=binomial) 
# base.neo$bm2 <- glmer(neonatal ~prs + (1|yr), data = neonatal.data, family=binomial) 
# base.neo$bm3 <- glmer(neonatal ~fem_z + (1|yr), data = neonatal.data, family=binomial) 
# base.neo$bm4 <- glmer(neonatal ~fem_tm1_z + (1|yr), data = neonatal.data, family=binomial)
# #base.neo$bm5 <- glmer(neonatal ~pred + (1|yr), data = neonatal.data, family=binomial) 
# base.neo$bm6 <- bglmer(neonatal ~pred_tm1 + (1|yr), data = neonatal.data, family=binomial, fixef.prior = normal) # singular fit in glmer 
# base.neo$bm7 <- glmer(neonatal ~prs+fem_z + (1|yr), data = neonatal.data, family=binomial)
# #$bm8 <- glmer(neonatal ~fem_z+ pred + (1|yr), data = neonatal.data, family=binomial) 
# base.neo$bm9 <- bglmer(neonatal ~prs+ pred_tm1 + (1|yr), data = neonatal.data, family=binomial, fixef.prior = normal) 

bbmle::AICctab(base.neo)
# dAICc df
# bm9  0.0  5 
# bm6  5.4  3 
# bm2  6.8  4 
# bm7  7.6  5 
# bm1  9.8  2 
# bm4 10.7  3 
# bm3 10.9  3 

# When specifying standard deviations, a vector of length less than the number of fixed effects will have its tail repeated, 
# while the first element is assumed to apply only to the intercept term. So in the default of ‘c(10, 2.5)’, the intercept receives a standard deviation of 10 and the various slopes are all given a standard deviation of 2.5


# # done with Original variables ONLY AND SCALED 
# colnames(neonatal.data)
# mod.neonatal<-list()
# mod.neonatal$m1 <- glmer(neonatal ~1 + (1|yr), data = neonatal.data, family=binomial)  # this one in glmer only otherwise does not converge
# mod.neonatal$m2 <- bglmer(neonatal ~prs+ pred_tm1 + (1|yr), data = neonatal.data, family=binomial, fixef.prior = normal)  # this one in glmer only otherwise does not converge
# mod.neonatal$m3 <- bglmer(neonatal ~ prs+ pred_tm1 +evi_im_z + evi_pm_z + (1|yr), data = neonatal.data, family=binomial, fixef.prior = normal)
# mod.neonatal$m4<- bglmer(neonatal ~ prs+ pred_tm1 +evi_im_z + I(evi_im_z^2)+ evi_pm_z + (1|yr), data = neonatal.data, family=binomial, fixef.prior = normal) #(sd = c(10,2.5)))
# mod.neonatal$m5<- bglmer(neonatal ~ prs+ pred_tm1 +evi_pm_z*(evi_im_z + I(evi_im_z^2)) + (1|yr), data = neonatal.data, family=binomial, fixef.prior = normal(sd = c(2.5)))# the intercept is a factor 
# summary(mod.neonatal$m5)
# mod.neonatal$m6<- bglmer(neonatal ~ prs+ pred_tm1 +bd_pop_z + bd_ind_z + (1|yr), data = neonatal.data, family=binomial, fixef.prior = normal) # (sd = c(10,2.5)))
# mod.neonatal$m7<- bglmer(neonatal ~ prs+ pred_tm1 + evi_pm_z + fem_tm1_z*(evi_im_z + I(evi_im_z^2)) + (1|yr), data = neonatal.data, family=binomial, fixef.prior = normal(sd = c(2.5)))# test the buffering hypothesis 
# 

# done with Original variables ONLY AND SCALED 
colnames(neonatal.data)
mod.neonatal<-list()
mod.neonatal$m1 <- glmer(neonatal ~1 + (1|year), data = neonatal.data, family=binomial)  # this one in glmer only otherwise does not converge
mod.neonatal$m2 <- bglmer(neonatal ~ evi_im_z + evi_pm_z + (1|year), data = neonatal.data, family=binomial,  cov.prior = NULL, fixef.prior = normal(sd = c(10, 2.5)))
mod.neonatal$m3<- bglmer(neonatal ~  evi_im_z + I(evi_im_z^2)+ evi_pm_z + (1|year), data = neonatal.data, family=binomial, fixef.prior = normal) #(sd = c(10,2.5)))
mod.neonatal$m4<- bglmer(neonatal ~  evi_pm_z*evi_im_z + (1|year), data = neonatal.data, family=binomial, fixef.prior = normal)# default setting, quad removed from interaction
mod.neonatal$m5<- bglmer(neonatal ~ evi_pm_z + fallMass_tm1_z*evi_im_z + (1|year), data = neonatal.data, family=binomial, cov.prior = NULL, fixef.prior = normal)# test the buffering hypothesis 
mod.neonatal$m6<- bglmer(neonatal ~ bd_pop_z + bd_ind_z + (1|year), data = neonatal.data, family=binomial, fixef.prior = normal) # (sd = c(10,2.5)))
mod.neonatal$m7 <- bglmer(neonatal ~fallMass_tm1_z + (1|year), data = neonatal.data, family=binomial, cov.prior = NULL, fixef.prior = normal(sd = c(10,2.5)))  # this one in glmer only otherwise does not converge

summary(mod.neonatal$m1)
summary(mod.neonatal$m2)
summary(mod.neonatal$m3)
summary(mod.neonatal$m4)
summary(mod.neonatal$m5)
summary(mod.neonatal$m6)
summary(mod.neonatal$m7) # not even NS

bbmle::AICctab(mod.neonatal)
# dAICc df
# m1 0.0   2 
# m7 0.4   3 
# m2 4.0   4 
# m6 4.5   4 
# m4 5.2   5 
# m3 6.3   5 
# m5 6.4   6 


t=AICctab(mod.neonatal)
t=kable(t, digits=3) %>%
  kable_styling(font_size = 10) %>%
  kable_styling("bordered") #%>%
# save_kable(file = "output/table/APPENDIX3_Table_S1_modSelNeonatal.html", self_contained = T)


# verification 
test <- glm(neonatal ~ evi_im_z + evi_pm_z, data = neonatal.data, family=binomial)
testW<- glm(neonatal ~ evi_im_z + evi_pm_z, data = neonatal.data, family=binomial, weights=birthdate)
testlogistf <- logistf(neonatal ~ evi_im_z + evi_pm_z, data = neonatal.data, family=binomial, weights=birthdate)
summary(test)
summary(testW)
testbg<- bglmer(neonatal ~ evi_im_z + evi_pm_z + (1|year), data = neonatal.data, family=binomial, fixef.prior = normal) # failed to converge
testbgW<- bglmer(neonatal ~ evi_im_z + evi_pm_z + (1|year), data = neonatal.data, family=binomial, fixef.prior = normal, weights=birthdate)
summary(testbg)
summary(testbgW)
# end verification 

t=AICctab(mod.neonatal)

# t=kable(t) %>%
#   kable_styling(font_size = 10) %>%
#   kable_styling("bordered") %>%
#   save_kable(file = "Graphs/TableS1_neonatalAppendix3.html", self_contained = T)


round(summary(mod.neonatal$m7)$coef,3)
#                  Estimate Std. Error z value Pr(>|z|)
# (Intercept)       1.828      0.238   7.694    0.000
# fallMass_tm1_z    0.241      0.187   1.292    0.196

# Groups Name        Variance Std.Dev.
# yr     (Intercept) 0.2465   0.4964  
# Number of obs: 250, groups:  yr, 18

round(confint(mod.neonatal$m7),3)
lower = 1.828 -1.96*0.238  # 1.36152
upper = 1.828 +1.96*0.238 # 2.29448

lower = 0.241 -1.96*0.187  # -0.12552
upper = 0.241 +1.96*0.187 # 0.60752

# to get rsquared, had to retransform to glmer
mod.neonatal$m7.glmer <- glmer(neonatal ~fallMass_tm1_z + (1|year), data = neonatal.data, family=binomial)
summary(mod.neonatal$m7.glmer)
confint(mod.neonatal$m7.glmer)
MuMIn::r.squaredGLMM(mod.neonatal$m7.glmer) 
# R2m        R2c
# theoretical 0.016349240 0.08500124
# delta       0.007266508 0.03777926
piecewiseSEM::rsquared(mod.neonatal$m7.glmer)
# Response   family  link method    Marginal Conditional
# 1 neonatal binomial logit  delta 0.007266508  0.03777926



# weaning mass  - year-centered ---------------------------------------------------------
# HERE # for this need to keep newRS 
fitness.data$prs <- as.factor(fitness.data$prs)
fitness.data$newRS <- as.factor(fitness.data$newRS)

levels(fitness.data$prs)
levels(fitness.data$newRS)

dat.weaning <- fitness.data[fitness.data$neonatal == "1"& fitness.data$newRS=="weaned",] # n = 197 some lambs NOT weaned have an estimate for weaning mass.
as.data.frame(colnames(dat.weaning))

dat.weaning <- dat.weaning[, c("year","mom_id","lamb_id", "weanMass","weanMass_z", "lambSurv_t1", "newRS", "prs","sex","birthdate_z","fallMass_tm1_z",
                                  "evi_pm_z","evi_im_z", "bd_ind_z", "bd_pop_z", "evi_pm","evi_im", "bd_ind", "bd_pop")]
dat.weaning <- dat.weaning[!is.na(dat.weaning$birthdate_z),] #  281
dat.weaning <- dat.weaning[!is.na(dat.weaning$prs),] # makes it drop to 250 
#dat.weaning <- dat.weaning[!is.na(dat.weaning$fem_z),] #
dat.weaning <- dat.weaning[!is.na(dat.weaning$evi_pm_z),] #
dat.weaning <- dat.weaning[!is.na(dat.weaning$evi_im_z),] #
dat.weaning <- dat.weaning[!is.na(dat.weaning$sex),] #
dat.weaning=dat.weaning[!is.na(dat.weaning$fallMass_tm1_z),]
#dat.weaning=dat.weaning[!is.na(dat.weaning$fem_tm1_z),] # n = 199
dat.weaning=dat.weaning[!is.na(dat.weaning$bd_ind_z),] # 
dat.weaning=dat.weaning[!is.na(dat.weaning$bd_pop_z),] # 
dat.weaning$mom_id <- droplevels(dat.weaning$mom_id)


# wean mass : sample sizes  -------------------------------------------------
t =table(dat.weaning$mom_id, dat.weaning$year) # this is for year 2001 to 2017
t=rowSums(t)
hist(t)
mean(t) #  3.05
range(t) # 1 10
sd(t) #  2.235499

t = dat.weaning %>% 
  group_by(mom_id) %>% 
  droplevels () %>%
  filter(newRS == "weaned") %>% # some mass were estimated even if not weaned
  summarise(lamb = length(unique(lamb_id)))
unique(t$lamb)
mean(t$lamb) #3.05
sum(t$lamb) # 183

dat.weaning$mom_id %>%
  droplevels() %>%
  table %>%
  mean() # 3.05

dat.weaning$mom_id %>%
  droplevels() %>%
  nlevels() # 63


# weaning mass : model selection and LRT on random effect structure   --------------------------------------------------------
dat.weaning$dummy <-as.factor(1)
m1 <- lmer(weanMass~ prs+ (1|dummy), data = dat.weaning,
           lmerControl(optimizer = "bobyqa", check.nlev.gtr.1 = "ignore"), REML = T)  # does not converge so a bit useless
m2 <- lmer(weanMass ~ prs + (1|mom_id), data = dat.weaning, REML = T)
m3 <- lmer(weanMass ~ prs +  (1|year),data = dat.weaning, REML = T)
m4 <- lmer(weanMass ~ prs + (1|mom_id) + (1|year),data = dat.weaning, REML = T)

VarCorr(m1) # no converge
VarCorr(m2)
VarCorr(m3)
VarCorr(m4)

anova(m1, m4) # S
# Data: dat.weaning
# Models:
#   m1: weanMass ~ prs + (1 | dummy)
# m4: weanMass ~ prs + (1 | ID) + (1 | yr)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# m1    5 1070.5 1086.3 -530.26   1060.5                         
# m4    6 1021.6 1040.6 -504.79   1009.6 50.931  1  9.566e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(m2, m4) # S
# Data: dat.weaning
# Models:
#   m2: weanMass ~ prs + (1 | ID)
# m4: weanMass ~ prs + (1 | ID) + (1 | yr)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# m2    5 1050.7 1066.5 -520.36   1040.7                         
# m4    6 1021.6 1040.6 -504.79   1009.6 31.145  1  2.395e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(m3, m4) # S
# best to add both yr and ID 
# Data: dat.weaning
# Models:
#   m3: weanMass ~ prs + (1 | yr)
# m4: weanMass ~ prs + (1 | ID) + (1 | yr)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# m3    5 1048.0 1063.9 -519.02   1038.0                         
# m4    6 1021.6 1040.6 -504.79   1009.6 28.466  1  9.534e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m4)

# check residuals
plot(weanMass~year, data = dat.weaning)
plot(weanMass_z~year, data = dat.weaning)
# end check 

# retest without control 
wean <- list()
wean$m1 <- lmer(weanMass ~ 1 + (1|mom_id) + (1|year),data = dat.weaning, REML = F)
wean$m2 <- lmer(weanMass ~ sex + (1|mom_id) + (1|year),data = dat.weaning, REML = F)
wean$m3 <- lmer(weanMass ~ evi_im_z + evi_pm_z +sex +  (1|mom_id) + (1|year),data = dat.weaning, REML = F)
wean$m4 <- lmer(weanMass ~ evi_im_z*sex + evi_pm_z +  (1|mom_id) + (1|year),data = dat.weaning, REML = F)
wean$m5 <- lmer(weanMass ~ evi_im_z + evi_pm_z + I(evi_im_z^2) + sex + (1|mom_id) + (1|year),data = dat.weaning, REML = F)
wean$m6 <- lmer(weanMass ~ evi_pm_z*(evi_im_z + I(evi_im_z^2)) + sex + (1|mom_id) + (1|year),data = dat.weaning, REML = F)
wean$m7 <- lmer(weanMass ~ bd_pop_z + bd_ind_z + sex +  (1|mom_id) + (1|year),data = dat.weaning, REML = F) # only birthdate ? 
wean$m8 <- lmer(weanMass ~ bd_ind_z + bd_pop_z + I(bd_ind_z^2)+ sex +  (1|mom_id) + (1|year),data = dat.weaning, REML = F) # only bdate %quad? 
wean$m9 <- lmer(weanMass ~ fallMass_tm1_z*(evi_im_z + I(evi_im_z^2)) + evi_pm_z + sex +   (1|mom_id) + (1|year),data = dat.weaning, REML = F) # buffering hypothesis 
wean$m10 <-lmer(weanMass ~ fallMass_tm1_z*evi_im_z + evi_pm_z + sex + (1|mom_id) + (1|year),data = dat.weaning, REML = F) # buffering hypothesis 


bbmle::AICtab(wean)
# dAIC  df
# m10   0.0 9 
# m9    0.9 11
# m8    8.9 8 
# m7   10.6 7 
# m5   12.3 8 
# m3   12.8 7 
# m4   13.7 8 
# m6   16.0 10
# m2  125.4 5 
# m1  135.5 4 

t=AICtab(wean)
t=kable(t, digits=3) %>%
  kable_styling(font_size = 10) %>%
  kable_styling("bordered") #%>%
#save_kable(file = "output/table/APPENDIX3_Table_S3_modSelWeaning.html", self_contained = T)


AICcmodavg::aictab(wean) # slightly different from bblme
# K   AICc Delta_AICc AICcWt Cum.Wt      LL
# m10  9 839.82       0.00   0.67   0.67 -410.34
# m9  11 841.31       1.48   0.32   0.98 -408.80
# m8   8 848.47       8.64   0.01   0.99 -415.78
# m7   7 850.01      10.18   0.00   1.00 -417.65
# m5   8 851.91      12.08   0.00   1.00 -417.49
# m3   7 852.20      12.37   0.00   1.00 -418.74
# m4   8 853.28      13.46   0.00   1.00 -418.18
# m6  10 856.05      16.23   0.00   1.00 -417.32
# m2   5 964.41     124.59   0.00   1.00 -477.02
# m1   4 974.44     134.61   0.00   1.00 -483.09

summary(wean$m10)

est=round(summary(wean$m10)$coef, 3) # 
# Estimate Std. Error t value
# (Intercept)               26.413      0.606  43.588
# fallMass_tm1_z             1.164      0.267   4.356
# evi_im_z                  -3.158      0.237 -13.313
# evi_pm_z                  -0.505      0.572  -0.883
# sexmale                    2.035      0.406   5.013
# fallMass_tm1_z:evi_im_z   -0.146      0.186  -0.788
est=kable(est) %>%
  kable_styling(font_size = 10) %>%
  kable_styling("bordered") #%>%
# save_kable(file = "output/graph/APPENDIX3_Table_S4_estWeaning.html", self_contained = T)

round(confint(wean$m10),3)
# 2.5 % 97.5 %
#   .sig01                   0.702  1.962
# .sig02                   1.474  3.180
# .sigma                   2.028  2.682
# (Intercept)             25.167 27.632
# fallMass_tm1_z           0.624  1.699
# evi_im_z                -3.633 -2.684
# evi_pm_z                -1.693  0.681
# sexmale                  1.231  2.844
# fallMass_tm1_z:evi_im_z -0.517  0.221
# 
MuMIn::r.squaredGLMM(wean$m10) # the first 
# R2m       R2c
# 0.4806193 0.7584177

round(summary(wean$m9)$coef, 3) # NS interaction 
# Estimate Std. Error t value
# (Intercept)                    26.231      0.634  41.347
# fallMass_tm1_z                  1.286      0.296   4.345
# evi_im_z                       -3.508      0.380  -9.229
# I(evi_im_z^2)                   0.217      0.201   1.076
# evi_pm_z                       -0.635      0.594  -1.068
# sexmale                         2.053      0.403   5.094
# fallMass_tm1_z:evi_im_z         0.204      0.332   0.614
# fallMass_tm1_z:I(evi_im_z^2)   -0.151      0.140  -1.076

round(confint(wean$m9),3)
#  2.5 % 97.5 %
# .sig01                        0.646  1.905
# .sig02                        1.512  3.253
# .sigma                        2.013  2.659
# (Intercept)                  24.925 27.500
# fallMass_tm1_z                0.687  1.878
# evi_im_z                     -4.266 -2.748
# I(evi_im_z^2)                -0.183  0.615
# evi_pm_z                     -1.870  0.589
# sexmale                       1.256  2.853
# fallMass_tm1_z:evi_im_z      -0.452  0.860
# fallMass_tm1_z:I(evi_im_z^2) -0.429  0.126

MuMIn::r.squaredGLMM(wean$m9) # second 
# R2m      R2c
# 0.486957 0.7652174

# save for figure
# save(wean, dat.weaning, fitness.data, file = "cache/modSel_weaning.RData")

# overwinter : tidying, basic stats ---------------------------------------------------
# make filter HERE : those who were weaned only
dat.overwinter <- dat.weaning[!is.na(dat.weaning$weanMass_z),] 
dat.overwinter <- dat.weaning[dat.weaning$newRS == "weaned",] # n = 183 some lambs NOT weaned have an estimate for weaning mass.
dat.overwinter$sex <- as.factor(dat.overwinter$sex)# n = 171 thanks to updated data 

cat(colnames(dat.overwinter),sep="\n")

dat.overwinter <- dat.overwinter[, c("mom_id","year", "lamb_id", "lambSurv_t1", "weanMass_z","newRS","prs","sex", "birthdate_z","fallMass_tm1_z", 
                                     "evi_pm_z", "evi_im_z", "bd_ind_z", "bd_pop_z")]
# feder : S + BD + LMS + S:BD
dat.overwinter$mom_id <- droplevels(dat.overwinter$mom_id )

# basic stats
t =table(dat.overwinter$mom_id, dat.overwinter$year) # this is for year 2001 to 2017
t=rowSums(t)
hist(t)
mean(t) # 3.05
range(t) # 1 10
sd(t) # 2.235499

t = dat.overwinter %>% 
  group_by(mom_id) %>% 
  droplevels () %>%
  summarise(lamb = length(unique(lamb_id))) # n = 56 

mean(t$lamb) # 3.05

dat.overwinter$mom_id %>%
  droplevels() %>%
  table %>%
  mean()

# overwinter : LRT, base model + model selection -------------------------------------------
m1 <- glmer(lambSurv_t1 ~ sex + (1|mom_id), data = dat.overwinter, family=binomial)
m2 <- glmer(lambSurv_t1 ~ sex +  (1|year),data = dat.overwinter,family=binomial)
m3 <- glmer(lambSurv_t1 ~ sex + (1|mom_id) + (1|year),data = dat.overwinter,family=binomial) # singular 
m3_bglmer <- bglmer(lambSurv_t1~sex+(1|year) + (1|mom_id),
                    data=dat.overwinter,
                          family=binomial,
                          fixef.prior = normal(sd = c(10, 10)))
VarCorr(m1) # no converge
VarCorr(m2)
VarCorr(m3)
VarCorr(m3_bglmer)

anova(m1, m3) # 3 seems best 
anova(m2, m3) # NS # same df unable to proceed but M2 has a lower AIC
anova(m2, m3_bglmer) # NS # same df unable to proceed

AICtab(m1, m2, m3, m3_bglmer)# m2        0.0  3 # keep year only 

# base model selection 

ow <- list()
ow$bm1 <- glmer(lambSurv_t1~sex + weanMass_z + (1|year), data = dat.overwinter, family=binomial)
ow$bm2 <- glmer(lambSurv_t1~sex*weanMass_z + (1|year), data = dat.overwinter, family=binomial)
ow$bm3 <- glmer(lambSurv_t1~sex + (1|year), data = dat.overwinter, family=binomial)

#$bm3 <- glmer(lambSurv_t1~sex*weanMass_z + fem_z + (1|yr), data = dat.overwinter, family=binomial)
#ow$bm4 <- glmer(lambSurv_t1~sex + weanMass_z + fem_z + (1|yr), data = dat.overwinter, family=binomial)
AICtab(ow)# m2        0.0  3 # keep year only 

# dAIC df
# bm1  0.0 4 # only sex and wean mass since wean mass can be influenced by mismatch
# bm2  2.0 5 
# bm3 23.4 3 


ow <- list()
ow$m1 <- glmer(lambSurv_t1~1 + (1|year), data = dat.overwinter, family=binomial)
ow$m2 <- glmer(lambSurv_t1~weanMass_z  + (1|year), data = dat.overwinter, family=binomial)
ow$m3 <- glmer(lambSurv_t1~weanMass_z + evi_im_z + evi_pm_z +(1|year), data = dat.overwinter, family=binomial)
ow$m4 <- glmer(lambSurv_t1~weanMass_z + evi_im_z + evi_pm_z + I(evi_im_z^2)+ (1|year), data = dat.overwinter, family=binomial)
ow$m5 <- glmer(lambSurv_t1~weanMass_z  + evi_pm_z*(evi_im_z + I(evi_im_z^2)) + (1|year), data = dat.overwinter, family=binomial)
ow$m6 <- glmer(lambSurv_t1~weanMass_z +bd_pop_z + bd_ind_z +  (1|year), data = dat.overwinter, family=binomial) # bdate 
ow$m7 <- glmer(lambSurv_t1~weanMass_z + bd_ind_z + bd_pop_z + I(bd_ind_z^2)+ (1|year), data = dat.overwinter, family=binomial) # bdate 
ow$m8 <- glmer(lambSurv_t1~fallMass_tm1_z*(evi_im_z + I(evi_im_z^2)) + evi_pm_z +weanMass_z+ (1|year), data = dat.overwinter, family=binomial) # buffering  - quad 
ow$m9 <- glmer(lambSurv_t1~ fallMass_tm1_z*evi_im_z + evi_pm_z +weanMass_z+ (1|year), data = dat.overwinter, family=binomial) # buffering - lin + wean mass 

AICtab(ow)
# dAIC df
# m2  0.0 3 
# m6  3.5 5 
# m3  3.8 5 
# m7  5.5 6 
# m4  5.8 6 
# m9  7.3 7 
# m5  8.7 8 
# m8 11.2 9 
# m1 20.7 2 
# AICcmodavg::aictab(ow)
# K   AICc Delta_AICc AICcWt Cum.Wt      LL
# m2 3 213.72       0.00   0.70   0.70 -103.78
# m6 5 217.43       3.71   0.11   0.81 -103.53
# m3 5 217.76       4.04   0.09   0.91 -103.69
# m7 6 219.55       5.83   0.04   0.95 -103.51
# m4 6 219.91       6.19   0.03   0.98 -103.69
# m9 7 221.56       7.85   0.01   0.99 -103.43
# m5 8 223.16       9.45   0.01   1.00 -103.12
# m8 9 225.88      12.16   0.00   1.00 -103.36
# m1 2 234.35      20.64   0.00   1.00 -115.14

t = AICcmodavg::aictab(ow)

t=kable(t, digits=3) %>%
  kable_styling(font_size = 10) %>%
  kable_styling("bordered") #%>%
  #save_kable(file = "output/table/APPENDIX3_Table_S5_OWSurv.html", self_contained = T)

ow$m2 <- glmer(lambSurv_t1~weanMass_z  + (1|yr), data = dat.overwinter, family=binomial)

round(summary(ow$m2)$coefficients, 3) # best model but not significant - only wean mass is marginal 
# Estimate Std. Error z value Pr(>|z|)
# (Intercept)    0.682      0.203   3.354    0.001
# weanMass_z     0.216      0.192   1.123    0.261
lower = 0.216 -1.96*0.192  # -0.16032
upper = 0.216+1.96*0.192  # 0.59232


# Random effects:
#   Groups Name        Variance Std.Dev.
# yr     (Intercept) 0.2126   0.4611  
# Number of obs: 166, groups:  yr, 18

round(confint(ow$m2),3)
# 2.5 % 97.5 %
#   .sig01       0.000  1.056
# (Intercept)  0.256  1.117
# weanMass_z  -0.171  0.592

MuMIn::r.squaredGLMM(ow$m2)
# R2m        R2c
# theoretical 0.011785975 0.07177110
# delta       0.008959471 0.05455901

round(summary(ow$m6)$coefficients, 3) # best model but not significant - only wean mass is marginal 
# Estimate Std. Error z value Pr(>|z|)
# (Intercept)    0.682      0.198   3.446    0.001
# weanMass_z     0.311      0.250   1.242    0.214
# bd_pop_z       0.145      0.232   0.628    0.530
# bd_ind_z       0.135      0.258   0.523    0.601
round(confint(ow$m6),3)
# 2.5 % 97.5 %
#   .sig01       0.000  1.024
# (Intercept)  0.266  1.107
# weanMass_z  -0.196  0.797
# bd_pop_z    -0.337  0.653
# bd_ind_z    -0.372  0.655

round(summary(ow$m3)$coefficients, 3) # best model but not significant - only wean mass is marginal 
# Estimate Std. Error z value Pr(>|z|)
# (Intercept)    0.683      0.201   3.399    0.001
# weanMass_z     0.272      0.244   1.116    0.265
# evi_im_z       0.090      0.249   0.362    0.717
# evi_pm_z       0.065      0.213   0.306    0.760
round(confint(ow$m3),3)
# .5 % 97.5 %
#   .sig01       0.000  1.044
# (Intercept)  0.261  1.115
# weanMass_z  -0.222  0.748
# evi_im_z    -0.403  0.591
# evi_pm_z    -0.384  0.531

round(summary(ow$m7)$coefficients, 3) # best model but not significant - only wean mass is marginal 
# Estimate Std. Error z value Pr(>|z|)
# (Intercept)      0.706      0.236   2.986    0.003
# weanMass_z       0.317      0.252   1.258    0.209
# bd_ind_z         0.186      0.380   0.489    0.625
# bd_pop_z         0.169      0.265   0.638    0.523
# I(bd_ind_z^2)   -0.030      0.163  -0.183    0.855
round(confint(ow$m7),3)
# 2.5 % 97.5 %
#   .sig01         0.000  1.021
# (Intercept)    0.216  1.191
# weanMass_z    -0.194  0.807
# bd_ind_z      -0.571  0.934
# bd_pop_z      -0.379  0.728
# I(bd_ind_z^2) -0.351  0.300
