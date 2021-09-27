# script modified 2020-01-21 by L. Renaud to show that temporal trends in the fall and spring are not the same 

library(plyr)
library (dplyr)
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
dat.trend <- read.csv2("data/mine/trends_df.csv")

# create dat.trend.yr for analyses at the annual scale --------------------
dat.trend.yr <- dat.trend[, c("year", "spring_temp", "fall_temp", "evi_up",  "evi_tm1","birthdate")] %>%
    group_by(year) %>%
    summarise_all(mean)


# fit models of temporal trends  ----------------------------------------------
mod.ts.l <- lm(spring_temp~year,data=unique(dat.trend[, c("year", "spring_temp")]) )
mod.ts.gam <- gam(spring_temp~s(year), gamma = 1.4, data=unique(dat.trend[, c("year", "spring_temp")]) )
# parametric coefficients : (Intercept)   4.3396     0.3197   13.57 7.87e-10 ***
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(year)   1      1 2.308   0.149

mod.ts.gam <- gam(spring_temp~s(year), data=unique(dat.trend[, c("year", "spring_temp")]) ) #7.146  8.185 3.226  0.0502 .
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   4.3396     0.2136   20.32 9.82e-09 ***
#   edf Ref.df     F p-value  
# s(year) 7.146  8.185 3.226  0.0502 .
summary(mod.ts.gam)
summary(mod.ts.l) # year           0.09913    0.06525   1.519    0.149
round(confint(mod.ts.l),3) # year          -0.040  0.238


mod.tf.l <- lm(fall_temp~year,data=unique(dat.trend[, c("year", "fall_temp")]) )
mod.tf.gam <- gam(fall_temp~s(year),gamma = 1.4, data=unique(dat.trend[, c("year", "fall_temp")]) ) # Woods 
summary(mod.tf.gam)# Parametric coefficients:  #4.3023     0.1397    30.8 5.62e-15 ***       #   s(year)   1      1 6.682  0.0207 *
summary(mod.tf.l) # 0.07369    0.02851   2.585   0.0207 *
round(confint(mod.tf.l),3)  #0.013   0.134
mod.tf.gam <- gam(fall_temp~s(year),data=unique(dat.trend[, c("year", "fall_temp")]) ) # param 4.3023     0.1397    30.8 5.62e-15 *** edf 1      1 6.682  0.0207 *

dat.trend.yr$ts.l=predict(mod.ts.l, data=data.trend.yr)
dat.trend.yr$tf.l=predict(mod.tf.l, data=data.trend.yr)
tf.l.ci.l=predict(mod.tf.l, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
tf.l.ci=data.frame(year = seq(2001, 2017, by=0.1), y=tf.l.ci.l$fit, ci.h=tf.l.ci.l$se.fit*1.96+tf.l.ci.l$fit, ci.l=-tf.l.ci.l$se.fit*1.96+tf.l.ci.l$fit)

ts.l.ci.l=predict(mod.ts.l, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
ts.l.ci=data.frame(year = seq(2001, 2017, by=0.1), y=ts.l.ci.l$fit, ci.h=ts.l.ci.l$se.fit*1.96+ts.l.ci.l$fit, ci.l=-ts.l.ci.l$se.fit*1.96+ts.l.ci.l$fit)

g.ts <- ggplot(dat.trend.yr, aes(x=year))+
  geom_point(aes(y=spring_temp), alpha=0.3, size =2) + geom_line(aes(y=spring_temp), alpha=0.3)+
 # geom_path(data = ts.l.ci, aes(y=y), color="black") +
  #geom_ribbon(data=ts.l.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) + 
  labs(x="Year", y="Spring temperature \n (°C)")+ 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


g.tf <- ggplot(dat.trend.yr, aes(x=year))+
  geom_point(aes(y=fall_temp), alpha=0.3, size =2) +
  geom_path(data = tf.l.ci, aes(y=y)) +
  geom_ribbon(data=tf.l.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) + 
  labs(x="Year", y="Autumn temperature \n (°C)")+ 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  
cowplot::plot_grid(g.tf, g.ts)

# we tried gams but were signifiant gamma 1.4 and the outcome was a linear regression significant (fall only)

# here gams 
# removed penalty gamma 1.4
mod.ts.gam <- gam(spring_temp~s(year), data=unique(dat.trend[, c("year", "spring_temp")]) )
mod.tf.gam <- gam(fall_temp~s(year),data=unique(dat.trend[, c("year", "fall_temp")]) ) # Woods 

dat.trend.yr$ts.g=predict(mod.ts.gam, data=data.trend.yr)
dat.trend.yr$tf.g=predict(mod.tf.gam, data=data.trend.yr)

tf.g.ci.l=predict(mod.tf.gam, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
tf.g.ci=data.frame(year = seq(2001, 2017, by=0.1), y=tf.g.ci.l$fit, ci.h=tf.g.ci.l$se.fit*1.96+tf.g.ci.l$fit, ci.l=-tf.g.ci.l$se.fit*1.96+tf.g.ci.l$fit)
g.tf.gam <- ggplot(dat.trend.yr, aes(x=year))+
  geom_point(aes(y=fall_temp), alpha=0.3, size =2) +
  geom_path(data = tf.g.ci, aes(y=y)) +
  geom_ribbon(data=tf.g.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) + 
  labs(x="Year", y="Autumn temperature \n (°C)") + 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


ts.g.ci.l=predict(mod.ts.gam, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
ts.g.ci=data.frame(year = seq(2001, 2017, by=0.1), y=ts.g.ci.l$fit, ci.h=ts.g.ci.l$se.fit*1.96+ts.g.ci.l$fit, ci.l=-ts.g.ci.l$se.fit*1.96+ts.g.ci.l$fit)

# KEEP FOR SUPP
g.ts.gam <- ggplot(dat.trend.yr, aes(x=year))+
  geom_point(aes(y=spring_temp), alpha=0.3, size =2) +
  geom_path(data = ts.g.ci, aes(y=y), color="black") +
  geom_ribbon(data=ts.g.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) + 
  labs(x="Year", y="Spring temperature \n (°C)")+ 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

cowplot::plot_grid(g.tf.gam, g.ts.gam)

mod.gu.l <- lm(evi_log_up_jul~year,data=unique(dat.trend[, c("year", "evi_log_up_jul")]) )
mod.gu.gam <- gam(evi_log_up_jul~s(year),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_log_up_jul")]) ) # Woods 
mod.gu.gam2 <- gam(evi_log_up_jul~s(spring_temp),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_log_up_jul", "spring_temp")]) ) # Woods 
mod.gu.lm2 <- lm(evi_log_up_jul~spring_temp,data=unique(dat.trend[, c("year", "evi_log_up_jul", "spring_temp")]) ) # Woods 

dat.trend.yr$pred.gu <- predict(mod.gu.gam2, newdata = dat.trend.yr) # this is to make figure 1 version 2

summary(mod.gu.gam) #s(year)   1      1 5.987  0.0272 *
summary(mod.gu.gam2) ##s(spring_temp) 1.56  1.934 18.26 0.000225 ***
summary(mod.gu.l) # year          -1.0049     0.4107  -2.447   0.0272 *
round(confint(mod.gu.l), 3) # year         -1.880   -0.130

summary(mod.gu.lm2) 
#spring_temp   -5.588      1.058   -5.28 9.26e-05 ***
round(confint(mod.gu.lm2),3) # spring_temp  -7.844  -3.332

summary(mod.gu.l)
up.l.ci.l=predict(mod.gu.l, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
up.l.ci=data.frame(year = seq(2001, 2017, by=0.1), y=up.l.ci.l$fit, ci.h=up.l.ci.l$se.fit*1.96+up.l.ci.l$fit, ci.l=-up.l.ci.l$se.fit*1.96+up.l.ci.l$fit)

g.up.yr <- ggplot(dat.trend.yr, aes(x=year))+
  geom_point(aes(y=evi_log_up_jul), alpha=0.3, size =2) +
  #geom_point(aes(y=pred.gu), colour="turquoise") + #for predeicted gu
  geom_path(data = up.l.ci, aes(y=y), color="black") +
  geom_ribbon(data=up.l.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
  labs(x="Year",y="Green-up date (Julian day)")+
  ylim(123,170)+ 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

summary(mod.gu.lm2)

# revisions - add brown down

mod.gd.l <- lm(evi_tm1~year,data=unique(dat.trend[, c("year", "evi_tm1")]) )
mod.gd.gam <- gam(evi_tm1~s(year),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_tm1")]) ) # Woods 
mod.gd.gam2 <- gam(evi_tm1~s(year),data=unique(dat.trend[, c("year", "evi_tm1")]) ) # Woods 
#mod.gd.gam2 <- gam(evi_tm1~s(spring_temp),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_tm1", "spring_temp")]) ) # Woods 
#mod.gd.lm2 <- lm(evi_tm1~spring_temp,data=unique(dat.trend[, c("year", "evi_tm1", "spring_temp")]) ) # Woods 

dat.trend.yr$pred.gd <- predict(mod.gd.gam2, newdata = dat.trend.yr) # this is to make figure 1 version 2

summary(mod.gd.gam) #s(year)   s(year)   1      1 0.084   0.776
summary(mod.gd.l) # year      year          0.1225     0.4228   0.290    0.776
summary(mod.gd.gam2) # year      year          0.1225     0.4228   0.290    0.776
round(confint(mod.gd.l), 3) # year         -0.779    1.024

do.l.ci.l=predict(mod.gd.gam, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
do.l.ci=data.frame(year = seq(2001, 2017, by=0.1), y=do.l.ci.l$fit, ci.h=do.l.ci.l$se.fit*1.96+do.l.ci.l$fit, ci.l=-do.l.ci.l$se.fit*1.96+do.l.ci.l$fit)

g.gd.yr <- ggplot(dat.trend.yr, aes(x=year))+
	geom_point(aes(y=evi_tm1), alpha=0.3, size =2) + geom_line(aes(y=evi_tm1), alpha=0.3)+
	#geom_point(aes(y=pred.gu), colour="turquoise") + #for predeicted gu
	#geom_path(data = do.l.ci, aes(y=y), color="black", linetype = "dashed") +
#	geom_ribbon(data=do.l.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
	labs(x="Year",y="Brown-down date \n (Julian day)")+
	theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

# end revisions



# SUPP MAT
up.l.ci.l=predict(mod.gu.lm2,newdata =data.frame(spring_temp = seq(1.4, 7.5, by=0.1)), se.fit=T)
up.l.ci=data.frame(spring_temp = seq(1.4, 7.5, by=0.1), y=up.l.ci.l$fit, ci.h=up.l.ci.l$se.fit*1.96+up.l.ci.l$fit, ci.l=-up.l.ci.l$se.fit*1.96+up.l.ci.l$fit)

g.up <- ggplot(dat.trend.yr, aes(x=spring_temp))+
  geom_point(aes(y=evi_log_up_jul), alpha=0.3, size =2) +
  #geom_point(aes(y=pred.gu), colour="turquoise") + #for predeicted gu 
  geom_path(data = up.l.ci, aes(y=y)) +
  geom_ribbon(data=up.l.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
  ylim(120,180) + 
  labs(x='Spring temperature (°C) ',
       y='Green-up date \n (Julian day)')+ 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


# gams of GU over time # wihtout penalty
mod.gu.gam <- gam(evi_log_up_jul~s(year),data=unique(dat.trend[, c("year", "evi_log_up_jul")]) ) # Woods 
summary(mod.gu.gam)#edf Ref.df     F p-value - s(year) 3.133   3.89 2.894  0.0744 

up.g.ci.l=predict(mod.gu.gam, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
up.g.ci=data.frame(year = seq(2001, 2017, by=0.1), y=up.g.ci.l$fit, ci.h=up.g.ci.l$se.fit*1.96+up.g.ci.l$fit, ci.l=-up.g.ci.l$se.fit*1.96+up.g.ci.l$fit)

g.up.yr <- ggplot(dat.trend.yr, aes(x=year))+
  geom_point(aes(y=evi_log_up_jul), alpha=0.3, size =2) +
  #geom_point(aes(y=pred.gu), colour="turquoise") + #for predeicted gu
  geom_path(data = up.g.ci, aes(y=y)) +
  geom_ribbon(data=up.g.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
  labs(x="Year",y="Green-up date \n (Julian day)")+
  ylim(123,175)+ 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


mod.bd.gam<- gamm4(birthdate ~ s(year), random=~(1|mom_id), REML=T, data=dat.trend)
mod.bd.lmer<- lmer(birthdate ~ year + (1|mom_id), data=dat.trend)

#mod.bd.gam2 <- gamm4(birthdate~s(fall_temp), random=~(1|mom_id) + (1|year), REML=T, data=dat.trend)
#mod.bd.lmer2 <- lmer(birthdate~fall_temp +(1|mom_id) + (1|year), REML=T, data=dat.trend)

summary(mod.bd.gam$gam) # s(year) 5.898  5.898 7.466 5.33e-07 ***
summary(mod.bd.lmer) # year          -0.5932     0.2145  -2.765
summary(mod.bd.gam2$gam) #s(fall_temp) 1.274  1.274 3.568  0.0762 .
summary(mod.bd.lmer2) #fall_temp     -4.528      2.355  -1.922
round(confint(mod.bd.lmer2),3) #fall_temp    -9.223   0.107

plot(mod.bd.gam2$gam)
plot(mod.bd.gam$gam)

pred=predict(mod.bd.gam$gam, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
pred.ci=data.frame(year = seq(2001, 2017, by=0.1), y=pred$fit, ci.h=pred$se.fit*1.96+pred$fit, ci.l=-pred$se.fit*1.96+pred$fit)

pred=predict(mod.bd.gam2$gam, newdata =dat.trend.yr, se.fit=T)
dat.trend.yr=dat.trend.yr %>% mutate(pred.bd=pred$fit, ci.h=pred$se.fit*1.96+pred$fit, ci.l=-pred$se.fit*1.96+pred$fit) # this is for fig 1, version 2


pred=predict(mod.bd.lmer, newdata =data.frame(year = seq(2001, 2017, by=0.1)),re.form=NA)
fun=function(x){predict(x, newdata=data.frame(year=seq(2001, 2017, by=0.1)), re.form=NA)}
pred=bootMer(mod.bd.lmer, FUN=fun, nsim = 1000)


pred2.ci=data.frame(year = seq(2001, 2017, by=0.1),
                    y=colMeans(pred$t),
                    ci.h=apply(pred$t, 2, quantile, 0.025),
                    ci.l=apply(pred$t, 2, quantile, 0.975))


g.bd.yr <- ggplot(dat.trend.yr, aes(x=year))+
  geom_point(aes(y=birthdate), alpha=0.3, size =2) +
  #=geom_point(data=dat.trend.yr, aes(y=pred.bd), color="turquoise") +
  #geom_point(aes(y=pred.gu), colour="turquoise") + #for predeicted gu
 # geom_path(data = pred2.ci, aes(y=y), linetype = "dotted") +
  #geom_ribbon(data=pred2.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
  geom_path(data = pred.ci, aes(y=y)) +
  geom_ribbon(data=pred.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
  ylim(123,175) + 
  labs(x='Year ',
       y='Parturition date \n (Julian day)')+ 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

# REVISIONS - ADD TEMPORAL TRENDS IN GUP AND BDOWN
cowplot::plot_grid(g.ts, g.tf.gam, g.up.yr, g.gd.yr, g.bd.yr,ncol=2, nrow=3,labels=c("a)", "b)", "c)", "d)", "e)"), align = 'v')
#ggsave("FIG_S3_gamsBDGU.png", width = 140, height = 140, units = 'mm')



pred=predict(mod.bd.gam2$gam, newdata =data.frame(fall_temp = seq(2.5, 5.5, by=0.1)), se.fit=T)
pred.ci=data.frame(fall_temp = seq(2.75, 5.5, by=0.1), y=pred$fit, ci.h=pred$se.fit*1.96+pred$fit, ci.l=-pred$se.fit*1.96+pred$fit)
range(dat.trend$fall_temp)

pred=predict(mod.bd.lmer2, newdata =data.frame(fall_temp = seq(2.88, 5.3, by=0.1)),re.form=NA)
fun=function(x){predict(x, newdata=data.frame(fall_temp = seq(2.88, 5.3, by=0.1)), re.form=NA)}
pred=bootMer(mod.bd.lmer2, FUN=fun, nsim = 1000)


pred2.ci=data.frame(fall_temp = seq(2.88, 5.3, by=0.1), 
                    y=colMeans(pred$t), 
                    ci.h=apply(pred$t, 2, quantile, 0.025), 
                    ci.l=apply(pred$t, 2, quantile, 0.975))


g.bd <- ggplot(dat.trend, aes(x=fall_temp))+
  geom_point(aes(y=birthdate), alpha=0.3, size =2) + 
  geom_path(data = pred2.ci, aes(y=y)) +
  geom_ribbon(data=pred2.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) + 
 # geom_path(data = pred.ci, aes(y=y), color="orange") +
  #geom_ribbon(data=pred.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
  ylim(120,180) + 
  labs(x='Autumn temperature (°C) ',
       y='Parturition date \n (Julian day)')+ 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

 

# p=cowplot::plot_grid(g.tf, g.ts,g.bd, g.up, ncol=2, 
#                      labels=c("a)", "b)", "c)", "d)"),
#                      align = "vh")


# SWING TO SUPP 
SUPP=cowplot::plot_grid(g.bd, g.up, ncol=2, 
										 labels=c("a)", "b)"),
										 align = "vh")
#ggsave("Graphs/revisions2/FIGS2_A2_driversBD_GU.png")
#("Graphs/revisions2/FIGS2_A2_driversBD_GU.png", width = 140, height = 140, units = 'mm')

# panel 
library(cowplot)
FIGS2_A2= cowplot::plot_grid(g.bd,g.up,
																	ncol=2,
																	align = "vh",
																	labels = c("a)", "b)"),
																	rel_widths = c(1,1))
#cowplot::save_plot("Graphs/revisions2/FIGS2_A2_revised_driversBD_GU.png", FIGS2_A2,
									 ncol = 2, #
									 nrow = 1, #
									 # each individual subplot should have an aspect ratio of 1.3
									 base_aspect_ratio = 1)


# REVISIONS - ADD TEMPORAL TRENDS IN GUP AND BDOWN
p <- cowplot::plot_grid(g.ts, g.tf.gam, g.up.yr, g.gd.yr, g.bd.yr,ncol=2, nrow=3,labels=c("a)", "b)", "c)", "d)", "e)"), align = 'v')
#ggsave("FIG_S3_gamsBDGU.png", width = 140, height = 140, units = 'mm')


mod.mis.gam2 <- gamm4(mismatch~s(year), random=~(1|mom_id), REML=T, data=dat.trend)
plot(mod.mis.gam2$gam)
mod.mis.l <- lmer(mismatch~year +I(year^2) + (1|mom_id), REML=T, data=dat.trend)
summary(mod.mis.l)

mod.mis.l2 <- lmer(mismatch~year+ (1|mom_id), REML=T, data=dat.trend)
summary(mod.mis.l2)


pred=predict(mod.mis.gam2$gam, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
pred.ci.mismatch=data.frame(year = seq(2001, 2017, by=0.1), y=pred$fit, ci.h=pred$se.fit*1.96+pred$fit, ci.l=-pred$se.fit*1.96+pred$fit)


pred=predict(mod.mis.l, newdata =data.frame(year = seq(2001, 2017, by=0.1)),re.form=NA)
fun=function(x){predict(x, newdata=data.frame(year=seq(2001, 2017, by=0.1)), re.form=NA)}
pred=bootMer(mod.mis.l, FUN=fun, nsim = 10)


pred.ci.mis=data.frame(year = seq(2001, 2017, by=0.1), 
                       y=colMeans(pred$t), 
                       ci.h=apply(pred$t, 2, quantile, 0.025), 
                       ci.l=apply(pred$t, 2, quantile, 0.975))


g.mis <- ggplot(dat.trend, aes(x=year))+
  geom_point(aes(y=mismatch), alpha =0.2) + 
  #geom_point(data=dat.trend.yr, aes(y=pred.bd), color="turquoise") +
  #geom_point(aes(y=pred.gu), colour="turquoise") + #for predeicted gu 
  geom_path(data = pred.ci.mis, aes(y=y), color="turquoise") +
  geom_ribbon(data=pred.ci.mis, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2, fill = "turquoise") + 
  geom_path(data = pred.ci.mismatch, aes(y=y)) +
  geom_ribbon(data=pred.ci.mismatch, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
  labs(x="Year", y ="Mismatch \n (number of days)")+ 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

#q = plot_grid(p, g.mis, ncol=1,labels=c("", "e)"), rel_widths = 1)


# REVISIONS - ADD TEMPORAL TRENDS IN GUP AND BDOWN
FIG1_multi= cowplot::plot_grid(g.ts, g.tf.gam, g.up.yr, g.gd.yr, g.bd.yr,g.mis, 
									 ncol=2, nrow=3,
									 labels=c("a)", "b)", "c)", "d)", "e)", 'f)'),
									 align = 'v', rel_widths = c(1,1))


cowplot::save_plot("Graphs/revisions2/FIG1_multi_revised.png", FIG1_multi,
ncol = 2, #
nrow = 3, #
# each individual subplot should have an aspect ratio of 1.3
base_aspect_ratio = 1.3)


#ggsave("FIG_S3_gamsBDGU.png", width = 140, height = 140, units = 'mm')

# p=cowplot::plot_grid(g.ts,g.tf,g.bd,  ncol=2, 
# 										 labels=c("a)", "b)", "c)", "d)"),
# 										 align = "vh")
# #cowplot::plot_grid(g.ts.gam, g.tf.gam,, g.bd.yr, g.up.yr, ncol=2, nrow=2,labels=c("a)", "b)", "c)", "d)"), align = 'v')
# cowplot::plot_grid(g.ts.gam, g.tf.gam,g.gd.yr, g.up.yr, g.bd.yr,  ncol=2, nrow=2,labels=c("a)", "b)", "c)", "d)"), align = 'v')
# 
# q = plot_grid(p, g.mis, ncol=1,labels=c("", "e)"), rel_widths = 1)




# figure of predictions
pred.mis=ggplot(dat.trend.yr,aes(y=birthdate,x=evi_log_up_jul))+
  geom_point(aes(color=year), alpha = 0.5)+
  #geom_text(aes(color=year,label=year))+
  geom_abline(intercept = 0, slope=1, linetype = "dotted")+
  geom_point(data=dat.trend.yr,aes(x=pred.gu,y=pred.bd,color=year),shape=15, size = 3)+
  #geom_path(data=dat.trend.yr,aes(x=pred.gu,y=pred.bd,color=year))+
  scale_y_continuous(limits = c(130, 170)) + 
  scale_x_continuous(limits = c(130, 170)) +
  labs(x='Green-up date (Julian day)',
       y='Parturition date \n (driven by conception date, in Julian day)')+
  coord_equal()

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

plot_grid(p, pred.mis, ncol=1,labels=c("", "e)"), rel_widths = 1)

#ggsave("Graphs/FIG1_trendsV2.png", width = 140, height = 140, units = "mm" )
ggsave("Graphs/FIGS_trendsV2.png", width = 140, height = 140, units = "mm")



# variance in BD ----------------------------------------------------------


t2=group_by(dat.trend, year) %>% 
  summarise(GroupVariance=var(mismatch), 
            mean=mean(mismatch)) 
p2=ggplot(t2,aes(x = year, y = GroupVariance)) +
  geom_point() + geom_path() + geom_point(aes(y=mean), shape=2)



load("Cache/20210426fitnessData_centeredDetrended.RData")

t3 = group_by(fitness.data, yr) %>% 
  filter(neonatal=="1") %>% 
  filter(birthdate<183) %>% 
  summarise(GroupVariance=var(evi_im), 
            mean=mean(evi_im)) 
p3 <- ggplot(t3, aes(x = yr, y = GroupVariance)) +
  geom_point()+ geom_path() + geom_point(aes(y=mean), shape=2)
plot_grid(p2, p3, 
          labels=c("a)", "b)"))




















# relative distributions - conceptual firgure for one year  -------------------------------------------------
load("Cache/20210412fitnessData_centeredDetrended.RData")
pheno <- read.csv2("Data/pheno_ram.csv", sep = ",")
pheno <- pheno[,c("year", "evi_log_up_jul")]
names(pheno) <- c("yr", "eviUP")
pheno$yr <- as.factor(as.character(pheno$yr))
# check out distribution of response variable 
fitness.data$ID<-droplevels(fitness.data$ID) # important because large amount of 0s

fitness.data <- fitness.data %>% select(yr, ID, birthdate, evi_pm, evi_im) %>% left_join(pheno)
boxplot(fitness.data$birthdate~fitness.data$yr)

# retidy in case 
colnames(fitness.data)
fitness.data$yr <- as.factor(fitness.data$yr)
fitness.data$ID <- as.factor(fitness.data$ID)

# choose year with no pop mismatch 
dt.noPM <- fitness.data %>% filter(birthdate<200) %>% filter(evi_pm>-1&evi_pm<1) %>% select(yr, evi_pm, evi_im, birthdate, eviUP)
mean(dt.noPM$birthdate)
sd(dt.noPM$birthdate)
# that parses the data
date_info <- with(dt.noPM, paste(yr, birthdate))
# Parse that character vector
dt.noPM <- cbind(dt.noPM, strptime(date_info, "%Y %j"))
names(dt.noPM)[6] <- "date_info" 
dt.noPM$newdate <- strptime(dt.noPM$date_info, "%Y-%m-%d") # this is the actual format
dt.noPM$false.evi <- rnorm(18, mean=151, sd=7)

dt.noPM$false.bd<- rnorm(18, mean=151.6667, sd=7)

#dt.noPM$newdate <- lubridate::ymd(dt.noPM$date_info)# for ggplot 

subset=dt.noPM %>% filter(evi_im>0)

 ggplot(dt.noPM) + 
  geom_histogram(aes(x=birthdate), fill = "pink", colour= "pink", alpha=0.4) + 
  #geom_density(aes(x=false.evi), fill = "turquoise", colour= "turquoise", alpha=0.4) +
  geom_vline(xintercept=151, linetype = "dotted") + # the value of EVI this year 
  geom_vline(xintercept=151.6667, linetype = "solid") + # the value of EVI this year 
  geom_histogram(data = subset,aes(x=birthdate), fill = "grey", colour= "grey")
 
 # positive value for mismatch 
dt.posPM <- fitness.data %>% filter(birthdate<200) %>% filter(evi_pm>1) %>% select(yr, evi_pm, evi_im, birthdate, eviUP) %>% filter (yr=="2016")

median(dt.posPM$birthdate)
sd(dt.posPM$birthdate)
 # that parses the data
 date_info <- with(dt.posPM, paste(yr, birthdate))
 # Parse that character vector
 dt.posPM <- cbind(dt.posPM, strptime(date_info, "%Y %j"))
 names(dt.posPM)[6] <- "date_info" 
 dt.posPM$newdate <- strptime(dt.posPM$date_info, "%Y-%m-%d") # this is the actual format
 
mean(dt.posPM$evi_pm)

df <- data.frame(matrix(ncol = 4, nrow = 100))
x <- c("pm", "im", "evi", "bd")

colnames(df) <- x
df$pm <- rnorm(100, mean=38, sd=0)
df$pm <- rnorm(100, mean=38, sd=0)
df$bd <- rnorm(19, mean=149, sd=5)

mean(dt.posPM$evi_pm)
sd(dt.posPM$evi_pm)


dt.posPM$false.bd <- rnorm(19, mean=149, sd=5)
 
 
 subset=dt.posPM %>% filter(evi_im>0)
 
 ggplot(dt.posPM) + 
   geom_histogram(aes(x=false.bd), fill = "pink", colour= "pink", alpha=0.4) + 
   geom_histogram(aes(x=false.evi), fill = "turquoise", colour= "turquoise", alpha=0.4) +
   geom_vline(xintercept=149, linetype = "dotted") + # the value of EVI this year 
   geom_vline(xintercept=123, linetype = "solid") + # the value of EVI this year 
   geom_histogram(data = subset,aes(x=birthdate), fill = "grey", colour= "grey")
 



# relative distributions supp figure - revisions  --------------------------------------------------

#date <- read.csv2("Data/tidybdates_noneonatal.csv", sep = ",")
load("Cache/20210426fitnessData_centeredDetrended.RData")
pheno <- read.csv2("Data/pheno_ram.csv", sep = ",")

date <- fitness.data[, c("yr", "birthdate", "mismatch")]

colnames(pheno)
# select columns of interest
keep <- c("year","evi_log_up_jul")
data.pheno <- pheno[keep]
colnames(data.pheno) <- c("year","greenup")

tmp2 <- merge(date, 
              data.pheno, 
              by.x = "yr", 
              by.y = "year")
colnames(tmp2)
# melt data to long format 
colnames(tmp2)
data.long <- melt(tmp2[, c("yr", "birthdate", "greenup")], id.vars="yr", variable.name="category")

all <- ggplot(tmp2, aes(x = birthdate, fill = category)) +
  geom_density(fill = "grey", alpha = .6, adjust = 2) +
  geom_density(data = data.long, aes(x = value), colour = "black", alpha = 0.8, adjust = 2) +
  scale_x_continuous(limits = c(100, 250), breaks = seq(100, 250,50)) +
  scale_y_continuous(limits = c(0, 0.075)) +
  guides(fill = FALSE) + 
  theme_pander(12) +
 labs(x= "Julian day", y = "Density")  

all = all + scale_fill_manual(values=c("grey", "#00BA38"))



#  pos mismatch  SELECT EARLY GREEN-UP < MEAN 
colnames(tmp2)
mean(tmp2$greenup, na.rm = T) # 147.5658

tmp3 <- tmp2[tmp2$mismatch>0&tmp2$greenup<147,]
data.long <- melt(tmp3[, c("yr", "birthdate", "greenup")], id.vars="yr", variable.name="category")

pos <- ggplot(tmp3, aes(x = birthdate, fill = category)) +
  geom_density(data = tmp3, fill = "grey", alpha = .6,  adjust = 2) +
  geom_density(data = data.long, aes(x = value), colour = "black", alpha = 0.8,  adjust = 2) +
  scale_x_continuous(limits = c(100, 250), breaks = seq(100, 250,50)) +
  scale_y_continuous(limits = c(0, 0.075)) +
    guides(fill = FALSE) + 
  theme_pander(12) +
  labs(x= "Julian day", y = "Density")  

pos=pos + scale_fill_manual(values=c("grey", "#00BA38"))


#  neg mismatch 
colnames(tmp2)
tmp4 <- tmp2[tmp2$mismatch<0&tmp2$greenup>147,]
data.long <- melt(tmp4[, c("yr", "birthdate", "greenup")], id.vars="yr", variable.name="category")

neg <- ggplot(tmp4, aes(x = birthdate, fill = category)) +
  geom_density(data = tmp4, fill = "grey", alpha = .6,  adjust = 2) +
  geom_density(data = data.long, aes(x = value), colour = "black", alpha = 0.8,  adjust = 2) +
  scale_x_continuous(limits = c(100, 250), breaks = seq(100, 250,50)) +
  scale_y_continuous(limits = c(0, 0.075)) +
    guides(fill = FALSE) + 
  theme_pander(12) +
  labs(x= "Julian day", y = "Density")  

neg = neg + scale_fill_manual(values=c("grey", "#00BA38"))

plot_grid(all, pos, neg, ncol=3, nrow=1, labels = c("a)", "b)", "c)"))

#ggsave("Graphs/FIGSrelativeDistributions.png", width = 200, height = 100, units = "mm")
##ggsave("Graphs/FIG2relativeDistributions.tiff")


# 
# 
# 
# 
# 
# 
# # make only 3 most important over time 
# tmp2 <- merge(date, 
#               data.pheno, 
#               by.x = "year_birth", 
#               by.y = "year")
# 
# colnames(tmp2)
# 
# data.long <- melt(tmp2[, c(1, 7, 9, 13)], id.vars="year_birth", variable.name="category")
# 
# p <- ggdensity(data.long, x = "value", fill = "category", 
#                palette = "jco", 
#                ggtheme = theme_light(), legend = "top")
# p
# facet(p, facet.by = "category")
# 
# # New facet label names for supp variable
# supp.labs <- c("NDVI", "EVI",  "GPP")
# names(supp.labs) <- c("ndviUP", "eviUP", "gppUP")
# 
# ggplot(tmp2, aes(x = birthdate, fill = category)) +
#   geom_density(data = tmp2, fill = "grey", alpha = .6) +
#   geom_density(data = data.long, aes(x = value), colour = "black", alpha = 0.8) +
#   facet_wrap(~ category, labeller = labeller(category = supp.labs)) + 
#  # scale_x_continuous(limits = c(100, 300), breaks = seq(100, 300,50)) +
#   guides(fill = FALSE) + 
#   theme_pander(10) + 
#   theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
#  labs(x= "Julian day", y= "Density") 
# # 
# #ggsave("Graphs/FIG2relativeDistributionsThree.png", width = 140, height = 110, units = "mm" )
# #ggsave("Graphs/FIG2relativeDistributionsThree.tiff", width = 140, height = 110, units = "mm")
# 


# quartz()
# ggplot(tmp2,aes(x=birthdate)) +
#   geom_histogram(aes(x =birthdate, y=..density.., fill = 'birthdate'), alpha = 0.2) +
#   geom_density(aes(x =birthdate, y=..density.., fill = "birthdate", col = 'orange', alpha = 0.2, linetype = 'solid')) +
#   geom_histogram(aes(x =ndviUP, y=..density.., fill = 'ndviUP'), alpha = 0.2) +
#   geom_density(aes(x =ndviUP, y=..density.., fill = "ndviUP", col = "blue",
#                    alpha = 0.2, linetype = 'dotted')) +
#   labs(x= "Mismatch (in number of days) between two spring events and parturition",
#        y= "Density") +
#   scale_fill_manual(values = c('birthdate' = 'orange', 'ndviUP' = 'blue')) +
#   scale_color_manual(values = c('birthdate' = 'orange', 'ndviUP' = 'blue')) +
#   # scale_linetype_manual(values = c('index3_snow' = 'solid','index3_gpp' = 'dotted')) +
#    theme_pander() +
#   theme(legend.title=element_blank())
#   dev.off()



# trend in parturition date  -------------------------------------------------------
load("Cache/scaledMismatch.RData")

# just add median birthdate into this ! v
tmp <- surv_bdate[, c("yr", "ID", "birthdate")] %>% 
  group_by(yr)%>%
  mutate(median_part = median(birthdate))

tmp$yr <- as.numeric(as.character(tmp$yr))


# trend in mediant part - do not take duplicated 
temp.trend = unique(tmp[, c("yr", "median_part")])
round(summary(lm(median_part ~ yr + I(yr^2), data=temp.trend))$coefficients,3)
 
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept) 441369.083 258524.987   1.707    0.110
# yr            -438.741    257.368  -1.705    0.110
# I(yr^2)          0.109      0.064   1.703    0.111

summary(lm(median_part ~ yr, data=tmp))$coefficients
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept) 1161.886    595.378   1.952     0.07
# yr            -0.504      0.296  -1.700     0.11
median(surv_bdate$birthdate[surv_bdate$neonatal==0]) # 153
sd(surv_bdate$birthdate[surv_bdate$neonatal==0]) # 13.65907

median(surv_bdate$birthdate[surv_bdate$neonatal==1]) # 148
sd(surv_bdate$birthdate[surv_bdate$neonatal==1]) # 11.89667

t = surv_bdate[surv_bdate$neonatal==0,] # 33 
t2 = surv_bdate[surv_bdate$neonatal==1,] # 197



# rates of change over year for phenology -------------------------------------------
tmp1 <- merge(unique(tmp[, c("yr", "median_part")]), 
              data.pheno, 
              by.x = "yr", 
              by.y = "year")
tmp1$yr <- as.numeric(as.character(tmp1$yr))

# outcome
out_start=2
out_end= 16
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

# exposure (here = year)
exp_start=1
exp_end=1
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)

number=1

dat  = tmp1[!is.na(tmp1$gppUP), ]

for (i in out_start:out_end){
  outcome = colnames(dat)[i]
  for (j in exp_start:exp_end){
    exposure = colnames(dat)[j]
    model <- lm(get(outcome) ~ get(exposure),
                na.action = na.exclude,
                data=dat)
    
    Vcov <- vcov(model, useScale = FALSE)
    beta <- coef(model)
    se <- sqrt(diag(Vcov))
    zval <- beta / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    
    out_beta[number] = as.numeric(beta[2])
    out_se[number] = as.numeric(se[2])
    out_pvalue[number] = as.numeric(pval[2])
    out_variable[number] = outcome
    number = number + 1
    
    exp_beta[number] = as.numeric(beta[2])
    exp_se[number] = as.numeric(se[2])
    exp_pvalue[number] = as.numeric(pval[2])
    exp_variable[number] = exposure
    number = number + 1
  }
}
#Create a dataframe with results:

outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)
#exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue)
outcome = outcome[!is.na(outcome$out_variable),]
rownames(outcome)<- 1:nrow(outcome)

outcome[, 2:4] = round(outcome[, 2:4], 3)

# verification on raw data
summary(lm(gppUP~as.numeric(yr), data= dat)) # earlierNS
summary(lm(ndviUP~as.numeric(yr), data= dat)) # earlier S 

summary(lm(pheno$ndvi_log_up_jul~as.numeric(year), data= pheno)) # earlier 
# end verification

#save(outcome, data.pheno, dat, file = "Graphs/rate of change.RData")

springDates <- outcome[c(2,4,6,8,10,12,14),]
springDates <- springDates %>%
  arrange(out_beta)

fallDates <- outcome[c(3,5,7,9,11,13,15),]
fallDates <- fallDates %>%
  arrange(out_beta)

# you might want to show a couple even if NS 
x.labels <- c("EVI","FPAR", "GPP", "LAI", "NDVI", "PSSNET", "SNOW")

p = springDates %>%
  mutate(name = fct_reorder(out_variable, out_beta)) %>%
  ggplot( aes(x=name, y=out_beta, ymin = out_beta-out_se, ymax = out_beta+out_se)) +
  geom_pointrange(size = 1, linetype = 1, color = "navyblue") + 
  labs(x='Spring phenologies') + 
  labs(y="Rate of change (days/year)") +  
  scale_y_continuous(limits = c(-2.5, 1)) +
  scale_x_discrete(labels= c("NDVI", "EVI", "GPP", "FPAR", "PSNNET", "LAI", "SNOW")) + 
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size =12))+
  theme_pander() +
  theme(plot.margin = unit(c(0,0,0,1), "cm"))


p
getwd()
#ggsave("Graphs/FIGS1RateSpring.pdf", width = 140, height = 140, units = "mm")
#  positive values indicate arrival (orange) or green-up (green) are getting later, negative values indicate they are getting earlier.


q = fallDates %>%
  mutate(name = fct_reorder(out_variable, out_beta)) %>%
  ggplot( aes(x=name, y=out_beta, ymin = out_beta-out_se, ymax = out_beta+out_se)) +
  geom_pointrange(size = 1, linetype = 1, color = "orange") + 
  labs(x='Autumn phenologies') + 
  labs(y="") +  
  scale_y_continuous(limits = c(-2.5, 1)) +
  scale_x_discrete(labels= c("LAI", "FPAR", "NDVI", "EVI", "PSNNET", "GPP", "SNOW")) + 
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size =12)) +
  theme_pander() +
  theme(plot.margin = unit(c(0,0,0,1), "cm"))


q

plot <- plot_grid(p, q, 
                  ncol = 2, 
                  labels = c("a)", "b)"))
plot

# save_plot("Graphs/supp/FIGS1RateBoth.png", plot,
#           base_aspect_ratio = 1.5,
#           ncol = 2,
#           nrow = 1)
# save_plot("Graphs/supp/FIGS1RateBoth.tiff", plot,
#           base_aspect_ratio = 1.5,
#           ncol = 2,
#           nrow = 1)







# TOP 3 WITH PART 
CutspringDates <- outcome[c(1,2,4,8),] # ADDED PART DATE 
CutspringDates <- CutspringDates %>%
  arrange(out_beta)

CutfallDates <- outcome[c(1,3,5,9),]
CutfallDates <- CutfallDates %>%
  arrange(out_beta)

p = CutspringDates %>%
  mutate(name = fct_reorder(out_variable, out_beta)) %>%
  ggplot( aes(x=name, y=out_beta, ymin = out_beta-out_se, ymax = out_beta+out_se)) +
  geom_pointrange(size = 1, linetype = 1, color = "navyblue") + 
  labs(x='Spring phenologies') + 
  labs(y="Rate of change (days/year)") +  
  scale_y_continuous(limits = c(-2.5, 1)) +
  scale_x_discrete(labels= c("NDVI", "EVI", "GPP", "PARTURITION")) + # doouble check in outcome ! 
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size =12)) +
  theme_pander() +
  theme(plot.margin = unit(c(0,0,0,1), "cm"))

p

q = CutfallDates %>%
  mutate(name = fct_reorder(out_variable, out_beta)) %>%
  ggplot( aes(x=name, y=out_beta, ymin = out_beta-out_se, ymax = out_beta+out_se)) +
  geom_pointrange(size = 1, linetype = 1, color = "orange") + 
  labs(x='Autumn phenologies') + 
  labs(y="") +  
  scale_y_continuous(limits = c(-2.5, 1)) +
  scale_x_discrete(labels= c("PARTURITION", "NDVI", "EVI", "GPP")) + 
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size =12)) +
  theme_pander() +
  theme(plot.margin = unit(c(0,0,0,1), "cm"))


plot <- plot_grid(p, q, 
                  ncol = 2, 
                  labels = c("a)", "b)"))
plot
# 
# save_plot("Graphs/FIG3RateTop3.png", plot,
#           base_aspect_ratio = 1.5,
#           ncol = 2,
#           nrow = 1)
# 
# save_plot("Graphs/FIG3RateTop3.tiff", plot,
#           base_aspect_ratio = 1.5,
#           ncol = 2,
#           nrow = 1)





# trend in mismatch -------------------------------------------------------

rm(list = ls ())
load("Cache/updatedMismatch.Rdata")


summary(lm(index1_gpp_z~ as.numeric(as.character(yr)) +I(as.numeric(as.character(yr))^2), data=surv_bdate))$coefficients
#                                       Estimate   Std. Error   t value     Pr(>|t|)
# (Intercept)                        5.783989e+04 9.266029e+03  6.242144 2.096879e-09
# as.numeric(as.character(yr))      -5.760336e+01 9.223806e+00 -6.245076 2.063598e-09
# I(as.numeric(as.character(yr))^2)  1.434185e-02 2.295432e-03  6.247994 2.030976e-09


# The x2 sign is positive when the model is convex and negative when the curve is concave.
# The coefficient b2 tells both the direction and steepness of the curvature (a positive value indicates the curvature is upwards while a 
#  negative value indicates the curvature is downwards).

yr = c(2001:2016)
yr* -5.760336e+01 + yr^2*1.434185e-02
yr[which.min(yr* -5.760336e+01 + yr^2*1.434185e-02)] #  2008
yr[which.max(yr* -5.760336e+01 + yr^2*1.434185e-02)] #  2016

5.783989e+04 + (2008*-5.760336e+01) + (2008^2*1.434185e-02) = -0.399


# THIS IS NOT SCALED - THIS IS FOR EFFECT SIZE 
summary(lm(index1_evi~ as.numeric(as.character(yr)) +I(as.numeric(as.character(yr))^2), data=surv_bdate))$coefficients
#                                       Estimate   Std. Error   t value     Pr(>|t|)
# (Intercept)                        7.899631e+05 1.610786e+05  4.904208 1.787734e-06
# as.numeric(as.character(yr))      -7.867124e+02 1.603446e+02 -4.906385 1.769902e-06
# I(as.numeric(as.character(yr))^2)  1.958691e-01 3.990329e-02  4.908596 1.751966e-06

yr = c(2001:2016)
yr* -7.867124e+02 + yr^2*1.958691e-01 # BASIC VALUES 

yr[which.min(yr* -7.867124e+02 + yr^2*1.958691e-01)]
yr[which.max(yr* -7.867124e+02 + yr^2*1.958691e-01)]

c <- 7.899631e+05 + (yr*-7.867124e+02) + (yr^2*1.958691e-01) # ALL VALUES FOR ALL YEAR

# [1] 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016
# [1] 11.659869  9.011476  6.754822  4.889906  3.416727  2.335288  1.645586  1.347622  1.441397  1.926910  2.804161  4.073150  5.733878
# [14]  7.786344 10.230547 13.066490

t.evi = data.frame(yr, c)








surv_bdate$yr <- as.numeric(as.character(surv_bdate$yr))

# quadratic  
colnames(surv_bdate)
ml=which(colnames(surv_bdate) %in% 
           c("index1_gpp","index1_ndvi","index1_evi", "index1_snow")) 

res2 <- ldply(ml,function(i){ # prepare dataframe of results
  
  # autumn mass and age are correlated 
  mod1 <- lm(surv_bdate[, i] ~ yr + I(yr^2), data=surv_bdate)
  r1=data.frame(colnames(surv_bdate)[i],
                linEstimate= coef(summary(mod1))[2,1], 
                linSE= coef(summary(mod1))[2,2], 
                linPvalue= coef(summary(mod1))[2,4],
                QuadEstimate= coef(summary(mod1))[3,1],
                QuadSE= coef(summary(mod1))[3,2],
                QuadPvalue= coef(summary(mod1))[3,4],
                adjR = summary(mod1)$adj.r.squared# pheno is the 2rd term in model output
  )
  
  return(r1)
})

summary(lm(index1_ndvi ~ yr + I(yr^2), data=surv_bdate))

QmismatchTrends <- xtable(res2)
QmismatchTrends[, 2:8] <- round(QmismatchTrends[, 2:8], digits = 3)
QmismatchTrends$adjR <-sort(QmismatchTrends$adjR, decreasing = F)




# linear trend over year 
colnames(surv_bdate)
ml=which(colnames(surv_bdate) %in% 
           c("index1_gpp","index1_ndvi","index1_evi", "index1_snow")) 

res3 <- ldply(ml,function(i){ # prepare dataframe of results
  
  # autumn mass and age are correlated 
  mod2 <- lm(surv_bdate[, i] ~ yr, data=surv_bdate)
  r1=data.frame(colnames(surv_bdate)[i],
                linEstimate= coef(summary(mod2))[2,1], 
                linSE= coef(summary(mod2))[2,2], 
                linPvalue= coef(summary(mod2))[2,4],
                adjR = summary(mod2)$adj.r.squared# pheno is the 2rd term in model output
  )
  
  return(r1)
})

summary(lm(index1_ndvi ~ yr, data=surv_bdate))

LinearmismatchTrends <- xtable(res3)
LinearmismatchTrends[, 2:5] <- round(LinearmismatchTrends[, 2:5], digits = 3)
LinearmismatchTrends$adjR <-sort(LinearmismatchTrends$adjR, decreasing = F)


# select columns of interest
colnames(surv_bdate)
keep <- c("yr", "index1_ndvi","index1_evi", "index1_gpp")
data.pheno <- na.omit(surv_bdate[keep])

# melt data to long format 
data.pheno.long <- melt(data.pheno[, c(1:4)], id.vars="yr", variable.name="category")
data.pheno.long$yr <- as.numeric(as.character(data.pheno.long$yr))

names <- as_labeller(c(`index1_ndvi` = "NDVI", `index1_evi` = "EVI",`index1_gpp` = "GPP")) #Necesarry to put RH% into the facet labels

ggplot(data.pheno.long, aes(yr, value)) + 
  geom_point(size=2, alpha = 0.5, colour = "navyblue") +
  geom_smooth(method = "lm", colour = "navyblue") + 
  facet_wrap( ~  category, ncol=3, labeller=names) +
  theme_pander() +
  theme(legend.position="none") + 
  scale_x_continuous("Year") + 
  labs(y = "Mismatch in number of days")


#ggsave("Graphs/FIG3mismatchTrends.pdf", width = 234, height = 174, units = "mm", pointsize = 8)
#ggsave("Graphs/FIG3mismatchTrends.png", width = 234, height = 174, units = "mm", pointsize = 8)
#ggsave("Graphs/FIG3mismatchTrends.tiff", width = 234, height = 174, units = "mm", pointsize = 8)


# prediction figure temporal trends quadratic  ------------------------------------------------------
ndvi <- lm(index1_ndvi ~ yr + I(yr^2), data=surv_bdate)
summary(ndvi)

pred.x <- expand.grid(yr=seq(min(surv_bdate$yr), max(surv_bdate$yr), length = 100))
pred.y <- predict(ndvi, newdata=pred.x, se.fit=T, interval="confidence", level=0.95)
pred.all <- data.frame(pred.x, pred.y$fit)

ndvi <- ggplot(data=pred.all, aes(x=yr, y=fit)) + 
  geom_smooth(se=F,size = 1.5, colour = "blue")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "blue") +
  geom_point(data=surv_bdate, aes(x=yr, y=index1_ndvi), size=2, alpha = 0.5, color = 'navyblue') + # removed the getCI function
  labs(x="Year",y="Mismatch based on NDVI (days)") +  
  theme(legend.position ="none") +
  scale_x_continuous(limits = c(2001, 2017),breaks=seq(2001, 2017,5)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size =12)) + 
  theme_pander(12) +
  theme(plot.margin = unit(c(0,0,0,1), "cm")) 

ndvi


evi <- lm(index1_evi ~ yr + I(yr^2), data=surv_bdate)
summary(evi)

pred.x <- expand.grid(yr=seq(min(surv_bdate$yr), max(surv_bdate$yr), length = 100))
pred.y <- predict(evi, newdata=pred.x, se.fit=T, interval="confidence", level=0.95)
pred.all <- data.frame(pred.x, pred.y$fit)

evi <- ggplot(data=pred.all, aes(x=yr, y=fit)) + 
  geom_smooth(se=F,size = 1.5, colour = "blue")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "blue") +
  geom_point(data=surv_bdate, aes(x=yr, y=index1_evi), size=2, alpha = 0.5, color = 'navyblue') + # removed the getCI function
  labs(x="Year",y="Mismatch based on EVI (days)") +  
  theme(legend.position ="none") +
  scale_x_continuous(limits = c(2001, 2017),breaks=seq(2001, 2017,5)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size =12)) + 
  theme_pander(12) +
  theme(plot.margin = unit(c(0,0,0,1), "cm")) 

evi

gpp <- lm(index1_gpp ~ yr + I(yr^2), data=surv_bdate)
summary(gpp)

pred.x <- expand.grid(yr=seq(min(surv_bdate$yr), max(surv_bdate$yr), length = 100))
pred.y <- predict(gpp, newdata=pred.x, se.fit=T, interval="confidence", level=0.95)
pred.all <- data.frame(pred.x, pred.y$fit)

gpp <- ggplot(data=pred.all, aes(x=yr, y=fit)) + 
  geom_smooth(se=F,size = 1.5, colour = "blue")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "blue") +
  geom_point(data=surv_bdate, aes(x=yr, y=index1_gpp), size=2, alpha = 0.5, color = 'navyblue') + # removed the getCI function
  labs(x="Year",y="Mismatch based on GPP (days)") +  
  theme(legend.position ="none") +
  scale_x_continuous(limits = c(2001, 2017),breaks=seq(2001, 2017,5)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size =12)) + 
  theme_pander(12) +
  theme(plot.margin = unit(c(0,0,0,1), "cm")) 

gpp



snow <- lm(index1_snow ~ yr + I(yr^2), data=surv_bdate)
summary(snow)

pred.x <- expand.grid(yr=seq(min(surv_bdate$yr), max(surv_bdate$yr), length = 100))
pred.y <- predict(snow, newdata=pred.x, se.fit=T, interval="confidence", level=0.95)
pred.all <- data.frame(pred.x, pred.y$fit)

snow <- ggplot(data=pred.all, aes(x=yr, y=fit)) + 
  geom_smooth(se=F,size = 1.5, colour = "blue")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "blue") +
  geom_point(data=surv_bdate, aes(x=yr, y=index1_snow), size=2, alpha = 0.5, color = 'navyblue') + # removed the getCI function
  labs(x="Year",y="Mismatch based on SNOW (days)") +  
  theme(legend.position ="none") +
  scale_x_continuous(limits = c(2001, 2017),breaks=seq(2001, 2017,5)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size =12)) + 
  theme_pander(12) +
  theme(plot.margin = unit(c(0,0,0,1), "cm")) 
snow




plot <- plot_grid(ndvi, evi, gpp, snow, 
                  ncol = 2, 
                  nrow =2,
                  labels = c("a)", "b)","c)", "d)"))
plot

# save_plot("Graphs/FIGS6QuadraticTrends.png", plot, 
#           base_aspect_ratio = 1.5,
#           ncol = 2, 
#           nrow = 2)
# save_plot("Graphs/FIGS6QuadraticTrends.tiff", plot, 
#           base_aspect_ratio = 1.5,
#           ncol = 2, 
#           nrow = 2)



# prediction figure temporal trends linear   ------------------------------------------------------
ndvi <- lm(index1_ndvi ~ yr, data=surv_bdate)
summary(ndvi)

pred.x <- expand.grid(yr=seq(min(surv_bdate$yr), max(surv_bdate$yr), length = 100))
pred.y <- predict(ndvi, newdata=pred.x, se.fit=T, interval="confidence", level=0.95)
pred.all <- data.frame(pred.x, pred.y$fit)

ndvi <- ggplot(data=pred.all, aes(x=yr, y=fit)) + 
  geom_smooth(se=F,size = 1.5, colour = "blue")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "blue") +
  geom_point(data=surv_bdate, aes(x=yr, y=index1_ndvi), size=2, alpha = 0.5, color = 'navyblue') + # removed the getCI function
  labs(x="Year",y="Mismatch based on NDVI (days)") +  
  theme(legend.position ="none") +
  scale_x_continuous(limits = c(2001, 2017),breaks=seq(2001, 2017,5)) +
  scale_y_continuous(limits = c(-30, 100),breaks=seq(-30, 100,20)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size =12)) + 
  theme_pander() +
  theme(plot.margin = unit(c(0,0,0,1), "cm")) 

ndvi


evi <- lm(index1_evi ~ yr, data=surv_bdate)
summary(evi)

pred.x <- expand.grid(yr=seq(min(surv_bdate$yr), max(surv_bdate$yr), length = 100))
pred.y <- predict(evi, newdata=pred.x, se.fit=T, interval="confidence", level=0.95)
pred.all <- data.frame(pred.x, pred.y$fit)

evi <- ggplot(data=pred.all, aes(x=yr, y=fit)) + 
  geom_smooth(se=F,size = 1.5, colour = "blue")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "blue") +
  geom_point(data=surv_bdate, aes(x=yr, y=index1_evi), size=2, alpha = 0.5, color = 'navyblue') + # removed the getCI function
  labs(x="Year",y="Mismatch based on EVI (days)") +  
  theme(legend.position ="none") +
  scale_x_continuous(limits = c(2001, 2017),breaks=seq(2001, 2017,5)) +
  scale_y_continuous(limits = c(-30, 100),breaks=seq(-30, 100,20)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size =12)) + 
  theme_pander() +
  theme(plot.margin = unit(c(0,0,0,1), "cm")) 

evi

gpp <- lm(index1_gpp ~ yr, data=surv_bdate)
summary(gpp)

pred.x <- expand.grid(yr=seq(min(surv_bdate$yr), max(surv_bdate$yr), length = 100))
pred.y <- predict(gpp, newdata=pred.x, se.fit=T, interval="confidence", level=0.95)
pred.all <- data.frame(pred.x, pred.y$fit)

gpp <- ggplot(data=pred.all, aes(x=yr, y=fit)) + 
  geom_smooth(se=F,size = 1.5, colour = "blue")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "blue") +
  geom_point(data=surv_bdate, aes(x=yr, y=index1_gpp), size=2, alpha = 0.5, color = 'navyblue') + # removed the getCI function
  labs(x="Year",y="Mismatch based on GPP (days)") +  
  theme(legend.position ="none") +
  scale_x_continuous(limits = c(2001, 2017),breaks=seq(2001, 2017,5)) +
  scale_y_continuous(limits = c(-30, 100),breaks=seq(-30, 100,20)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size =12)) + 
  theme_pander() +
  theme(plot.margin = unit(c(0,0,0,1), "cm")) 

gpp


# 
# snow <- lm(index1_snow ~ yr, data=surv_bdate)
# summary(snow)
# 
# pred.x <- expand.grid(yr=seq(min(surv_bdate$yr), max(surv_bdate$yr), length = 100))
# pred.y <- predict(snow, newdata=pred.x, se.fit=T, interval="confidence", level=0.95)
# pred.all <- data.frame(pred.x, pred.y$fit)
# 
# snow <- ggplot(data=pred.all, aes(x=yr, y=fit)) + 
#   geom_smooth(se=F,size = 1.5, colour = "blue")+
#   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "blue") +
#   geom_point(data=surv_bdate, aes(x=yr, y=index1_snow), size=2, alpha = 0.5, color = 'navyblue') + # removed the getCI function
#   labs(x="Year",y="Mismatch based on SNOW (days)") +  
#   theme(legend.position ="none") +
#   scale_x_continuous(limits = c(2001, 2017),breaks=seq(2001, 2017,5)) +
#   theme(axis.text=element_text(size=12)) +
#   theme(axis.title = element_text(size =12)) + 
#   theme_pander(12) +
#   theme(plot.margin = unit(c(0,0,0,1), "cm")) 
# snow




plot <- plot_grid(ndvi,gpp,evi,
                  ncol = 3, 
                  nrow =1,
                  labels = c("a)", "b)","c)"))
plot

# save_plot("Graphs/FIG4LinearTrends.png", plot, 
#           base_aspect_ratio = 1,
#           ncol = 3, 
#           nrow = 1)
# 
# save_plot("Graphs/FIG4LinearTrends.tiff", plot, 
#           base_aspect_ratio = 1,
#           ncol = 3, 
#           nrow = 1)







# SUPPLEMENTARY -----------------------------------------------------------


# temporal trends and rates of change  -----------------------------------------------------


# 7 proxies 
date <- read.csv2("Data/tidybdates_noneonatal.csv", sep = ",")
pheno <- read.csv2("Data/pheno_ram.csv", sep = ",")

colnames(pheno)
# select columns of interest
keep <- c("year","ndvi_log_up_jul","ndvi_log_do_jul","evi_log_up_jul","evi_log_do_jul","lai_log_up_jul", "lai_log_do_jul","gpp_log_up_jul",
          "gpp_log_do_jul","snow_log_up_jul","snow_log_do_jul","psnnet_log_up_jul", "psnnet_log_do_jul", "fpar_log_up_jul","fpar_log_do_jul")
data.pheno <- pheno[keep]
colnames(data.pheno) <- c("year","ndviUP","ndviDO","eviUP","eviDO","laiUP", "laiDO","gppUP",
                          "gppDO","snowUP","snowDO","psnnetUP", "psnnetDO", "fparUP","fparDO")

# melt data to long format 
#data.pheno.long <- melt(data.pheno[, c(1, 2, 4, 6, 8, 10, 12, 14)], id.vars="year", variable.name="category")
data.pheno.long <- melt(data.pheno[, c(1, 4)], id.vars="year", variable.name="category")


# temporal trends  in UP
up=ggplot(data.pheno.long, aes(year, value, colour = 'green')) + 
  geom_line(size=2, colour = 'green') +
  facet_wrap( ~  category, ncol=2) +
  theme_pander() +
  theme(legend.position="none") + 
  # scale_y_continuous("County Population within Age Categories (thousands)", 
  #                 limits=c(0, max(df$value[df$categorie==df[i]]))) +
  scale_x_continuous("Year") + 
  ylab("Julian Day")

#+
# ggtitle('Phenological indices used to calculate date of green-up at Ram Mountain, Alberta, \n\ Canada, 2000-2017')
# ggsave("Graphs/supp/FIGS2pheno_trends_UP.pdf", width = 8.5, height = 11, units = "in", pointsize = 8)
# ggsave("Graphs/supp/FIGS2pheno_trends_UP.png", width = 8.5, height = 11, units = "in", pointsize = 8)
# ggsave("Graphs/supp/FIGS2pheno_trends_UP.tiff", width = 8.5, height = 11, units = "in", pointsize = 8)

# temporal trends in DO

#data.pheno.long <- melt(data.pheno[, c(1, 3, 5, 7, 9, 11, 13, 15)], id.vars="year", variable.name="category")
data.pheno.long <- melt(data.pheno[, c(1, 5)], id.vars="year", variable.name="category")

do=ggplot(data.pheno.long, aes(year, value, colour = 'green')) + 
  geom_line(size=2, colour = 'brown') +
  facet_wrap( ~  category, ncol=2) +
  theme_pander() +
  theme(legend.position="none") + 
  scale_x_continuous("Year") +
  ylab("Julian Day")#+
# ggtitle('Phenological indices used to calculate date of green-down at Ram Mountain, Alberta, \n\ Canada, 2000-2016')


# ggsave("Graphs/supp/FIGS3pheno_trends_DO.png", width = 8.5, height = 11, units = "in", pointsize = 8)
# ggsave("Graphs/supp/FIGS3pheno_trends_DO.tiff", width = 8.5, height = 11, units = "in", pointsize = 8)

plot_grid(up, do, labels =c("a)", "b)"), ncol = 2, nrow=1)
#ggsave("Graphs/supp/FIGS1pheno_trends_UPDO.png", width = 200, height = 140, units = "mm", pointsize = 10)

# linear trend fall vs spring (for fp)---------------------------------------------

ggplot(data.pheno, aes(year, ndviUP)) + 
  geom_line(size=2, color = "green") +
  geom_line(data = data.pheno, aes(year, ndviDO), size = 2, color = "orange") + 
  theme_pander() +
  theme(legend.position="none") + 
  scale_x_continuous("Year") + 
  ylab("NDVI")#+
#ggsave("Graphs/supp/trendsNDVI.png", width = 8.5, height = 11, units = "in", pointsize = 8)

ggplot(data.pheno, aes(year, eviUP)) + 
  geom_line(size=2, color = "green") +
  geom_line(data = data.pheno, aes(year, eviDO), size = 2, color = "orange") + 
  theme_pander() +
  theme(legend.position="none") + 
  scale_x_continuous("Year") + 
  ylab("EVI")#+
#ggsave("Graphs/supp/trendsEVI.png", width = 8.5, height = 11, units = "in", pointsize = 8)








# rate of change local weather variables per month/year -------------------------

# get weather variables # determine which period 
# (16 May) earliest parturition 
# renaud et al 2019
# Precipitation was thus averaged from
# 21 October to 15 November (AIC lower than model without
# precipitation by 22.16) while mean daily temperatures were
# averaged from 30 August to 19 November (AIC lower than
# model without temperature by 21.52) in the conception year.

getwd()

rm(list = ls())
load("Data/weatherdataFR.RData")

temp<-aggregate(mean_temp~day+month+year+date,data=d,mean)
temp = temp%>%filter(year>= 2000) # added later for fig for fp
temp$day <- as.numeric(as.character(temp$day))
temp$month<-as.integer(as.character(temp$month))
#select good period
temp<-temp[temp$month>10&temp$month<=11,] # should do it by date 
temp<-temp%>%
  group_by(year) %>%
  summarise(mean(mean_temp, na.rm =T)) 
names(temp) <- c("year", "fallTemp")

# 
# # test 
# weather <- read.csv("/Users/LimoilouARenaud/Documents/PhD/Analyses/Plasticity/Data/2018-10-31_daily_myenv.csv", sep = ",") 
# 
# weather<-weather[weather$month>10&weather$month<=11,] # should do it by date 
# weather<-weather%>%
#   group_by(year) %>%
#   summarise(mean(mean_temp, na.rm =T)) 
# end test 

precip<-aggregate(total_precip~day+month+year+date,data=d,sum)
precip= precip%>%filter(year>=2000)
precip$day <- as.numeric(as.character(precip$day))
precip$month<-as.integer(as.character(precip$month))
precip<-precip[precip$month>=10&precip$month<=11,]
precip<-precip%>%
  group_by(year) %>%
  summarise(mean(total_precip, na.rm =T)) 

names(precip) <- c("year", "fallPrecip")
# 
# # test 
# weather <- read.csv("/Users/LimoilouARenaud/Documents/PhD/Analyses/Plasticity/Data/2018-10-31_daily_myenv.csv", sep = ",") 
# 
# weather<-weather[weather$month>=10&weather$month<=11,]
# weather<-weather%>%
#   group_by(year) %>%
#   summarise(mean(total_precip_mm, na.rm =T)) 
# end test 


mTf <- lm(fallTemp ~ year,
          na.action = na.exclude,
          data=temp)
mPf <- lm(fallPrecip ~ year, 
          na.action = na.exclude, 
          data = precip)

# this is for below correlation
fallprecip = precip
falltemp = temp

#select SPRING
temp<-aggregate(mean_temp~day+month+year+date,data=d,mean)
temp = temp%>%filter(year>= 2000) # added later for fig for fp

temp$day <- as.numeric(as.character(temp$day))
temp$month<-as.integer(as.character(temp$month))
temp<-temp[temp$month>04&temp$month<=05,] # should do it by date 

temp<-temp%>%
  group_by(year) %>%
  summarise(mean(mean_temp, na.rm =T)) 

names(temp) <- c("year", "springTemp")

precip<-aggregate(total_precip~day+month+year+date,data=d,sum)
precip = precip%>%filter(year>= 2000) # added later for fig for fp

precip$day <- as.numeric(as.character(precip$day))
precip$month<-as.integer(as.character(precip$month))
precip<-precip[precip$month>=04&precip$month<=05,]
precip<-precip%>%
  group_by(year) %>%
  summarise(mean(total_precip, na.rm =T)) # SHOULD I TAKE SUM ??? 

names(precip) <- c("year", "springPrecip")


mTs <- lm(springTemp ~ year,
          na.action = na.exclude,
          data=temp)
mPs <- lm(springPrecip ~ year, 
          na.action = na.exclude, 
          data = precip)
summary(mPs)


springtemp = temp 
springprecip = precip

#Create a dataframe with results:
outcome <- data.frame()
outcome<- data.frame(matrix(ncol = 4, nrow = 4))
names(outcome) <- c("variable","beta", "se", "p")

outcome[1,1] = "fallTemp"
outcome[1,2] = summary(mTf)$coefficients[2,1]
outcome[1,3] = summary(mTf)$coefficients[2,2]
outcome[1,4] = summary(mTf)$coefficients[2,4]

outcome[2,1] = "fallPrecip"
outcome[2,2] = summary(mPf)$coefficients[2,1]
outcome[2,3] = summary(mPf)$coefficients[2,2]
outcome[2,4] = summary(mPf)$coefficients[2,4]


outcome[3,1] = "springTemp"
outcome[3,2] = summary(mTs)$coefficients[2,1]
outcome[3,3] = summary(mTs)$coefficients[2,2]
outcome[3,4] = summary(mTs)$coefficients[2,4]

outcome[4,1] = "springPrecip"
outcome[4,2] = summary(mPs)$coefficients[2,1]
outcome[4,3] = summary(mPs)$coefficients[2,2]
outcome[4,4] = summary(mPs)$coefficients[2,4]
outcome[, 2:4] = round(outcome[, 2:4], 3)


# you might want to show a couple even if NS 
outcome %>%
  mutate(name = fct_reorder(variable, beta)) %>%
  ggplot(aes(x=name, y=beta, ymin = beta-se, ymax = beta+se, colour = name)) +
  geom_pointrange(size = 1, linetype = 1) + 
  labs(x='') + 
  labs(y="Rate of change (°C or mm/year) ") +  
  theme_pander() +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size =12)) +
  scale_color_manual(breaks = c("fallPrecip","springPrecip", "fallTemp", "springTemp"),
                     values = c('orange','navyblue', "orange", "navyblue")) +
  scale_x_discrete(labels=c("Spring temperature","Fall precipitation", "Fall temperature","Spring precipitation")) +
  theme(legend.position = "none")



# ggsave("Graphs/supp/FIGS4RatesTP.png", width = 160, height = 140, units = "mm")
# ggsave("Graphs/supp/FIGS4RatesTP.tiff", width = 160, height = 140, units = "mm")


# spring vs fall temp et precip pour FP -----------------------------------------------------------------

ggplot(springprecip, aes(year, springPrecip)) + 
  geom_line(size=2, color = "green") +
  geom_line(data = fallprecip, aes(year, fallPrecip), size = 2, color = "orange") + 
  theme_pander() +
  theme(legend.position="none") + 
  scale_x_continuous("Year") + 
  ylab("Precipitation (mm)")#+
#ggsave("Graphs/supp/trendsSpringVSFallPrecip.png", width = 8.5, height = 11, units = "in", pointsize = 8)


ggplot(springtemp, aes(year, springTemp)) + 
  geom_line(size=2, color = "green") +
  geom_line(data = falltemp, aes(year, fallTemp), size = 2, color = "orange") + 
  theme_pander() +
  theme(legend.position="none") + 
  scale_x_continuous("Year") + 
  ylab("Temperature (°C)")#
#ggsave("Graphs/supp/trendsSpringVSFallTemp.png", width = 8.5, height = 11, units = "in", pointsize = 8)


# fall prec
load("/Users/LimoilouARenaud/Documents/PhD/Analyses/Plasticity/Scripts/2018-12-06_dataframe_windowed_absolute.RData")
ggplot(tidy_dates, aes(as.numeric(as.character(year_birth)), prec_Windowed2)) + 
  geom_line(size=2, color = "green") +
  geom_line(data = fallprecip, aes(year, fallPrecip), size = 2, color = "orange") + 
  theme_pander() +
  theme(legend.position="none") + 
  scale_x_continuous("Year") + 
  ylab("Precipitation (mm)")#+
#ggsave("Graphs/supp/trendsSpringVSFallPrecip.png", width = 8.5, height = 11, units = "in", pointsize = 8)




# SUPP correlation table --------------------------------------------------

pheno <- read.csv2("Data/pheno_ram.csv", sep = ",")

colnames(pheno)
# select columns of interest
keep <- c("year","ndvi_log_up_jul","ndvi_log_do_jul","evi_log_up_jul","evi_log_do_jul","lai_log_up_jul", "lai_log_do_jul","gpp_log_up_jul",
          "gpp_log_do_jul","snow_log_up_jul","snow_log_do_jul","psnnet_log_up_jul", "psnnet_log_do_jul", "fpar_log_up_jul","fpar_log_do_jul")
data.pheno <- pheno[keep]
colnames(data.pheno) <- c("year","ndviUP","eviUP","laiUP", "gppUP","snowUP","psnnetUP", "fparUP", 
                         "ndviDO","eviDO","laiDO", "gppDO","snowDO","psnnetDO","fparDO")
data.env <- data.pheno %>% left_join(springprecip) %>% left_join(springtemp) %>% left_join(fallprecip) %>% left_join(falltemp) %>% arrange(year)
colnames (data.env)
cor.test(data.env$ndviUP, data.env$eviUP)

corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex")){
  #Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .0001, "****", ifelse(p < .0001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
  
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 3))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html", file = "dataEnvCorrelations.html") # mod to get output
    else print(xtable(Rnew), type="latex", file = "dataEnvCorrelations.tex")  # mod to get output
  }
} 
data.env <- na.omit(data.env[, -1])

corstars(data.env, result="html")
str(data.env)



