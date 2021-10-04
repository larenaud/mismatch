library(ggplot2)
library(ggthemes)
library(pander)
library(cowplot)
library(tidyverse)
library(forcats)
library(ggpubr)
library(mgcv)
library(gamm4)
library(lme4)


# ggplot custom  ----------------------------------------------------------
theme_set(theme_bw())
scale_colour_discrete <- function(...,palette="Set1") {
	scale_colour_brewer(...,palette=palette)
}
scale_colour_orig <- ggplot2::scale_colour_discrete
scale_fill_discrete <- function(...,palette="Set1") {
	scale_fill_brewer(...,palette=palette)
}


# load R object of models -------------------------------------------------
load("cache/models_trends.RData")

# load data ---------------------------------------------------------------
# dat.trend <- read.csv2("data/mine/trends_data_noID.csv")
# dat.trend.yr <- dat.trend[, c("year", "spring_temp", "fall_temp", "evi_up",  "evi_tm1","birthdate")] %>%
#     group_by(year) %>%
#     summarise_all(mean)

# Figure 1 main ms, panels a-f ---------------------------------------------------------
mod.ts.l <- lm(spring_temp~year,data=unique(dat.trend[, c("year", "spring_temp")]) )
mod.tf.l <- lm(fall_temp~year,data=unique(dat.trend[, c("year", "fall_temp")]) )

dat.trend.yr$ts.l=predict(mod.ts.l, data=data.trend.yr)
dat.trend.yr$tf.l=predict(mod.tf.l, data=data.trend.yr)

tf.l.ci.l=predict(mod.tf.l, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
tf.l.ci=data.frame(year = seq(2001, 2017, by=0.1), y=tf.l.ci.l$fit, ci.h=tf.l.ci.l$se.fit*1.96+tf.l.ci.l$fit, ci.l=-tf.l.ci.l$se.fit*1.96+tf.l.ci.l$fit)

ts.l.ci.l=predict(mod.ts.l, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
ts.l.ci=data.frame(year = seq(2001, 2017, by=0.1), y=ts.l.ci.l$fit, ci.h=ts.l.ci.l$se.fit*1.96+ts.l.ci.l$fit, ci.l=-ts.l.ci.l$se.fit*1.96+ts.l.ci.l$fit)

# A
g.ts <- ggplot(dat.trend.yr, aes(x=year))+
	geom_point(aes(y=spring_temp), alpha =0.4, size =2) + geom_line(aes(y=spring_temp), alpha=0.3)+
	# geom_path(data = ts.l.ci, aes(y=y), color="black") +
	#geom_ribbon(data=ts.l.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) + 
	labs(x="Year", y="Spring temperature \n (°C)")+ 
	theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

# B 
g.tf <- ggplot(dat.trend.yr, aes(x=year))+
	geom_point(aes(y=fall_temp), alpha =0.4, size =2) +
	geom_path(data = tf.l.ci, aes(y=y)) +
	geom_ribbon(data=tf.l.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) + 
	labs(x="Year", y="Autumn temperature \n (°C)")+ 
	theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

# green-up/down
mod.gu.l <- lm(evi_up~year,data=unique(dat.trend[, c("year", "evi_up")]) )
dat.trend.yr$gu.l=predict(mod.gu.l, data=data.trend.yr)

up.l.ci.l=predict(mod.gu.l, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
up.l.ci=data.frame(year = seq(2001, 2017, by=0.1), y=up.l.ci.l$fit, ci.h=up.l.ci.l$se.fit*1.96+up.l.ci.l$fit, ci.l=-up.l.ci.l$se.fit*1.96+up.l.ci.l$fit)

# C
g.up.yr <- ggplot(dat.trend.yr, aes(x=year))+
	geom_point(aes(y=evi_up), alpha =0.4, size =2) +
	geom_path(data = up.l.ci, aes(y=y), color="black") +
	geom_ribbon(data=up.l.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
	labs(x="Year",y="Green-up date \n (Julian day)")+
	ylim(123,170)+ 
	theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

# D 
mod.gd.gam <- gam(evi_tm1~s(year),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_tm1")]) ) # Woods 
do.l.ci.l=predict(mod.gd.gam, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
do.l.ci=data.frame(year = seq(2001, 2017, by=0.1), y=do.l.ci.l$fit, ci.h=do.l.ci.l$se.fit*1.96+do.l.ci.l$fit, ci.l=-do.l.ci.l$se.fit*1.96+do.l.ci.l$fit)

g.gd.yr <- ggplot(dat.trend.yr, aes(x=year))+
	geom_point(aes(y=evi_tm1), alpha =0.4, size =2) + geom_line(aes(y=evi_tm1), alpha=0.3)+
	#geom_path(data = do.l.ci, aes(y=y), color="black", linetype = "dashed") +
	#geom_ribbon(data=do.l.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
	labs(x="Year",y="Brown-down date \n (Julian day)")+
	theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

# E 
mod.bd.gam<- gamm4(birthdate ~ s(year), random=~(1|mom_id), REML=T, data=dat.trend)

pred=predict(mod.bd.gam$gam, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
pred.ci=data.frame(year = seq(2001, 2017, by=0.1), y=pred$fit, ci.h=pred$se.fit*1.96+pred$fit, ci.l=-pred$se.fit*1.96+pred$fit)

# another way to predict from lmer # would have to use pred2 in ggplot
pred=predict(mod.bd.lmer, newdata =data.frame(year = seq(2001, 2017, by=0.1)),re.form=NA)
fun=function(x){predict(x, newdata=data.frame(year=seq(2001, 2017, by=0.1)), re.form=NA)}
pred=bootMer(mod.bd.lmer, FUN=fun, nsim = 1000)

pred2.ci=data.frame(year = seq(2001, 2017, by=0.1),
                    y=colMeans(pred$t),
                    ci.h=apply(pred$t, 2, quantile, 0.025),
                    ci.l=apply(pred$t, 2, quantile, 0.975))


g.bd.yr <- ggplot(dat.trend.yr, aes(x=year))+
	geom_point(aes(y=birthdate), alpha =0.4, size =2) +
	# geom_path(data = pred2.ci, aes(y=y), linetype = "dotted") +
	# geom_ribbon(data=pred2.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
	geom_path(data = pred.ci, aes(y=y)) +
	geom_ribbon(data=pred.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
	ylim(123,175) + 
	labs(x='Year ',
			 y='Parturition date \n (Julian day)')+ 
	theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


# F mismatch- both gamm + poly on figure 
mod.mis.gam2 <- gamm4(mismatch~s(year), random=~(1|mom_id), REML=T, data=dat.trend)
pred=predict(mod.mis.gam2$gam, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
pred.ci.mismatch=data.frame(year = seq(2001, 2017, by=0.1), y=pred$fit, ci.h=pred$se.fit*1.96+pred$fit, ci.l=-pred$se.fit*1.96+pred$fit)

mod.mis.l <- lmer(mismatch~year +I(year^2) + (1|mom_id), REML=T, data=dat.trend) # convergence issue
pred=predict(mod.mis.l, newdata =data.frame(year = seq(2001, 2017, by=0.1)),re.form=NA)
fun=function(x){predict(x, newdata=data.frame(year=seq(2001, 2017, by=0.1)), re.form=NA)}
pred=bootMer(mod.mis.l, FUN=fun, nsim = 100)

pred.ci.mis=data.frame(year = seq(2001, 2017, by=0.1), 
											 y=colMeans(pred$t), 
											 ci.h=apply(pred$t, 2, quantile, 0.025), 
											 ci.l=apply(pred$t, 2, quantile, 0.975))


g.mis <- ggplot(dat.trend, aes(x=year))+
	geom_point(aes(y=mismatch), alpha =0.4, size = 2) + 
	geom_path(data = pred.ci.mis, aes(y=y), color="turquoise", linetype = 'dotted') + # this is the poly term
	geom_ribbon(data=pred.ci.mis, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2, fill = "turquoise") + 
	geom_path(data = pred.ci.mismatch, aes(y=y)) + # this is the gamm
	geom_ribbon(data=pred.ci.mismatch, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
	labs(x="Year", y ="Mismatch \n (number of days)")+ 
	theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

# REVISED
FIG1_multi= cowplot::plot_grid(g.ts, g.tf, g.up.yr, g.gd.yr, g.bd.yr,g.mis, 
															 ncol=2, nrow=3,
															 labels=c("(a)", "(b)", "(c)", "(d)", "(e)", '(f)'),
															 align = 'v', rel_widths = c(1,1))

#cowplot::save_plot("output/graph/FIG1_multi_revised_20210928.png", FIG1_multi,
#									 ncol = 2, #
#									 nrow = 3, #
#									 base_aspect_ratio = 1.3) # each individual subplot should have an aspect ratio of 1.3

# saving with minimum 300 dpi in tiff # MAKE SURE OUTPUT IS OUT OF GIT HISTORY!! 
# ggsave("output/graph/FIG1_multi_revised_20210928.tiff", units="mm", device='tiff', dpi=300)


# Figure 2  ---------------------------------------------------------------
rm(list = ls())
load("Cache/finalDF.RData")

m10 <- gamm4(mismatch ~ s(spring_temp, k=10) + s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
summary(m10$gam)
summary(m10$mer)

newdat <- data.frame(spring_temp=seq(min(dt4$spring_temp), max(dt4$spring_temp), length = 100),
										 fall_prec=mean(dt4$fall_prec, na.rm = T), 
										 mother_id="A44",
										 year="2010")
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

FIG2_mismatch= cowplot::plot_grid(p1,p2,
																	ncol=2,
																	align = "vh",
																	labels = c("a)", "b)"),
																	rel_widths = c(1,1))

# cowplot::save_plot("Graphs/revisions2/FIG2_mismatch_revised_2.png", FIG2_mismatch,
# 									 ncol = 2, #
# 									 nrow = 1, #
# 									 base_aspect_ratio = 1)


# Figure 3 Lamb weaning mass ----------------------------------------------------------------
library(partR2)
wean$m9 <- lmer(weanMass ~ fallMass_tm1_z*(evi_im_z + I(evi_im_z^2)) + evi_pm_z + sex +   (1|mom_id) + (1|yr),data = dat.weaning, REML = T) # buffering hypothesis 
wean$m10 <-lmer(weanMass ~ fallMass_tm1_z*evi_im_z + evi_pm_z + sex + (1|mom_id) + (1|yr),data = dat.weaning, REML = T) # buffering hypothesis 
# res <- partR2(wean$m9, partvars = c("fallMass_tm1_z", "evi_im_z", "evi_pm_z", "sex", "fallMass_tm1_z:evi_im_z"), nboot=100)
# res <- partR2(wean$m3, partvars = c("evi_im_z", "evi_pm_z", "sex"), nboot=300)

summary(wean$m10)
# mom_id   (Intercept) 1.777    1.333   
# yr       (Intercept) 5.114    2.262   
# Residual             5.533    2.352  
# 
# Fixed effects:
# 	Estimate Std. Error t value
# (Intercept)              26.4056     0.6371  41.450
# fallMass_tm1_z            1.1527     0.2717   4.243
# evi_im_z                 -3.1605     0.2405 -13.143
# evi_pm_z                 -0.5049     0.6071  -0.832
# sexmale                   2.0341     0.4116   4.942
# fallMass_tm1_z:evi_im_z  -0.1446     0.1882  -0.768

t=round(summary(wean$m10)$coef, 3) # 
t=kable(t) %>%
	kable_styling(font_size = 10) %>%
	kable_styling("bordered") %>%
 #save_kable(file = "Graphs/TableS5_weaning_A3.html", self_contained = T)



MuMIn::r.squaredGLMM(wean$m10) # the first 
# R2m       R2c
# 0.4628682 0.7607955

MuMIn::r.squaredGLMM(wean$m9)
# R2m       R2c
# [1,] 0.4670302 0.7651987


str(dat.weaning)
newd <- data.frame()
nd <- expand.grid(evi_im_z=seq(-1.7, 4.2, length = 100), 
									fallMass_tm1_z = mean(dat.weaning$fallMass_tm1_z, na.rm =T), # without interaction since NS
									evi_pm_z= mean(dat.weaning$evi_pm_z, na.rm=T),
									sex = "male")
newdWM <- rbind(newd,nd)
# compared t0 
min(dat.weaning$evi_im) #  -18.7619
max(dat.weaning$evi_im) #  43.33333
min(dat.weaning$evi_im_z) #  -1.514306
max(dat.weaning$evi_im_z) #   3.49751
sd(dat.weaning$evi_im)


myfun <- function(x) predict(x,newdata=newdWM,type="response",re.form=NA) # removed allow.new.levels=T

boo <- bootMer(m10, myfun, 1000, seed = NULL, use.u = FALSE, re.form=NA,
							 verbose = T)

boo <- data.frame(boo)
str(boo)

newdWM$predi <- apply(boo, 2, mean) # removed boo$t because not in seq of x anymore
newdWM$max <- apply(boo, 2, function(x) quantile(x, 0.975))
newdWM$min <- apply(boo, 2, function(x) quantile(x, 0.025))

# # pour 'déscaler'
newdWM$evi_im<-(newdWM[,"evi_im_z"]*q.SD.IM) + q.MEAN.IM

# prediction
IM<- ggplot(newdWM, aes(y=predi, x = evi_im)) +
	geom_ribbon(aes(ymin = min, ymax = max),  colour=NA, alpha=0.2) + 
	geom_line(aes(y = predi)) + 
	geom_point(data=dat.weaning, aes(y=weanMass, x=evi_im), alpha = .5) + 
	labs(x=expression('Individual mismatch (days)')) + 
	labs(y=expression('Lamb weaning mass (kg)')) +  
	#scale_x_continuous(limits = c(-2, 4),breaks=seq(-2, 4,1)) +
	# theme_pander(12) +
	theme(legend.title=element_blank(), legend.position = c(0.3, 0.4)) +
	theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

# population mismatch

str(dat.weaning)
newd <- data.frame()
nd <- expand.grid(evi_pm_z=seq(min(dat.weaning$evi_pm_z, na.rm = T), 
															 max(dat.weaning$evi_pm_z, na.rm =T), length = 100), 
									fallMass_tm1_z = mean(dat.weaning$fallMass_tm1_z, na.rm =T), 
									evi_im_z= mean(dat.weaning$evi_im_z, na.rm=T),
									#fem_z = mean(dat.weaning$fem_z, na.rm = T), 
									sex = "female")
newdWM <- rbind(newd,nd)

myfun <- function(x) predict(x,newdata=newdWM,type="response",re.form=NA) # removed allow.new.levels=T

boo <- bootMer(m10, myfun, 1000, seed = NULL, use.u = FALSE, re.form=NA,
							 verbose = T)

boo <- data.frame(boo)
str(boo)

newdWM$predi <- apply(boo, 2, mean) # removed boo$t because not in seq of x anymore
newdWM$max <- apply(boo, 2, function(x) quantile(x, 0.975))
newdWM$min <- apply(boo, 2, function(x) quantile(x, 0.025))

# # pour 'déscaler'
newdWM$evi_pm<-(newdWM[,"evi_pm_z"]*q.SD.PM) + q.MEAN.PM

# prediction
PM<- ggplot(newdWM, aes(y=predi, x = evi_pm)) +
	geom_ribbon(aes(ymin = min, ymax = max),  colour=NA, alpha=0.2) + 
	geom_line(aes(y = predi)) + 
	geom_point(data=fitness.data, aes(y=weanMass, x=evi_pm), alpha = .5) + # perhaps take dat.weaning??? 
	labs(x=expression('Population mismatch (days)')) + 
	labs(y=expression('Lamb weaning mass (kg)')) 
#scale_x_continuous(limits = c(-2, 4),breaks=seq(-2, 4,1)) +
# theme_pander(12) +
#theme(legend.title=element_blank(), legend.position = c(0.3, 0.4)) +
#theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

# panel 
library(cowplot)
plot <- plot_grid(IM, PM, labels = c("a)", "b)"), align = 'vh')

ggsave("Graphs/revisions2/FIG3PanelWeanMass_1.png", width = 150, height = 100, units="mm")
save_plot("Graphs/revisions2FIG3PanelWeanMass_2.png", plot,
					ncol = 2, #
					nrow = 1, #
					# each individual subplot should have an aspect ratio of 1.3
					base_aspect_ratio = 1
)


# Appendix 1 --------------------------------------------------------------

# Figure S1 gams without penalty for temporal trends 
mod.ts.gam <- gam(spring_temp~s(year), data=unique(dat.trend[, c("year", "spring_temp")]))
mod.tf.gam <- gam(fall_temp~s(year),data=unique(dat.trend[, c("year", "fall_temp")]))

dat.trend.yr$ts.g=predict(mod.ts.gam, data=data.trend.yr)
dat.trend.yr$tf.g=predict(mod.tf.gam, data=data.trend.yr)

ts.g.ci.l=predict(mod.ts.gam, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
ts.g.ci=data.frame(year = seq(2001, 2017, by=0.1), y=ts.g.ci.l$fit, ci.h=ts.g.ci.l$se.fit*1.96+ts.g.ci.l$fit, ci.l=-ts.g.ci.l$se.fit*1.96+ts.g.ci.l$fit)

tf.g.ci.l=predict(mod.tf.gam, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
tf.g.ci=data.frame(year = seq(2001, 2017, by=0.1), y=tf.g.ci.l$fit, ci.h=tf.g.ci.l$se.fit*1.96+tf.g.ci.l$fit, ci.l=-tf.g.ci.l$se.fit*1.96+tf.g.ci.l$fit)

# fall panel S1 A
g.tf.gam <- ggplot(dat.trend.yr, aes(x=year))+
    geom_point(aes(y=fall_temp), alpha=0.3, size =2) +
    geom_path(data = tf.g.ci, aes(y=y)) +
    geom_ribbon(data=tf.g.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) + 
    labs(x="Year", y="Autumn temperature \n (°C)") + 
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


# spring panel S1 B
g.ts.gam <- ggplot(dat.trend.yr, aes(x=year))+
    geom_point(aes(y=spring_temp), alpha=0.3, size =2) +
    geom_path(data = ts.g.ci, aes(y=y), color="black") +
    geom_ribbon(data=ts.g.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) + 
    labs(x="Year", y="Spring temperature \n (°C)")+ 
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

cowplot::plot_grid(g.tf.gam, g.ts.gam)


# panel S1 C - Parturition date 





# panel S1 D - GU
# gams of GU over time # wihtout penalty
mod.gu.gam <- gam(evi_up~s(year),data=unique(dat.trend[, c("year", "evi_up")]) ) # Woods 
summary(mod.gu.gam)#edf Ref.df     F p-value - s(year) 3.133   3.89 2.894  0.0744 

up.g.ci.l=predict(mod.gu.gam, newdata =data.frame(year = seq(2001, 2017, by=0.1)), se.fit=T)
up.g.ci=data.frame(year = seq(2001, 2017, by=0.1), y=up.g.ci.l$fit, ci.h=up.g.ci.l$se.fit*1.96+up.g.ci.l$fit, ci.l=-up.g.ci.l$se.fit*1.96+up.g.ci.l$fit)

g.up.yr <- ggplot(dat.trend.yr, aes(x=year))+
    geom_point(aes(y=evi_up), alpha=0.3, size =2) +
    #geom_point(aes(y=pred.gu), colour="turquoise") + #for predeicted gu
    geom_path(data = up.g.ci, aes(y=y)) +
    geom_ribbon(data=up.g.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
    labs(x="Year",y="Green-up date \n (Julian day)")+
    ylim(123,175)+ 
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))



#ggsave("FIG_S3_gamsBDGU.png", width = 140, height = 140, units = 'mm')

# p=cowplot::plot_grid(g.ts,g.tf,g.bd,  ncol=2, 
# 										 labels=c("a)", "b)", "c)", "d)"),
# 										 align = "vh")
# #cowplot::plot_grid(g.ts.gam, g.tf.gam,, g.bd.yr, g.up.yr, ncol=2, nrow=2,labels=c("a)", "b)", "c)", "d)"), align = 'v')
# cowplot::plot_grid(g.ts.gam, g.tf.gam,g.gd.yr, g.up.yr, g.bd.yr,  ncol=2, nrow=2,labels=c("a)", "b)", "c)", "d)"), align = 'v')
# 
# q = plot_grid(p, g.mis, ncol=1,labels=c("", "e)"), rel_widths = 1)


# Figure S2 - Density distributions for illustrating mismatch  --------------------------------------------------

#date <- read.csv2("Data/tidybdates_noneonatal.csv", sep = ",")
load("cache/20210426fitnessData_centeredDetrended.RData")
pheno <- read.csv2("data/mine/pheno_ram.csv", sep = ",")

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




# Appendix 2 --------------------------------------------------------------
# Table S1 & S2 - Model selection of parturition date  and mismatch

# Figure S1 - correlations, detrended 
detrended_df <- read.csv2("data/mine/corDf_detrended.csv")
# correlations - remporal trends removed either linear or quad (for 2 var)
detrended_df %>% ggpairs(columns = c("bd.r","springP.r", "springT.r", "evi_GU.r", "fP.r","fT.r", "fEVI.r"),
                      upper = list(continuous = wrap('cor', size = 3)),
                      columnLabels = c("Parturition date","Spring precipitation", "Spring temperature","Green-up","Fall precipitation" , "Fall temperature" , "Green-down"))
getwd()

#ggsave("Graphs/FIGSpairPlot_detrended20210426.png", width = 11, height = 8.5, units = "in")

# verif 
cor.test(detrended_df$springP.r,detrended_df$springT.r, method = "pearson")
cor.test(detrended_df$fEVI.r,detrended_df$springT.r, method = "pearson")
cor.test(detrended_df$bd.r,detrended_df$fT.r, method = "pearson")
# end verif 


# Figure S2 - Drivers of mismatch 
# Figure S3 - Predicted mismatch based on changes in part date and green-up date 

# predict greenup~spring temp
mod.gu.gam <- gam(evi_up~s(year),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_up")]) ) # Woods 
mod.gu.gam2 <- gam(evi_up~s(spring_temp),gamma = 1.4, data=unique(dat.trend[, c("year", "evi_up", "spring_temp")]) ) # Woods 
mod.gu.lm2 <- lm(evi_up~spring_temp,data=unique(dat.trend[, c("year", "evi_up", "spring_temp")]) ) # Woods 

dat.trend.yr$pred.gu <- predict(mod.gu.gam2, newdata = dat.trend.yr) # this is to make figure 1 version 2


# predict bdate ~ fall temp
mod.bd.gam2 <- gamm4(birthdate~s(fall_temp), random=~(1|mom_id) + (1|year), REML=T, data=dat.trend)
mod.bd.lmer2 <- lmer(birthdate~fall_temp +(1|mom_id) + (1|year), REML=T, data=dat.trend)

pred=predict(mod.bd.gam2$gam, newdata =dat.trend.yr, se.fit=T)
dat.trend.yr=dat.trend.yr %>% mutate(pred.bd=pred$fit, ci.h=pred$se.fit*1.96+pred$fit, ci.l=-pred$se.fit*1.96+pred$fit) 

# make figure S3
pred.mis=ggplot(dat.trend.yr,aes(y=birthdate,x=evi_up))+
    geom_point(aes(color=year), alpha = 0.5)+
    #geom_text(aes(color=year,label=year))+
    geom_abline(intercept = 0, slope=1, linetype = "dotted")+
    geom_point(data=dat.trend.yr,aes(x=pred.gu,y=pred.bd,color=year),shape=15, size = 3)+
    #geom_path(data=dat.trend.yr,aes(x=pred.gu,y=pred.bd,color=year))+
    scale_y_continuous(limits = c(130, 170)) + 
    scale_x_continuous(limits = c(130, 170)) +
    labs(x='Green-up date (Julian day) \n (driven by spring temperature, °C',
         y='Parturition date \n (driven by conception date, in Julian day)')+
    coord_equal()

# ggsave("output/graph/APPENDIX3_FIG_S3_predictedChanges.png", width = 140, height = 140, units = "mm" )


# Appendix 3 --------------------------------------------------------------
# Table S1 Model selection of neonatal survival
# Table S2 Effect of mass on neonatal survival
# Table S3 Model selection of lamb weaning mass 
# Table S4 Effect of mismatch on lamb weaning mass 
# Table S5 Model selection of lamb overwinter survival 
# Table S6 Effect of lamb weaning mass on overwinter survival 



# extras - removed from ms ------------------------------------------------
# the effect of spring temps on GU

up.l.ci.l=predict(mod.gu.lm2,newdata =data.frame(spring_temp = seq(1.4, 7.5, by=0.1)), se.fit=T)
up.l.ci=data.frame(spring_temp = seq(1.4, 7.5, by=0.1), y=up.l.ci.l$fit, ci.h=up.l.ci.l$se.fit*1.96+up.l.ci.l$fit, ci.l=-up.l.ci.l$se.fit*1.96+up.l.ci.l$fit)

g.up <- ggplot(dat.trend.yr, aes(x=spring_temp))+
    geom_point(aes(y=evi_up), alpha=0.3, size =2) +
    #geom_point(aes(y=pred.gu), colour="turquoise") + #for predeicted gu 
    geom_path(data = up.l.ci, aes(y=y)) +
    geom_ribbon(data=up.l.ci, aes(ymin = ci.l, ymax = ci.h), alpha = 0.2) +
    ylim(120,180) + 
    labs(x='Spring temperature (°C) ',
         y='Green-up date \n (Julian day)')+ 
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

# the effect of fall temp on birthdate 
mod.bd.gam2 <- gamm4(birthdate~s(fall_temp), random=~(1|mom_id) + (1|year), REML=T, data=dat.trend)
mod.bd.lmer2 <- lmer(birthdate~fall_temp +(1|mom_id) + (1|year), REML=T, data=dat.trend)
summary(mod.bd.gam2$gam) #s(fall_temp) 1.274  1.274 3.568  0.0762 .
summary(mod.bd.lmer2) #fall_temp     -4.528      2.355  -1.922

# pred=predict(mod.bd.gam2$gam, newdata =dat.trend.yr, se.fit=T)
# dat.trend.yr=dat.trend.yr %>% mutate(pred.bd=pred$fit, ci.h=pred$se.fit*1.96+pred$fit, ci.l=-pred$se.fit*1.96+pred$fit) # this is for fig 1, version 2

pred=predict(mod.bd.gam2$gam, newdata =data.frame(fall_temp = seq(2.75, 5.5, by=0.1)), se.fit=T) # was 2.5 for min
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

FIGS2_A2= cowplot::plot_grid(g.bd,g.up,
                             ncol=2,
                             align = "vh",
                             labels = c("a)", "b)"),
                             rel_widths = c(1,1))

