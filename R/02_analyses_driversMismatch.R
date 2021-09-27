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


#load("Cache/scaledMismatch.RData")

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

# set where and load
setwd("/Users/LimoilouARenaud/Documents/Scolarite/PostDocI/Mismatch")

# get data on birthdate ---------------------------------------------------------------
load("Cache/20210426fitnessData_centeredDetrended.RData")
fitness.data$date <- as.Date(fitness.data$date, format = "%Y-%m-%d")
fitness.data <- fitness.data %>% rename(year = yr)

# tidy_dates <- read.csv("~/Documents/PhD/Analyses/Plasticity/Data/tidybdates.csv")
# str(tidy_dates)
# # add TRUE date from julian day (including leap years)
# date_info <- with(tidy_dates, paste(year_birth, birthdate))
# # Parse that character vector
# tidy_dates <- cbind(tidy_dates, strptime(date_info, "%Y %j"))
# names(tidy_dates)[7] <- "date_info" # this is correct with leap years
# #tidy_dates$date_info <- as.Date(tidy_dates$date_info, format = "%Y-%m-%d")
# 
# # that parses the data
# tidy_dates$newdate <- strptime(tidy_dates$date_info, "%Y-%m-%d") # this is the actual format
# # then format it like you want
# tidy_dates$newdate <- format(tidy_dates$newdate, "%d/%m/%Y")
# 
# tidy_dates <- tidy_dates[!is.na(tidy_dates$birthdate),]
# tidy_dates <- tidy_dates[!is.na(tidy_dates$newdate),]
# tidy_dates$newdate2 <- lubridate::dmy(tidy_dates$newdate)
# 
# #tidy_dates <-tidy_dates[tidy_dates$birthdate<240,]#remove outlier birthdates 
# tidy_dates$birthdate <- as.numeric(tidy_dates$birthdate)
# tidy_dates$log_bd=log1p(tidy_dates$birthdate-min(tidy_dates$birthdate, na.rm = T)) 
# hist(tidy_dates$log_bd)
# tidy_dates <- tidy_dates[, -c(1,5:8)]
# 
# colnames(tidy_dates)
# names(tidy_dates) <- c("year", "mother_id", "birthdateJJ", "birthdateYMD","log_birthdate")


# get weather variables - update data-------------------------------------

### all nordegg stations
s<-stations_search("nordegg", interval = "day")
s

s<-as.data.frame(s)
# coordinates(s)<-~lon+lat
# proj4string(s)<-"+init=epsg:4326"
# 
# ### location of RS seems approximate or rounded. One is 12m higjher than the second.
# mapview(s,zcol="station_name")

### download nordegg cs
s<-stations_search("NORDEGG CS", interval = "day")
cs<-weather_dl(station_ids = s$station_id, start = "2000-01-01", end = "2020-01-01",interval="day")
cs<-as.data.frame(cs)

### download nordegg rs
s<-stations_search("NORDEGG RS", interval = "day")
rs<-weather_dl(station_ids = s$station_id, start = "1970-01-01", end = "2020-01-01",interval="day")
rs<-as.data.frame(rs)

### show both
plot(cs$date,cs$mean_temp,xlim=range(c(cs$date,rs$date)))
points(rs$date,rs$mean_temp,col="red")

### get overlapping dates
od<-intersect(cs$date,rs$date)

### show both for the overlap period
plot(cs$date,cs$mean_temp,xlim=range(od))
points(rs$date,rs$mean_temp,col="red")

cor.test(cs$mean_temp, rs$mean_temp)

### calculate temperature difference between the two and the mean difference
diff<-cs$mean_temp[cs$date%in%od]-rs$mean_temp[rs$date%in%od]
hist(diff,breaks=100)
offset<-mean(diff,na.rm=TRUE)
offset
cor.test(cs$mean_temp[cs$date%in%od],rs$mean_temp[rs$date%in%od])

### calculate temperature difference between the two and the mean difference
diff<-cs$total_precip[cs$date%in%od]-rs$total_precip[rs$date%in%od]
hist(diff,breaks=100)
offset<-mean(diff,na.rm=TRUE)
offset

### number of cases where rs has NA and cs has data (only 2)
table(is.na(cs$mean_temp[cs$date%in%od]) & !is.na(rs$mean_temp[rs$date%in%od]))

### merge both and keep CS over RS (do we add the offset for each variable?)
cs$station_historic<-"CS"
rs$station_historic<-"RS"
d<-rbind(cs,rs)
d<-d[order(d$date,d$station_historic),]
d<-d[!duplicated(d$date),]

### plot NA values along time
vars<-c("mean_temp","total_precip")
plot(d$date,is.na(d$mean_temp),ylim=c(0,1),yaxt="n",ylab="")
lapply(seq_along(vars),function(i){
  v<-vars[i]
  h<-seq(0,1,length.out=length(vars)+1)[i+1]
  points(d$date,ifelse(is.na(d[,v]),h,0))
  text(min(d$date),h,paste(v,"NAs"),adj=c(0,2),font=2,xpd=TRUE,col="red")
})

### lengths of runs of NA values for each variable
l<-lapply(vars,function(i){
  r<-rle(is.na(d[,i]))
  table(r[[1]][r[[2]]])  
})
names(l)<-vars
l

### add some values/columns
d$jul<-as.integer(format(d$date,"%j"))
d$day<-as.integer(d$day) # is there an error here?? replace date by day
d$year<-as.integer(d$year)

#######################################################
### linear interpolation of daily temps 
plot(mean_temp~date,d)
k<-is.na(d$mean_temp)
d$mean_temp<-zoo::na.approx(d$mean_temp,na.rm=FALSE,maxgap=Inf)
plot(d$date,d$mean_temp, xlab = "Year", ylab = "Mean daily temperature (°C)")
points(d$date,d$mean_temp,col=ifelse(k,"red","black"),pch=ifelse(k,16,1))

########################################################
### interpolate precip using seasonal and yearly trend using a cyclic gam 
m<-gam(total_precip~s(jul,bs="cc")+year,data=d,family=tw())
g<-ggpredict(m)
plot(g,facet=TRUE,raw=TRUE)

### replace NA values
wNA<-is.na(d$total_precip)
dNA<-d[wNA,]
p<-predict(m,dNA,type="response")
d$total_precip[wNA]<-p

### and show them
par(mfrow=c(2,1))
plot(d$jul,d$total_precip,ylim=c(0,5), xlab = "Date in Julian Day (Day 1 = January 1)", ylab = "Total precipitation (in mm)")
points(dNA$jul,p,col="red",pch=16)
plot(d$date,d$total_precip, ylab = "Total precipitation (in mm)", xlab = "Year")
lines(d$date,predict(m,d,type="response"),col="blue",lwd=2)
points(dNA$date,p,col="red",pch=16)

### check precip model with mgcViz as we use a tweedie distribution
sims<-t(simulate(m,100))
dsims<-createDHARMa(simulatedResponse = sims, observedResponse = model.frame(m)[,1],fittedPredictedResponse = predict(m),integerResponse=FALSE)
plot(dsims,quantreg=FALSE)


##############################################################
### compute daily time-series over which apply "temporal" window for fall and spring 

temp<-aggregate(mean_temp~month+year+day,data=d,mean) # ok day was replaced above as integer 
precip<-aggregate(total_precip~month+year+day,data=d,sum)
dates<-as.Date(paste(temp$year,temp$month,"15",sep="-")) # assigne mid dates for visualization
temp$date<-dates 
precip$date<-dates

precip$date2=as.Date(paste(precip$year,precip$month,precip$day,sep="-")) 
temp$date2=as.Date(paste(temp$year,temp$month,temp$day,sep="-")) 


### show monthly time-series with na?ve trends
par(mfrow=c(2,1))
plot(temp$date,temp$mean_temp,type="b",xaxt="n")
m<-lm(mean_temp~date,data=temp)
abline(m)
axis.Date(1,at=seq(min(temp$date),max(temp$date),by="quarter"),las=2,format="%b-%y")
plot(precip$date,precip$total_precip,type="b",xaxt="n")
m<-lm(total_precip~date,data=precip)
abline(m)
axis.Date(1,at=seq(min(temp$date),max(temp$date),by="quarter"),las=2,format="%b-%y")



#select SPRING
spring.temp<-aggregate(mean_temp~day+month+year+date,data=d,mean)
spring.temp$day <- as.numeric(as.character(spring.temp$day))
spring.temp$month<-as.integer(as.character(spring.temp$month))
spring.temp<-spring.temp[spring.temp$month>=04&spring.temp$month<=05,] # here april AND may = check douhard 

spring.temp<-spring.temp%>%
  group_by(year) %>%
  dplyr::summarise(mean(mean_temp, na.rm =T)) 


spring.precip<-aggregate(total_precip~day+month+year+date,data=d,sum)
spring.precip$day <- as.numeric(as.character(spring.precip$day))
spring.precip$month<-as.integer(as.character(spring.precip$month))
spring.precip<-precip[spring.precip$month>=04&spring.precip$month<=05,]
spring.precip<-spring.precip%>%
  group_by(year) %>%
  dplyr::summarise(mean(total_precip, na.rm =T)) # SHOULD I TAKE SUM ??? 

# mTs <- lm(temp$`mean(mean_temp, na.rm = T)` ~ year,
#           na.action = na.exclude,
#           data=temp)
# mPs <- lm(precip$`mean(total_precip, na.rm = T)` ~ year, 
#           na.action = na.exclude, 
#           data = precip)
# summary(mPs)
# summary(mTs)


spring.pheno <- cbind(spring.precip, spring.temp[,2])
names(spring.pheno) <- c("year", "spring_prec", "spring_temp")




# include spring phenology (same year) ------------------------------------


# get spring variables 
pheno_ram <- read.csv("~/Documents/PhD/Analyses/Mismatch/Data/pheno_ram.csv")
spring.pheno$year <- as.numeric(as.character(spring.pheno$year))
fitness.data$year <- as.numeric(as.character(fitness.data$year))

medBD <- fitness.data %>% group_by(year) %>% summarise(medBD = median(birthdate)) # this now includes neonatal dates 

# create pop and ind datasets 
spring.pheno.pop <- spring.pheno%>%left_join(medBD) %>% left_join(pheno_ram[, c(1, 18)]) %>% arrange(year) %>% na.omit()

spring.pheno.ind <- fitness.data[c("year", "mom_id", "birthdate", "mismatch")] %>% left_join(spring.pheno) %>% left_join(pheno_ram[, c(1, 18)]) %>% arrange(year, mom_id) %>% na.omit()


# detrend spring ----------------------------------------------------------

colnames(spring.pheno.pop)
spring.pheno.pop<-spring.pheno.pop[, c("year","spring_prec","spring_temp","medBD","evi_log_up_jul")]
spring.pheno.pop <- na.omit(spring.pheno.pop) # tt part à 2000
summary(lm(spring_prec~year+ I(year^2), data=spring.pheno.pop))
summary(lm(spring_temp~year+ I(year^2), data=spring.pheno.pop))
summary(lm(evi_log_up_jul~year+ I(year^2), data=spring.pheno.pop))

detrended.spring.pheno <- spring.pheno.pop %>% 
  mutate(springP.r = resid(lm(spring_prec ~ year)),
         springT.r = resid(lm(spring_temp ~ year)),
         evi_GU.r = resid(lm(evi_log_up_jul ~ year)))


# clean

rm(spring.precip, spring.pheno, spring.temp, rs, s)


# re-extract fall weather based on sliding window ---------------------------------

fitness.data$ym1 <- as.numeric(as.character(fitness.data$year)) - 1 # here fitness.data replaced tidy dates 

fitness.data$prec_open <- paste(fitness.data$ym1, "10-21", sep = "-")
fitness.data$prec_close <- paste(fitness.data$ym1, "11-15", sep = "-")
fitness.data$prec_open <- ymd(fitness.data$prec_open)
fitness.data$prec_close<- ymd(fitness.data$prec_close)

#class(precip$date2) == “Date”
for(i in 1:nrow(fitness.data)){
  mywindow <- precip[precip$date2 >= fitness.data$prec_open[i] & precip$date2 <= fitness.data$prec_close[i],]
  fitness.data$prec_Windowed2[i] <-mean(mywindow$total_precip,na.rm=T)
}

# temperature 
fitness.data$temp_open <- paste(fitness.data$ym1, "08-30", sep = "-") # attention to leap years. 
fitness.data$temp_close <-paste(fitness.data$ym1, "11-19", sep = "-")
fitness.data$temp_open <- ymd(fitness.data$temp_open)
fitness.data$temp_close<- ymd(fitness.data$temp_close)

for(i in 1:nrow(fitness.data)){
  mywindow <- temp$mean_temp[temp$date2 >= fitness.data$temp_open[i] & temp$date2 <= fitness.data$temp_close[i]]
  fitness.data$temp_Windowed2[i] <-mean(mywindow,na.rm=T)
}

fitness.data%>% 
  group_by(year) %>%
  dplyr::summarise(mean(temp_Windowed2, na.rm =T),
                   mean(prec_Windowed2, na.rm=T))
# 2000                                 3.55                             1.10 
# 2  2001                              2.88                             0.492
# 3  2002                              5.15                             0.731
# 4  2003                              3.21                             0.942
# 5  2004                              3.93                             1.14 
# 6  2005                              3.83                             0.477
# 7  2006                              4.09                             0.192
# 8  2007                              3.86                             0.668
# 9  2008                              4.20                             0.381
# 10  2009                              5.25                             0.115
# 11  2010                              4.55                             0.223
# 12  2011                              4.83                             0.242
# 13  2012                              4.31                             0.482
# 14  2013                              3.75                             0.602
# 15  2014                              4.82                             1.16 
# 16  2015                              4.75                             0.812
# 17  2016                              4.95                             0.748
# 18  2017                              4.77                             0.131


# compare with plast df 
load("/Users/LimoilouARenaud/Documents/PhD/Analyses/Plasticity/Cache/2018-12-06_asshole_curDT_curDT2_revisions.RData")

curDT%>% 
  group_by(year) %>%
  filter(year>1999) %>% 
  dplyr::summarise(mean(temp_Windowed2, na.rm =T),
                   mean(prec_Windowed2, na.rm=T))
# 2000                                 3.55                             1.10 
# 2  2001                              2.48                             0.5  
# 3  2002                              4.68                             0.919
# 4  2003                              2.90                             0.988
# 5  2004                              3.86                             1.54 
# 6  2005                              3.79                             0.427
# 7  2006                              4.09                             0.181
# 8  2007                              3.86                             0.664
# 9  2008                              4.20                             0.381
# 10  2009                              5.25                             0.115
# 11  2010                              4.55                             0.183
# 12  2011                              5.34                             0.173
# 13  2012                              4.31                             0.472
# 14  2013                              3.65                             0.6  
# 15  2014                              4.73                             1.33 
# 16  2015                              4.70                             0.812
# 17  2016                              5.02                             0.754
# 18  2017                              4.68                             0.131

str(fitness.data)
# fitness.data <- fitness.data[, c(1:6, 9,12)]
# colnames(fitness.data)
# names(fitness.data) <- c("year", "mother_id", "birthdateJJ", "birthdateYMD","log_birthdate", "fall_prec", "fall_temp")
# end verificaiton 


# create two df : one pop and one ind -------------------------------------


# merge data to birthdates 
fall.pheno.pop <- fitness.data %>% group_by(year) %>% summarise(medBD = median(birthdate), 
                                                          fall_prec = unique(prec_Windowed2),
                                                          fall_temp = unique(temp_Windowed2))


# include fall pheno with time lag and detrend fall -----------------------------------------------------------
pheno_ram <- read.csv("~/Documents/PhD/Analyses/Mismatch/Data/pheno_ram.csv")

fall.pheno.pop <-fall.pheno.pop %>% left_join(pheno_ram[, c(1, 19)]) %>% arrange(year) 

# add time lag 
dt2 <- fall.pheno.pop[, c("year", "evi_log_do_jul")] # jeu de données temporaire 
colnames(dt2)

dt2$ym1<- dt2$year + 1 # on enlève une année au jeu de donnée temporaire pour reculer dans le temps 
colnames(dt2)
dt2 <- dt2 %>% 
  dplyr ::rename(evi_tm1=evi_log_do_jul)

head(dt2)
head(fall.pheno.pop)

# on recombine les jeux de données avec le pas de temps 
dt3 <- merge(fall.pheno.pop, # these are ok for time steps - done in sliding window 
             dt2, # notre jeu de donées avec pas de temps en arrière 
             all.x = T, 
             by.x = c( "year"), 
             by.y = c("ym1"), 
             na.rm = T)


colnames(dt3)
fall.pheno.pop <- dt3[, c("year","medBD","fall_prec","fall_temp","evi_tm1")] %>% na.omit()


dt4 <- fall.pheno.pop %>% # see Grosbois p. 374 for linear regression technique (Graham 2003)
  mutate(bd.r = resid(lm(medBD~year,na.action = na.exclude)),
         fEVI.r = resid(lm(evi_tm1~year,na.action = na.exclude)),
         fP.r = resid(lm(fall_prec ~ year,na.action = na.exclude)),
         fT.r = resid(lm(fall_temp ~ year,na.action = na.exclude))
  )


detrended.fall.pheno <- dt4

# individual data 
fall.pheno.ind <- fitness.data[, c("year", "mom_id", "mismatch", "birthdate")] %>% left_join(fall.pheno.pop) %>% arrange(year, mom_id) %>% na.omit()

# clean up 
rm(sims, d,dNA, dt2, dt3, dt4,pheno_ram, precip, temp, tmp1j, curDT, curDT2, Birthdates, s, cs, rs, g, l, dsims, tmp2, tmp3, tmp, m, dna, dt, m, medBD, tmp1,tmp2)

# save(detrended.fall.pheno, detrended.spring.pheno, fall.pheno.pop, fall.pheno.ind, spring.pheno.pop, spring.pheno.ind, fitness.data, 
# file = "Cache/20210426revisedDf.RData")


# correlations between detrended variables  --------------------------------------

rm(list = ls())
load("Cache/20210426revisedDf.RData")


# spring phenology detrending - pop level  --------

# join fall and spring 
colnames(detrended.fall.pheno)

detrended.fall <- detrended.fall.pheno[, c("year", "bd.r","fEVI.r","fP.r","fT.r" )]
detrended.fall <- na.omit(detrended.fall)

detrended <- detrended.spring.pheno %>% left_join(detrended.fall)

# correlations - remporal trends removed either linear or quad (for 2 var)
detrended %>% ggpairs(columns = c("bd.r","springP.r", "springT.r", "evi_GU.r", "fP.r","fT.r", "fEVI.r"),
                      upper = list(continuous = wrap('cor', size = 3)),
                      columnLabels = c("Parturition date","Spring precipitation", "Spring temperature","Green-up","Fall precipitation" , "Fall temperature" , "Green-down"))
getwd()

#ggsave("Graphs/FIGSpairPlot_detrended20210426.png", width = 11, height = 8.5, units = "in")

# verif 
cor.test(detrended$springP.r,detrended$springT.r, method = "pearson")
cor.test(detrended$fEVI.r,detrended$springT.r, method = "pearson")
cor.test(detrended$bd.r,detrended$fT.r, method = "pearson")
# end verif 

rm(fallprecip, falltemp, mPf, mTf, mPs, mTs, pheno_ram, springprecip, springtemp, d)



# what are the predictors of parturition date and mismatch ----------------------------

#dataSurvUnscld<-na.omit(dataSurvUnscld) %>% arrange(yr,ID)
rm(list = ls())
load("Cache/20210426revisedDf.RData")

# take spring phenology to calculate mismatch 
colnames(spring.pheno.ind)

names(spring.pheno.ind) <-  c("year","mom_id","birthdate","mismatch","spring_prec","spring_temp","evi_up")# merge with birthdate data 

# dt4 <- merge(tidy_dates[, 1:3], 
#              spring.pheno[, c("year","spring_prec","spring_temp","medBD","ndvi_up", "evi_up")], 
#              by.x="year", 
#              by.y = "year")
# colnames(dt4)

dt4 <- na.omit(spring.pheno.ind) %>% arrange(year, mom_id) 

head(fall.pheno.ind)
head(dt4)

# join fall 
dt4 <- dt4 %>% left_join(fall.pheno.ind) %>% na.omit()
dt4$mom_id <- as.factor(as.character(dt4$mom_id))
dt4$year <- as.factor(as.character(dt4$year))
table(dt4$mom_id)
nlevels(dt4$mom_id) # 75?? 
dt4$mom_id <- droplevels(dt4$mom_id)

#save(dt4, file="Cache/finalDF.RData")


m1.g <- gamm4(birthdate ~ s(fall_temp, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m2.g <- gamm4(birthdate ~ s(fall_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m3.g <- gamm4(birthdate ~ s(evi_tm1, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m4.g <- gamm4(birthdate ~ s(spring_prec, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)
m5.g <- gamm4(birthdate ~ s(spring_temp, k=10), random = ~ (1|mom_id) + (1|year), data = dt4)

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
