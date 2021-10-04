# codes for tidying data  for analyses in 01_analysesTempTrends 
# the end result of this dataframe is on Dryad 
# extraction of raw weather and pheno variables initially by F. Rousseu - modified by L.Renaud 

# load libraries
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
library(ggeffects)
library(visreg)
library(scales)
library(zoo)
library(weathercan)
library(sp)
library(mapview)
library(mgcViz)
library(DHARMa)
library(cowplot)
library(lubridate)
library(readxl)

getwd()
rm(list = ls())


# ADD TIDYING STEPS FOR THIS DF - TO BE ADDED ON THE README.TXT 
# CAN BE FOUND IN REV 1 FILE


# create dataframe from raw data ------------------------------------------

load("Cache/neonatalData.RData")
table(surv_ind$yr, surv_ind$neonatal)

# include birthdates to individual data -----------------------------------------------------
tidybdates<- read.csv("Data/tidybdates_noneonatal.csv") # 
head(tidybdates)

summary(tidybdates$birthdate) # med = 150

tidybdates = tidybdates %>%
    group_by(year_birth) %>% 
    mutate(median = median(birthdate, na.rm = T))
# summary(lm(median~year_birth, data = tidybdates)) # annual median 
# #year_birth    -0.52025    0.04072  -12.78   <2e-16 ***
# summary(lm(birthdate~year_birth, data = tidybdates))
# # year_birth   -0.38683    0.08545  -4.527 7.93e-06 *** # all individual dates

# merge birthdates of survivors to surv_ind file 
tmp1<-merge(surv_ind, 
            tidybdates[,-1], 
            by.x =c("yr", "mom_id"), 
            by.y = c("year_birth", "mother_id"), 
            all.x=T) # garde toutes celes ci 
names(tmp1) # n = 397 # NOW 332



# check who's nas
tmp1[is.na(tmp1$birthdate),]

tmp1 <- tmp1 %>% 
    dplyr::rename(live_birthdate = birthdate)
head(tmp1)


# merge birthdates of neonatal to this tmp1 file 

neonatal_dates <- read.csv("Data/neonatal_dates.csv")
head(neonatal_dates)
summary(neonatal_dates$birthdate)


names(neonatal_dates)
neonatal_dates <- neonatal_dates %>% 
    rename(dead_birthdate = dateJJ)

tmp2 <- merge(surv_ind, 
              neonatal_dates[, c("year","mother_id","dead_birthdate" )], 
              by.x =c("yr", "mom_id"), 
              by.y = c("year", "mother_id"), 
              all.x=T) # n=125

names(tmp1)
names(tmp2)
droplevels(tmp1$lamb_id)
droplevels(tmp2$lamb_id)
tmp1$lamb_id= as.factor(tmp1$lamb_id)
tmp2$lamb_id = as.factor(tmp2$lamb_id) # pas supposés d'avoir D'ID


names(tmp1)
tmp1= tmp1[, c("yr","mom_id","wt12","fallMass","newRS", "prs", "sex","fallMass_tm1","pred","fem","pred_tm1","fem_tm1","total","neonatal", "weanMass", "lamb.surv.t1", "lamb_id","live_birthdate", "eviUP")]
names(tmp2)
tmp2= tmp2[, c("yr","mom_id","wt12","fallMass","newRS","prs", "sex","fallMass_tm1","pred","fem","pred_tm1","fem_tm1","total","neonatal", "weanMass", "lamb.surv.t1", "lamb_id","dead_birthdate","eviUP")]
names(tmp1)= c("yr","mom_id","wt12","fallMass","newRS","prs", "sex","fallMass_tm1","pred","fem","pred_tm1","fem_tm1","total","neonatal", "weanMass", "lamb.surv.t1", "lamb_id","birthdate","eviUP")
names(tmp2)= c("yr","mom_id","wt12","fallMass","newRS","prs","sex", "fallMass_tm1","pred","fem","pred_tm1","fem_tm1","total","neonatal", "weanMass", "lamb.surv.t1", "lamb_id","birthdate","eviUP")
tmp2 = rbind(tmp1, tmp2) # n = 522 # now 454
tmp2<-tmp2[!is.na(tmp2$birthdate),] %>% arrange(yr, mom_id) # n= 600+# NO NOW n=281 2020-08-06


# verif
tmp1 %>%
    filter(yr %in% c(2002, 2014))
tmp2%>%
    filter(yr %in% c(2002, 2014))
# end verif


# unite two columns each having birthdates 
# library(tidyr)
# names(tmp2)
# tmp2$birthdate <- coalesce(tmp2$live_birthdate, tmp2$dead_birthdate)


surv_bdate <- tmp2


# convert dates to nice format  ----------------------------------------------------------
# add TRUE date from julian day (including leap years)
date_info <- with(surv_bdate, paste(yr, birthdate))

tmp = strptime(date_info, "%Y %j")

# Parse that character vector
surv_bdate$date <- tmp
surv_bdate$conception <- as.Date(surv_bdate$date)-172

colnames(surv_bdate)
surv_bdate <- surv_bdate[, c("yr","mom_id","wt12","fallMass","newRS","prs","sex","fallMass_tm1","pred",
                             "fem","pred_tm1","fem_tm1","total", "lamb_id","birthdate","date","conception", "weanMass", 
                             "lamb.surv.t1","neonatal","eviUP")]


colnames(surv_bdate) <- c("yr","mom_id","springMass","fallMass","newRS","prs","sex","fallMass_tm1","pred",
                          "fem","pred_tm1","fem_tm1","total", "lamb_id","birthdate","date","conception", "weanMass",
                          "lambSurv_t1","neonatal","eviUP")

# verification and basic stat for ms  ------------------------------------------------------
# verification 
Ram_mass <- read.csv("~/Documents/PhD/Analyses/OWPC/OWPC/data/Ram_mass.csv")
colnames(surv_ind)
tmp <-  surv_ind[, c("yr","mom_id","age","wt12","fallMass", "neonatal","lamb_id","weanMass")]

tmp_lb <- merge(tmp,
                Ram_mass, 
                by.x = c("yr", "lamb_id"), 
                by.y=c("yr", "ID"))

# ok for some reasons lots of lambs which survived the neonatal stage could not have their mass estimated.


# verification 
table(surv_bdate$yr, surv_bdate$neonatal)

median(surv_bdate$birthdate, na.rm = T) # 150 total 

median(surv_bdate$birthdate[surv_bdate$neonatal == "0"],na.rm = T) # 153 - June 2 julian calendar  
median(surv_bdate$birthdate[surv_bdate$neonatal == "1"],na.rm = T) # 149 - 29 may


# end verification 

# calculate mismatch ------------------------------------------------------
# index available for year 2000-2001

# took birthdate, not conception and reversed the subtraction # reed et al. 2013

surv_bdate$mismatch <- as.Date(surv_bdate$date) - as.Date(surv_bdate$eviUP)


# calculate individual timing  ----------------------------------------------

surv_bdate <- surv_bdate %>%
    group_by(yr) %>%
    mutate(median = median(birthdate, na.rm =T), 
           timing = birthdate - median)


# clean up file for analyses 
colnames(surv_bdate)

surv_bdate <- surv_bdate[, c("yr","mom_id","neonatal","weanMass", "lambSurv_t1","newRS", "prs","sex","fallMass","fallMass_tm1","pred",
                             "fem","pred_tm1","fem_tm1","lamb_id","birthdate","date",
                             "mismatch", "timing")] # removed indices 2 and 3 may 20 2020

rm(tmp, tmp1, tmp2)
# prepare for analyses 

colnames(surv_bdate)
head(surv_bdate) # with birthdates and conception dates 
#surv_bdate$index1_gpp= as.integer(surv_bdate$index1_gpp)
#surv_bdate$index1_snow= as.integer(surv_bdate$index1_snow)
#surv_bdate$index1_ndvi= as.integer(surv_bdate$index1_ndvi)
surv_bdate$mismatch= as.integer(surv_bdate$mismatch)
surv_bdate$timing= as.integer(surv_bdate$timing)

# keep timing as number of days

surv_bdate$yr <- as.factor(as.character(surv_bdate$yr))

str(surv_bdate)
surv_bdate<- as.data.frame(surv_bdate)
surv_bdate$pred_tm1 <- as.factor(surv_bdate$pred_tm1)
surv_bdate$pred <- as.factor(surv_bdate$pred)
surv_bdate$neonatal <- as.factor(surv_bdate$neonatal)
surv_bdate$prs <- as.factor(surv_bdate$prs)
surv_bdate$prs<-droplevels(surv_bdate$prs)

surv_bdate$newRS <- as.factor(surv_bdate$newRS)
surv_bdate$lambSurv_t1 <- as.factor(surv_bdate$lambSurv_t1)
surv_bdate$newRS <- as.factor(surv_bdate$newRS)

# surv_bdate <- surv_bdate[!is.na(surv_bdate$birthdate),] #  264
# surv_bdate <- surv_bdate[!is.na(surv_bdate$fallMass_tm1),] # DROPS TO -
# surv_bdate <- surv_bdate[!is.na(surv_bdate$prs),] # 
# surv_bdate <- surv_bdate[!is.na(surv_bdate$fem),] #   
# surv_bdate=surv_bdate[!is.na(surv_bdate$fem_tm1),] # 250
#dataSurvUnscld = surv_bdate
# save unscaled 
#save(dataSurvUnscld, file = "Cache/unscldMismatch.RData")


# save data for differential selection analysis (no need of control variable without NAS)  ------------------------------

# data up to 2016 (complete)
surv_bdate <- surv_bdate[!is.na(surv_bdate$yr),] #  284
surv_bdate$yr<-droplevels(surv_bdate$yr)
surv_bdate$mom_id<-droplevels(surv_bdate$mom_id)

colnames(surv_bdate)

# # modified apr 8 2021 in response to rev#2 comment to do sel analyses on bd only, not mismatch 
# selection.data <- surv_bdate[, c("yr","ID","neonatal", "lambSurv_t1", "birthdate" )]
# #save(selection.data, file = "Cache/20210408selectionData.RData")
# # this goes to scrip 04_selection_analysis
# rm(dataSurvUnscld, females, neonatal_dates, pheno, popPheno, selection.data, neosurv)


# mean centering - Reed et al. 2013 --------------------------------------------
colnames(surv_bdate)

# stats for table S1 0 and for descaling later # this is a reference for raw unique data per year and dt2
tmp_avg <- surv_bdate %>% # not inflated
    group_by(yr) %>% 
    summarise(mean_evi = mean(mismatch, na.rm = T), 
              sd_evi = sd(mismatch, na.rm = T),
              mean_fem = mean(fem, na.rm = T),
              sd_fem = sd(fem, na.rm = T),
              mean_femtm1 = mean(fem_tm1, na.rm = T),
              sd_femtm1 = sd(fem_tm1, na.rm = T),
              mean_masstm1= mean(fallMass_tm1, na.rm = T),
              sd_masstm1 = sd(fallMass_tm1, na.rm = T),
              mean_bdate = mean(birthdate, na.rm = T), 
              sd_bdate = sd(birthdate, na.rm = T), 
              mean_weanmass = mean(weanMass, na.rm = T), 
              sd_weanmass = sd(weanMass, na.rm = T), 
              medBD = median(birthdate, na.rm = T))


# here we define annual population mismatch (PM) as the arithmeic average of IM values each year
dt <- surv_bdate %>%
    dplyr::group_by(yr) %>% # la moyenne annuelle et on centre à cette moyenne annuelle
    dplyr::mutate(evi_pm= mean(mismatch, na.rm = T), evi_im = mismatch - evi_pm, # then IM is standardized relative to annual average... a within-year IM 
                  bd_pop= mean(timing, na.rm = T), bd_ind=timing-bd_pop) # the ind is exactly the same as evi_im - mismatch cancels out but this is a timing to pop rather than a mismatch to GU

plot(evi_pm~yr, data=dt)
plot(evi_im~yr, data=dt)



# detrend -----------------------------------------------------------------
# detrending based on lm # ATTENTION SOME ARE ANNUAL VALUES THUS DUPLICATED

dt1 = dt %>%
    mutate(bd.r = resid(lm(birthdate~as.numeric(as.character(yr)),na.action = na.exclude)),
           bd_ind.r = resid(lm(bd_ind~as.numeric(as.character(yr)),na.action = na.exclude)),
           # evi_pm.r = resid(lm(evi_pm~as.numeric(as.character(yr)),na.action = na.exclude)),
           evi_im.r = resid(lm(evi_im~as.numeric(as.character(yr)),na.action = na.exclude)),
           #  fem.r = resid(lm(fem ~ as.numeric(as.character(yr)),na.action = na.exclude)),
           # fem_tm1.r = resid(lm(fem_tm1 ~ as.numeric(as.character(yr)),na.action = na.exclude)),
           weanMass.r =resid(lm(weanMass ~ as.numeric(as.character(yr)),na.action = na.exclude)),
           mass_tm1.r =resid(lm(fallMass_tm1 ~ as.numeric(as.character(yr)),na.action = na.exclude))) %>% arrange(yr)

dt2 = dt %>%
    group_by(yr) %>% 
    summarise(evi_pm = unique(evi_pm), 
              fem = unique(fem), 
              fem_tm1 = unique(fem_tm1),
              bd_pop = unique(bd_pop)) %>% 
    na.omit() 


# detrend 

evi_pm.r= resid(lm(evi_pm~as.numeric(as.character(yr)),data=dt2)) # this is equivalent to bdate since annual effect is removed - mismatch cancels out 
fem.r=resid(lm(fem~as.numeric(as.character(yr)),data=dt2)) # this is equivalent to bdate since annual effect is removed - mismatch cancels out 
fem_tm1.r=resid(lm(fem_tm1~as.numeric(as.character(yr)),data=dt2)) # this is equivalent to bdate since annual effect is removed - mismatch cancels out 
bd_pop.r=resid(lm(bd_pop~as.numeric(as.character(yr)),data=dt2)) # this is equivalent to bdate since annual effect is removed - mismatch cancels out 

dt2 <- cbind(dt2, evi_pm.r, fem.r, fem_tm1.r, bd_pop.r)
colnames(dt2) <- c( "yr","evi_pm","fem","fem_tm1", "bd_pop", "evi_pm.r","fem.r","fem_tm1.r", "bd_pop.r")

dt1 = merge(dt1, dt2[, c("yr", "evi_pm.r","fem.r","fem_tm1.r", "bd_pop.r")], by.x="yr", by.y = "yr", all.x = T)
# clean dataset 
rm(male, RamMtnPop, surv, tmp, tmp2, tmp3, females, neonatal_dates,tmp_lb, surv_pop,popPheno, tidybdates, surv_pop, Ram_mass, pheno, neosurv, dt)
View(t1<-table(dt1$mom_id)) # 80 females total (NA removed in environment)
dt1$mom_id <-  droplevels(dt1$mom_id) # 

# keep existing females since droplevels doesn't work well 
dt1 <- dt1[dt1$mom_id%in% names(t1)[t1>=1],] # keep only those with 2 and more captures for mean centering ?? NO NEED ADDED =1
dt1$mom_id<- droplevels(dt1$mom_id) # 
View(mom_id<-table(dt1$mom_id))  # now 80 but some are NAs

# DO NOT remove NA OR LOOSE ALL NEONATAL
# dt1 <- dt1[!is.na(dt1$bd.r),]
# dt1 <- dt1[!is.na(dt1$EVI.r),]
# dt1 <- dt1[!is.na(dt1$fem.r),]
# dt1 <- dt1[!is.na(dt1$fem_tm1.r),]
# dt1 <- dt1[!is.na(dt1$weanMass.r),]
# dt1 <- dt1[!is.na(dt1$mass_tm1.r),]

# applies to all models
dt1$prs <- as.factor(dt1$prs)
dt1$yr<- as.factor(as.character(dt1$yr))


colnames(dt1)
cols.num <-c("bd.r", "evi_pm.r", "evi_im.r","fem.r",
             "fem_tm1.r","weanMass.r","mass_tm1.r")
dt1[cols.num] <- sapply(dt1[cols.num],as.numeric)
sapply(dt1, class)

fitness.data <- dt1

#save(fitness.data, file = "Cache/20210426fitnessData_centeredDetrended.RData")



# basic stats and verification  -------------------------------------------------------------
colnames(dt1)

hist(dt1$evi_pm) # OK
hist(dt1$evi_pm.r) # OK

mean(dt1$evi_pm.r) # there's a problem - too small value ?!! OK now =probably caused by regressing over duplicated years 
mean(dt1$evi_pm)
mean(dt1$evi_im)
mean(dt1$evi_im.r)

dt1$ID %>% 
    droplevels() %>% 
    table %>% 
    mean() # on average 3.5125 per female  - did not exclude any female , starting 2000
str(dt1$ID) # 69 niveaux 



View(t1<-table(dt1$ID)) # de 1 à 12 captures 
dt1$ID<- droplevels(dt1$ID) # n = 69 females


# compare to "Raw"
SURV20200608 <- read.csv("~/Documents/PhD/Analyses/Mismatch/Data/SURV20200608.csv")
raw_df <- SURV20200608

raw_df <- raw_df[raw_df$yr>1999,] # n= 714

lambs <- raw_df[raw_df$age==0,]
code0 <- raw_df[raw_df$code.sr %in% "0",] #n = 148
code1 <- raw_df[raw_df$code.sr %in% "1",] #n = 61

code2 <- raw_df[raw_df$code.sr %in% "2",] %>% select(c("yr", "ID", "wt12", "wt114")) %>%arrange(yr)#n = 38 
print(code2)
# yr  ID     wt12    wt114
# 1  2000 37T 50.89435 70.80708
# 2  2000 51T 55.54690 71.92865
# 3  2000 28U 57.87515 74.95261
# 4  2001 28W 59.46933 69.94026
# 5  2001 28T 57.06784 70.64924
# 6  2001 A44 60.99539 73.43043
# 7  2001 42U 67.47750 74.47740
# 8  2001 A50 58.87118 68.97804
# 9  2002 E13 42.09051 61.90885
# 10 2002 A35 56.19937 63.53771
# 11 2003 A44 63.96367 76.31414
# 12 2003 A35 53.58069 66.20738
# 13 2004 A44 64.32383 81.62945
# 14 2004 G14 60.48006 74.70594
# 15 2004  F1 58.36308 73.79375
# 16 2005  D8 56.83168 69.80055
# 17 2005 A35 60.07749 69.34564
# 18 2007 C14 51.33252 63.35377
# 19 2008 C14 48.59432 69.44156
# 20 2008 E12 61.38846 82.41512
# 21 2009  K1 53.18254 68.79174
# 22 2010 M28 42.02237 62.34650
# 23 2010  J5 59.40132 79.03256
# 24 2012  N4 50.87282 66.02334
# 25 2012  J5 63.86907 76.28546
# 26 2012  N6 64.33342       NA
# 27 2013  M7 58.86456 70.82768
# 28 2013 M15 57.45078 70.52497
# 29 2013  J5 62.76116 75.18145
# 30 2013  L2 62.26739 74.43463
# 31 2013 P13 57.54786 69.43493
# 32 2015 P11 52.47062 69.55153
# 33 2015  O6 57.53812 74.88706
# 34 2016 S10 46.24042 60.61860
# 35 2017  U7 51.27235 72.77973
# 36 2018 M28 52.61488 65.15623
# 37 2018  L2 62.84172 77.41216
# 38 2018  P3 52.61533 68.18099

colnames(fitness.data)
na.fit <- fitness.data[fitness.data$weanMass %in% NA & fitness.data$neonatal %in% "1",] %>% select(c("yr", "mom_id", "weanMass", "fallMass", "fallMass_tm1")) %>% arrange(yr)#n = 38 
print(na.fit)
# yr    ID    weanMass fallMass fallMass_tm1
# <fct> <fct>    <dbl>    <dbl>        <dbl>
# 1 2000  28U         NA     75.0         74.7
# 2 2000  37T         NA     70.8         68.5
# 3 2000  51T         NA     71.9         69.7
# 4 2001  28T         NA     70.6         71.8
# 5 2001  A33         NA     70.4         73.2
# 6 2001  A44         NA     73.4         76.3
# 7 2002  A35         NA     63.5         67.9
# 8 2003  A35         NA     66.2         63.5
# 9 2003  A44         NA     76.3         74.6
# 10 2004  E12         NA     74.5         72.0
# # … with 27 more rows










# load fitness data and parturition date for later --------------------------------------
load("cache/20210426fitnessData_centeredDetrended.RData")
fitness.data$date <- as.Date(fitness.data$date, format = "%Y-%m-%d")
fitness.data <- fitness.data %>% rename(year = yr)


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
pheno_ram <- read_csv("data/mine/pheno_ram.csv")
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

# do some quick verification and compare to plasticity data here 
# looks fine
# end verification 


# create two df : one pop and one ind -------------------------------------
# merge data to birthdates 
fall.pheno.pop <- fitness.data %>% group_by(year) %>% summarise(medBD = median(birthdate), 
                                                                fall_prec = unique(prec_Windowed2),
                                                                fall_temp = unique(temp_Windowed2))


# include fall pheno with time lag and detrend fall -----------------------------------------------------------
pheno_ram <- read.csv("data/mine/pheno_ram.csv")

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

#  re-merge with time step 
dt3 <- merge(fall.pheno.pop, # these are ok for time steps - done in sliding window 
             dt2, # data with time lag of 1 yr
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
# file = "cache/20210426revisedDf.RData")




# prepare export for Dryad ------------------------------------------------
# load("cache/20210426revisedDf.RData") 

# create dat.trend dataframe for temporal trends from raw R object
dat.trend <- fall.pheno.ind[, c(1:8)]%>% left_join(spring.pheno.ind)
dat.trend$mom_id <- droplevels(dat.trend$mom_id) # 75 levels

# REVISIONS ADDED EVI_TM1

# rename columns
colnames(dat.trend) <- c("year","mom_id","mismatch","birthdate","medBD","fall_prec","fall_temp",  
                         "evi_tm1","spring_prec","spring_temp","evi_up")

#  remove ID for Dryad storage
dat.trend$mom_id <- as.numeric(dat.trend$mom_id)
dat.trend$mom_id  <- as.factor(as.character(dat.trend$mom_id))
head(dat.trend)
#   year mom_id mismatch birthdate medBD fall_prec fall_temp evi_tm1 spring_prec spring_temp evi_up
# 1 2001      1        2       160 160.5 0.4923077  2.879268     279      1.5056    4.913115    158
# 2 2001      2       39       197 160.5 0.4923077  2.879268     279      1.5056    4.913115    158
# 3 2001      3       -8       150 160.5 0.4923077  2.879268     279      1.5056    4.913115    158
# 4 2001      4        3       161 160.5 0.4923077  2.879268     279      1.5056    4.913115    158
# 5 2001      5      -10       148 160.5 0.4923077  2.879268     279      1.5056    4.913115    158
# 6 2001      6        9       167 160.5 0.4923077  2.879268     279      1.5056    4.913115    158


fitness.data$mom_id <- droplevels(fitness.data$mom_id) # 80 levels
fitness.data$mom_id <- as.numeric(fitness.data$mom_id)
fitness.data$mom_id  <- as.factor(as.character(fitness.data$mom_id))
head(fitness.data)
#   year mom_id neonatal weanMass lambSurv_t1     newRS       prs  sex fallMass fallMass_tm1 pred fem pred_tm1 fem_tm1 lamb_id birthdate
# 1 2000      1        1 36.13143           1    weaned notWeaned male 71.92339     72.57620    1  37        1      45      F3       146
# 2 2000      2        1 27.50417           1    weaned    noLamb male 59.66906     58.35663    1  37        1      45     F12       160
# 3 2000      3        0       NA           0 notWeaned    weaned <NA> 71.97534     70.10714    1  37        1      45    <NA>       150
# 4 2000      4        1 30.48450           0    weaned notWeaned male 71.81758     67.93633    1  37        1      45      F2       154
# 5 2000      5        1       NA           0 notWeaned notWeaned <NA> 74.95261     74.70597    1  37        1      45    <NA>       150
# 6 2000      6        0       NA           0 notWeaned notWeaned <NA> 73.14261     74.43051    1  37        1      45    <NA>       160


fitness.data$lamb_id <- droplevels(fitness.data$lamb_id) # 226 levels
fitness.data$lamb_id <- as.numeric(fitness.data$lamb_id)
fitness.data$lamb_id  <- as.factor(as.character(fitness.data$lamb_id))


# save as csv for Dryad and analyses

# CHANGE TO FITNESS DATA

# write.csv2(dat.trend, file = "data/mine/trends_data_noID.csv", row.names = FALSE)
# write.csv2(fitness.data, file = "data/mine/mismatch_data_noID.csv", row.names = FALSE)

