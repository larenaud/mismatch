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

getwd()
rm(list = ls())

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
# write.csv2(dat.trend, file = "data/mine/trends_data_noID.csv", row.names = FALSE)
# write.csv2(fitness.data, file = "data/mine/mismatch_data_noID.csv", row.names = FALSE)

