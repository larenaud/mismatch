# tidying data for script 01_analysesTempTrends # not on Dryad but on github


# load libraries
library(plyr)
library (dplyr) # very important
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


rm(list = ls())

# load R object with main files to be tidyied
load("~/Documents/PostDocI/Projects/mismatch/cache/20210421revisedDf.RData") 

# create dat.trend dataframe for temporal trends from raw R object
dat.trend <- fall.pheno.ind%>% left_join(spring.pheno.ind)
#dat.trend$mismatch <- dat.trend$birthdateJJ - dat.trend$evi_log_up_jul
# REVISIONS ADDED EVI_TM1


# rename columns
colnames(dat.trend) <- c("year","mother_id","birthdate","median_bd",
                         "fall_prec","fall_temp","evi_tm1","birthdateYMD","log_birthdate",
                         "spring_prec","spring_temp","evi_up")

# save as csv for Dryad and other analyses
write.csv2(dat.trend, file = "data/mine/trends_df.csv", row.names = FALSE)


