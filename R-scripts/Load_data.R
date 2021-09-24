library('INLA')
library("MASS")
library('JSM')
library('nlme')
library(tidyr)
library(Rcpp)
library(reshape2)
library(lme4)
library(nlme)
library(MCMCglmm)
library(lme4)
library(ggplot2)
library(reshape2)
library(MCMCvis)
theme_set(theme_bw())
library(xtable)

################################################################################
#Loading Pediatric Pain data
#https://robweiss.faculty.biostat.ucla.edu/book-data-sets
################################################################################

data<-read.csv('pain.txt', sep='')
names(data)<-c('id', 'ses', 'age_years', 'age_months', 'sex', 'treatment', 'coping_style', 
               ' helpful', 'easy', 'tol11', 'tol12', 'rat11', 'rat12', 'tol13', 'tol14',
               'rat13', 'rat14', 'senscat1', 'cpantic1')
data[data==-999]<-NA
data$sex<-factor(data$sex,levels=c(1,2), labels=c('male', 'female'))
data$treatment<-factor(data$treatment,levels=c(1,2,3), labels=c('attend', 'distract', 'none'))
data$coping_style<-factor(data$coping_style,levels=c(1,2), labels=c('attend', 'distract'))
data_long<-reshape(data, varying=c(cbind('rat11','tol11', 'rat12','tol12', 'rat13', 'tol13','rat14', 'tol14'  )), 
                   timevar='Period', times=c('1', '2', '3', '4'), idvar='id' ,direction='long', 
                   v.names=c('pain_tolerance', 'pain_rating'))
data_long$pain_tolerance<-log(data_long$pain_tolerance)
data_long$Period<-as.factor(data_long$Period)
data_long<-data_long[with(data_long,order(id)),]
data_long_melted<-melt(data_long, measure.vars=c('pain_tolerance', 'pain_rating'))
data_long_melted<-data_long_melted[with(data_long_melted,order(id)),]
data_INLA<-data_long[, c(1,2,5, 12)]
y<-data_long[,c(13,14)]
nLong<-nrow(data_INLA)