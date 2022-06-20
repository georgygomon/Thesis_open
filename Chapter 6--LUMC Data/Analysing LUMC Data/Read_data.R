if(!require("xtable")){install.packages("xtable"); library("xtable")}
if(!require("reshape2")){install.packages("reshape2"); library("reshape2")}
if(!require("ggplot2")){install.packages("ggplot2"); library("ggplot2")}
if(!require("nlme")){install.packages("nlme"); library("nlme")}
if(!require("egg")){install.packages("egg"); library("egg")}
if(!require("mda")){install.packages("mda"); library("mda")}
if(!require("tidyverse")){install.packages("tidyverse"); library("tidyverse")}
#To install the INLA-package in R, you have to manually add the r-inla repository as it is not on CRAN.
#https://www.r-inla.org/download-install
library('INLA')


################################################################################
#Due to privacy concerns the data will not be publicly shared.
#Thus, it will not be possible to run the code presented in this section.
################################################################################


create_data<-function(){
  load('Data_Severity.rda')
  data<-data_share
  #Log transformation best for all cytokines
  data$Cy_1<-log(data$Cytokine_1)
  data$Cy_2<-log(data$Cytokine_2)
  data$Cy_3<-log(data$Cytokine_3)
  data<-data[!is.na(data$day_no),]
  length(unique(data$ID))
  data$ID<-rep(1:length(unique(data$ID)), table(data$ID))
  return(data)
}








  
