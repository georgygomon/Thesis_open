if(!require("ggplot2")){install.packages("ggplot2"); library("ggplot2")}
if(!require("xtable")){install.packages("xtable"); library("xtable")}
if(!require("reshape2")){install.packages("reshape2"); library("reshape2")}
if(!require("ggpubr")){install.packages("ggpubr"); library("ggpubr")}

################################################################################
#Reads all different simulation result files and combines them into 1 large list
################################################################################


read_data<-function(){
  load('All_simulation_data.RData')
  return(All_simulation_data)
}
