if(!require("xtable")){install.packages("xtable"); library("xtable")}
if(!require("ggpubr")){install.packages("ggpubr"); library("ggpubr")}
if(!require("ggplot2")){install.packages("ggplot2"); library("ggplot2")}

load(file="comparison_all.RData")

################################################################################
#Extracting GOF values: Tables 6.2 & 6.4 in thesis
################################################################################
#Producing Table 6.2
GOF_1<-t(comparison_all$Cytokine1$CV1$results[,-c(2,8,9)])
colnames(GOF_1)<-c('LMM', 'JMM', 'JSM')
mdat <- matrix(c(rep(0,(3*3)),rep(4,3), rep(0, 3), rep(2,3)),
               nrow = 6, ncol=3, byrow=TRUE)
xtable(GOF_1, digits=cbind(1,mdat))

#Producing Table 6.4
GOF<-t(rbind(comparison_all$Cytokine2$CV1$results[,-c(2,8,9)], 
             comparison_all$Cytokine3$CV1$results[,-c(2,8,9)]))
colnames(GOF)<-rep(c('LMM', 'JMM', 'JSM'), 2)
mdat <- matrix(c(rep(0,(3*6)),rep(4,6), rep(0, 6), rep(2,6)),
               nrow = 6, ncol=6, byrow=TRUE)
xtable(GOF, digits=cbind(1,mdat))
