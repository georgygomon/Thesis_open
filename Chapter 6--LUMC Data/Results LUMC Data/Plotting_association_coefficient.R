if(!require("xtable")){install.packages("xtable"); library("xtable")}
if(!require("ggpubr")){install.packages("ggpubr"); library("ggpubr")}
if(!require("ggplot2")){install.packages("ggplot2"); library("ggplot2")}

################################################################################
#This file creates figures 6.3 & 6.4 in the thesis (pages 47-48)
################################################################################

#Function which introduces breakpoint
b1<-function(x, bp) ifelse(x<bp, 0, x-bp)

#Function that calculates the association coefficient for the JMM in time
beta_jmm<-function(t, cov_mat){
  (cov_mat$Cov_0$Cov[1,2]+t^2*cov_mat$Cov_t$Cov[1,2]+b1(t, 28)^2*cov_mat$Cov_t28$Cov[1,2])/
    (cov_mat$Cov_0$Cov[1,1]+t^2*cov_mat$Cov_t$Cov[1,1]+b1(t, 28)^2*cov_mat$Cov_t28$Cov[1,1]+cov_mat$Cov_error[1,1])
}

#Function that calculates the association coefficient for the JSM in time
beta_jsm<-function(t, cov_mat, gamma){
  gamma$mean*(1-cov_mat$Cov_error[1,1]/(cov_mat$Cov_C$Cov[1,1]+t^2*cov_mat$Cov_C$Cov[2,2]+b1(t, 28)^2*cov_mat$Cov_C$Cov[3,3]+cov_mat$Cov_error[1,1]+
                                          2*t*cov_mat$Cov_C$Cov[1,2]+2*b1(t, 28)*cov_mat$Cov_C$Cov[1,3]+2*t*b1(t, 28)*cov_mat$Cov_C$Cov[2,3]))
}



plot_association<-function(list, t){
  coef_time<-data.frame(time=rep(t,3),
                        beta=c(beta_jmm(t, list$`Model 2`$model_specific$cov_mat),
                               beta_jsm(t, list$`Model 3_0`$model_specific$cov_mat, list$`Model 3_0`$model_specific$gamma),
                               rep(list$`Model 0`$coefficients[2,1], length(t))),
                        l=c(list$`Model 2`$model_specific$association_time[,1],
                            list$`Model 3_0`$model_specific$coef_time[,1],
                            rep(list$`Model 0`$coefficients[2,2], length(t))),
                        u=c(list$`Model 2`$model_specific$association_time[,3],
                            list$`Model 3_0`$model_specific$coef_time[,3],
                            rep(list$`Model 0`$coefficients[2,3], length(t))),
                        Model=rep(c('JMM', 'JSM', 'LMM'), each=length(t)))
  ggplot(data=coef_time, aes(x=time, y=beta, color=Model, linetype=Model))+geom_line(aes(color=Model), size=1)+
    geom_ribbon(aes(ymin=l, ymax=u, fill=Model), alpha=0.15,colour = NA)+theme_bw()+
    labs(x = "time", y='Association coefficient')+theme(legend.text=element_text(size=16), text = element_text(size=14))
}

load(file="comparison_all.RData")

p_0<-plot_association(comparison_all$Cytokine1$CV1, t=c(0:50))
p_0<-p_0+ggtitle("Association coefficient for Cytokine 1")
p_0
p_2<-plot_association(comparison_all$Cytokine2$CV1, t=c(0:50))
p_2<-p_2+ggtitle("Association coefficient for Cytokine 2")+theme(legend.title=element_blank())
p_3<-plot_association(comparison_all$Cytokine3$CV1, t=c(0:50))
p_3<-p_3+ggtitle("Association coefficient for Cytokine 3")
ggarrange(p_2, p_3)
