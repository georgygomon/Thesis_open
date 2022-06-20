source('Read_data.R') #Reads simulation data

#################################################################################
#Figures 5.1-5.2 as well as table 5.4 are produced here (pages 32-33)
#################################################################################

obtain_coefficient<-function(list){
  CV<-250
  N_s<-c(10, 25, 50, 75, 100, 150, 250, 500, 750)
  for (N in N_s){
    if (N %in% c(250, 500, 750)){
      CV<-100
    }
    b_x<-rep(0, CV)
    b_x_jmm<-rep(0, CV)
    b_x_jsm<-rep(0, CV)
    for (cv in 1:CV){
      b_x[cv]<-list[[paste0('N', N)]][[paste0('CV',cv)]][['Coefficients']][4,1]
      b_x_jsm[cv]<-list[[paste0('N', N)]][[paste0('CV',cv)]]$`Model 3`$model_specific$gamma$gamma
      b_x_jmm[cv]<-tryCatch(list[[paste0('N', N)]][[paste0('CV',cv)]]$`Model 2`$model_specific$cov_both$Cov[3,4]/
                              list[[paste0('N', N)]][[paste0('CV',cv)]]$`Model 2`$model_specific$cov_both$Cov[3,3], 
                            error=function(error_message) {
                              return(NA)
                            })
      
    }
    list[[paste0('Ass_coef', N)]]<-matrix(c(mean(b_x, na.rm=TRUE),sd(b_x, na.rm=TRUE),
                                            mean(b_x_jmm, na.rm=TRUE),sd(b_x_jmm, na.rm=TRUE),
                                            mean(b_x_jsm, na.rm=TRUE),sd(b_x_jsm, na.rm=TRUE)), nrow=3,byrow=TRUE)
    list[[paste0('Asso_coef_NA', N)]]<-matrix(c(sum(is.na((b_x))),sum(is.na((b_x_jmm))), sum(is.na((b_x_jsm)))))
    
  }
  return(list)
}

plot_coefficient<-function(list, actual){
  N_s<-c(10, 25, 50, 75, 100, 150, 250, 500, 750)
  plotting_simulation<-as.data.frame(do.call(rbind,list[grepl('Ass_coef', names(list))]))
  plotting_simulation$N<-rep(N_s, each=3)
  plotting_simulation$Model<-factor(rep(c('LMM', 'JMM', 'JSM'), length(N_s)), levels=c('LMM', 'JMM', 'JSM', 'Actual'))
  
  plotting_simulation<-rbind(plotting_simulation, data.frame(V1=rep(actual, 9), V2=NA, N=N_s, Model=c("Actual")))
  ggplot(data=plotting_simulation, aes(x=N, y=V1, group=Model, linetype=Model))+geom_line(aes(color=Model), size=1)+  
    geom_ribbon(aes(ymin=V1-V2, ymax=V1+V2, fill=Model), alpha=0.15)+
    labs(x = "N", y='Association coefficient')+theme_bw()+theme(legend.title=element_blank())
}


################################################################################
#Read data
################################################################################
data<-read_data()

################################################################################
#Plot association coefficient with time: Figure 5.1 in thesis (page 32)
################################################################################
load('association_coefficients_time.RData')
p0<-ggplot(data=association_coefficients_time$LMM, aes(x=t, y=beta, color=Model, linetype=Model))+geom_line(size=1.25)+
  labs(y= "Association coefficient", x = "Time")+ theme_bw()+ggtitle("Data simulated according to LMM")+theme(legend.title=element_blank())
p2<-ggplot(data=association_coefficients_time$JMM, aes(x=t, y=beta, color=Model, linetype=Model))+geom_line(size=1.25)+
  labs(y= "Association coefficient", x = "Time")+ theme_bw()+ggtitle("Data simulated according to JMM")+theme(legend.title=element_blank())
p3<-ggplot(data=association_coefficients_time$JSM, aes(x=t, y=beta, color=Model, linetype=Model))+geom_line(size=1.25)+labs(y= "Association coefficient", x = "Time")+
  theme_bw()+ggtitle("Data simulated according to JSM")+theme(legend.title=element_blank())


ggarrange(p0,p2, p3,  common.legend = TRUE, nrow=1, legend='bottom')

################################################################################
#Plot association coefficient: Figure 5.2 in thesis (page 32)
################################################################################
# # LMM & JSM
# Actual association coefficient as time to infinity=1.2
# # JMM
# Actual association coefficient as time to infinity=0.87
p_0<-plot_coefficient(obtain_coefficient(data$Model_0_all), 1.2)
p_0<-p_0+ggtitle("Data simulated from LMM")
p_2<-plot_coefficient(obtain_coefficient(data$Model_2_all), 0.87)
p_2<-p_2+ggtitle("Data simulated from JMM")
p_3<-plot_coefficient(obtain_coefficient(data$Model_3_all), 1.2)
p_3<-p_3+ggtitle("Data simulated from JSM")
ggarrange(p_0,p_2, p_3,  as_ggplot(cowplot::get_legend(p_2)), common.legend = TRUE, legend="none")

################################################################################
#Construct table 5.4 (page 33): Number of times JMM elements are not fit
################################################################################
failed_cov<-rbind(as.data.frame(do.call(cbind,obtain_coefficient(data$Model_0_all)
                                        [grepl('Asso_coef', names(obtain_coefficient(data$Model_0_all)))])),
                  as.data.frame(do.call(cbind,obtain_coefficient(data$Model_2_all)
                                        [grepl('Asso_coef', names(obtain_coefficient(data$Model_2_all)))])),
                  as.data.frame(do.call(cbind,obtain_coefficient(data$Model_3_all)
                                        [grepl('Asso_coef', names(obtain_coefficient(data$Model_3_all)))])))

failed_cov<-t(t(failed_cov)/c(rep(250, 6), rep(100,3)))[c(2,5,8),]*100
colnames(failed_cov)<-c('N=10', 'N=25', 'N=50', 'N=75', 'N=100', 'N=150', 'N=250', 'N=500', 'N=750')
rownames(failed_cov)<-c('LMM', 'JMM', 'JSM')
xtable(failed_cov, digit=0)


