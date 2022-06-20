source('Read_data.R') #Read all simulation results

#################################################################################
#Figures 5.3-5.7 are produced here (pages 34-37)
#################################################################################

data<-read_data()

obtain_GOF<-function(list){
  CV<-250
  N_s<-c(10, 25, 50, 75, 100, 150, 250, 500, 750)
  for (N in N_s){
    if (N %in% c(250, 500, 750)){
      CV<-100
    }
    temp<-0
    counter<-0
    for (cv in 1:CV){
      counter<-counter+!is.na(list[[paste0('N', N)]][[paste0('CV',cv)]][['GOF']][1,])
      temp=temp+replace(list[[paste0('N', N)]][[paste0('CV',cv)]][['GOF']],
                        is.na(list[[paste0('N', N)]][[paste0('CV',cv)]][['GOF']]),
                        0)
    }
    list[[paste0('Total GOF:', N)]]<-t(t(temp)/counter)
  }
  return(list)
}

#################################################################################
#Plot MLIK results: Figure 5.3 on page 34 of thesis
#################################################################################
plot_Mlik<-function(list,indexes, name){
  N_s<-c(10, 25, 50, 75, 100, 150, 250, 500, 750)
  plotting_simulation<-as.data.frame(t(do.call(cbind,list[grepl('GOF', names(list))])))
  names(plotting_simulation)<-c('Mlik', 'dont care', 'DIC', 'WAIC', 'PIT', 'CPO', 'MSE on training set',
                                'MSE on subsequent measurements of subjects in training set', 'MSE test')
  plotting_simulation$N<-rep(N_s, each=3)
  plotting_simulation$Model<-factor(rep(c('LMM', 'JMM', 'JSM'), length(N_s)))
  plotting_simulation<-plotting_simulation[,c(1,10,11)]
  plotting_simulation<-rbind(plotting_simulation,data.frame(Mlik=100*(c(plotting_simulation$Mlik[seq(2,27, by=3)]-plotting_simulation$Mlik[seq(3,28, by=3)])-100),
                                                            N=N_s, Model=rep(c('Bayes Factor (JMM vs JSM)'), 9)))
  ggplot(data=plotting_simulation, aes(x=N, y=Mlik, group=Model, linetype=Model))+geom_line(aes(color=Model), size=1)+
    labs(x = "N", y='Marginal Likelihood')+theme_bw()+ theme(legend.text=element_text(size=16), text = element_text(size=12))+
    scale_y_continuous(name='Marginal Likelihood', sec.axis = sec_axis(~.*(1/100)+100, name='Bayes Factor K'))
}
p_0<-plot_Mlik(obtain_GOF(data$Model_0_all), c(1), c('Mlik'))
p_0<-p_0+ggtitle("Data simulated according to LMM")+theme(legend.title=element_blank())
p_2<-plot_Mlik(obtain_GOF(data$Model_2_all), c(1), c('Information Criteria'))
p_2<-p_2+ggtitle("Data simulated according to JMM")+theme(legend.title=element_blank())
p_3<-plot_Mlik(obtain_GOF(data$Model_3_all), c(1), c('Information Criteria'))
p_3<-p_3+ggtitle("Data simulated according to JSM")+theme(legend.title=element_blank())
ggarrange(p_0,p_2, p_3, ncol=3, common.legend = TRUE, legend='bottom')





#################################################################################
#Plot PIT+CPO results: Figure 5.4 page 35 thesis
#################################################################################
plot_CPO_PIT<-function(list,indexes, name){
  N_s<-c(10, 25, 50, 75, 100, 150, 250, 500, 750)
  plotting_simulation<-as.data.frame(t(do.call(cbind,list[grepl('GOF', names(list))])))
  names(plotting_simulation)<-c('Mlik', 'dont care', 'DIC', 'WAIC', 'PIT', 'CPO', 'MSE train', 'MSE same', 'MSE test')
  plotting_simulation$N<-rep(N_s, each=3)
  plotting_simulation$PIT<-plotting_simulation$PIT*100000
  plotting_simulation$Model<-factor(rep(c('LMM', 'JMM', 'JSM'), length(N_s)))
  plotting_MSE<-melt(plotting_simulation[,c(indexes,10,11)],  id = c("N","Model"), variable.name='GOF')
  plotting_MSE$group<-rep(seq(0,30, by=3)[1:length(indexes)], each=27)+rep(rep(c(1:3), 9),length(indexes))
  ggplot(data=plotting_MSE, aes(x=N, y=value, group=group, linetype=GOF))+geom_line(aes(color=Model), size=1)+
    geom_point(aes(shape=Model, color=Model), size=2)+
    scale_y_continuous(name='CPO', sec.axis = sec_axis(~.*(1/100000), name='PIT'))+theme_bw()+
    theme(legend.text=element_text(size=16), text = element_text(size=12))
}

p_0<-plot_CPO_PIT(obtain_GOF(data$Model_0_all), c(5,6), c('Information Criteria'))
p_0<-p_0+ggtitle("Data simulated according to LMM")
p_2<-plot_CPO_PIT(obtain_GOF(data$Model_2_all), c(5,6), c('Information Criteria'))
p_2<-p_2+ggtitle("Data simulated according to JMM")
p_3<-plot_CPO_PIT(obtain_GOF(data$Model_3_all), c(5,6), c('Information Criteria'))
p_3<-p_3+ggtitle("Data simulated according to JSM")
ggarrange(p_0,p_2, p_3,  as_ggplot(cowplot::get_legend(p_2)), common.legend = TRUE, legend="none")




#################################################################################
##plotting MSE train+same results: Figure 5.5 page 36 in thesis
#################################################################################
plot_MSE_train<-function(list,indexes, name){
  N_s<-c(10, 25, 50, 75, 100, 150, 250, 500, 750)
  plotting_simulation<-as.data.frame(t(do.call(cbind,list[grepl('GOF', names(list))])))
  names(plotting_simulation)<-c('Mlik', 'dont care', 'DIC', 'WAIC', 'PIT', 'CPO', 'MSE on training set', 'MSE on subsequent measurements of subjects in training set', 'MSE test')
  plotting_simulation$N<-rep(N_s, each=3)
  plotting_simulation$Model<-factor(rep(c('LMM', 'JMM', 'JSM'), length(N_s)))
  plotting_MSE<-melt(plotting_simulation[,c(indexes,10,11)],  id = c("N","Model"), variable.name='MSE')
  plotting_MSE$group<-rep(seq(0,30, by=3)[1:length(indexes)], each=27)+rep(rep(c(1:3), 9),length(indexes))
  ggplot(data=plotting_MSE, aes(x=N, y=value, group=group, linetype=MSE))+geom_line(aes(color=Model), size=1)+
    labs(x = "N", y='MSE')+theme_bw()
}

p_0<-plot_MSE_train(obtain_GOF(data$Model_0_all), c(7,8), c('MSE'))
p_0<-p_0+ggtitle("Data simulated according to LMM")
p_2<-plot_MSE_train(obtain_GOF(data$Model_2_all), c(7,8), c('MSE'))
p_2<-p_2+ggtitle("Data simulated according to JMM")
p_3<-plot_MSE_train(obtain_GOF(data$Model_3_all), c(7,8), c('MSE'))
p_3<-p_3+ggtitle("Data simulated according to JSM")
ggarrange(p_0,p_2, p_3,  as_ggplot(cowplot::get_legend(p_2)), common.legend = TRUE, legend="none")



#################################################################################
#Plot MSE test results: Figure 5.6 page 36 in thesis
#################################################################################
plot_MSE_test<-function(list,indexes, name){
  N_s<-c(10, 25, 50, 75, 100, 150, 250, 500, 750)
  plotting_simulation<-as.data.frame(t(do.call(cbind,list[grepl('GOF', names(list))])))
  names(plotting_simulation)<-c('Mlik', 'dont care', 'DIC', 'WAIC', 'PIT', 'CPO', 'MSE on training set', 'MSE on subsequent measurements of subjects in training set', 'MSE test')
  plotting_simulation$N<-rep(N_s, each=3)
  plotting_simulation$Model<-factor(rep(c('LMM', 'JMM', 'JSM'), length(N_s)))
  plotting_MSE<-melt(plotting_simulation[,c(indexes,10,11)],  id = c("N","Model"), variable.name='MSE')
  plotting_MSE$group<-rep(seq(0,30, by=3)[1:length(indexes)], each=27)+rep(rep(c(1:3), 9),length(indexes))
  ggplot(data=plotting_MSE, aes(x=N, y=value, group=group, linetype=Model))+geom_line(aes(color=Model), size=1)+
    labs(x = "N", y='MSE on test set')+theme_bw()+geom_point(aes(shape=Model, color=Model), size=2)+
    theme(legend.text=element_text(size=16), text = element_text(size=12))
}

p_0<-plot_MSE_test(obtain_GOF(data$Model_0_all), c(9), c('Mlik'))
p_0<-p_0+ggtitle("Data simulated according to LMM")+theme(legend.title=element_blank())
p_2<-plot_MSE_test(obtain_GOF(data$Model_2_all), c(9), c('Information Criteria'))
p_2<-p_2+ggtitle("Data simulated according to JMM")+theme(legend.title=element_blank())
p_3<-plot_MSE_test(obtain_GOF(data$Model_3_all), c(9), c('Information Criteria'))
p_3<-p_3+ggtitle("Data simulated according to JSM")+theme(legend.title=element_blank())
ggarrange(p_0,p_2, p_3, ncol=3, common.legend = TRUE, legend='bottom')


#################################################################################
#Plot DIC results: Figure 5.7 page 37 in thesis
#################################################################################
plot_DIC<-function(list,indexes, name){
  N_s<-c(10, 25, 50, 75, 100, 150, 250, 500, 750)
  plotting_simulation<-as.data.frame(t(do.call(cbind,list[grepl('GOF', names(list))])))
  names(plotting_simulation)<-c('Mlik', 'dont care', 'DIC', 'WAIC', 'PIT', 'CPO', 'MSE on training set',
                                'MSE on subsequent measurements of subjects in training set', 'MSE test')
  plotting_simulation$N<-rep(N_s, each=3)
  plotting_simulation$Model<-factor(rep(c('LMM', 'JMM', 'JSM'), length(N_s)))
  plotting_MSE<-melt(plotting_simulation[,c(indexes,10,11)],  id = c("N","Model"), variable.name='MSE')
  plotting_MSE$group<-rep(seq(0,30, by=3)[1:length(indexes)], each=27)+rep(rep(c(1:3), 9),length(indexes))
  ggplot(data=plotting_MSE, aes(x=N, y=log(value), group=group, linetype=Model))+geom_line(aes(color=Model), size=1)+
    labs(x = "N", y='log(DIC)')+theme_bw()+ theme(legend.text=element_text(size=16), text = element_text(size=12))
}
p_0<-plot_DIC(obtain_GOF(data$Model_0_all), c(3), c('Information Criteria'))
p_0<-p_0+ggtitle("Data simulated according to LMM")+theme(legend.title=element_blank())
p_2<-plot_DIC(obtain_GOF(data$Model_2_all), c(3), c('Information Criteria'))
p_2<-p_2+ggtitle("Data simulated according to JMM")+theme(legend.title=element_blank())
p_3<-plot_DIC(obtain_GOF(data$Model_3_all), c(3), c('Information Criteria'))
p_3<-p_3+ggtitle("Data simulated according to JSM")+theme(legend.title=element_blank())
ggarrange(p_0,p_2, p_3, ncol=3, common.legend = TRUE, legend='bottom')