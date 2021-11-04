source('Make_models.R')
source('Fit_models.R')

Create_CV_lists<-function(df_standard){
  CV_list<-list()
  CV_all<-list()
  CV_all$Model_0<-lapply(c(1:CV), function(x) CV_list[[x]]<-list(fitted=df_standard,
                                                                 simulation_parameters=list(U1=matrix(c(2,1,1,3),2, byrow=TRUE),
                                                                                            U2=matrix(c(3,1.5,1.5,4),2, byrow=TRUE))))
  CV_all$Model_1A<-lapply(c(1:CV), function(x) CV_list[[x]]<-list(fitted=df_standard,
                                                             simulation_parameters=list(Sigma=rbind(cbind(matrix(c(5,2,2,4), 2), 
                                                                                                          matrix(c(-1, -0.5, -0.5, -1.5), 2)),
                                                                                                    cbind(t(matrix(c(-1, -0.5, -0.5, -1.5), 2)), 
                                                                                                          matrix(c(4,1,1,4),2))))))
  CV_all$Model_2A<-lapply(c(1:CV), function(x) CV_list[[x]]<-list(fitted=df_standard,simulation_parameters=list(U=matrix(c(2,-1,
                                                                                                                         -1,3),2, byrow=TRUE))))
  CV_all$Model_2B<-lapply(c(1:CV), function(x) CV_list[[x]]<-list(fitted=df_standard,
                                                             simulation_parameters=list(U=matrix(c(2,-1,-1,3),2, byrow=TRUE),
                                                                                        Sigma=rbind(cbind(matrix(c(5,2,2,4), 2), 
                                                                                                          matrix(c(-1, -0.5, -0.5, -1.5), 2)),
                                                                                                    cbind(t(matrix(c(-1, -0.5, -0.5, -1.5), 2)), 
                                                                                                          matrix(c(4,1,1,4),2))))))
  CV_all$Model_2C1<-lapply(c(1:CV), function(x) CV_list[[x]]<-list(fitted=df_standard,
                                                              simulation_parameters=list(U=rbind(cbind(matrix(c(2,-1,-1,3),2, byrow=TRUE), 
                                                                                                       matrix(c(0,0,0,0),2, byrow=TRUE)),
                                                                                                 cbind(t(matrix(c(0,0,0,0),2, byrow=TRUE)), 
                                                                                                       matrix(c(3,1.5,1.5,4),2, byrow=TRUE))))))
  CV_all$Model_2C2<-lapply(c(1:CV), function(x) CV_list[[x]]<-list(fitted=df_standard,
                                                              simulation_parameters=list(U=rbind(cbind(matrix(c(2,-1,-1,3),2, byrow=TRUE), 
                                                                                                       matrix(c(1.25,1,-0.5,0.3),2, byrow=TRUE)),
                                                                                                 cbind(t(matrix(c(1.25,1,-0.5,0.3),2, byrow=TRUE)), 
                                                                                                       matrix(c(3,1.5,1.5,4),2, byrow=TRUE))))))
  CV_all$Model_3A1<-lapply(c(1:CV), function(x) CV_list[[x]]<-list(fitted=df_standard,
                                                             simulation_parameters=list(U1=matrix(c(2,0,0,3),2, byrow=TRUE),
                                                                                        U2=matrix(c(3,1.5,1.5,4),2, byrow=TRUE),
                                                                                        gamma=1.2)))
  CV_all$Model_3A2<-lapply(c(1:CV), function(x) CV_list[[x]]<-list(fitted=df_standard,
                                                                   simulation_parameters=list(U1=matrix(c(2,1.5,1.5,3),2, byrow=TRUE),
                                                                                              U2=matrix(c(3,1.5,1.5,4),2, byrow=TRUE),
                                                                                              gamma=1.2)))
  CV_all$Model_3B1<-lapply(c(1:CV), function(x) CV_list[[x]]<-list(fitted=df_standard,
                                                              simulation_parameters=list(U1=matrix(c(2,0,0,3),2, byrow=TRUE),
                                                                                         U2=matrix(c(3,1.5,1.5,4),2, byrow=TRUE),
                                                                                         gamma1=1.2,gamma2=-1.3)))
  return(CV_all)
}


#Overall parameters
CV<-5
true_betas<-c(2,4,2.5,3, 1.5, 3.5)
N<-75
n<-7
### Probability p that a measurement is observed
p1<-0.7
### probability p2 that covariate x is measured
p2<-0.9

Fit_function_list<-list(Model_0=fit_model_0, Model_1A=fit_model_1A, Model_2A=fit_model_2A, Model_2B=fit_model_2B,
                        Model_2C1=fit_model_2C1, Model_2C2=fit_model_2C2, Model_3A1=fit_model_3A1, Model_3A2=fit_model_3A2, Model_3B1=fit_model_3B1)

Make_function_list<-list(Model_0=make_model_0, Model_1A=make_model_1A, Model_2A=make_model_2A, Model_2B=make_model_2B,
                         Model_2C1=make_model_2C, Model_2C2=make_model_2C, Model_3A1=make_model_3A, Model_3A2=make_model_3A, Model_3B1=make_model_3B)

outcomes<-data.frame(true=c(true_betas, rep(NA, 6)),
                     Model_0=c(rep(NA, 12)),
                     Model_1A=c(rep(NA, 12)),
                     Model_2A=c(rep(NA, 12)),
                     Model_2B=c(rep(NA, 12)),
                     Model_2C1=c(rep(NA, 12)),
                     Model_2C2=c(rep(NA, 12)),
                     Model_3A1=c(rep(NA, 12)),
                     Model_3A2=c(rep(NA, 12)),
                     Model_3B1=c(rep(NA, 12)))
rownames(outcomes)<-c( 'beta_0^1', 'beta_x^1','beta_t^1', 'beta_0^2', 'beta_x^2','beta_t^2', 'MLIK','DIC_approx', 'DIC', 'WAIC', 'PIT', 'CPO')
df_standard=list(outcomes=outcomes,
                 Model_0_parameters=list(),
                 Model_1A_parameters=list(),
                 Model_2A_parameters=list(),
                 Model_2B_parameters=list(),
                 Model_2C1_parameters=list(),
                 Model_2C2_parameters=list(),
                 Model_3A1_parameters=list(),
                 Model_3B1_parameters=list())
CV_all<-Create_CV_lists(df_standard)
Model_comparison_big<-list(Model_0=list(outcomes=outcomes,
                                          simulation_parameters=CV_all$Model_0[[1]]$simulation_parameters),
                           Model_1A=list(outcomes=outcomes,
                                          simulation_parameters=CV_all$Model_1A[[1]]$simulation_parameters),
                           Model_2A=list(outcomes=outcomes,
                                          simulation_parameters=CV_all$Model_2A[[1]]$simulation_parameters),
                           Model_2B=list(outcomes=outcomes,
                                          simulation_parameters=CV_all$Model_2B[[1]]$simulation_parameters),
                           Model_2C1=list(outcomes=outcomes,
                                          simulation_parameters=CV_all$Model_2C1[[1]]$simulation_parameters),
                           Model_2C2=list(outcomes=outcomes,
                                          simulation_parameters=CV_all$Model_2C2[[1]]$simulation_parameters),
                           Model_3A1=list(outcomes=outcomes,
                                        simulation_parameters=CV_all$Model_3A1[[1]]$simulation_parameters),
                           Model_3A2=list(outcomes=outcomes,
                                          simulation_parameters=CV_all$Model_3A2[[1]]$simulation_parameters),
                           Model_3B1=list(outcomes=outcomes,
                                          simulation_parameters=CV_all$Model_3B1[[1]]$simulation_parameters))
  
# data<-Make_function_list[[1]](true_betas,CV_all[[1]][[5]]$simulation_parameters , N, n, p1, p2, 5)
# temp<-Fit_function_list[[3]](data$data, N, n)
for (model in c(1,3,5:9)){
  cat('\n Generated Model', model, '\n')
  for (cv in 1:CV){
    cat('In CV', cv, '/', CV, '\n')
    data<-Make_function_list[[model]](true_betas,CV_all[[model]][[cv]]$simulation_parameters , N, n, p1, p2, cv)
    cat('Have made data! \n')
    CV_all[[model]][[cv]]$simulation_parameters<-data[-1]
    for(fit in c(1,3,5:9)){
      cat('Fitting Model', fit, '\n')
      temp<-Fit_function_list[[fit]](data$data, N, n)
      CV_all[[model]][[cv]]$fitted$outcomes[,(fit+1)]<-unlist(temp[c(1:7)])
      CV_all[[model]][[cv]]$fitted[[(fit+1)]]<-temp[-c(1:7)]
    }
  }
}
# for (model in c(9)){
#   temp<-0
#   for (cv in c(1,2,3,5)){
#     temp=temp+CV_all[[model]][[cv]]$fitted$outcomes$Model_3A2
#   }
#   Model_comparison_big[[model]]$outcomes$Model_3A2<-temp/CV
# }
# CV_all[[1]][[5]]$fitted$outcomes
# Model_comparison_big$Model_3A1$outcomes
# Model_comparison_big$Model_3B1$outcomes
# 
# save(Model_comparison_big, file='04_11_big.RData')
# save(CV_all, file='04_11_all.RData')
# load('04_11_big.RData')
# Model_comparison_big
