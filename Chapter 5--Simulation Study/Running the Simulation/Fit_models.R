calculate_MSE<-function(y, y_hat){
  return(mean((y-y_hat)^2))
}

calculate_all_MSE<-function(fitted_model, data, indexes, N, n){
  MSE_train<-calculate_MSE(fitted_model$summary.fitted.values[c(indexes$train+N*n),1],
                 data$y[indexes$train])
  MSE_test_same<-calculate_MSE(fitted_model$summary.fitted.values[c(indexes$test_same+N*n),1],
                       data$y[indexes$test_same])
  MSE_test_others<-calculate_MSE(fitted_model$summary.fitted.values[c(indexes$test_others+N*n),1],
                     data$y[indexes$test_others])
  return(list(MSE_train=MSE_train, MSE_test_same=MSE_test_same, 
              MSE_test_others=MSE_test_others))
}

ready_data<-function(data, indexes){
  data$y[-indexes$train]<-data$x[-indexes$train]<-NA
  data$x[data$observed_x==0]<-NA
  data$y[data$observed_y==0]<-NA
  return(data)
}


#Function which obtains Variance Covariance matrices D of the random effects
obtainVarCov<-function(model, index, n1){
  SD_s<-diag(n1)
  CIs<-matrix(nrow=length(index), ncol=2)
  tryCatch(CIs[c(1:n1),]<-t(sapply(index[1:n1],function(i) unlist(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                                                                model$internal.marginals.hyperpar[[i]]), silent=TRUE)[c(3,7)]))), 
           error = function(e) print("Some covariance elements CI not well defined"))
  tryCatch(SD_s<-sqrt(sapply(index[1:n1],function(i) inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                                                   model$internal.marginals.hyperpar[[i]]), silent=TRUE)$mean)), 
           error = function(e) print("Some covariance elements not well defined"))
  cor<-matrix(1, nrow=n1, ncol=n1)
  cor[lower.tri(cor)]<-model$summary.hyperpar[index[-c(1:n1)],1]
  CIs[(n1+1):length(index),]<-unlist(model$summary.hyperpar[index[-c(1:n1)],c(3,5)])
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  answer<-diag(SD_s) %*% cor %*% diag(SD_s)
  return(list(Cov=answer, CI=CIs))
}

Fit_LMM<-function(data, indexes, failed){
  INLA_data<-ready_data(data, indexes)
  N_obs<-length(data$x)
  INLA_data$t_random<-INLA_data$id+indexes$N
  
  #Imputing values of the time-varying covariate that are not available
  tryCatch(imputing_model<-inla(x~f(t, model = "linear")+f(w, model = "linear")+
                                 f(id, model = "iid2d", n=2*indexes$N, hyper = list(theta=list(param=c(4, 1,1,0))))+
                                  f(t_random, t, copy="id"),
                               data = INLA_data, family = "gaussian", verbose=FALSE, 
                               control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE)),
           error = function(e) {print("INLA failed: Exo Model 0");assign('failed',1,envir=globalenv()); failed<<-1})
  if (failed==1){
    print("Returning NA's")
    return(list(coefficients=rep(NA,6), mlik=NA, dic_approx=NA, dic=NA, 
                waic=NA, pit=NA,
                cpo=NA, MSE=c(NA, NA, NA), model_specific=list(cov_random_1=NA, cov_random_2=NA)))
  }

  #Incorporating the imputed values
  INLA_data$x[-indexes$train]<-imputing_model$summary.fitted.values[-indexes$train,1]
  tryCatch(LMM_INLA<-inla(y~f(x, model = "linear")+f(w, model = "linear")+f(t, model = "linear")+
                               f(id, model = "iid2d", n=2*indexes$N, hyper = list(theta=list(param=c(4, 1,1,0))))+
                                 f(t_random, t, copy="id"),
                             data = INLA_data, family = "gaussian", verbose=TRUE, 
                             control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE)),
           error = function(e) {print("INLA failed: Model 0");assign('failed',1,envir=globalenv()); failed<<-1})
  if (failed==1){
    print("Returning NA's")
    return(list(coefficients=rep(NA,6), mlik=NA, dic_approx=NA, dic=NA, 
                waic=NA, pit=NA,
                cpo=NA, MSE=c(NA, NA, NA), model_specific=list(cov_random_1=NA, cov_random_2=NA)))
  }
  Cov1_0_t<-obtainVarCov(LMM_INLA, c(2:4), 2)
  MSE<-list(MSE_train=calculate_MSE(LMM_INLA$summary.fitted.values[indexes$train,1], data$y[indexes$train]), 
            MSE_test_same=calculate_MSE(LMM_INLA$summary.fitted.values[indexes$test_same ,1], 
                                                     data$y[indexes$test_same]), 
            MSE_test_others=calculate_MSE(LMM_INLA$summary.fitted.values[indexes$test_others,1], 
                                          data$y[indexes$test_others]))
  
  return(list(coefficients=LMM_INLA$summary.fixed[c(1,3,4,2),c(1,2,3,5)],
              mlik=LMM_INLA$mlik[1], dic_approx=-2*LMM_INLA$mlik[1], dic=LMM_INLA$dic$dic, 
              waic=LMM_INLA$waic$waic, pit=ks.test(LMM_INLA$cpo$pit,"punif",0,1)$statistic,
              cpo=-sum(log(LMM_INLA$cpo$cpo), na.rm=TRUE), MSE=MSE, model_specific=list(cov_random_1=Cov1_0_t)))
}

Fit_JMM<-function(data, indexes, failed){
  JMM_data<-ready_data(data, indexes)
  N_obs<-dim(JMM_data)[1]
  fixed.effects<-list(Intercept_x=c(rep(1, N_obs), rep(NA, N_obs)),
                      Intercept_y=c(rep(NA, N_obs), rep(1, N_obs)),
                      w_x=c(JMM_data$w, rep(NA, N_obs)),
                      w_y=c(rep(NA, N_obs), JMM_data$w), 
                      t_x=c(JMM_data$t, rep(NA, N_obs)),
                      t_y=c(rep(NA, N_obs), JMM_data$t),
                      t=c(JMM_data$t, JMM_data$t))
  #The 4 random effects should be jointly distributed, so they are indexed from 1 till 4*N
  #With N being the total number of subjects (as both intercept and slope are unique to a subject)
  random.effects<-list(Random_Intercept=c(JMM_data$id, JMM_data$id+indexes$N),
                       Random_Slope=c(JMM_data$id+2*indexes$N, JMM_data$id+3*indexes$N))
  
  INLA_data<-c(fixed.effects, random.effects)
  INLA_data$Y<-list(c(JMM_data$x, rep(NA, N_obs)),
                    c(rep(NA, N_obs),JMM_data$y))
  
  INLA_formula=Y~-1+f(Intercept_x, model='linear')+  f(Intercept_y, model='linear')+
    f(w_x, model='linear')+  f(w_y, model='linear')+f(t_x, model='linear')+f(t_y, model='linear')+
    f(Random_Intercept, model="iid4d", n=4*indexes$N)+  f(Random_Slope, t,copy="Random_Intercept")
  tryCatch(JMM_INLA<-inla(INLA_formula, family = c("gaussian","gaussian"),
                        data = INLA_data, verbose=FALSE, control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE)),
           error = function(e) {print("INLA failed: Model 2");assign('failed',1,envir=globalenv()); failed<<-1})
  if (failed==1){
    print("Returning NA's")
    return(list(coefficients=rep(NA,6), mlik=NA, dic_approx=NA, dic=NA, 
                waic=NA, pit=NA,
                cpo=NA, MSE=c(NA, NA, NA), model_specific=list(cov_both=NA)))
  }
  Cov_0_t<-obtainVarCov(JMM_INLA, c(3:12), 4)
  
  #Computing the association coefficient
  #For exact formula see equation 4.8 in thesis (page 24)
  temp<-inla.hyperpar.sample(100000, JMM_INLA)[,c(1:12)]
  var<-mean(rowSums(1/temp[,c(1,3, 5)])+2*temp[,8]*sqrt(1/(temp[,3]*temp[,5])))
  cov<-temp[,7]*sqrt(1/(temp[,3]*temp[,4]))+temp[,9]*sqrt(1/(temp[,3]*temp[,6]))+
    temp[,10]*sqrt(1/(temp[,4]*temp[,5]))+temp[,12]*sqrt(1/(temp[,5]*temp[,6]))
  coefficient_sampled<-cov/var

  MSE<-calculate_all_MSE(JMM_INLA, data, indexes, indexes$N, indexes$n)
  Cov_error<-diag(2)
  diag(Cov_error)<-NA
  tryCatch(Cov_error<-diag(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),JMM_INLA$internal.marginals.hyperpar[[1]]), silent=TRUE)$mean,
                                inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),JMM_INLA$internal.marginals.hyperpar[[2]]), silent=TRUE)$mean)))
  rownames(Cov_error)<-colnames(Cov_error)<-c('e_x', 'e_y')
  
  
  return(list(coefficients=JMM_INLA$summary.fixed[c(1,3,5,2,4,6),c(1,2,3,5)], 
              mlik=JMM_INLA$mlik[1], 
              dic_approx=-2*JMM_INLA$mlik[1], 
              dic=sum(JMM_INLA$dic$local.dic[(length(data$id)+indexes$train)]), 
              waic=sum(JMM_INLA$waic$local.waic[(length(data$id)+indexes$train)]), 
              pit=ks.test(JMM_INLA$cpo$pit[(length(data$id)+indexes$train)],"punif",0,1)$statistic,
              cpo=-sum(log(JMM_INLA$cpo$cpo[(length(data$id)+indexes$train)]), na.rm=TRUE), 
              MSE=MSE, model_specific=list(cov_both=Cov_0_t,Cov_error=Cov_error),
              hyperparameters=summary(JMM_INLA)$hyperpar,
         exo_coefficient=list(coefficient_sampled=quantile(coefficient_sampled, probs=c(0.025, 0.5, 0.975), na.rm=TRUE))))
}

Fit_JSM<-function(data, indexes, failed){
  JSM_data<-ready_data(data, indexes)
  #Modelling as 2 separate likelihoods
  N_obs<-dim(JSM_data)[1]
  fixed.effects<-list(w=c(JSM_data$w, JSM_data$w), 
                         t=c(JSM_data$t, JSM_data$t),
                         Intercept_y=c(rep(NA, N_obs), rep(1, N_obs)),
                         t_y=c(rep(NA, N_obs),JSM_data$t),
                         w_y=c(rep(NA, N_obs),JSM_data$w))
  random.effects<-list(Intercept_x=c(rep(1, N_obs), rep(NA, N_obs)),
                          Intercept_x_scaled=c(rep(NA, N_obs), rep(1, N_obs)),
                          w_x=c(rep(1, N_obs),rep(NA, N_obs)),
                          w_x_scaled=c(rep(NA, N_obs),rep(1, N_obs)),
                          t_x=c(rep(1, N_obs),rep(NA, N_obs)),
                          t_x_scaled=c(rep(NA, N_obs),rep(1, N_obs)),
                          Random_intercept_x=c(JSM_data$id, rep(NA, N_obs)),
                          Random_intercept_x_scaled=c(rep(NA, N_obs), JSM_data$id),
                          Random_intercept_y=c(rep(NA, N_obs), JSM_data$id),
                          Random_slope_x=c(JSM_data$id+indexes$N, rep(NA, N_obs)),
                          Random_slope_x_scaled=c(rep(NA, N_obs), JSM_data$id+indexes$N),
                          Random_slope_y=c(rep(NA, N_obs), JSM_data$id+indexes$N))
  
  INLA_data<-c(fixed.effects, random.effects)
  INLA_data$Y<-list(c(JSM_data$x, rep(NA, N_obs)),
                    c(rep(NA, N_obs),JSM_data$y))
  
  INLA_formula=Y~-1+
    f(Intercept_x)+ 
    #Copying and scaling fixed effects
    f(Intercept_x_scaled, copy="Intercept_x",
      hyper = list(beta = list(fixed=FALSE)))+ 
    f(w_x, w)+
    #same.as = 'Intercept_x_scaled' is added to ensure that all factors are scaled with the same gamma
    f(w_x_scaled, w, copy="w_x", same.as = 'Intercept_x_scaled',
      hyper = list(beta = list(fixed=FALSE)))+ 
    f(t_x, t)+
    f(t_x_scaled, t, copy="t_x", same.as = 'Intercept_x_scaled',
      hyper = list(beta = list(fixed=FALSE)))+
    f(Intercept_y, model='linear')+f(w_y, model='linear')+f(t_y, model='linear')+ #Model 2 fixed effects
    f(Random_intercept_x, model="iid2d", n=2*indexes$N)+f(Random_slope_x, t, copy="Random_intercept_x")+
    #Copying and scaling random effects
    f(Random_intercept_x_scaled, copy='Random_intercept_x', same.as='Intercept_x_scaled', fixed=FALSE)+
    f(Random_slope_x_scaled, t, copy='Random_intercept_x', same.as='Intercept_x_scaled', fixed=FALSE)+
    f(Random_intercept_y, model="iid2d", n=2*indexes$N)+f(Random_slope_y, t, copy="Random_intercept_y")
  
  tryCatch(JSM_INLA<-inla(INLA_formula, family = c("gaussian","gaussian"),
                                data = INLA_data,verbose=FALSE,
                                control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE)),
           error = function(e) {print("INLA failed: Model 3");assign('failed',1,envir=globalenv()); failed<<-1})
  if (failed==1){
    print("Returning NA's")
    return(list(coefficients=rep(NA,6), mlik=NA, dic_approx=NA, dic=NA, 
                waic=NA, pit=NA, 
                cpo=NA, MSE=c(NA, NA, NA), model_specific=list(cov_both_1=NA, cov_both_2=NA, 
                                                               gamma_s=NA)))
  }
  gamma<-summary(JSM_INLA)$hyperpar[12,1]


  x_coefficients=rbind(JSM_INLA$summary.random$Intercept1[c(2,3,4,6)],
                       JSM_INLA$summary.random$idv1[c(2,3,4,6)],
                       JSM_INLA$summary.random$t1[c(2,3,4,6)])
  rownames(x_coefficients)<-c('Intercept1', 'v1', 't1')
  Cov1_0_t<-obtainVarCov(JSM_INLA, c(6:8), 2)
  rownames(Cov1_0_t$Cov)<-colnames(Cov1_0_t$Cov)<-c('u_x0','u_xt')
  Cov2_0_t<-obtainVarCov(JSM_INLA, c(9:11), 2)
  rownames(Cov2_0_t$Cov)<-colnames(Cov2_0_t$Cov)<-c('u_y0','u_yt')
  Cov_error<-diag(2)
  diag(Cov_error)<-NA
  tryCatch(Cov_error<-diag(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),JSM_INLA$internal.marginals.hyperpar[[1]]), silent=TRUE)$mean,
                             inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),JSM_INLA$internal.marginals.hyperpar[[2]]), silent=TRUE)$mean)))
  rownames(Cov_error)<-colnames(Cov_error)<-c('e_x', 'e_y')
  Cov2_0_t_scaled<-gamma^2*Cov1_0_t$Cov+Cov2_0_t$Cov
  rownames(Cov2_0_t_scaled)<-colnames(Cov2_0_t_scaled)<-c('u_y0_combined','u_yt_combined')
  MSE<-calculate_all_MSE(JSM_INLA, data, indexes, indexes$N, indexes$n)
  return(list(actual_coefficients=x_coefficients,
              mlik=JSM_INLA$mlik[1], dic_approx=-2*JSM_INLA$mlik[1], 
              dic=sum(JSM_INLA$dic$local.dic[(length(data$id)+indexes$train)]), 
              waic=sum(JSM_INLA$waic$local.waic[(length(data$id)+indexes$train)]), 
              pit=ks.test(JSM_INLA$cpo$pit[(length(data$id)+indexes$train)],"punif",0,1)$statistic,
              cpo=-sum(log(JSM_INLA$cpo$cpo[(length(data$id)+indexes$train)]), na.rm=TRUE), 
              MSE=MSE, hyperparameters=summary(JSM_INLA)$hyperpar,
              model_specific=list(cov_both_1=Cov1_0_t, cov_both_2_actual=Cov2_0_t,cov_both_2_scaled=Cov2_0_t_scaled, Cov_error=Cov_error,
                                  gamma=list(gamma=gamma, CI=summary(JSM_INLA)$hyperpar[12,c(3,5)]),
                                  model_coefficients=rbind(x_coefficients, JSM_INLA$summary.fixed[,c(1,2,3,5)]))))
}


