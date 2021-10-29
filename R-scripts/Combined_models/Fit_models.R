fit_model_0<-function(data, N, n){
  INLA_data_0<-data
  #Modelling as 2 separate likelihoods
  fixed.effects_0<-list(Intercept1=c(rep(1, n*N), rep(NA, n*N)),
                        Period1=c(INLA_data_0$Period, rep(NA, n*N)),
                        x1=c(INLA_data_0$x, rep(NA, n*N)),
                        Intercept2=c(rep(NA, n*N), rep(1, n*N)),
                        Period2=c(rep(NA, n*N),INLA_data_0$Period),
                        x2=c(rep(NA, n*N),INLA_data_0$x),
                        Period=c(INLA_data_0$Period, INLA_data_0$Period))
  random.effects_0<-list(Intercept_random1=c(INLA_data_0$id, rep(NA, n*N)),
                         Period_random1=c(INLA_data_0$id+N, rep(NA, n*N)),
                         Intercept_random2=c(rep(NA, n*N), INLA_data_0$id),
                         Period_random2=c(rep(NA, n*N), INLA_data_0$id+N))
  y_final_0<-list(c(INLA_data_0$y1, rep(NA, n*N)),
                  c(rep(NA, n*N),INLA_data_0$y2))
  final_data_0<-c(fixed.effects_0, random.effects_0)
  final_data_0$Y<-y_final_0
  
  formula.model_0=Y~-1+
    f(Intercept1, model='linear')+f(x1, model='linear')+f(Period1, model='linear')+ #Model 1 fixed effects
    f(Intercept_random1, model="iid2d", n=2*N)+f(Period_random1, Period, copy="Intercept_random1")+
    f(Intercept2, model='linear')+f(x2, model='linear')+f(Period2, model='linear')+ #Model 2 fixed effects
    f(Intercept_random2, model="iid2d", n=2*N)+f(Period_random2, Period, copy="Intercept_random2")
  
  final_model_0<-inla(formula.model_0, family = c("gaussian","gaussian"),
                      data = final_data_0,verbose=FALSE,
                      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))
  
  coefficients_0<-final_model_0$summary.fixed[,1]
  
  #Random effects covariance Y1
  SD_s<-diag(2)
  tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_0$internal.marginals.hyperpar[[3]], silent=TRUE))$mean,
                        inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_0$internal.marginals.hyperpar[[4]]))$mean)), 
           error = function(e) print("Model 0: Some covariance elements of Y1 not well defined"))
  cor<-diag(2)
  cor[lower.tri(cor)]<-summary(final_model_0)$hyperpar[5, 1]
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  Cov1_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)
  
  #Random effects covariance Y2
  SD_s<-diag(2)
  tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_0$internal.marginals.hyperpar[[6]], silent=TRUE))$mean,
                        inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_0$internal.marginals.hyperpar[[7]]))$mean)), 
           error = function(e) print("Model 0: Some covariance elements of Y2 not well defined"))
  cor<-diag(2)
  cor[lower.tri(cor)]<-summary(final_model_0)$hyperpar[8, 1]
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  Cov2_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)
  return(list(coefficients=coefficients_0, mlik=final_model_0$mlik[1], dic=final_model_0$dic$dic, waic=final_model_0$waic$waic,
              cpo=-sum(log(final_model_0$cpo$cpo)), model_specific=list(cov_random_1=Cov1_0_t, cov_random_2=Cov2_0_t)))
}

fit_model_1A<-function(data, N, n){
  data<-rbind(data[data$Period==0,], data[data$Period==1,])
  fixed.effects<-list(Intercept1=c(rep(1, n*N), rep(NA, n*N)),
                      Intercept2=c(rep(NA, n*N), rep(1, n*N)),
                      x1=c(data$x, rep(NA, n*N)),
                      x2=c(rep(NA, n*N), data$x),
                      Period1=c(data$Period, rep(NA, n*N)),
                      Period2=c(rep(NA, n*N), data$Period))
  INLA_y<-list(c(data$y1, rep(NA, n*N)),
               c(rep(NA, n*N),data$y2))
  INLA_data<-c(fixed.effects)
  INLA_data$Y<-INLA_y
  INLA_data$i<-c(1:(2*n*N))
  control_family = list(list(initial=15, fixed=T), list(initial=15, fixed=T))
  INLA_formula_1=Y~-1+f(Intercept1, model='linear')+f(Intercept2, model='linear')+
    f(x1, model='linear')+f(x2, model='linear')+
    f(Period1, model='linear')+f(Period2, model='linear')+
    f(i, model="iid4d", n=2*n*N, hyper = list(theta1 = list(prior = "wishart4d")))
  INLA_model_1<-inla(INLA_formula_1, family = c("gaussian","gaussian"),
                     data = INLA_data, verbose=FALSE, control.family=control_family, 
                     control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))
  
  coefficients_1<-INLA_model_1$summary.fixed[c(1,3,5, 2,4, 6),1]
  SD_s<-diag(4)
  for (i in 1:4){
    tryCatch(SD_s[i,i]<-sqrt(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),INLA_model_1$internal.marginals.hyperpar[[i]], silent=TRUE))$quant0.5),
             error = function(e) print("Model 1A: Some error covariance elements not well defined"))  
  }
  cor<-matrix(1, nrow=4, ncol=4)
  cor[lower.tri(cor)]<-INLA_model_1$summary.hyperpar[5:10,1]
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  Errors<-SD_s %*% cor %*% SD_s
  return(list(coefficients=coefficients_1, mlik=INLA_model_1$mlik[1], dic=INLA_model_1$dic$dic, waic=INLA_model_1$waic$waic,
              cpo=-sum(log(INLA_model_1$cpo$cpo)), model_specific=list(cov_errors=Errors)))
}

fit_model_2A<-function(data, N, n){
  INLA_data_2a<-data
  #Modelling as 2 separate likelihoods
  fixed.effects_2a<-list(Intercept1=c(rep(1, n*N), rep(NA, n*N)),
                         Intercept2=c(rep(NA, n*N), rep(1, n*N)),
                         x1=c(INLA_data_2a$x, rep(NA, n*N)),
                         x2=c(rep(NA, n*N), INLA_data_2a$x), 
                         Period1=c(INLA_data_2a$Period, rep(NA, n*N)),
                         Period2=c(rep(NA, n*N), INLA_data_2a$Period))
  
  random.effects_2a<-list(ID_all=c(INLA_data_2a$id, INLA_data_2a$id+N))
  y_final_2a<-list(c(INLA_data_2a$y1, rep(NA, n*N)),
                   c(rep(NA, n*N),INLA_data_2a$y2))
  final_data_2a<-c(fixed.effects_2a, random.effects_2a)
  final_data_2a$Y<-y_final_2a
  
  formula.model_2a=Y~-1+f(Intercept1, model='linear')+  f(Intercept2, model='linear')+
    f(x1, model='linear')+  f(x2, model='linear')+
    f(Period1, model='linear')+f(Period2, model='linear')+
    f(ID_all, model="iid2d", n=2*N,hyper = list(theta1 = list(prior = "wishart2d", param =c(4,1,1, 0))))
  final_model_2a<-inla(formula.model_2a, family = c("gaussian","gaussian"),
                       data = final_data_2a,verbose=FALSE,
                       control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))
  
  
  coefficients_2a<-final_model_2a$summary.fixed[c(1,3,5,2,4,6),1]
  SD_s<-diag(2)
  tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_2a$internal.marginals.hyperpar[[3]], silent=TRUE))$mean,
                        inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_2a$internal.marginals.hyperpar[[4]]))$mean)), 
           error = function(e)
             print("Model 2A: Some random effect covariance elements not well defined"))
  
  
  cor<-matrix(1, nrow=2, ncol=2)
  cor[lower.tri(cor)]<-summary(final_model_2a)$hyperpar[5, 1]
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  U_cov<-diag(SD_s) %*% cor %*% diag(SD_s)
  
  return(list(coefficients=coefficients_2a, mlik=final_model_2a$mlik[1], dic=final_model_2a$dic$dic, waic=final_model_2a$waic$waic,
              cpo=-sum(log(final_model_2a$cpo$cpo)), model_specific=list(cov_intercept=U_cov)))
}

fit_model_2B<-function(data, N, n){
  data<-rbind(data[data$Period==0,], data[data$Period==1,])
  INLA_data_2b<-data
  fixed.effects_2b<-list(Intercept1=c(rep(1, n*N), rep(NA, n*N)),
                         Intercept2=c(rep(NA, n*N), rep(1, n*N)),
                         x1=c(INLA_data_2b$x, rep(NA, n*N)),
                         x2=c(rep(NA, n*N), INLA_data_2b$x),
                         Period1=c(INLA_data_2b$Period, rep(NA, n*N)),
                         Period2=c(rep(NA, n*N), INLA_data_2b$Period))
  random.effects_2b<-list(ID_all=c(INLA_data_2b$id, INLA_data_2b$id+N))
  y_final_2b<-list(c(INLA_data_2b$y1, rep(NA, n*N)),
                   c(rep(NA, n*N),INLA_data_2b$y2))
  final_data_2b<-c(fixed.effects_2b, random.effects_2b)
  final_data_2b$Y<-y_final_2b
  final_data_2b$i<-c(1:(2*n*N))
  
  #Set-up Model
  control_family = list(list(initial=15, fixed=T), list(initial=15, fixed=T))
  formula.model_2b=Y~-1+f(Intercept1, model='linear')+  f(Intercept2, model='linear')+
    f(x1, model='linear')+  f(x2, model='linear')+
    f(Period1, model='linear')+f(Period2, model='linear')+
    f(ID_all, model="iid2d", n=2*N,hyper = list(theta1 = list(prior = "wishart2d")))+
    f(i, model="iid4d", n=2*n*N, hyper = list(theta1 = list(prior = "wishart4d")))
  
  #Run Model
  final_model_2b<-inla(formula.model_2b, family = c("gaussian","gaussian"),
                       data = final_data_2b,verbose=FALSE, control.family=control_family,
                       control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))
  
  
  #Set random effect coefficients
  coefficients_2b<-final_model_2b$summary.fixed[c(1,3,5,2,4,6),1]
  SD_s<-diag(2)
  
  tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                        final_model_2b$internal.marginals.hyperpar[[1]], silent=TRUE))$mean,
          inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                        final_model_2b$internal.marginals.hyperpar[[2]]))$mean)), 
          error = function(e)
            print("Model 2B: Some random effect covariance elements not well defined"))
  cor<-matrix(1, nrow=2, ncol=2)
  cor[lower.tri(cor)]<-summary(final_model_2b)$hyperpar[3, 1]
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  U_cov<-diag(SD_s) %*% cor %*% diag(SD_s)

  
  #Set correlated errors
  SD_s<-diag(4)
  for (i in 1:4){
    tryCatch(SD_s[i,i]<-sqrt(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),final_model_2b$internal.marginals.hyperpar[[(i+3)]], silent=TRUE))$quant0.5),
             error = function(e)   print("Model 2B: Some residual error covariance elements not well defined"))
  }
  cor<-matrix(1, nrow=4, ncol=4)
  cor[lower.tri(cor)]<-final_model_2b$summary.hyperpar[8:13,1]
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  Errors<-SD_s %*% cor %*% SD_s
  
  
  return(list(coefficients=coefficients_2b, mlik=final_model_2b$mlik[1], dic=final_model_2b$dic$dic, waic=final_model_2b$waic$waic,
              cpo=-sum(log(final_model_2b$cpo$cpo)), model_specific=list(cov_intercept=U_cov, cov_errors=Errors)))

}

fit_model_2C1<-function(data, N, n){
  INLA_data_2c1<-data
  #Modelling as 2 separate likelihoods
  Period_id<-rep(NA, times=n*N)
  Period_id[INLA_data_2c1$Period2!=0]<-rep(1:N, each=(n-1))
  fixed.effects_2c1<-list(Intercept1=c(rep(1, n*N), rep(NA, n*N)),
                          Intercept2=c(rep(NA, n*N), rep(1, n*N)),
                          x1=c(INLA_data_2c1$x, rep(NA, n*N)),
                          x2=c(rep(NA, n*N), INLA_data_2c1$x), 
                          Period1=c(INLA_data_2c1$Period, rep(NA, n*N)),
                          Period2=c(rep(NA, n*N), INLA_data_2c1$Period),
                          Period_all=c(INLA_data_2c1$Period, INLA_data_2c1$Period))
  random.effects_2c1<-list(ID_all=c(INLA_data_2c1$id, INLA_data_2c1$id+N),
                           Period=c(Period_id, Period_id+N))
  
  y_final_2c1<-list(c(INLA_data_2c1$y1, rep(NA, n*N)),
                    c(rep(NA, n*N),INLA_data_2c1$y2))
  final_data_2c1<-c(fixed.effects_2c1, random.effects_2c1)
  final_data_2c1$Y<-y_final_2c1
  
  formula.model_2c1=Y~-1+f(Intercept1, model='linear')+  f(Intercept2, model='linear')+
    f(x1, model='linear')+  f(x2, model='linear')+f(Period1, model='linear')+f(Period2, model='linear')+
    f(ID_all, model="iid2d", n=2*N)+f(Period, Period_all, model="iid2d", n=2*N)
  #  f(ID_all_2,Period2, model="iid2d", n=2*N,hyper = list(theta1 = list(prior = "wishart2d", param =c(4,1,1, 0))))
  final_model_2c1<-inla(formula.model_2c1, family = c("gaussian","gaussian"),
                        data = final_data_2c1,verbose=FALSE,
                        control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))
  #Fixed effects
  coefficients_2c<-final_model_2c1$summary.fixed[c(1,3,5,2,4,6),1]
  #Intercept covariance matrix
  SD_s<-diag(2)
  tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_2c1$internal.marginals.hyperpar[[3]], silent=TRUE))$mean, 
                        inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_2c1$internal.marginals.hyperpar[[4]]))$mean)), 
           error = function(e)
             print("Model 2C1: Some covariance elements not well defined"))
  
  cor<-matrix(1, nrow=2, ncol=2)
  cor[lower.tri(cor)]<-summary(final_model_2c1)$hyperpar[5, 1]
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  U_cov_0<-diag(SD_s) %*% cor %*% diag(SD_s)
  #Slope covariance matrix
  SD_s<-diag(2)
  tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                        final_model_2c1$internal.marginals.hyperpar[[6]], silent=TRUE))$mean, 
          inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                        final_model_2c1$internal.marginals.hyperpar[[7]]))$mean)),
          error = function(e)
            print("Model 2C1: Some covariance elements not well defined"))
  cor<-matrix(1, nrow=2, ncol=2)
  cor[lower.tri(cor)]<-summary(final_model_2c1)$hyperpar[8, 1]
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  U_cov_t<-diag(SD_s) %*% cor %*% diag(SD_s)
  
  return(list(coefficients=coefficients_2c, mlik=final_model_2c1$mlik[1], dic=final_model_2c1$dic$dic, waic=final_model_2c1$waic$waic,
              cpo=-sum(log(final_model_2c1$cpo$cpo)), model_specific=list(cov_intercept=U_cov_0, cov_slope=U_cov_t)))
  
  
}

fit_model_2C2<-function(data, N, n){
  INLA_data_2c2<-data
  #Modelling as 2 separate likelihoods
  #Use copy feature to create dependence between Intercept and Slope
  #https://arxiv.org/pdf/1210.0333.pdf
  Period_id<-rep(NA, times=n*N)
  Period_id[INLA_data_2c2$Period2!=0]<-rep(1:N, each=(n-1))
  fixed.effects_2c2<-list(Intercept1=c(rep(1, n*N), rep(NA, n*N)),
                          Intercept2=c(rep(NA, n*N), rep(1, n*N)),
                          x1=c(INLA_data_2c2$x, rep(NA, n*N)),
                          x2=c(rep(NA, n*N), INLA_data_2c2$x), 
                          Period1=c(INLA_data_2c2$Period, rep(NA, n*N)),
                          Period2=c(rep(NA, n*N), INLA_data_2c2$Period),
                          Period_all=c(INLA_data_2c2$Period, INLA_data_2c2$Period))
  random.effects_2c2<-list(ID_all=c(INLA_data_2c2$id, INLA_data_2c2$id+N),
                           Period_random=c(Period_id+2*N, Period_id+3*N))
  
  y_final_2c2<-list(c(INLA_data_2c2$y1, rep(NA, n*N)),
                    c(rep(NA, n*N),INLA_data_2c2$y2))
  final_data_2c2<-c(fixed.effects_2c2, random.effects_2c2)
  final_data_2c2$Y<-y_final_2c2
  
  formula.model_2c2=Y~-1+f(Intercept1, model='linear')+  f(Intercept2, model='linear')+
    f(x1, model='linear')+  f(x2, model='linear')+f(Period1, model='linear')+f(Period2, model='linear')+
    f(ID_all, model="iid4d", n=4*N)+  f(Period_random, Period_all,copy="ID_all")
  #  f(ID_all_2,Period2, model="iid2d", n=2*N,hyper = list(theta1 = list(prior = "wishart2d", param =c(4,1,1, 0))))
  final_model_2c2<-inla(formula.model_2c2, family = c("gaussian","gaussian"),
                        data = final_data_2c2,verbose=FALSE, control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))
  
  
  coefficients_2c2<-final_model_2c2$summary.fixed[c(1,3,5,2,4,6),1]
  SD_s<-diag(4)
  tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                             final_model_2c2$internal.marginals.hyperpar[[3]], silent=TRUE))$mean,
               inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                             final_model_2c2$internal.marginals.hyperpar[[4]]))$mean,
               inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                             final_model_2c2$internal.marginals.hyperpar[[5]]))$mean,
               inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                             final_model_2c2$internal.marginals.hyperpar[[6]]))$mean)), 
           error = function(e)
             print("Model 2C2: Some covariance elements not well defined"))
  cor<-matrix(1, nrow=4, ncol=4)
  cor[lower.tri(cor)]<-final_model_2c2$summary.hyperpar[7:12,1]
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  Cov_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)
  
  return(list(coefficients=coefficients_2c2, mlik=final_model_2c2$mlik[1], dic=final_model_2c2$dic$dic, waic=final_model_2c2$waic$waic,
              cpo=-sum(log(final_model_2c2$cpo$cpo)), model_specific=list(cov_both=Cov_0_t)))
  
  
}

fit_model_3A1<-function(data, N, n){
  INLA_data_3a<-data
  #Modelling as 2 separate likelihoods
  fixed.effects_3a<-list(x=c(INLA_data_3a$x, INLA_data_3a$x), 
                         Period=c(INLA_data_3a$Period, INLA_data_3a$Period),
                         Intercept2=c(rep(NA, n*N), rep(1, n*N)),
                         Period2=c(rep(NA, n*N),INLA_data_3a$Period),
                         x2=c(rep(NA, n*N),INLA_data_3a$x))
  random.effects_3a<-list(Intercept1=c(rep(1, n*N), rep(NA, n*N)),
                          Intercept12=c(rep(NA, n*N), rep(1, n*N)),
                          idx1=c(rep(1, n*N),rep(NA, n*N)),
                          idx12=c(rep(NA, n*N),rep(1, n*N)),
                          Period1=c(rep(1, n*N),rep(NA, n*N)),
                          Period12=c(rep(NA, n*N),rep(1, n*N)),
                          Intercept_random1=c(INLA_data_3a$id, rep(NA, n*N)),
                          Intercept_random12=c(rep(NA, n*N), INLA_data_3a$id),
                          Intercept_random2=c(rep(NA, n*N), INLA_data_3a$id),
                          Period_random1=c(INLA_data_3a$id, rep(NA, n*N)),
                          Period_random12=c(rep(NA, n*N), INLA_data_3a$id),
                          Period_random2=c(rep(NA, n*N), INLA_data_3a$id+N))
  y_final_3a<-list(c(INLA_data_3a$y1, rep(NA, n*N)),
                   c(rep(NA, n*N),INLA_data_3a$y2))
  final_data_3a<-c(fixed.effects_3a, random.effects_3a)
  final_data_3a$Y<-y_final_3a
  
  formula.model_3a=Y~-1+
    f(Intercept1)+ #Copying fixed effects
    f(Intercept12, copy="Intercept1",
      hyper = list(beta = list(fixed=FALSE)))+
    f(idx1, x)+
    f(idx12, x, copy="idx1", same.as = 'Intercept12',
      hyper = list(beta = list(fixed=FALSE)))+
    f(Period1, Period)+
    f(Period12, Period, copy="Period1", same.as = 'Intercept12',
      hyper = list(beta = list(fixed=FALSE)))+
    f(Intercept_random1, model="iid")+ #Copying random effects
    f(Intercept_random12, copy="Intercept_random1", same.as = 'Intercept12',
      hyper = list(beta = list(fixed=FALSE)))+
    f(Period_random1, Period, model="iid")+
    f(Period_random12, copy="Period_random1", same.as = 'Intercept12',
      hyper = list(beta = list(fixed=FALSE)))+
    f(Intercept2, model='linear')+f(x2, model='linear')+f(Period2, model='linear')+ #Model 2 fixed effects
    f(Intercept_random2, model="iid2d", n=2*N)+f(Period_random2, Period, copy="Intercept_random2")
  
  final_model_3a<-inla(formula.model_3a, family = c("gaussian","gaussian"),
                       data = final_data_3a,verbose=FALSE, control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))
  
  model_coefficients<-c(final_model_3a$summary.random$Intercept1$mean, 
                                       final_model_3a$summary.random$idx1$mean,
                                       final_model_3a$summary.random$Period1$mean,
                                       final_model_3a$summary.fixed[c(1,2,3),1])
  gamma<-summary(final_model_3a)$hyperpar[11,1]
  actual_coefficients<-c(model_coefficients[1:3], model_coefficients[1:3]*gamma+model_coefficients[4:6])
  
  
  
  
  
  SD_s<-diag(2)
  tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_3a$internal.marginals.hyperpar[[6]], silent=TRUE))$mean,
                        inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_3a$internal.marginals.hyperpar[[7]]))$mean)), 
           error = function(e) print("Model 3A1: Some covariance elements of Y1 not well defined"))
  cor<-diag(2)
  Cov1_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)
  
  SD_s<-diag(2)
  tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_3a$internal.marginals.hyperpar[[8]], silent=TRUE))$mean,
                        inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_3a$internal.marginals.hyperpar[[9]]))$mean)), 
           error = function(e) print("Model 3A1: Some covariance elements of Y2 not well defined"))
  cor<-diag(2)
  cor[lower.tri(cor)]<-summary(final_model_3a)$hyperpar[10, 1]
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  Cov2_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)
  
  
  return(list(coefficients=actual_coefficients, mlik=final_model_3a$mlik[1], dic=final_model_3a$dic$dic, waic=final_model_3a$waic$waic,
              cpo=-sum(log(final_model_3a$cpo$cpo)), model_specific=list(cov_both_1=Cov1_0_t, cov_both_2=Cov2_0_t, 
                                                                          model_coefficients=model_coefficients, gamma=gamma)))
}

fit_model_3B1<-function(data, N, n){
  INLA_data_3b1<-data
  #Modelling as 2 separate likelihoods
  fixed.effects_3b1<-list(x=c(INLA_data_3b1$x, INLA_data_3b1$x), 
                          Period=c(INLA_data_3b1$Period, INLA_data_3b1$Period),
                          Intercept2=c(rep(NA, n*N), rep(1, n*N)),
                          Intercept1=c(rep(1, n*N), rep(NA, n*N)),
                          Period2=c(rep(NA, n*N),INLA_data_3b1$Period),
                          Period1=c(INLA_data_3b1$Period, rep(NA, n*N)),
                          x2=c(rep(NA, n*N),INLA_data_3b1$x),
                          x1=c(INLA_data_3b1$x, rep(NA, n*N)))
  random.effects_3b1<-list(Intercept_random1=c(INLA_data_3b1$id, rep(NA, n*N)),
                           Intercept_random12=c(rep(NA, n*N), INLA_data_3b1$id),
                           Intercept_random2=c(rep(NA, n*N), INLA_data_3b1$id),
                           Period_random1=c(INLA_data_3b1$id, rep(NA, n*N)),
                           Period_random12=c(rep(NA, n*N), INLA_data_3b1$id),
                           Period_random2=c(rep(NA, n*N), INLA_data_3b1$id+N))
  y_final_3b1<-list(c(INLA_data_3b1$y1, rep(NA, n*N)),
                    c(rep(NA, n*N),INLA_data_3b1$y2))
  final_data_3b1<-c(fixed.effects_3b1, random.effects_3b1)
  final_data_3b1$Y<-y_final_3b1
  
  formula.model_3b1=Y~-1+
    f(Intercept1, model='linear')+f(x1, model='linear')+f(Period1, model='linear')+ #Model 1 fixed effects
    f(Intercept_random1, model="iid")+ #Copying random effects
    f(Intercept_random12, copy="Intercept_random1",
      hyper = list(beta = list(fixed=FALSE)))+
    f(Period_random1, Period, model="iid")+
    f(Period_random12, copy="Period_random1",
      hyper = list(beta = list(fixed=FALSE)))+
    f(Intercept2, model='linear')+f(x2, model='linear')+f(Period2, model='linear')+ 
    f(Intercept_random2, model="iid2d", n=2*N)+f(Period_random2, Period, copy="Intercept_random2")
  
  final_model_3b1<-inla(formula.model_3b1, family = c("gaussian","gaussian"),
                        data = final_data_3b1,verbose=FALSE, control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))
  
  #Coefficients and gamma_s
  coefficients_3b1<-c(final_model_3b1$summary.fixed[1:6,1])
  gamma_s<-summary(final_model_3b1)$hyperpar[8:9,1]
  
  #Y1 random effects covariance structure
  SD_s<-diag(2)
  tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_3b1$internal.marginals.hyperpar[[3]], silent=TRUE))$mean,
                        inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_3b1$internal.marginals.hyperpar[[4]]))$mean)), 
           error = function(e) print("Model 3B1: Some covariance elements of Y1 not well defined"))
  cor<-diag(2)
  Cov1_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)
  #Y2 random effects covariance structure
  SD_s<-diag(2)
  tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_3b1$internal.marginals.hyperpar[[5]], silent=TRUE))$mean,
                        inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                      final_model_3b1$internal.marginals.hyperpar[[6]]))$mean)), 
           error = function(e) print("Model 3A1: Some covariance elements of Y2 not well defined"))
  cor<-diag(2)
  cor[lower.tri(cor)]<-summary(final_model_3b1)$hyperpar[7, 1]
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  Cov2_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)

  
  return(list(coefficients=coefficients_3b1, mlik=final_model_3b1$mlik[1], dic=final_model_3b1$dic$dic, waic=final_model_3b1$waic$waic,
              cpo=-sum(log(final_model_3b1$cpo$cpo)), model_specific=list(cov_both_1=Cov1_0_t, cov_both_2=Cov2_0_t, 
                                                                          gamma_s=gamma_s)))
}

fit_all_models<-function(data, N, n, df){
  print('Fitting model 1/8')
  df$Model_0<-unlist(fit_model_0(data, N, n)[c(1:5)])
  print('Fitting model 2/8')
  df$Model_1A<-unlist(fit_model_1A(data, N, n)[c(1:5)])
  print('Fitting model 3/8')
  df$Model_2A<-unlist(fit_model_2A(data, N, n)[c(1:5)])
  print('Fitting model 4/8')
  df$Model_2B<-unlist(fit_model_2B(data, N, n)[c(1:5)])
  print('Fitting model 5/8')
  df$Model_2C1<-unlist(fit_model_2C1(data, N, n)[c(1:5)])
  print('Fitting model 6/8')
  df$Model_2C2<-unlist(fit_model_2C2(data, N, n)[c(1:5)])
  print('Fitting model 7/8')
  df$Model_3A1<-unlist(fit_model_3A1(data, N, n)[c(1:5)])
  print('Fitting model 8/8')
  df$Model_3B1<-unlist(fit_model_3B1(data, N, n)[c(1:5)])
  return(df)
}
