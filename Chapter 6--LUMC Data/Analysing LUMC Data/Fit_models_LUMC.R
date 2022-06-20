################################################################################
################################################################################
#
################################################################################
################################################################################

#Function that introduces breakpoints
b1<-function(x, bp) ifelse(x<bp, 0, x-bp)

#Gives indexes for training set, test set and subsequent measurements of subjects in test set
#Based on number of subjects that need to be included in test set
#Subjects in test set are chosen at random
test_indexes<-function(data, number){
  if (number==0){
    return(indexes=list(train=c(1:length(data$ID)), test_same=NA, test=NA))
  }
  test_same_subjects<-sort(sample(unique(data$ID), number, replace=FALSE))
  test_subject_indexes<-which(data$ID %in% test_same_subjects)
  indexes_per_subject<-split(test_subject_indexes,data$ID[data$ID %in% test_same_subjects])
  number_to_keep<-lapply(indexes_per_subject,function(i) c((1:(floor(length(i)/2)))))
  train<-unlist(sapply(c(1:number), function(i) indexes_per_subject[[i]][number_to_keep[[i]]]))
  
  test_subjects<-sort(sample(setdiff(unique(data$ID), test_same_subjects), number, replace=FALSE))
  test<-which(data$ID %in% test_subjects)
  
  test_same<-setdiff(test_subject_indexes, train)
  train<-setdiff(c(1:length(data$ID)), c(test_same, test))
  return(indexes=list(train=train, test_same=test_same, test=test))
}

calculate_MSE<-function(y, y_hat){
  return(mean((y-y_hat)^2))
}

calculate_all_MSE<-function(fitted_model, data, indexes){
  total_length<-length(data$ID)
  MSE_train<-calculate_MSE(fitted_model$summary.fitted.values[c(indexes$train+total_length),1],
                           data$y[indexes$train])
  MSE_test_same<-calculate_MSE(fitted_model$summary.fitted.values[c(indexes$test_same+total_length),1],
                               data$y[indexes$test_same])
  MSE_test_others<-calculate_MSE(fitted_model$summary.fitted.values[c(indexes$test+total_length),1],
                                 data$y[indexes$test])
  return(list(MSE_train=MSE_train, MSE_test_same=MSE_test_same, 
              MSE_test_others=MSE_test_others))
}

#Obtain random effects variance-covariance matrix of the fitted model
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

#Delete elements in the test set + subsequent measurements of subjects in training set
ready_data<-function(data, indexes){
  data$y[-indexes$train]<-data$exo[-indexes$train]<-NA
  return(data)
}

Fit_LMM<-function(data, indexes){
  INLA_data<-ready_data(data, indexes)
  N<-length(unique(data$ID))
  INLA_data$time_random<-data$ID+N
  INLA_data$bp1_random<-data$ID+2*N
  exo_predict_INLA<-inla(exo~f(day_no, model = "linear")+
                           f(b1(day_no, 28), model = "linear")+
                           f(ID, model = "iid3d", n=3*N)+f(time_random, day_no, copy="ID")+
                           f(bp1_random, b1(day_no, 28), copy="ID"), 
                         data = INLA_data, family = "gaussian", verbose=FALSE)

  INLA_data$exo[-indexes$train]<-exo_predict_INLA$summary.fitted.values[-indexes$train,1]
  LMM_INLA<-inla(y~f(exo, model = "linear")+f(day_no, model = "linear")+
                  f(b1(day_no, 28), model = "linear")+
                  f(ID, model = "iid3d", n=3*N)+f(time_random, day_no, copy="ID")+
                  f(bp1_random, b1(day_no, 28), copy="ID"), 
                data = INLA_data, family = "gaussian", verbose=FALSE, 
                control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))
  
  coefficients<-LMM_INLA$summary.fixed[,c(1,3,5)]
  Cov_0_t<-obtainVarCov(LMM_INLA, c(2:7), 3)
  MSE<-list(MSE_train=calculate_MSE(LMM_INLA$summary.fitted.values[indexes$train,1], data$y[indexes$train]), 
            MSE_test_same=calculate_MSE(LMM_INLA$summary.fitted.values[indexes$test_same ,1], 
                                        data$y[indexes$test_same]), 
            MSE_test_others=calculate_MSE(LMM_INLA$summary.fitted.values[indexes$test,1], 
                                          data$y[indexes$test]))
  
  return(list(coefficients=coefficients, mlik=LMM_INLA$mlik[1], dic_approx=-2*LMM_INLA$mlik[1], dic=LMM_INLA$dic$dic, 
              waic=LMM_INLA$waic$waic, pit=ks.test(LMM_INLA$cpo$pit,"punif",0,1)$statistic,
              cpo=-sum(log(LMM_INLA$cpo$cpo), na.rm=TRUE), MSE=MSE, model_specific=list(cov_both=Cov_0_t)
              ))
}

Fit_JMM<-function(data, indexes){
  JMM_data<-ready_data(data, indexes)
  N_obs<-length(JMM_data$ID)
  N<-length(unique(JMM_data$ID))
  #Abbreviations: c=cytokine, s=severity score
  fixed.effects<-list(Intercept_c=c(rep(1, N_obs), rep(NA, N_obs)),
                      t_c=c(JMM_data$day_no, rep(NA, N_obs)),
                      t28_c=c(b1(JMM_data$day_no, 28), rep(NA, N_obs)),
                      Intercept_s=c(rep(NA, N_obs), rep(1, N_obs)),
                      t_s=c(rep(NA, N_obs), JMM_data$day_no),
                      t28_s=c(rep(NA, N_obs), b1(JMM_data$day_no, 28)),
                      t=c(JMM_data$day_no, JMM_data$day_no),
                      t28=c(b1(JMM_data$day_no, 28), b1(JMM_data$day_no, 28)))
  random.effects<-list(Intercept_random_c=c(JMM_data$ID, rep(NA, N_obs)),
                       t_random_c=c(JMM_data$ID, rep(NA, N_obs)),
                       t28_random_c=c(JMM_data$ID, rep(NA, N_obs)),
                       Intercept_random_s=c(rep(NA, N_obs), JMM_data$ID+N),
                       t_random_s=c(rep(NA, N_obs), JMM_data$ID+N),
                       t28_random_s=c(rep(NA, N_obs), JMM_data$ID+N))
  
  INLA_data<-c(fixed.effects, random.effects)
  INLA_data$Y<-list(c(JMM_data$exo, rep(NA, N_obs)),
                    c(rep(NA, N_obs),JMM_data$y))
  
  INLA_formula=Y~-1+
    #Fixed effects
    f(Intercept_c, model='linear')+  f(t_c, model='linear')+
    f(t28_c, model='linear')+ f(Intercept_s, model='linear')+f(t_s, model='linear')+
    f(t28_s, model='linear')+
    #Random effects, pairwise jointly distributed
    f(Intercept_random_c, model="iid2d", n=2*N)+ f(Intercept_random_s, copy="Intercept_random_c")+
    f(t_random_c,t, model="iid2d", n=2*N)+ f(t_random_s, t, copy="t_random_c")+
    f(t28_random_c, t28,  model="iid2d", n=2*N)+f(t28_random_s, t28, copy="t28_random_c")
  
  JMM_INLA<-inla(INLA_formula, family = c("gaussian","gaussian"),
                data = INLA_data,verbose=FALSE,
                control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config=TRUE))
  
  #Sampling from the joint distribution of the variance-covariance elements to obtain the association coefficient
  temp<-inla.hyperpar.sample(1000000, JMM_INLA)[,c(1:11)]
  temp_mu_2<-rowSums(1/temp[,c(1,3, 6, 9)])
  temp_mu_1<-rowSums(1/temp[,c(2,4,7,10)])
  cov_temp<-temp[,5]*sqrt(1/temp[,3])*sqrt(1/temp[,4])+temp[,8]*sqrt(1/temp[,6])*sqrt(1/temp[,7])+
    temp[,11]*sqrt(1/temp[,10])*sqrt(1/temp[,9])
  coefficient_sampled<-cov_temp/temp_mu_2
  cor_temp<-cov_temp/(temp_mu_2*temp_mu_1)
  
  #Obtaining all random effects variance-covariance matrices
  Cov_0<-obtainVarCov(JMM_INLA, c(3:5), 2)
  rownames(Cov_0$Cov)<-colnames(Cov_0$Cov)<-c('C_0', 'S_0')
  Cov_t<-obtainVarCov(JMM_INLA, c(6:8), 2)
  rownames(Cov_t$Cov)<-colnames(Cov_t$Cov)<-c('C_t', 'S_t')
  Cov_t28<-obtainVarCov(JMM_INLA, c(9:11), 2)
  rownames(Cov_t28$Cov)<-colnames(Cov_t28$Cov)<-c('C_t28', 'S_t28')
  Cov_error<-diag(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),JMM_INLA$internal.marginals.hyperpar[[1]]), silent=TRUE)$mean,
                    inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),JMM_INLA$internal.marginals.hyperpar[[2]]), silent=TRUE)$mean))
  MSE<-calculate_all_MSE(JMM_INLA, data, indexes)
  
  return(list(coefficients=JMM_INLA$summary.fixed[,c(1,3,5)], mlik=JMM_INLA$mlik[1], dic_approx=-2*JMM_INLA$mlik[1], 
              dic=sum(JMM_INLA$dic$local.dic[c((length(data$ID)+1): (2*length(data$ID)))], na.rm=TRUE), 
              waic=sum(JMM_INLA$waic$local.waic[c((length(data$ID)+1): (2*length(data$ID)))], na.rm=TRUE), 
              pit=ks.test(JMM_INLA$cpo$pit[c((length(data$ID)+1): (2*length(data$ID)))],"punif",0,1)$statistic,
              cpo=-sum(log(JMM_INLA$cpo$cpo[c((length(data$ID)+1): (2*length(data$ID)))]), na.rm=TRUE), MSE=MSE, 
              model_specific=list(cov_mat=list(Cov_0=Cov_0, Cov_t=Cov_t, Cov_t28=Cov_t28, Cov_error=Cov_error),
                                  Association_coefficient=list(coefficient_sampled=quantile(coefficient_sampled, probs=c(0.025, 0.5, 0.975), na.rm=TRUE),
                                                       cor_sampled=quantile(cor_temp, probs=c(0.025,0.5, 0.975), na.rm=TRUE)))
              ))
}

Fit_JSM<-function(data, lagg, indexes){
  JSM_data<-ready_data(data, indexes)
  
  #Lagged values up to lagg of order 5
  JSM_data$t_m1<-NA
  JSM_data$t_m2<-NA
  JSM_data$t_m3<-NA
  JSM_data$t_m4<-NA
  JSM_data$t_m5<-NA
  #Determine the lagged values
  for (i in 2:length(JSM_data$ID)){
    for (j in (i-5):(i-1)){
      if (j>0){
        if(!is.na(JSM_data$day_no[i]) & !is.na(JSM_data$day_no[j])) {
          temp<-JSM_data$day_no[i]-JSM_data$day_no[j]
          if (temp<=5 & temp>0 & JSM_data$ID[i]==JSM_data$ID[j]){
            JSM_data[i,11+temp]<-JSM_data$day_no[j]
          }  
        }  
      }
    }
  }
  
  N_obs<-length(JSM_data$ID)
  N<-length(unique(JSM_data$ID))
  
  #Abbreviations: c=cytokine, s=severity score
  fixed.effects<-list(#Overall values
    t=c(JSM_data$day_no, JSM_data$day_no),
    t28=c(b1(JSM_data$day_no, 28), b1(JSM_data$day_no, 28)),
    
    #Lagged values
    t_m1=c(JSM_data$t_m1, JSM_data$t_m1),
    t28_m1=c(b1(JSM_data$t_m1,28), b1(JSM_data$t_m1,28)),
    
    t_m2=c(JSM_data$t_m2, JSM_data$t_m2),
    t28_m2=c(b1(JSM_data$t_m2,28), b1(JSM_data$t_m2,28)),
    
    t_m3=c(JSM_data$t_m3, JSM_data$t_m3),
    t28_m3=c(b1(JSM_data$t_m3,28), b1(JSM_data$t_m3,28)),
    
    t_m4=c(JSM_data$t_m4, JSM_data$t_m4),
    t28_m4=c(b1(JSM_data$t_m4,28), b1(JSM_data$t_m4,28)),
    
    t_m5=c(JSM_data$t_m5, JSM_data$t_m5),
    t28_m5=c(b1(JSM_data$t_m5,28), b1(JSM_data$t_m5,28)),
    
    #Severity Score fixed effects
    Intercept_s=c(rep(NA, N_obs), rep(1, N_obs)),
    t_s=c(rep(NA, N_obs),JSM_data$day_no),
    
    t28_s=c(rep(NA, N_obs), b1(JSM_data$day_no, 28)))
  
  random.effects<-list(#Fixed effects for cytokine
    Intercept_c=c(rep(1, N_obs), rep(NA, N_obs)),
    t_c=c(rep(1, N_obs),rep(NA, N_obs)),
    t28_c=c(rep(1, N_obs),rep(NA, N_obs)),
    
    #Scaled Cytokine fixed effects
    Intercept_c_scaled=c(rep(NA, N_obs), rep(1, N_obs)),
    t_c_scaled=c(rep(NA, N_obs), rep(1, N_obs)),
    t28_c_scaled=c(rep(NA, N_obs),rep(1, N_obs)),
    
    
    #FIXED EFFECTS FOR LAGGED TERMS!
    #Fixed effects for lagg order 1
    Intercept_1_m1=c(rep(NA, N_obs), rep(1, N_obs)),
    t_1_m1=c(rep(NA, N_obs),rep(1, N_obs)),
    t28_1_m1=c(rep(NA, N_obs),rep(1, N_obs)),
    
    #Fixed effects for lagg order 2
    Intercept_1_m2=c(rep(NA, N_obs), rep(1, N_obs)),
    t_1_m2=c(rep(NA, N_obs),rep(1, N_obs)),
    t28_1_m2=c(rep(NA, N_obs),rep(1, N_obs)),
    
    #Fixed effects for lagg order 3
    Intercept_1_m3=c(rep(NA, N_obs), rep(1, N_obs)),
    t_1_m3=c(rep(NA, N_obs),rep(1, N_obs)),
    t28_1_m3=c(rep(NA, N_obs),rep(1, N_obs)),
    
    #Fixed effects for lagg order 4
    Intercept_1_m4=c(rep(NA, N_obs), rep(1, N_obs)),
    t_1_m4=c(rep(NA, N_obs),rep(1, N_obs)),
    t28_1_m4=c(rep(NA, N_obs),rep(1, N_obs)),
    
    #Fixed effects for lagg order 5
    Intercept_1_m5=c(rep(NA, N_obs), rep(1, N_obs)),
    t_1_m5=c(rep(NA, N_obs),rep(1, N_obs)),
    t28_1_m5=c(rep(NA, N_obs),rep(1, N_obs)),
    
    
    #Random effects for cytokine
    Intercept_random_c=c(JSM_data$ID, rep(NA, N_obs)),
    t_random_c=c(JSM_data$ID+N, rep(NA, N_obs)),
    t28_random_c=c(JSM_data$ID+2*N, rep(NA, N_obs)),
    
    #Scaled random effects for cytokine
    Intercept_random_c_scaled=c(rep(NA, N_obs), JSM_data$ID),
    t_random_c_scaled=c(rep(NA, N_obs), JSM_data$ID+N),
    t28_random_c_scaled=c(rep(NA, N_obs), JSM_data$ID+2*N),
    
    
    #RANDOM EFFECTS FOR LAGGED M
    #Random effects for lagg order 1
    Intercept_random_1_m1=c(rep(NA, N_obs), JSM_data$ID),
    t_random_1_m1=c(rep(NA, N_obs), JSM_data$ID+N),
    t28_random_1_m1=c(rep(NA, N_obs), JSM_data$ID+2*N),
    
    #Random effects for lagg order 2
    Intercept_random_1_m2=c(rep(NA, N_obs), JSM_data$ID),
    t_random_1_m2=c(rep(NA, N_obs), JSM_data$ID+N),
    t28_random_1_m2=c(rep(NA, N_obs), JSM_data$ID+2*N),
    
    #Random effects for lagg order 3
    Intercept_random_1_m3=c(rep(NA, N_obs), JSM_data$ID),
    t_random_1_m3=c(rep(NA, N_obs), JSM_data$ID+N),
    t28_random_1_m3=c(rep(NA, N_obs), JSM_data$ID+2*N),
    
    #Random effects for lagg order 4
    Intercept_random_1_m4=c(rep(NA, N_obs), JSM_data$ID),
    t_random_1_m4=c(rep(NA, N_obs), JSM_data$ID+N),
    t28_random_1_m4=c(rep(NA, N_obs), JSM_data$ID+2*N),
    
    #Random effects for lagg order 5
    Intercept_random_1_m5=c(rep(NA, N_obs), JSM_data$ID),
    t_random_1_m5=c(rep(NA, N_obs), JSM_data$ID+N),
    t28_random_1_m5=c(rep(NA, N_obs), JSM_data$ID+2*N),
    
    
    #Random effects for severity score
    Intercept_random_s=c(rep(NA, N_obs), JSM_data$ID),
    t_random_s=c(rep(NA, N_obs), JSM_data$ID+N),
    t28_random_s=c(rep(NA, N_obs), JSM_data$ID+2*N))
  
  
  INLA_data<-c(fixed.effects, random.effects)
  INLA_data$Y<-list(c(JSM_data$exo, rep(NA, N_obs)),
                    c(rep(NA, N_obs),JSM_data$y)) 
  
  INLA_formula=Y~-1+
    #Fixed effects for cytokine
    f(Intercept_c)+f(t_c, t)+f(t28_c, t28)+
    #Scaled fixed effects for cytokine
    f(Intercept_c_scaled, copy="Intercept_c", model='iid', fixed=FALSE, param=c(0,0.01), initial = 0.1)+
    f(t_c_scaled, t, copy="t_c", same.as = 'Intercept_c_scaled',fixed=FALSE)+
    f(t28_c_scaled, t28, copy="t28_c", same.as = 'Intercept_c_scaled',fixed=FALSE)+
    #RANDOM EFFECTS
    #Random effects for cytokine
    f(Intercept_random_c, model="iid3d", n=3*N)+f(t_random_c, t, copy="Intercept_random_c")+
    f(t28_random_c, t28, copy="Intercept_random_c")+
    #Scaled random effects for cytokine
    f(Intercept_random_c_scaled, copy='Intercept_random_c', same.as='Intercept_c_scaled', fixed=FALSE)+
    f(t_random_c_scaled, t, copy='Intercept_random_c', same.as='Intercept_c_scaled', fixed=FALSE)+
    f(t28_random_c_scaled, t28, copy='Intercept_random_c', same.as='Intercept_c_scaled', fixed=FALSE)+
    #y fixed effects
    f(Intercept_s, model='linear')+f(t_s, model='linear')+f(t28_s, model='linear')+
    #y random effects
    f(Intercept_random_s, model="iid3d", n=3*N)+f(t_random_s, t, copy="Intercept_random_s")+
    f(t28_random_s, t28, copy="Intercept_random_s")
  
  #Elaborate formula if lagg is needed
  if (lagg>=1){
    INLA_formula=update(INLA_formula, ~.+
                           #Copied fixed effects m1
                           f(Intercept_1_m1, copy="Intercept_1", model='iid', fixed=FALSE, param=c(0,0.01), initial = 0.1)+
                           f(t_1_m1, t_m1, copy="t_1", same.as = 'Intercept_1_m1',fixed=FALSE)+
                           f(t28_1_m1, t28_m1, copy="t28_1", same.as = 'Intercept_1_m1',fixed=FALSE)+
                           #Random effects for m1
                           f(Intercept_random_1_m1, copy='Intercept_random_1', same.as='Intercept_1_m1', fixed=FALSE)+
                           f(t_random_1_m1, t_m1, copy='Intercept_random_1', same.as='Intercept_1_m1', fixed=FALSE)+
                           f(t28_random_1_m1, t28_m1, copy='Intercept_random_1', same.as='Intercept_1_m1', fixed=FALSE))
  }
  if (lagg>=2){
    INLA_formula=update(INLA_formula, ~.+
                           #Copied fixed effects m2
                           f(Intercept_1_m2, copy="Intercept_1", model='iid', fixed=FALSE, param=c(0,0.01), initial = 0.1)+
                           f(t_1_m2, t_m2, copy="t_1", same.as = 'Intercept_1_m2',fixed=FALSE)+
                           f(t28_1_m2, t28_m2, copy="t28_1", same.as = 'Intercept_1_m2',fixed=FALSE)+
                           #Random effects for m2
                           f(Intercept_random_1_m2, copy='Intercept_random_1', same.as='Intercept_1_m2', fixed=FALSE)+
                           f(t_random_1_m2, t_m2, copy='Intercept_random_1', same.as='Intercept_1_m2', fixed=FALSE)+
                           f(t28_random_1_m2, t28_m2, copy='Intercept_random_1', same.as='Intercept_1_m2', fixed=FALSE))
  }
  if (lagg>=3){
    INLA_formula=update(INLA_formula, ~.+
                           #Copied fixed effects m3
                           f(Intercept_1_m3, copy="Intercept_1", model='iid', fixed=FALSE, param=c(0,0.01), initial = 0.1)+
                           f(t_1_m3, t_m3, copy="t_1", same.as = 'Intercept_1_m3',fixed=FALSE)+
                           f(t28_1_m3, t28_m3, copy="t28_1", same.as = 'Intercept_1_m3',fixed=FALSE)+
                           #Random effects for m3
                           f(Intercept_random_1_m3, copy='Intercept_random_1', same.as='Intercept_1_m3', fixed=FALSE)+
                           f(t_random_1_m3, t_m3, copy='Intercept_random_1', same.as='Intercept_1_m3', fixed=FALSE)+
                           f(t28_random_1_m3, t28_m3, copy='Intercept_random_1', same.as='Intercept_1_m3', fixed=FALSE))
  }
  if (lagg>=4){
    INLA_formula=update(INLA_formula, ~.+
                           #Copied fixed effects m4
                           f(Intercept_1_m4, copy="Intercept_1", model='iid', fixed=FALSE, param=c(0,0.01), initial = 0.1)+
                           f(t_1_m4, t_m4, copy="t_1", same.as = 'Intercept_1_m4',fixed=FALSE)+
                           f(t28_1_m4, t28_m4, copy="t28_1", same.as = 'Intercept_1_m4',fixed=FALSE)+
                           #Random effects for m4
                           f(Intercept_random_1_m4, copy='Intercept_random_1', same.as='Intercept_1_m4', fixed=FALSE)+
                           f(t_random_1_m4, t_m4, copy='Intercept_random_1', same.as='Intercept_1_m4', fixed=FALSE)+
                           f(t28_random_1_m4, t28_m4, copy='Intercept_random_1', same.as='Intercept_1_m4', fixed=FALSE))
  }
  if (lagg>=5){
    INLA_formula=update(INLA_formula, ~.+
                           #Copied fixed effects m5
                           f(Intercept_1_m5, copy="Intercept_1", model='iid', fixed=FALSE, param=c(0,0.01), initial = 0.1)+
                           f(t_1_m5, t_m5, copy="t_1", same.as = 'Intercept_1_m5',fixed=FALSE)+
                           f(t28_1_m5, t28_m5, copy="t28_1", same.as = 'Intercept_1_m5',fixed=FALSE)+
                           #Random effects for m5
                           f(Intercept_random_1_m5, copy='Intercept_random_1', same.as='Intercept_1_m5', fixed=FALSE)+
                           f(t_random_1_m5, t_m5, copy='Intercept_random_1', same.as='Intercept_1_m5', fixed=FALSE)+
                           f(t28_random_1_m5, t28_m5, copy='Intercept_random_1', same.as='Intercept_1_m5', fixed=FALSE))
  }
    
    
  JSM_INLA<-inla(INLA_formula, family = c("gaussian","gaussian"),
                    data = INLA_data, verbose=FALSE, control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))
  
  coefficients<-rbind(JSM_INLA$summary.random$Intercept_1[c(2,4,6)],
                        JSM_INLA$summary.random$t_1[c(2,4,6)],
                        JSM_INLA$summary.random$t28_1[c(2,4,6)],
                        JSM_INLA$summary.fixed[c(1,2,3),c(1,3,5)])
  gamma<-JSM_INLA$summary.hyperpar[c(18:length(JSM_INLA$summary.hyperpar$mean)), c(1,3,5)]
  Cov_C<-obtainVarCov(JSM_INLA, c(6:11), 3)
  Cov_S<-obtainVarCov(JSM_INLA, c(12:17), 3)
  Cov_error<-diag(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),JSM_INLA$internal.marginals.hyperpar[[1]]), silent=TRUE)$mean,
                    inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),JSM_INLA$internal.marginals.hyperpar[[2]]), silent=TRUE)$mean))
  MSE<-calculate_all_MSE(JSM_INLA, data, indexes)
  
  return(list(coefficients=coefficients, mlik=JSM_INLA$mlik[1], dic_approx=-2*JSM_INLA$mlik[1], 
              dic=sum(JSM_INLA$dic$local.dic[c((length(JSM_data$ID)+1): (2*length(JSM_data$ID)))], na.rm=TRUE), 
              waic=sum(JSM_INLA$waic$local.waic[c((length(JSM_data$ID)+1): (2*length(JSM_data$ID)))], na.rm=TRUE), 
              pit=ks.test(JSM_INLA$cpo$pit[c((length(JSM_data$ID)+1): (2*length(JSM_data$ID)))],"punif",0,1)$statistic,
              cpo=-sum(log(JSM_INLA$cpo$cpo[c((length(JSM_data$ID)+1): (2*length(JSM_data$ID)))]), na.rm=TRUE), MSE=MSE,
              model_specific=list(cov_mat=list(Cov_C=Cov_C, Cov_S=Cov_S, Cov_error=Cov_error), 
                                  gamma=gamma)
              ))
}