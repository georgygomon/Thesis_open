if(!require("reshape2")){install.packages("reshape2"); library("reshape2")}
if(!require("lme4")){install.packages("lme4"); library("lme4")}
if(!require("nlme")){install.packages("nlme"); library("nlme")}
if(!require("MCMCglmm")){install.packages("MCMCglmm"); library("MCMCglmm")}
if(!require("MCMCvis")){install.packages("MCMCvis"); library("MCMCvis")}
if(!require("optimx")){install.packages("optimx"); library("optimx")}
if(!require("inlabru")){install.packages("inlabru"); library("inlabru")}
if(!require("Rlab")){install.packages("Rlab"); library("Rlab")}
if(!require("dgof")){install.packages("dgof"); library("dgof")}
#To install the INLA-package in R, you have to manually add the r-inla repository as it is not on CRAN.
#https://www.r-inla.org/download-install
library('INLA')

#Creates general data-structure
create_data<-function(N, n, p){
  data <- data.frame(
    id=rep(1:N, each=n),
    w=rep(rnorm(N), each=n),
    t=rep(c(0:(n-1)), N), 
    observed_x=rbern(N*n, p[1]),
    observed_y=rbern(N*n, p[2]),
    x=NA, y=NA)
  return(data)
}

Create_multivariate_data<-function(true_betas, parameters,N,n, p){
  Multivariate_data <-create_data(N, n, p)
  errors1 <- matrix(rnorm(N*n, sd=0.5), N)
  errors2 <- matrix(rnorm(N*n, sd=0.5), N)
  u<-matrix(rnorm(N*2*2), N) %*% chol(parameters$true)
  Multivariate_data$x <- true_betas[1] +rep(u[,1], each=n)+ 
    true_betas[2]*Multivariate_data$w + true_betas[3]*Multivariate_data$t+
    rep(u[,3], each=n)*Multivariate_data$t+as.vector(t(errors1))
  Multivariate_data$y<- true_betas[4]+rep(u[,2], each=n)+
    true_betas[5]*Multivariate_data$w+true_betas[6]*Multivariate_data$t+
    rep(u[,4], each=n)*Multivariate_data$t+as.vector(t(errors2))
  return(list(data=Multivariate_data))
}

Create_JSM_data<-function(true_betas, parameters,N,n, p){
  JSM_data <- create_data(N, n, p)
  errors1 <- matrix(rnorm(N*n, sd=0.5), N)
  errors2 <- matrix(rnorm(N*n, sd=0.5), N)
  u1<-matrix(rnorm(N*2), N) %*% chol(parameters$D1)
  u2<-matrix(rnorm(N*2), N) %*% chol(parameters$D2)
  m_3a<-true_betas[1]+rep(u1[,1], each=n)+true_betas[2]*JSM_data$w+(true_betas[3]+rep(u1[,2], each=n))*JSM_data$t
  JSM_data$x <-  m_3a+as.vector(t(errors1))
  JSM_data$y<- parameters$gamma*(m_3a)+true_betas[4]+rep(u2[,1], each=n)+true_betas[5]*JSM_data$w+
    (true_betas[6]+rep(u2[,2], each=n))*JSM_data$t+as.vector(t(errors2))
  return(list(data=JSM_data, actual_cov_u1=cov(u1), modelled_cov_u2=cov(u2),
              actual_cov_u2=parameters$gamma^2*cov(u1)+cov(u2)))
}


Create_JMM_data<-function(true_betas, parameters,N,n, p){
  JMM_data <- create_data(N, n, p)
  errors1 <- matrix(rnorm(N*n, sd=0.5), N)
  errors2 <- matrix(rnorm(N*n, sd=0.5), N)
  u<-matrix(rnorm(N*2*2), N) %*% chol(parameters$D)
  JMM_data$x <- true_betas[1] +rep(u[,1], each=n)+ 
    true_betas[2]*JMM_data$w + true_betas[3]*JMM_data$t+
    rep(u[,3], each=n)*JMM_data$t+as.vector(t(errors1))
  JMM_data$y<- true_betas[4]+rep(u[,2], each=n)+
    true_betas[5]*JMM_data$w+true_betas[6]*JMM_data$t+
    rep(u[,4], each=n)*JMM_data$t+as.vector(t(errors2))
  return(list(data=JMM_data))
}

Create_LMM_data<-function(true_betas, parameters,N,n, p){
  LMM_data <- create_data(N, n, p)
  errors1 <- matrix(rnorm(N*n, sd=0.5), N)
  errors2 <- matrix(rnorm(N*n, sd=0.5), N)
  d1<-matrix(rnorm(N*2), N) %*% chol(parameters$D1)
  d2<-matrix(rnorm(N*2), N) %*% chol(parameters$D2)
  LMM_data$x <-  true_betas[4]+rep(d1 [,1], each=n)+true_betas[5]*LMM_data$w+
    (true_betas[6]+rep(d1[,2], each=n))*LMM_data$t+as.vector(t(errors1))
  LMM_data$y<-parameters$beta_v*LMM_data$x +true_betas[1]+rep(d2[,1], each=n)+true_betas[2]*LMM_data$w+
    (true_betas[3]+rep(d2[,2], each=n))*LMM_data$t+as.vector(t(errors2))
  return(list(data=LMM_data))
}

#Obtain random effects variance-covariance matrices
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

