if(!require("xtable")){install.packages("xtable"); library("xtable")}
if(!require("mvtnorm")){install.packages("mvtnorm"); library("mvtnorm")}
if(!require("MASS")){install.packages("MASS"); library("MASS")}
if(!require("ggplot2")){install.packages("ggplot2"); library("ggplot2")}
if(!require("tidyr")){install.packages("tidyr"); library("tidyr")}
if(!require("reshape2")){install.packages("reshape2"); library("reshape2")}
if(!require("optimx")){install.packages("optimx"); library("optimx")}
if(!require("inlabru")){install.packages("inlabru"); library("inlabru")}
if(!require("Rlab")){install.packages("Rlab"); library("Rlab")}
if(!require("dgof")){install.packages("dgof"); library("dgof")}
#To install the INLA-package in R, you have to manually add the r-inla repository as it is not on CRAN.
#https://www.r-inla.org/download-install
library('INLA')

get_indexes<-function(n, n_share, N, N_extra){
  N_index<-c(N, N_extra)
  n_index<-c(n, round(n_share*n))
  n<-sum(n_index)
  N<-sum(N_index)
  train<-c(1:n_index[1])+ rep(seq(0, N_index[1]*n-1, by=n), each=n_index[1])
  test_same<-c((n_index[1]+1):(sum(n_index[1:2])))+ rep(seq(0, N_index[1]*n-1, by=n), each=n_index[2])
  test_others<-c((N_index[1]*n+1):(N*n))
  indexes<-list(train=train, test_same=test_same, 
                test_others=test_others, 
                N_index=N_index, n_index=n_index, 
                n=n, N=N)
  return(indexes)
}

create_data<-function(N, n, p, n_index){
  data <- data.frame(
    id=rep(1:N, each=n),
    w=rep(rnorm(N), each=n),
    t=rep(c(0:(n-1)), N), 
    observed_x=rbern(N*n, p[1]),
    observed_y=rbern(N*n, p[2]),
    x=NA, y=NA)
  return(data)
}

Create_LMM_data<-function(true_betas, parameters, p, indexes, cv){
  LMM_data <- create_data(indexes$N, indexes$n, p, indexes$n_index)
  errors1 <- matrix(rnorm(indexes$N*indexes$n, sd=0.5), indexes$N)
  errors2 <- matrix(rnorm(indexes$N*indexes$n, sd=0.5), indexes$N)
  u1<-matrix(rnorm(indexes$N*2), indexes$N) %*% chol(parameters$D1)
  u2<-matrix(rnorm(indexes$N*2), indexes$N) %*% chol(parameters$D2)
  LMM_data$x <-  true_betas[1]+rep(u1[,1], each=indexes$n)+true_betas[2]*LMM_data$w+
    (true_betas[3]+rep(u1[,2], each=indexes$n))*LMM_data$t+as.vector(t(errors1))
  LMM_data$y<-parameters$beta_x*LMM_data$x +true_betas[4]+rep(u2[,1], each=indexes$n)+true_betas[5]*LMM_data$w+
    (true_betas[6]+rep(u2[,2], each=indexes$n))*LMM_data$t+as.vector(t(errors2))
  return(list(data=LMM_data))
}

Create_JMM_data<-function(true_betas, parameters, p, indexes, cv){
  JMM_data <-create_data(indexes$N, indexes$n, p, indexes$n_index)
  errors1 <- matrix(rnorm(indexes$N*indexes$n, sd=0.5), indexes$N)
  errors2 <- matrix(rnorm(indexes$N*indexes$n, sd=0.5), indexes$N)
  u<-matrix(rnorm(indexes$N*2*2), indexes$N) %*% chol(parameters$D)
  JMM_data$x <- true_betas[1] +rep(u[,1], each=sum(indexes$n_index[1:2]))+ 
    true_betas[2]*JMM_data$w + true_betas[3]*JMM_data$t+
    rep(u[,3], each=indexes$n)*JMM_data$t+as.vector(t(errors1))
  JMM_data$y<- true_betas[4]+rep(u[,2], each=indexes$n)+
    true_betas[5]*JMM_data$w+true_betas[6]*JMM_data$t+
    rep(u[,4], each=indexes$n)*JMM_data$t+as.vector(t(errors2))
  return(list(data=JMM_data, true_cov=cov(u[1:indexes$N_index[1],])))
}

Create_JSM_data<-function(true_betas, parameters, p, indexes, CV){
  true_betas[4:6]<-true_betas[4:6]-parameters$gamma*true_betas[1:3]
  JSM_data <- create_data(indexes$N, indexes$n, p, indexes$n_index)
  errors1 <- matrix(rnorm(indexes$N*indexes$n, sd=0.5), indexes$N)
  errors2 <- matrix(rnorm(indexes$N*indexes$n, sd=0.5), indexes$N)
  u1<-matrix(rnorm(indexes$N*2), indexes$N) %*% chol(parameters$D1)
  u2<-matrix(rnorm(indexes$N*2), indexes$N) %*% chol(parameters$D2)
  m_3a<-true_betas[1]+rep(u1[,1], each=indexes$n)+true_betas[2]*JSM_data$w+(true_betas[3]+rep(u1[,2], each=indexes$n))*JSM_data$t
  JSM_data$x <-  m_3a+as.vector(t(errors1))
  JSM_data$y<- parameters$gamma*(m_3a)+true_betas[4]+rep(u2[,1], each=indexes$n)+true_betas[5]*JSM_data$w+
    (true_betas[6]+rep(u2[,2], each=indexes$n))*JSM_data$t+as.vector(t(errors2))
  return(list(data=JSM_data, actual_cov_u1=cov(u1[1:indexes$N_index[1],]), modelled_cov_u2=cov(u2[1:indexes$N_index[1],]), 
              actual_cov_u2=parameters$gamma^2*cov(u1[1:indexes$N_index[1],])+cov(u2[1:indexes$N_index[1],])))
}

