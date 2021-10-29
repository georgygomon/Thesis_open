library('INLA')
library('nlme')
library(tidyr)
library(reshape2)
library(lme4)
library(nlme)
library(MCMCglmm)
library(lme4)
library(ggplot2)
library(reshape2)
library(MCMCvis)
library(MASS)
library(mvtnorm)
library(xtable)

make_model_0<-function(true_betas, U1, U2, N, n){
  data_0 <- data.frame(
    id=rep(1:N, each=n),
    x=rnorm(N*n),
    Period=rep(0:(n-1), times=N), 
    Period2=factor(rep(0:(n-1), times=N)))
  errors1 <- matrix(rnorm(N*n), N)
  errors2 <- matrix(rnorm(N*n), N)
  u1<-matrix(rnorm(N*2), N) %*% chol(U1)
  u2<-matrix(rnorm(N*2), N) %*% chol(U2)
  data_0$y1 <-  true_betas[1]+rep(u1[,1], each=n)+true_betas[2]*data_0$x+(true_betas[3]+rep(u1[,2], each=n))*data_0$Period+as.vector(t(errors1))
  data_0$y2 <- true_betas[4]+rep(u2[,1], each=n)+true_betas[5]*data_0$x+(true_betas[6]+rep(u2[,2], each=n))*data_0$Period+as.vector(t(errors2))
  return(list(data=data_0, true_cov_1=cov(u1), true_cov_2=cov(u2)))
}

make_model_1A<-function(true_betas, Sigma, N, n){
  ################################################################################
  #Only Errors dependent
  ################################################################################
  set.seed(2022)
  data<- data.frame(
    id=rep(1:N, each=n),
    x=rnorm(N*n), 
    Period=rep(0:(n-1), times=N), 
    Period2=factor(rep(0:(n-1), times=N)))
  
  errors <- matrix(rnorm(N*2*n), N) %*% chol(Sigma)
  data$y1 <- true_betas[1] + true_betas[2]*data$x + true_betas[3]*data$Period+ as.vector(t(errors[,1:n]))
  data$y2<- true_betas[4]+true_betas[5]*data$x+ true_betas[6]*data$Period+ as.vector(t(errors[,(n+1):(2*n)]))
  return(list(data=data, true_errors=cov(errors)))
}

make_model_2A<-function(true_betas, U, N, n){
  ################################################################################
  #Only Random intercept
  ################################################################################
  set.seed(2022)
  data_2a<- data.frame(
    id=rep(1:N, each=n),
    x=rnorm(N*n), 
    Period=rep(0:(n-1), times=N), 
    Period2=factor(rep(0:(n-1), times=N)))
  
  errors1 <- matrix(rnorm(N*n), N)
  errors2 <- matrix(rnorm(N*n), N)
  u<-matrix(rnorm(N*2), N) %*% chol(U)
  data_2a$y1 <- true_betas[1] + true_betas[2]*data_2a$x + true_betas[3]*data_2a$Period+rep(u[,1], each=n)+as.vector(t(errors1))
  data_2a$y2<- true_betas[4]+true_betas[5]*data_2a$x+true_betas[6]*data_2a$Period+rep(u[,2], each=n)+as.vector(t(errors2))
  return(list(data=data_2a, true_cov=cov(u)))
}

make_model_2B<-function(true_betas, U, Sigma, N, n){
  ################################################################################
  #Random intercept & dependent errors
  ################################################################################
  set.seed(2022)
  data_2b<- data.frame(
    id=rep(1:N, each=n),
    x=rnorm(N*n), 
    Period=rep(0:(n-1), times=N), 
    Period2=factor(rep(0:(n-1), times=N)))
  errors <- matrix(rnorm(N*2*n), N) %*% chol(Sigma)
  u<-matrix(rnorm(N*2), N) %*% chol(U)
  data_2b$y1 <- true_betas[1] + true_betas[2]*data_2b$x +true_betas[3]*data_2b$Period+rep(u[,1], each=n)+as.vector(t(errors[,1:n]))
  data_2b$y2<- true_betas[4]+true_betas[5]*data_2b$x+true_betas[6]*data_2b$Period+rep(u[,2], each=n)+as.vector(t(errors[,(n+1):(2*n)]))
  return(list(data=data_2b, true_cov=cov(u), true_errors=cov(errors)))
}

make_model_2C<-function(true_betas, U, N, n){
  ################################################################################
  #Random slope & intercept: Dependent
  ################################################################################
  set.seed(2022)
  data_2c <- data.frame(
    id=rep(1:N, each=n),
    x=rnorm(N*n),  
    Period=rep(0:(n-1), times=N), 
    Period2=factor(rep(0:(n-1), times=N)))
  
  errors1 <- matrix(rnorm(N*n), N)
  errors2 <- matrix(rnorm(N*n), N)
  u<-matrix(rnorm(N*2*2), N) %*% chol(U)
  data_2c$y1 <- true_betas[1] +rep(u[,1], each=n)+ true_betas[2]*data_2c$x + true_betas[3]*data_2c$Period+
    rep(u[,3], each=n)*data_2c$Period+as.vector(t(errors1))
  data_2c$y2<- true_betas[4]+rep(u[,2], each=n)+true_betas[5]*data_2c$x+true_betas[6]*data_2c$Period+
    rep(u[,4], each=n)*data_2c$Period+as.vector(t(errors2))
  return(list(data=data_2c, true_cov=cov(u)))
}

make_model_3A<-function(true_betas, U1, U2, gamma, N, n){
  ################################################################################
  #Random slope & intercept: Dependent
  ################################################################################
  set.seed(2022)
  true_betas[4:6]<-true_betas[4:6]-gamma*true_betas[1:3]
  data_3a <- data.frame(
    id=rep(1:N, each=n),
    x=rnorm(N*n),
    Period=rep(0:(n-1), times=N), 
    Period2=factor(rep(0:(n-1), times=N)))
  
  errors1 <- matrix(rnorm(N*n), N)
  errors2 <- matrix(rnorm(N*n), N)
  u1<-matrix(rnorm(N*2), N) %*% chol(U1)
  u2<-matrix(rnorm(N*2), N) %*% chol(U2)
  m_3a<-true_betas[1]+rep(u1[,1], each=n)+true_betas[2]*data_3a$x+(true_betas[3]+rep(u1[,2], each=n))*data_3a$Period
  data_3a$y1 <-  m_3a+as.vector(t(errors1))
  data_3a$y2<- gamma*(m_3a)+true_betas[4]+rep(u2[,1], each=n)+true_betas[5]*data_3a$x+(true_betas[6]+rep(u2[,2], each=n))*data_3a$Period+as.vector(t(errors2))
  return(list(data=data_3a, true_cov_u1=cov(u1), true_cov_u2=cov(u2)))
}

make_model_3B<-function(true_betas, U1, U2, gamma1, gamma2, N, n){
  ################################################################################
  #Random slope & intercept: Dependent, cloned apart
  ################################################################################
  set.seed(2022)
  data_3b <- data.frame(
    id=rep(1:N, each=n),
    x=rnorm(N*n),
    Period=rep(0:(n-1), times=N), 
    Period2=factor(rep(0:(n-1), times=N)))
  errors1 <- matrix(rnorm(N*n), N)
  errors2 <- matrix(rnorm(N*n), N)
  u1<-matrix(rnorm(N*2), N) %*% chol(U1)
  u2<-matrix(rnorm(N*2), N) %*% chol(U2)
  data_3b$y1 <-  true_betas[1]+rep(u1[,1], each=n)+true_betas[2]*data_3b$x+(true_betas[3]+rep(u1[,2], each=n))*data_3b$Period+as.vector(t(errors1))
  data_3b$y2<- gamma1*rep(u1[,1], each=n)+true_betas[4]+rep(u2[,1], each=n)+true_betas[5]*data_3b$x+
    (true_betas[6]+rep(u2[,2], each=n)+gamma2*rep(u1[,2], each=n))*data_3b$Period+as.vector(t(errors2))
  return(list(data=data_3b, true_cov_u1=cov(u1), true_cov_u2=cov(u2)))
}


# #Same for all models
# true_betas<-c(2,4,2.5,3, 1.5, 3.5)
# N<-900
# n<-2
# #Model 1A special Parameters
# Sigma1 <- matrix(c(5,2,2,4), 2)
# Sigma2<-matrix(c(4,1,1,4),2)
# Sigma3<-matrix(c(-1, -0.5, -0.5, -1.5), 2)
# Sigma<-rbind(cbind(Sigma1, Sigma3), cbind(t(Sigma3), Sigma2))
# #Model 2A special parameters
# U_0<-matrix(c(2,-1,
#               -1,3),2, byrow=TRUE)
# #Model 2B special parameters
# #U_0 and Sigma
# #Model 2C special parameters
# U0<-matrix(c(2,-1,
#              -1,3),2, byrow=TRUE)
# Ut<-matrix(c(3,1.5,
#              1.5,4),2, byrow=TRUE)
# U0t<-matrix(c(0,0,
#               0,0),2, byrow=TRUE)
# U<-rbind(cbind(U0, U0t), cbind(t(U0t), Ut))
# #Model 3A Special parameters
# U1<-matrix(c(2,0,
#              0,3),2, byrow=TRUE)
# U2<-matrix(c(3,1.5,
#              1.5,4),2, byrow=TRUE)
# gamma<-0.45
# #Model 3B Special parameters
# U1<-matrix(c(2,0,
#              0,3),2, byrow=TRUE)
# U2<-matrix(c(3,1.5,
#              1.5,4),2, byrow=TRUE)
# gamma1<-0.45
# gamma2<-1.2

