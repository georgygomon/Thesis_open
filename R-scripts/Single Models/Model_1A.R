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

matrix2latex <- function(matr) {
  
  printmrow <- function(x) {
    
    cat(cat(x,sep=" & "),"\\\\ \n")
  }
  
  cat("\\begin{bmatrix}","\n")
  body <- apply(matr,1,printmrow)
  cat("\\end{bmatrix}")
}

#Parameter list
true_betas<-c(2,4,2.5,3, 1.5, 3.5)
Sigma1 <- matrix(c(5,2,2,4), 2)
Sigma2<-matrix(c(4,1,1,4),2)
Sigma3<-matrix(c(-1, -0.5, -0.5, -1.5), 2)
Sigma<-rbind(cbind(Sigma1, Sigma3), cbind(t(Sigma3), Sigma2))

coefficients_1<-list(coefficients=data.frame(true=true_betas, gls=rep(0,6), gls_sd=rep(0,6),
                                             MCMCglmm=rep(0,6), MCMCglmm_sd=rep(0,6),
                                             INLA=rep(0,6), INLA_sd=rep(0,6)),
                     VarCov=list(true=NA, gls=NA, MCMCglmm=NA, INLA=NA))
rownames(coefficients_1$coefficients)<-c( 'beta_0^1', 'beta_x^1', 'beta_t^1', 'beta_0^2', 'beta_x^2', 'beta_t^2')

################################################################################
#Model 1: Association via residual errors
################################################################################
set.seed(2022)
### Total of N individuals
N <- 1000
### n correlated measurements per individual
n <- 2

data<- data.frame(
  id=rep(1:N, each=n),
  x=rnorm(N*n), 
  Period=rep(0:(n-1), times=N), 
  Period2=factor(rep(0:(n-1), times=N)))

### sampling a covariate setting up data such that the n measurements on 1 individual are correlated
errors <- matrix(rnorm(N*2*n), N) %*% chol(Sigma)
coefficients_1$VarCov$true<-cov(errors)
#matrix2latex(round(coefficients_1$VarCov$true, digits=2))
data$y1 <- true_betas[1] + true_betas[2]*data$x + true_betas[3]*data$Period+ as.vector(t(errors[,1:n]))
data$y2<- true_betas[4]+true_betas[5]*data$x+ true_betas[6]*data$Period+ as.vector(t(errors[,(n+1):(2*n)]))



#This is very important!!!!!!!!!
data<-rbind(data[data$Period==0,], data[data$Period==1,])

################################################################################
#Model 1: gls
################################################################################
#First melt data
data_melted<-melt(data, measure.vars=c('y1', 'y2'))
data_melted<-data_melted[with(data_melted,order(id)),]
#Fit gls model
gls_unstructured<-gls(value~variable+variable:x+variable:Period-1,
                      correlation=corSymm(form=~1|id), 
                      weights= varIdent(form=~1|variable*Period),
                      data=data_melted, na.action=na.omit, method='REML')

coefficients_1$coefficients[,c(2,3)]<-summary(gls_unstructured)$tTable[c(1,3,5, 2,4,6),c(1:2)]
coefficients_1$VarCov$gls<-getVarCov(gls_unstructured, individual=1)

################################################################################
#Model 1: MCMCglmm
################################################################################
prior = list(B = list(mu = rep(0,6),V = 100*diag(6)), 
             R = list(V = diag(4), n = 1))

MCMC_1<-MCMCglmm(cbind(y1, y2) ~ trait+trait:x+trait:Period - 1,
                 rcov = ~ us(Period2:trait):units, 
                 family = rep("gaussian", 2), nitt = 10000, burnin = 1000,
                 prior=prior, thin=25, data = data)
coefficients_1$coefficients[,4:5]<-MCMCsummary(MCMC_1$Sol)[c(1,3,5, 2,4,6),c(1,2)]
coefficients_1$VarCov$MCMCglmm<-matrix(MCMCsummary(MCMC_1$VCV)[,1], nrow=4, byrow=FALSE)

################################################################################
#Model 1: INLA
################################################################################
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
#One can get rid of the Gaussian noise term by setting fixed=TRUE in the control family parameter
#https://www.r-inla.org/faq question 5
control_family = list(list(initial=15, fixed=T), list(initial=15, fixed=T))
INLA_formula_1=Y~-1+f(Intercept1, model='linear')+f(Intercept2, model='linear')+
  f(x1, model='linear')+f(x2, model='linear')+
  f(Period1, model='linear')+f(Period2, model='linear')+
  f(i, model="iid4d", n=2*n*N, hyper = list(theta1 = list(prior = "wishart4d")))
INLA_model_1<-inla(INLA_formula_1, family = c("gaussian","gaussian"),
                  data = INLA_data, verbose=TRUE, control.family=control_family,
                  control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))



coefficients_1$coefficients[,c(6:7)]<-INLA_model_1$summary.fixed[c(1,3,5, 2,4, 6),c(1:2)]
SD_s<-diag(4)
for (i in 1:4){
  SD_s[i,i]<-sqrt(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),INLA_model_1$internal.marginals.hyperpar[[i]]))$quant0.5)
}
cor<-matrix(1, nrow=4, ncol=4)
cor[lower.tri(cor)]<-INLA_model_1$summary.hyperpar[5:10,1]
cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
coefficients_1$VarCov$INLA<-SD_s %*% cor %*% SD_s



lapply(coefficients_1$VarCov,  round, digits=2)
coefficients_1$coefficients
 