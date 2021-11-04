library('INLA')
library(reshape2)
library(lme4)
library(nlme)
library(MCMCglmm)
library(MCMCvis)
library(optimx)


################################################################################
#Parameters
################################################################################
parameters<-list(U1=matrix(c(2,0,
                             0,3),2, byrow=TRUE),
                 U2=matrix(c(3,1.5,
                              1.5,4),2, byrow=TRUE),
                 gamma1=1.5 ,
                 gamma2=2.1)
true_betas<-c(2,4,2.5,3, 1.5, 3.5)

coefficients_3b1<-list(coefficients=data.frame(true= c(true_betas, parameters$gamma1, parameters$gamma2),  INLA=rep(0,8)),
                      U1=list(true=parameters$U1, INLA=diag(2)),
                      U2=list(true=parameters$U2, INLA=diag(2)))


rownames(coefficients_3b1$coefficients)<-c( 'beta_0^1', 'beta_x^1','beta_t^1', 'beta_0^2', 'beta_x^2','beta_t^2', 'gamma1', 'gamma2')
################################################################################
#Model 2 Make data: Only random intercept, uncorrelated errors
################################################################################
set.seed(2022)
### Total of N individuals
N <- 250
### Max # of measurements per individual
n <- 5
### Probability p that a measurement is observed
p1<-1
### probability p2 that covariate x is measured
p2<-1


data_3b1 <- data.frame(
  id=rep(1:N, each=n),
  x=rnorm(N*n),
  Period=rep(0:(n-1), times=N), 
  Period2=factor(rep(0:(n-1), times=N)),
  observed=rbern(N*n, p1),
  x_measured=rbern(N*n, p2))
errors1 <- matrix(rnorm(N*n), N)
errors2 <- matrix(rnorm(N*n), N)
u1<-matrix(rnorm(N*2), N) %*% chol(parameters$U1)
u2<-matrix(rnorm(N*2), N) %*% chol(parameters$U2)
coefficients_3b1$U1$true<-cov(u1)
coefficients_3b1$U2$true<-cov(u2)
data_3b1$y1 <-  true_betas[1]+rep(u1[,1], each=n)+true_betas[2]*data_3b1$x+(true_betas[3]+rep(u1[,2], each=n))*data_3b1$Period+as.vector(t(errors1))
data_3b1$y2<- parameters$gamma1*rep(u1[,1], each=n)+true_betas[4]+rep(u2[,1], each=n)+true_betas[5]*data_3b1$x+
  (true_betas[6]+rep(u2[,2], each=n)+parameters$gamma2*rep(u1[,2], each=n))*data_3b1$Period+as.vector(t(errors2))
data_3b1$x[data_3b1$x_measured==0]<-NA
data<-data_3b1[data_3b1$observed==1,]



################################################################################
#Model 2a using INLA 
################################################################################
#https://groups.google.com/g/r-inla-discussion-group/c/ClNVlx1lgwY
#https://arxiv.org/pdf/1210.0333.pdf
#One can copy only random effects. If one wants to copy fixed effect simply write it as random effect
#This one contains the answer!!!!!
#https://groups.google.com/g/r-inla-discussion-group/c/bodCLSuCzEw/m/hJMFbidpWAoJ
#https://groups.google.com/g/r-inla-discussion-group/c/ClNVlx1lgwY/m/C2aR-bF9O4YJ
#https://groups.google.com/g/r-inla-discussion-group/c/7I2wHxL6NBM/m/iJagfXdPBQAJ

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

formula.model_3b1=Y~-1+Intercept1+x1+Period1+Intercept2+x2+Period2+ #Model 1 fixed effects
  f(Intercept_random1, model="iid")+ #Copying random effects
  f(Intercept_random12, copy="Intercept_random1",hyper = list(beta = list(fixed=FALSE)))+
  f(Period_random1, Period, model="iid")+
  f(Period_random12, copy="Period_random1",hyper = list(beta = list(fixed=FALSE)))+
  f(Intercept_random2, model="iid2d", n=2*N)+f(Period_random2, Period, copy="Intercept_random2") #Random effects Y2

final_model_3b1<-inla(formula.model_3b1, family = c("gaussian","gaussian"),
                     data = final_data_3b1,verbose=TRUE)

coefficients_3b1$coefficients$INLA<-c(final_model_3b1$summary.fixed[1:6,1],
                                     summary(final_model_3b1)$hyperpar[8:9,1])

coefficients_3b1$coefficients$INLA[7:8]<-gamma_s<-summary(final_model_3b1)$hyperpar[8:9,1]
coefficients_3b1$coefficients$INLA[1:6]<-model_coefficients<-c(final_model_3b1$summary.fixed[1:6,1])


#Y1 random effects covariance structure
SD_s<-diag(2)
tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_3b1$internal.marginals.hyperpar[[3]]), silent=TRUE)$mean,
                      inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_3b1$internal.marginals.hyperpar[[4]]), silent=TRUE)$mean)), 
         error = function(e) print("Model 3B1: Some covariance elements of Y1 not well defined"))
cor<-diag(2)
coefficients_3b1$U1$INLA<-Cov1_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)
#Y2 random effects covariance structure
SD_s<-diag(2)
tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_3b1$internal.marginals.hyperpar[[5]]), silent=TRUE)$mean,
                      inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_3b1$internal.marginals.hyperpar[[6]]), silent=TRUE)$mean)), 
         error = function(e) print("Model 3A1: Some covariance elements of Y2 not well defined"))
cor<-diag(2)
cor[lower.tri(cor)]<-summary(final_model_3b1)$hyperpar[7, 1]
cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
coefficients_3b1$U2$INLA<-Cov2_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)
