library('INLA')
library(reshape2)
library(lme4)
library(nlme)
library(MCMCglmm)
library(MCMCvis)
library(optimx)
library(inlabru)


################################################################################
#Parameters
################################################################################
parameters<-list(U1=matrix(c(2,1,
                             1,3),2, byrow=TRUE),
                 U2=matrix(c(3,1.5,
                             1.5,4),2, byrow=TRUE))
true_betas<-c(2,4,2.5,3, 1.5, 3.5)

coefficients_0<-list(coefficients=data.frame(true_betas, INLA=rep(NA, 6)),
                      U1=list(true=parameters$U1, INLA=diag(2)),
                      U2=list(true=parameters$U2, INLA=diag(2)))


rownames(coefficients_0$coefficients)<-c( 'beta_0^1', 'beta_x^1','beta_t^1', 'beta_0^2', 'beta_x^2','beta_t^2')
################################################################################
#Model 0
################################################################################
set.seed(2021)
### Total of N individuals
N <- 120
### Max # of measurements per individual
n <- 7
### Probability p that a measurement is observed
p1<-0.7
### probability p2 that covariate x is measured
p2<-0.9

set.seed(2022)
data_0 <- data.frame(
  id=rep(1:N, each=n),
  x=rnorm(N*n),
  Period=rep(0:(n-1), times=N), 
  Period2=factor(rep(0:(n-1), times=N)),
  observed=rbern(N*n, p1),
  x_measured=rbern(N*n, p2))
errors1 <- matrix(rnorm(N*n, sd=0.1), N)
errors2 <- matrix(rnorm(N*n, sd=0.1), N)
u1<-matrix(rnorm(N*2), N) %*% chol(parameters$U1)
coefficients_0$U1$true<-cov(u1)
u2<-matrix(rnorm(N*2), N) %*% chol(parameters$U2)
coefficients_0$U2$true<-cov(u2)
data_0$y1 <-  true_betas[1]+rep(u1[,1], each=n)+true_betas[2]*data_0$x+(true_betas[3]+rep(u1[,2], each=n))*data_0$Period+as.vector(t(errors1))
data_0$y2 <- true_betas[4]+rep(u2[,1], each=n)+true_betas[5]*data_0$x+(true_betas[6]+rep(u2[,2], each=n))*data_0$Period+as.vector(t(errors2))
data_0$x[data_0$x_measured==0]<-NA
data<-data_0[data_0$observed==1,]


################################################################################
#Model 0 using INLA 
################################################################################
INLA_data_0<-data
N_obs<-length(data$x)
#Modelling as 2 separate likelihoods
fixed.effects_0<-list(Intercept1=c(rep(1, N_obs), rep(NA, N_obs)),
                      Period1=c(INLA_data_0$Period, rep(NA, N_obs)),
                      x1=c(INLA_data_0$x, rep(NA, N_obs)),
                      Intercept2=c(rep(NA, N_obs), rep(1, N_obs)),
                      Period2=c(rep(NA, N_obs),INLA_data_0$Period),
                      x2=c(rep(NA, N_obs),INLA_data_0$x),
                      Period=c(INLA_data_0$Period, INLA_data_0$Period))
random.effects_0<-list(Intercept_random1=c(INLA_data_0$id, rep(NA, N_obs)),
                       Period_random1=c(INLA_data_0$id+N, rep(NA, N_obs)),
                       Intercept_random2=c(rep(NA, N_obs), INLA_data_0$id),
                       Period_random2=c(rep(NA, N_obs), INLA_data_0$id+N))
y_final_0<-list(c(INLA_data_0$y1, rep(NA, N_obs)),
                c(rep(NA, N_obs),INLA_data_0$y2))
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

coefficients_0$coefficients$INLA<-final_model_0$summary.fixed[,1]

#Random effects covariance Y1
SD_s<-diag(2)
tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_0$internal.marginals.hyperpar[[3]]), silent=TRUE)$mean,
                      inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_0$internal.marginals.hyperpar[[4]]), silent=TRUE)$mean)), 
         error = function(e) print("Model 0: Some covariance elements of Y1 not well defined"))
cor<-diag(2)
cor[lower.tri(cor)]<-summary(final_model_0)$hyperpar[5, 1]
cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
coefficients_0$U1$INLA<-Cov1_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)

#Random effects covariance Y2
SD_s<-diag(2)
tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_0$internal.marginals.hyperpar[[6]]), silent=TRUE)$mean,
                      inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_0$internal.marginals.hyperpar[[7]]), silent=TRUE)$mean)), 
         error = function(e) print("Model 0: Some covariance elements of Y2 not well defined"))
cor<-diag(2)
cor[lower.tri(cor)]<-summary(final_model_0)$hyperpar[8, 1]
cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
coefficients_0$U2$INLA<-Cov2_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)
-2*final_model_0$mlik[1]
final_model_0$dic$dic
final_model_0$marginals.fixed
predict(final_model_0, data)

final_model_0$marginals.fitted.values$fitted.Predictor.0001
