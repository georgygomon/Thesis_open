library('INLA')
library(reshape2)
library(lme4)
library(nlme)
library(MCMCglmm)
library(MCMCvis)
library(optimx)
library(Rlab)
library(dgof)


################################################################################
#Parameters
################################################################################
U1<-matrix(c(2,0,
             0,3),2, byrow=TRUE)
U2<-matrix(c(3,1.5,
             1.5,4),2, byrow=TRUE)
gamma<-1.2
true_betas<-c(2,4,2.5,3, 1.5, 3.5)

coefficients_3a<-list(coefficients=data.frame(true= c(true_betas, gamma),  INLA=rep(0,7), 
                                              true_combi=c(true_betas[1:3], true_betas[1:3]*gamma+true_betas[4:6], NA), INLA_combi=c(rep(0,6), NA)),
                      U1=list(true=U1, INLA=diag(2)),
                      U2=list(true=U2, INLA=diag(2)))
  
  
rownames(coefficients_3a$coefficients)<-c( 'beta_0^1', 'beta_x^1','beta_t^1', 'beta_0^2', 'beta_x^2','beta_t^2', 'gamma')
################################################################################
#Model 2 Make data: Only random intercept, uncorrelated errors
################################################################################
set.seed(2022)
### Total of N individuals
N <- 500
### Max # of measurements per individual
n <- 7
### Probability p that a measurement is observed
p1<-1
### probability p2 that covariate x is measured
p2<-1

data_3a <- data.frame(
  id=rep(1:N, each=n),
  x=rnorm(N*n),
  Period=rep(0:(n-1), times=N), 
  Period2=factor(rep(0:(n-1), times=N)),
  observed=rbern(N*n, p1),
  x_measured=rbern(N*n, p2))
errors1 <- matrix(rnorm(N*n), N)
errors2 <- matrix(rnorm(N*n), N)
u1<-matrix(rnorm(N*2), N) %*% chol(U1)
u2<-matrix(rnorm(N*2), N) %*% chol(U2)
coefficients_3a$U1$true<-cov(u1)
coefficients_3a$U2$true<-cov(u2)
m_3a<-true_betas[1]+rep(u1[,1], each=n)+true_betas[2]*data_3a$x+(true_betas[3]+rep(u1[,2], each=n))*data_3a$Period
data_3a$y1 <-  m_3a+as.vector(t(errors1))
data_3a$y2<- gamma*(m_3a)+true_betas[4]+rep(u2[,1], each=n)+true_betas[5]*data_3a$x+(true_betas[6]+rep(u2[,2], each=n))*data_3a$Period+as.vector(t(errors2))
data_3a$x[data_3a$x_measured==0]<-NA
data<-data_3a[data_3a$observed==1,]

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

INLA_data_3a<-data
#Modelling as 2 separate likelihoods
N_obs<-length(INLA_data_3a$x)
fixed.effects_3a<-list(x=c(INLA_data_3a$x, INLA_data_3a$x), 
                       Period=c(INLA_data_3a$Period, INLA_data_3a$Period),
                       Intercept2=c(rep(NA, N_obs), rep(1, N_obs)),
                       Period2=c(rep(NA, N_obs),INLA_data_3a$Period),
                       x2=c(rep(NA, N_obs),INLA_data_3a$x))
random.effects_3a<-list(Intercept1=c(rep(1, N_obs), rep(NA, N_obs)),
                        Intercept12=c(rep(NA, N_obs), rep(1, N_obs)),
                        idx1=c(rep(1, N_obs),rep(NA, N_obs)),
                        idx12=c(rep(NA, N_obs),rep(1, N_obs)),
                        Period1=c(rep(1, N_obs),rep(NA, N_obs)),
                        Period12=c(rep(NA, N_obs),rep(1, N_obs)),
                        Intercept_random1=c(INLA_data_3a$id, rep(NA, N_obs)),
                        Intercept_random12=c(rep(NA, N_obs), INLA_data_3a$id),
                        Intercept_random2=c(rep(NA, N_obs), INLA_data_3a$id),
                        Period_random1=c(INLA_data_3a$id, rep(NA, N_obs)),
                        Period_random12=c(rep(NA, N_obs), INLA_data_3a$id),
                        Period_random2=c(rep(NA, N_obs), INLA_data_3a$id+N))
y_final_3a<-list(c(INLA_data_3a$y1, rep(NA, N_obs)),
                 c(rep(NA, N_obs),INLA_data_3a$y2))
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
summary(final_model_3a)
model_coefficients<-c(final_model_3a$summary.random$Intercept1$mean, 
                      final_model_3a$summary.random$idx1$mean,
                      final_model_3a$summary.random$Period1$mean,
                      final_model_3a$summary.fixed[c(1,2,3),1])
gamma<-summary(final_model_3a)$hyperpar[11,1]
actual_coefficients<-c(model_coefficients[1:3], model_coefficients[1:3]*gamma+model_coefficients[4:6])


summary(final_model_3a)
SD_s<-diag(2)
tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_3a$internal.marginals.hyperpar[[6]]), silent=TRUE)$mean,
                      inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_3a$internal.marginals.hyperpar[[7]]), silent=TRUE)$mean)), 
         error = function(e) print("Model 3A1: Some covariance elements of Y1 not well defined"))
cor<-diag(2)
Cov1_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)

SD_s<-diag(2)
tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_3a$internal.marginals.hyperpar[[8]]), silent=TRUE)$mean,
                      inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_3a$internal.marginals.hyperpar[[9]]), silent=TRUE)$mean)), 
         error = function(e) print("Model 3A1: Some covariance elements of Y2 not well defined"))
cor<-diag(2)
cor[lower.tri(cor)]<-summary(final_model_3a)$hyperpar[10, 1]
cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
Cov2_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)


coefficients_3a$coefficients$INLA<-c(model_coefficients, gamma)
coefficients_3a$coefficients$INLA_combi<-c(actual_coefficients, NA)
coefficients_3a$U1$INLA<-Cov1_0_t
coefficients_3a$U2$INLA<-Cov2_0_t

coefficients_3a
