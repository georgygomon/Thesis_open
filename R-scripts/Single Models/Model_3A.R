library('INLA')
library(reshape2)
library(lme4)
library(nlme)
library(MCMCglmm)
library(MCMCvis)


################################################################################
#Parameters
################################################################################
U1<-matrix(c(2,0,
             0,3),2, byrow=TRUE)
U2<-matrix(c(3,1.5,
             1.5,4),2, byrow=TRUE)
gamma<-0.45 
true_betas<-c(2,4,2.5,3, 1.5, 3.5)

coefficients_3a<-list(coefficients=data.frame(true= c(true_betas, gamma),  INLA=rep(0,7), 
                                              true_combi=c(true_betas[1:3], true_betas[1:3]*gamma+true_betas[4:6], NA), INLA_combi=c(rep(0,6), NA)),
                      U1=list(true=U1, INLA=diag(2)),
                      U2=list(true=U2, INLA=diag(2)))
  
  
rownames(coefficients_3a$coefficients)<-c( 'beta_0^1', 'beta_x^1','beta_t^1', 'beta_0^2', 'beta_x^2','beta_t^2', 'gamma')
################################################################################
#Model 3A: Y1 linear predictor scaled with independent Random Slope & Intercept
################################################################################
set.seed(2022)
### Total of N individuals
N <- 100
### n correlated measurements per individual
n <- 3

data_3a <- data.frame(
  id=rep(1:N, each=n),
  x=rnorm(N*n),
  Period=rep(0:(n-1), times=N), 
  Period2=factor(rep(0:(n-1), times=N)))
errors1 <- matrix(rnorm(N*n), N)
errors2 <- matrix(rnorm(N*n), N)
u1<-matrix(rnorm(N*2), N) %*% chol(U1)
u2<-matrix(rnorm(N*2), N) %*% chol(U2)
coefficients_3a$U1$true<-cov(u1)
coefficients_3a$U2$true<-cov(u2)
m_3a<-true_betas[1]+rep(u1[,1], each=n)+true_betas[2]*data_3a$x+(true_betas[3]+rep(u1[,2], each=n))*data_3a$Period
data_3a$y1 <-  m_3a+as.vector(t(errors1))
data_3a$y2<- gamma*(m_3a)+true_betas[4]+rep(u2[,1], each=n)+true_betas[5]*data_3a$x+(true_betas[6]+rep(u2[,2], each=n))*data_3a$Period+as.vector(t(errors2))



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

INLA_data_3a<-data_3a
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
                     data = final_data_3a,verbose=TRUE)

coefficients_3a$coefficients$INLA<-c(final_model_3a$summary.random$Intercept1$mode, 
                        final_model_3a$summary.random$idx1$mean,
                        final_model_3a$summary.random$Period1$mean,
                        final_model_3a$summary.fixed[c(1,2,3),1],
                        summary(final_model_3a)$hyperpar[11,1])
diag(coefficients_3a$U1$INLA)<-c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                               final_model_3a$internal.marginals.hyperpar[[6]]))$mean,
                                 inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                               final_model_3a$internal.marginals.hyperpar[[7]]))$mean)

diag(coefficients_3a$U2$INLA)<-c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                               final_model_3a$internal.marginals.hyperpar[[8]]))$mean,
                                 inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                               final_model_3a$internal.marginals.hyperpar[[9]]))$mean)
coefficients_3a$U2$INLA[1,2]<-coefficients_3a$U2$INLA[2,1]<- 
  summary(final_model_3a)$hyperpar[10, 1]*sqrt(coefficients_3a$U2$INLA[1,1]*coefficients_3a$U2$INLA[2,2])
coefficients_3a$coefficients$INLA_combi<-c(coefficients_3a$coefficients$INLA[1:3], 
                                           coefficients_3a$coefficients$INLA[1:3]*coefficients_3a$coefficients$INLA[7]+
                                             coefficients_3a$coefficients$INLA[4:6], NA)

coefficients_3a
