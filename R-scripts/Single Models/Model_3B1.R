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
gamma1<-0.45 
gamma2<-1.2
true_betas<-c(2,4,2.5,3, 1.5, 3.5)

coefficients_3b1<-list(coefficients=data.frame(true= c(true_betas, gamma1, gamma2),  INLA=rep(0,8)),
                      U1=list(true=U1, INLA=diag(2)),
                      U2=list(true=U2, INLA=diag(2)))


rownames(coefficients_3b1$coefficients)<-c( 'beta_0^1', 'beta_x^1','beta_t^1', 'beta_0^2', 'beta_x^2','beta_t^2', 'gamma1', 'gamma2')
################################################################################
#Model 3B1: Random Y1-effects scaled independently; independent Slope & Intercept
################################################################################
set.seed(2022)
### Total of N individuals
N <- 1000
### n correlated measurements per individual
n <- 3

data_3b1 <- data.frame(
  id=rep(1:N, each=n),
  x=rnorm(N*n),
  Period=rep(0:(n-1), times=N), 
  Period2=factor(rep(0:(n-1), times=N)))
errors1 <- matrix(rnorm(N*n), N)
errors2 <- matrix(rnorm(N*n), N)
u1<-matrix(rnorm(N*2), N) %*% chol(U1)
u2<-matrix(rnorm(N*2), N) %*% chol(U2)
coefficients_3b1$U1$true<-cov(u1)
coefficients_3b1$U2$true<-cov(u2)
data_3b1$y1 <-  true_betas[1]+rep(u1[,1], each=n)+true_betas[2]*data_3b1$x+(true_betas[3]+rep(u1[,2], each=n))*data_3b1$Period+as.vector(t(errors1))
data_3b1$y2<- gamma1*rep(u1[,1], each=n)+true_betas[4]+rep(u2[,1], each=n)+true_betas[5]*data_3b1$x+
  (true_betas[6]+rep(u2[,2], each=n)+gamma2*rep(u1[,2], each=n))*data_3b1$Period+as.vector(t(errors2))



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

INLA_data_3b1<-data_3b1
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
                     data = final_data_3b1,verbose=TRUE)
summary(final_model_3b1)
coefficients_3b1$coefficients$INLA<-c(final_model_3b1$summary.fixed[1:6,1],
                                     summary(final_model_3b1)$hyperpar[8:9,1])
names(final_model_3b1$internal.marginals.hyperpar)
diag(coefficients_3b1$U1$INLA)<-c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                               final_model_3b1$internal.marginals.hyperpar[[3]]))$mean,
                                 inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                               final_model_3b1$internal.marginals.hyperpar[[4]]))$mean)

diag(coefficients_3b1$U2$INLA)<-c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                               final_model_3b1$internal.marginals.hyperpar[[5]]))$mean,
                                 inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                               final_model_3b1$internal.marginals.hyperpar[[6]]))$mean)
coefficients_3b1$U2$INLA[1,2]<-coefficients_3b1$U2$INLA[2,1]<- 
  summary(final_model_3b1)$hyperpar[7, 1]*sqrt(coefficients_3b1$U2$INLA[1,1]*coefficients_3b1$U2$INLA[2,2])

coefficients_3b1
