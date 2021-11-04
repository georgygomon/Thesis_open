library('INLA')
library(reshape2)
library(lme4)
library(nlme)
library(MCMCglmm)
library(MCMCvis)


################################################################################
#Parameters
################################################################################
true_betas<-c(3,2,1.5, 4,3, 2.5)
U1<-matrix(c(2,-1,
             -1,3),2, byrow=TRUE)
U2<-matrix(c(3,1.5,
             1.5,4),2, byrow=TRUE)
U3<-matrix(c(0,0,
             0,0),2, byrow=TRUE)
U<-rbind(cbind(U1, U3), cbind(t(U3), U2))



coefficients_2c<-list(coefficients=data.frame(true= c(true_betas), 
                                              lmer=rep(0,6),nlme=rep(0,6),
                                              INLA=rep(0,6),  MCMCglmm=rep(0,6)),
                      VarCov=list(true=U, nlme=diag(4), lmer=diag(4), MCMCglmm=diag(4), INLA=diag(4)))
rownames(coefficients_2c$coefficients)<-c( 'beta_0^1', 'beta_x^1','beta_t^1', 'beta_0^2', 'beta_x^2','beta_t^2')

################################################################################
#Make data 2c1: Random slope & intercept, uncorrelated
################################################################################
set.seed(2022)
### Total of N individuals
N <- 750
### n correlated measurements per individual
n <- 3

data_2c <- data.frame(
  id=rep(1:N, each=n),
  x=rnorm(N*n),  
  Period=rep(0:(n-1), times=N), 
  Period2=factor(rep(0:(n-1), times=N)))

errors1 <- matrix(rnorm(N*n), N)
errors2 <- matrix(rnorm(N*n), N)
u<-matrix(rnorm(N*2*2), N) %*% chol(U)
coefficients_2c$VarCov$true<-cov(u)
data_2c$y1 <- true_betas[1] +rep(u[,1], each=n)+ true_betas[2]*data_2c$x + true_betas[3]*data_2c$Period+
  rep(u[,3], each=n)*data_2c$Period+as.vector(t(errors1))
data_2c$y2<- true_betas[4]+rep(u[,2], each=n)+true_betas[5]*data_2c$x+true_betas[6]*data_2c$Period+
  rep(u[,4], each=n)*data_2c$Period+as.vector(t(errors2))

################################################################################
#Model 2a using nlme
################################################################################
data_melted_2c<-melt(data_2c, measure.vars=c('y1', 'y2'))
data_melted_2c<-data_melted_2c[with(data_melted_2c,order(id)),]

#https://stats.stackexchange.com/questions/58669/specifying-multiple-separate-random-effects-in-lme
#list(year=~1, date=~time)
#Correlated intercept and slope
#random =~0 + variable + variable:Period2 | id
#uncorrelated intercept and slope
model_lme_2c<-lme(value~variable+variable:x+variable:Period-1,
                  random =list(id=~variable+0, id=~0+variable:Period), 
                  data=data_melted_2c, na.action=na.omit)
coefficients_2c$coefficients$nlme<-model_lme_2c$coefficients$fixed[c(1,3,5,2,4,6)]
coefficients_2c$VarCov$nlme[lower.tri(coefficients_2c$VarCov$true, diag=TRUE)]<-c(as.numeric(VarCorr(model_lme_2c)[2,1]),
                                                                                  as.numeric(VarCorr(model_lme_2c)[3,3])*as.numeric(VarCorr(model_lme_2c)[2,2])*as.numeric(VarCorr(model_lme_2c)[3,2]),
                                                                                  0,0,
                                                                                  as.numeric(VarCorr(model_lme_2c)[3,1]),
                                                                                  0,0,
                                                                                  as.numeric(VarCorr(model_lme_2c)[5,1]),
                                                                                  as.numeric(VarCorr(model_lme_2c)[6,3])*as.numeric(VarCorr(model_lme_2c)[5,2])*as.numeric(VarCorr(model_lme_2c)[6,2]),
                                                                                  as.numeric(VarCorr(model_lme_2c)[6,1]))
coefficients_2c$VarCov$nlme[upper.tri(coefficients_2c$VarCov$nlme)] <- t(coefficients_2c$VarCov$nlme)[upper.tri(coefficients_2c$VarCov$nlme)]


################################################################################
#Model 2a using lmer
################################################################################
model_lmer_2c<-lmer(value~ variable+variable:x+variable:Period-1+(variable+0+variable:Period||id), na.action=na.omit, 
                    data=data_melted_2c)
#LMER control options
#library(optimx)
#control = lmerControl(optimizer = "optimx",optCtrl = list(method = "L-BFGS-B", maxit = 1e6))
sum_lmer_2c<-summary(model_lmer_2c)

coefficients_2c$coefficients$lmer<-sum_lmer_2c$coefficients[c(1,3,5,2,4,6),1]
coefficients_2c$VarCov$lmer[c(1,2),c(1,2)]<-sum_lmer_2c$varcor$id[c(1,2),c(1,2)]
coefficients_2c$VarCov$lmer[c(3,4),c(3,4)]<-sum_lmer_2c$varcor$id.1[c(1,2),c(1,2)]

################################################################################
#Model 2a using MCMCglmm
################################################################################

#using an inverse Wishart distribution with parameter R11, R12, R21, R22, v
#Check out the Inverse Wishart priors!!!! Where does the variance come from?
prior = list(G = list(G1=list(V = coefficients_2c$VarCov$true[c(1,2), c(1,2)], n = 4),
                      G2=list(V = coefficients_2c$VarCov$true[c(3,4), c(3,4)], n = 4)),
             R = list(V = diag(2*n), n = 4))

MCMC_2C1<-MCMCglmm(cbind(y1, y2) ~ trait+trait:x +trait:Period- 1,
                   random = ~ us(trait):units+us(trait:Period):units,  rcov = ~ idh(trait:Period2):units, 
                   family = rep("gaussian", 2), prior=prior,nitt = 10000, burnin = 1000,
                   thin=25, data = data_2c)
coefficients_2c$VarCov$MCMCglmm[c(1,2),c(1,2)]<-matrix(MCMCsummary(MCMC_2C1$VCV)[c(1:4),1], nrow=2)
coefficients_2c$VarCov$MCMCglmm[c(3,4),c(3,4)]<-matrix(MCMCsummary(MCMC_2C1$VCV)[c(5:8),1], nrow=2)
coefficients_2c$coefficients$MCMCglmm<-c(MCMCsummary(MCMC_2C1$Sol)[c(1,3,5,2,4,6),1])

################################################################################
#Model 2a using INLA 
################################################################################
INLA_data_2c1<-data_2c
#Modelling as 2 separate likelihoods
Period_id<-rep(NA, times=n*N)
Period_id[INLA_data_2c1$Period2!=0]<-rep(1:N, each=(n-1))
fixed.effects_2c1<-list(Intercept1=c(rep(1, n*N), rep(NA, n*N)),
                        Intercept2=c(rep(NA, n*N), rep(1, n*N)),
                        x1=c(INLA_data_2c1$x, rep(NA, n*N)),
                        x2=c(rep(NA, n*N), INLA_data_2c1$x), 
                        Period1=c(INLA_data_2c1$Period, rep(NA, n*N)),
                        Period2=c(rep(NA, n*N), INLA_data_2c1$Period),
                        Period_all=c(INLA_data_2c1$Period, INLA_data_2c1$Period))
random.effects_2c1<-list(ID_all=c(INLA_data_2c1$id, INLA_data_2c1$id+N),
                         Period=c(Period_id, Period_id+N))

y_final_2c1<-list(c(INLA_data_2c1$y1, rep(NA, n*N)),
                  c(rep(NA, n*N),INLA_data_2c1$y2))
final_data_2c1<-c(fixed.effects_2c1, random.effects_2c1)
final_data_2c1$Y<-y_final_2c1

formula.model_2c1=Y~-1+f(Intercept1, model='linear')+  f(Intercept2, model='linear')+
  f(x1, model='linear')+  f(x2, model='linear')+f(Period1, model='linear')+f(Period2, model='linear')+
  f(ID_all, model="iid2d", n=2*N)+f(Period, Period_all, model="iid2d", n=2*N)
#  f(ID_all_2,Period2, model="iid2d", n=2*N,hyper = list(theta1 = list(prior = "wishart2d", param =c(4,1,1, 0))))
final_model_2c1<-inla(formula.model_2c1, family = c("gaussian","gaussian"),
                      data = final_data_2c1,verbose=TRUE)
summary(final_model_2c1)
coefficients_2c$coefficients$INLA<-final_model_2c1$summary.fixed[c(1,3,5,2,4,6),1]
diag(coefficients_2c$VarCov$INLA[c(1,2),c(1,2)])<-c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                                                  final_model_2c1$internal.marginals.hyperpar[[3]]))$mean, 
                                                    inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                                                  final_model_2c1$internal.marginals.hyperpar[[4]]))$mean)
coefficients_2c$VarCov$INLA[1,2]<-coefficients_2c$VarCov$INLA[2,1]<-
  summary(final_model_2c1)$hyperpar[5, 1]*sqrt(coefficients_2c$VarCov$INLA[1,1]*coefficients_2c$VarCov$INLA[2,2])

names(final_model_2c1$internal.marginals.hyperpar)
diag(coefficients_2c$VarCov$INLA[c(3,4),c(3,4)])<-c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                                                  final_model_2c1$internal.marginals.hyperpar[[6]]))$mean, 
                                                    inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                                                  final_model_2c1$internal.marginals.hyperpar[[7]]))$mean)
coefficients_2c$VarCov$INLA[3,4]<-coefficients_2c$VarCov$INLA[4,3]<-
  summary(final_model_2c1)$hyperpar[8, 1]*sqrt(coefficients_2c$VarCov$INLA[3,3]*coefficients_2c$VarCov$INLA[4,4])

coefficients_2c

