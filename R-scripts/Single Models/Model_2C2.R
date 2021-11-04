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
true_betas<-c(3,2,1.5, 4,3, 2.5)
U1<-matrix(c(2,-1,
             -1,3),2, byrow=TRUE)
U2<-matrix(c(3,1.5,
             1.5,4),2, byrow=TRUE)
U3<-matrix(c(1.25,1,-0.5,0.3),2, byrow=TRUE)
U<-rbind(cbind(U1, U3), cbind(t(U3), U2))



coefficients_2c2<-list(coefficients=data.frame(true= c(true_betas), 
                                              lmer=rep(0,6),nlme=rep(0,6),
                                              INLA=rep(0,6),  MCMCglmm=rep(0,6)),
                      VarCov=list(true=U, nlme=diag(4), lmer=diag(4), MCMCglmm=diag(4), INLA=diag(4)))
rownames(coefficients_2c2$coefficients)<-c( 'beta_0^1', 'beta_x^1','beta_t^1', 'beta_0^2', 'beta_x^2','beta_t^2')

################################################################################
#Model 2 Make data: Dependent Slope and Intercept!!
################################################################################
set.seed(2022)
### Total of N individuals
N <- 750
### n correlated measurements per individual
n <- 3

data_2c2 <- data.frame(
  id=rep(1:N, each=n),
  x=rnorm(N*n), 
  Period=rep(0:(n-1), times=N), 
  Period2=factor(rep(0:(n-1), times=N)))
errors1 <- matrix(rnorm(N*n), N)
errors2 <- matrix(rnorm(N*n), N)
u<-matrix(rnorm(N*2*2), N) %*% chol(U)
coefficients_2c2$VarCov$true<-cov(u)
data_2c2$y1 <- true_betas[1] +rep(u[,1], each=n)+ true_betas[2]*data_2c2$x + true_betas[3]*data_2c2$Period+
  rep(u[,3], each=n)*data_2c2$Period+as.vector(t(errors1))
data_2c2$y2<- true_betas[4]+rep(u[,2], each=n)+true_betas[5]*data_2c2$x+true_betas[6]*data_2c2$Period+
  rep(u[,4], each=n)*data_2c2$Period+as.vector(t(errors2))

################################################################################
#Model 2a using nlme
################################################################################
data_melted_2c2<-melt(data_2c2, measure.vars=c('y1', 'y2'))
data_melted_2c2<-data_melted_2c2[with(data_melted_2c2,order(id)),]

#https://stats.stackexchange.com/questions/58669/specifying-multiple-separate-random-effects-in-lme
#list(year=~1, date=~time)
#Correlated intercept and slope
#random =~0 + variable + variable:Period2 | id
#uncorrelated intercept and slope
model_lme_2c2<-lme(value~variable+variable:x+variable:Period-1,
                  random =list(id=~variable+variable:Period+0), 
                  data=data_melted_2c2, na.action=na.omit)

# model_lme_2c<-lme(value~variable+variable:x+variable:Period-1,
#                   random =list(id=~variable+0, id=~0+variable:Period2), 
#                   data=data_melted_2c, na.action=na.omit)

 


coefficients_2c2$coefficients$nlme<-model_lme_2c2$coefficients$fixed[c(1,3,5,2,4,6)]
coefficients_2c2$VarCov$nlme<-matrix(getVarCov(model_lme_2c2), nrow=4)
################################################################################
#Model 2a using lmer
################################################################################
model_lmer_2c2<-lmer(value~ variable+variable:x+variable:Period-1+(variable+0+variable:Period|id), na.action=na.omit, 
                    data=data_melted_2c2, 
                    control = lmerControl(
                      optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
#LMER control options
#library(optimx)
#control = lmerControl(optimizer = "optimx",optCtrl = list(method = "L-BFGS-B", maxit = 1e6))
sum_lmer_2c2<-summary(model_lmer_2c2)

coefficients_2c2$coefficients$lmer<-sum_lmer_2c2$coefficients[c(1,3,5,2,4,6),1]
coefficients_2c2$VarCov$lmer<-matrix(sum_lmer_2c2$varcor$id[c(1:4),c(1:4)], nrow=4)

################################################################################
#Model 2a using MCMCglmm
################################################################################

#using an inverse Wishart distribution with parameter R11, R12, R21, R22, v
#Check out the Inverse Wishart priors!!!! Where does the variance come from?
prior = list(G = list(G1=list(V = coefficients_2c2$VarCov$true, n = 4)),
             R = list(V = diag(2*n), n = 4))

MCMC_2c2<-MCMCglmm(cbind(y1, y2) ~ trait+trait:x +trait:Period- 1,
                   random = ~ us(trait:Period+trait):units,  rcov = ~ idh(trait:Period2):units, 
                   family = rep("gaussian", 2), prior=prior, nitt = 20000, burnin = 2000,
                   thin=25, data = data_2c2)

coefficients_2c2$VarCov$MCMCglmm<-matrix(MCMCsummary(MCMC_2c2$VCV)[c(1:16),1], nrow=4)
coefficients_2c2$coefficients$MCMCglmm<-c(MCMCsummary(MCMC_2c2$Sol)[c(1,3,5,2,4,6),1])

################################################################################
#Model 2a using INLA 
################################################################################
INLA_data_2c2<-data_2c2
#Modelling as 2 separate likelihoods
#Use copy feature to create dependence between Intercept and Slope
#https://arxiv.org/pdf/1210.0333.pdf
Period_id<-rep(NA, times=n*N)
Period_id[INLA_data_2c2$Period2!=0]<-rep(1:N, each=(n-1))
fixed.effects_2c2<-list(Intercept1=c(rep(1, n*N), rep(NA, n*N)),
                        Intercept2=c(rep(NA, n*N), rep(1, n*N)),
                        x1=c(INLA_data_2c2$x, rep(NA, n*N)),
                        x2=c(rep(NA, n*N), INLA_data_2c2$x), 
                        Period1=c(INLA_data_2c2$Period, rep(NA, n*N)),
                        Period2=c(rep(NA, n*N), INLA_data_2c2$Period),
                        Period_all=c(INLA_data_2c2$Period, INLA_data_2c2$Period))
random.effects_2c2<-list(ID_all=c(INLA_data_2c2$id, INLA_data_2c2$id+N),
                         Period_random=c(Period_id+2*N, Period_id+3*N))

y_final_2c2<-list(c(INLA_data_2c2$y1, rep(NA, n*N)),
                  c(rep(NA, n*N),INLA_data_2c2$y2))
final_data_2c2<-c(fixed.effects_2c2, random.effects_2c2)
final_data_2c2$Y<-y_final_2c2

formula.model_2c2=Y~-1+f(Intercept1, model='linear')+  f(Intercept2, model='linear')+
  f(x1, model='linear')+  f(x2, model='linear')+f(Period1, model='linear')+f(Period2, model='linear')+
  f(ID_all, model="iid4d", n=4*N)+  f(Period_random, Period_all,copy="ID_all")
#  f(ID_all_2,Period2, model="iid2d", n=2*N,hyper = list(theta1 = list(prior = "wishart2d", param =c(4,1,1, 0))))
final_model_2c2<-inla(formula.model_2c2, family = c("gaussian","gaussian"),
                      data = final_data_2c2,verbose=TRUE, control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))


 


coefficients_2c2$coefficients$INLA<-final_model_2c2$summary.fixed[c(1,3,5,2,4,6),1]
SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                           final_model_2c2$internal.marginals.hyperpar[[3]]))$mean,
             inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                           final_model_2c2$internal.marginals.hyperpar[[4]]))$mean,
             inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                           final_model_2c2$internal.marginals.hyperpar[[5]]))$mean,
             inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                           final_model_2c2$internal.marginals.hyperpar[[6]]))$mean))
cor<-matrix(1, nrow=4, ncol=4)
cor[lower.tri(cor)]<-final_model_2c2$summary.hyperpar[7:12,1]
cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
coefficients_2c2$VarCov$INLA<-diag(SD_s) %*% cor %*% diag(SD_s)



coefficients_2c2

# final_model_2c2$mlik
# final_model_2c2$dic$dic
# final_model_2c2$waic$waic
# final_model_2c2$cpo

