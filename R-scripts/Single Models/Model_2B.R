library('INLA')
library(reshape2)
library(lme4)
library(nlme)
library(MCMCglmm)
library(MCMCvis)


################################################################################
#Parameter List
################################################################################
true_betas<-c(2,4,2.5,3, 1.5, 3.5)
U<-matrix(c(2,-1,
            -1,3),2, byrow=TRUE)
Sigma1 <- matrix(c(5,2,2,4), 2)
Sigma2<-matrix(c(4,1,1,4),2)
Sigma3<-matrix(c(-1, -0.5, -0.5, -1.5), 2)
Sigma<-rbind(cbind(Sigma1, Sigma3), cbind(Sigma3, Sigma2))


coefficients_2b<-list(coefficients=data.frame(true= c(true_betas, diag(U), U[1,2]), 
                            lmer=rep(0,9),nlme=rep(0,9),
                            INLA=rep(0,9)),
                      VarCov=list(true=Sigma, nlme=NA, INLA=NA))
rownames(coefficients_2b$coefficients)<-c('beta_0^1', 'beta_x^1', 'beta_t^1', 'beta_0^2', 'beta_x^2', 'beta_t^2', 'u_1', 'u_2', 'u_12')




################################################################################
#Model 2b Make data: Both random intercept and correlated terms within 1 variable
################################################################################
set.seed(2021)
### Total of N individuals
N <- 700
### n correlated measurements per individual
n <- 2

data_2b<- data.frame(
  id=rep(1:N, each=n),
  x=rnorm(N*n), 
  Period=rep(0:(n-1), times=N), 
  Period2=factor(rep(0:(n-1), times=N)))

### sampling a covariate setting up data such that the n measurements on 1 individual are correlated
errors <- matrix(rnorm(N*2*n), N) %*% chol(Sigma)
u<-matrix(rnorm(N*2), N) %*% chol(U)
coefficients_2b$VarCov$true<-cov(errors)
coefficients_2b$coefficients$true[7:9]<-c(diag(cov(u)), cov(u)[1,2])

data_2b$y1 <- true_betas[1] + true_betas[2]*data_2b$x +true_betas[3]*data_2b$Period+rep(u[,1], each=n)+as.vector(t(errors[,1:n]))
data_2b$y2<- true_betas[4]+true_betas[5]*data_2b$x+true_betas[6]*data_2b$Period+rep(u[,2], each=n)+as.vector(t(errors[,(n+1):(2*n)]))



#This is very important!!!!!
data_2b<-rbind(data_2b[data_2b$Period==0,], data_2b[data_2b$Period==1,])



################################################################################
#Model 2b using nlme
################################################################################
data_melted_2b<-melt(data_2b, measure.vars=c('y1', 'y2'))
data_melted_2b<-data_melted_2b[with(data_melted_2b,order(id)),]

model_lme_2b<-lme(value~variable+variable:x+variable:Period-1,
                  random = ~0+variable | id, 
                  correlation=corSymm(form=~1|id), 
                  weights= varIdent(form=~1|variable*Period),
                  data=data_melted_2b, na.action=na.omit)

coefficients_2b$VarCov$nlme<-getVarCov(model_lme_2b, individual=1, type='conditional')$`1`
coefficients_2b$coefficients$nlme<-c(model_lme_2b$coefficients$fixed[c(1,3,5,2,4,6)],
                        diag(getVarCov(model_lme_2b, individual=1, type='random.effects')),
                        getVarCov(model_lme_2b, individual=1, type='random.effects')[1,2] )

################################################################################
#Model 2b using MCMCglmm
################################################################################
# #MCMCglmm can do either error covariance structure or random effects, however not both!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! :(
# prior = list(B = list(mu = rep(0,6),V = 10*diag(6)),
#              R = list(V = diag(4), n = 4))
# #??covu_missing not found error
# #G = list(G1=list(V = matrix(c(2,-1,-1, 3), nrow=2), n = 4)),
# 
# MCMC_2b<-MCMCglmm(cbind(y1, y2) ~ trait+trait:x+trait:Period - 1,
#                   rcov = ~ us(Period2:trait):units,
#                   family = rep("gaussian", 2),prior=prior, nitt = 20000, burnin = 2000,
#                   thin=25, data = data_2b)
# #random = ~ us(trait):units,
# coefficients_2b$MCMCglmm<-c(MCMCsummary(MCMC_2a$Sol)[c(1,3,2,4),1],MCMCsummary(MCMC_2a$VCV)[c(1,4,3),1])
# summary(MCMC_2b)
# # 
# # 
# # 
# # 
# coefficients_2b$coefficients$MCMCglmm<-c(MCMCsummary(MCMC_2b$Sol)[c(1,3,2,4),1], MCMCsummary(MCMC_2b$VCV)[c(1,4,3), 1])
# 
# MCMCsummary(MCMC_2b$VCV)
# # coefficients_1$VarCov$MCMCglmm<-matrix(MCMCsummary(MCMC_1$VCV)[,1], nrow=4, byrow=FALSE)

################################################################################
#Model 2b using INLA 
################################################################################
#Set-up data: Modelling as 2 separate likelihoods
INLA_data_2b<-data_2b
fixed.effects_2b<-list(Intercept1=c(rep(1, n*N), rep(NA, n*N)),
                       Intercept2=c(rep(NA, n*N), rep(1, n*N)),
                       x1=c(INLA_data_2b$x, rep(NA, n*N)),
                       x2=c(rep(NA, n*N), INLA_data_2b$x),
                       Period1=c(data_2b$Period, rep(NA, n*N)),
                       Period2=c(rep(NA, n*N), data_2b$Period))
random.effects_2b<-list(ID_all=c(INLA_data_2b$id, INLA_data_2b$id+N))
y_final_2b<-list(c(INLA_data_2b$y1, rep(NA, n*N)),
                 c(rep(NA, n*N),INLA_data_2b$y2))
final_data_2b<-c(fixed.effects_2b, random.effects_2b)
final_data_2b$Y<-y_final_2b
final_data_2b$i<-c(1:(2*n*N))

#Set-up Model
control_family = list(list(initial=15, fixed=T), list(initial=15, fixed=T))
formula.model_2b=Y~-1+f(Intercept1, model='linear')+  f(Intercept2, model='linear')+
  f(x1, model='linear')+  f(x2, model='linear')+
  f(Period1, model='linear')+f(Period2, model='linear')+
  f(ID_all, model="iid2d", n=2*N,hyper = list(theta1 = list(prior = "wishart2d")))+
  f(i, model="iid4d", n=2*n*N, hyper = list(theta1 = list(prior = "wishart4d")))

#Run Model
final_model_2b<-inla(formula.model_2b, family = c("gaussian","gaussian"),
                     data = final_data_2b,verbose=TRUE, control.family=control_family)

summary(final_model_2b)
#Set random effect coefficients
coefficients_2b$coefficients$INLA[1:8]<-c(final_model_2b$summary.fixed[c(1,3,5,2,4,6),1],
                             inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                           final_model_2b$internal.marginals.hyperpar[[1]]))$mean,
                             inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                           final_model_2b$internal.marginals.hyperpar[[2]]))$mean)
coefficients_2b$coefficients$INLA[9]<-summary(final_model_2b)$hyperpar[3, 1]*
  sqrt(coefficients_2b$coefficients$INLA[5]*coefficients_2b$coefficients$INLA[6])

#Set correlated errors
SD_s<-diag(4)
for (i in 1:4){
  SD_s[i,i]<-sqrt(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),final_model_2b$internal.marginals.hyperpar[[(i+3)]]))$quant0.5)
}
cor<-matrix(1, nrow=4, ncol=4)
cor[lower.tri(cor)]<-final_model_2b$summary.hyperpar[8:13,1]
cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
coefficients_2b$VarCov$INLA<-SD_s %*% cor %*% SD_s

#All results
round(coefficients_2b$coefficients[,c(1,3,4)], digits=2)
lapply(coefficients_2b$VarCov,  round, digits=2)
