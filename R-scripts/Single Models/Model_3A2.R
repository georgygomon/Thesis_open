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
U1<-matrix(c(2,1.5,
             1.5,3),2, byrow=TRUE)
U2<-matrix(c(3,1.5,
             1.5,4),2, byrow=TRUE)
gamma<-1.3
true_betas<-c(2,4,2.5,3, 1.5, 3.5)

coefficients_3a2<-list(coefficients=data.frame(true= c(true_betas, gamma),  INLA=rep(0,7), 
                                              true_combi=c(true_betas[1:3], true_betas[1:3]*gamma+true_betas[4:6], NA), INLA_combi=c(rep(0,6), NA)),
                      U1=list(true=U1, INLA=diag(2)),
                      U2=list(true=U2, INLA=diag(2)))


rownames(coefficients_3a2$coefficients)<-c( 'beta_0^1', 'beta_x^1','beta_t^1', 'beta_0^2', 'beta_x^2','beta_t^2', 'gamma')
################################################################################
#Model 2 Make data: Only random intercept, uncorrelated errors
################################################################################
set.seed(2022)
### Total of N individuals
N <- 175
### Max # of measurements per individual
n <- 7
### Probability p that a measurement is observed
p1<-1
### probability p2 that covariate x is measured
p2<-1

data_3a2 <- data.frame(
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
coefficients_3a2$U1$true<-cov(u1)
coefficients_3a2$U2$true<-cov(u2)
m_3a2<-true_betas[1]+rep(u1[,1], each=n)+true_betas[2]*data_3a2$x+(true_betas[3]+rep(u1[,2], each=n))*data_3a2$Period
data_3a2$y1 <-  m_3a2+as.vector(t(errors1))
data_3a2$y2<- gamma*(m_3a2)+true_betas[4]+rep(u2[,1], each=n)+true_betas[5]*data_3a2$x+(true_betas[6]+rep(u2[,2], each=n))*data_3a2$Period+as.vector(t(errors2))
data_3a2$x[data_3a2$x_measured==0]<-NA
data<-data_3a2[data_3a2$observed==1,]

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
#To copy an entire linear predictor see:
#https://www.r-inla.org/faq#h.pb4g1pwbmtli


INLA_data_3a2<-data
#Modelling as 2 separate likelihoods
N_obs<-length(INLA_data_3a2$x)
fixed.effects_3a2<-list(Intercept1=c(rep(1, 2*N_obs), rep(NA, N_obs)),
                       X1=c(INLA_data_3a2$x, INLA_data_3a2$x, rep(NA, N_obs)), 
                       Period1=c(INLA_data_3a2$Period, INLA_data_3a2$Period, rep(NA, N_obs)),
                       Intercept2=c(rep(NA, 2*N_obs), rep(1, N_obs)),
                       X2=c(rep(NA, N_obs), rep(NA, N_obs), INLA_data_3a2$x), 
                       Period2=c(rep(NA, N_obs), rep(NA, N_obs), INLA_data_3a2$Period),
                       Period=c(INLA_data_3a2$Period, INLA_data_3a2$Period, INLA_data_3a2$Period),
                       w=c(rep(NA, N_obs), rep(-1, N_obs), rep(NA, N_obs)))
random.effects_3a2<-list(Intercept_random1=c(INLA_data_3a2$id,INLA_data_3a2$id, rep(NA, N_obs)),
                        Period_random1=c(INLA_data_3a2$id+N, INLA_data_3a2$id+N, rep(NA, N_obs)),
                        Intercept_random2=c(rep(NA, 2*N_obs), INLA_data_3a2$id),
                        Period_random2=c(rep(NA, 2*N_obs), INLA_data_3a2$id+N),
                        u=c(rep(NA, N_obs), c(1:N_obs), rep(NA, N_obs)),
                        b.eta2=c(rep(NA, 2*N_obs),  c(1:N_obs)))
y_final_3a2<-list(c(INLA_data_3a2$y1, rep(NA, 2*N_obs)),
                 c(rep(NA, N_obs),rep(0, N_obs),  rep(NA, N_obs)),
                 c(rep(NA, 2*N_obs),  INLA_data_3a2$y2))
final_data_3a2<-c(fixed.effects_3a2, random.effects_3a2)
final_data_3a2$Y<-y_final_3a2

formula.model_3a2=Y~-1+Intercept1+X1+Period1+Intercept2+X2+Period2+
  f(Intercept_random1, model="iid2d", n=2*N_obs)+f(Period_random1, Period, copy="Intercept_random1")+
  f(Intercept_random2, model="iid2d", n=2*N_obs)+f(Period_random2, Period, copy="Intercept_random2")+
  f(u, w, model="iid", hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy="u", hyper = list(beta = list(fixed = FALSE)))

final_model_3a2<-inla(formula.model_3a2, family = c("gaussian","gaussian", "gaussian"),
                     data = final_data_3a2,verbose=TRUE, control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
                     control.family = list(list(),list(hyper = list(prec = list(initial = 10, fixed=TRUE))),list()))
summary(final_model_3a2)

gamma<-summary(final_model_3a2)$hyperpar[9,1]
coefficients_3a2$coefficients$INLA<-model_coefficients<-c(final_model_3a2$summary.fixed[,1], gamma)
coefficients_3a2$coefficients$INLA_combi[1:6]<-actual_coefficients<-c(model_coefficients[1:3], model_coefficients[1:3]*gamma+model_coefficients[4:6])

SD_s<-diag(2)
names(final_model_3a2$internal.marginals.hyperpar)
tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_3a2$internal.marginals.hyperpar[[3]]), silent=TRUE)$mean,
                      inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_3a2$internal.marginals.hyperpar[[4]]), silent=TRUE)$mean)),
         error = function(e) print("Model 3A2: Some covariance elements of Y1 not well defined"))
cor<-diag(2)
cor[lower.tri(cor)]<-summary(final_model_3a2)$hyperpar[5, 6]
cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
coefficients_3a2$U1$INLA<-Cov1_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)

SD_s<-diag(2)
tryCatch(SD_s<-sqrt(c(inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_3a2$internal.marginals.hyperpar[[6]]), silent=TRUE)$mean,
                      inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                    final_model_3a2$internal.marginals.hyperpar[[7]]), silent=TRUE)$mean)),
         error = function(e) print("Model 3A1: Some covariance elements of Y2 not well defined"))
cor<-diag(2)
cor[lower.tri(cor)]<-summary(final_model_3a2)$hyperpar[8, 6]
cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
coefficients_3a2$U2$INLA<-Cov2_0_t<-diag(SD_s) %*% cor %*% diag(SD_s)

coefficients_3a2
# 
# 
# coefficients_3a$coefficients$INLA<-c(model_coefficients, gamma)
# coefficients_3a$coefficients$INLA_combi<-c(actual_coefficients, NA)
# coefficients_3a$U1$INLA<-Cov1_0_t
# coefficients_3a$U2$INLA<-Cov2_0_t
# 
# coefficients_3a