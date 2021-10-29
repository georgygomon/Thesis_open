library('INLA')
library(reshape2)
library(lme4)
library(nlme)
library(MCMCglmm)
library(MCMCvis)

matrix2latex <- function(matr) {
  printmrow <- function(x) {
    cat(cat(x,sep=" & "),"\\\\ \n")
  }
  cat("\\begin{bmatrix}","\n")
  body <- apply(matr,1,printmrow)
  cat("\\end{bmatrix}")
}


################################################################################
#Parameters
################################################################################
true_betas<-c(2,4,2.5,3, 1.5, 3.5)
U<-matrix(c(2,-1,
            -1,3),2, byrow=TRUE)

coefficients_2a<-data.frame(true= c(true_betas, diag(U), U[1,2]), 
                            lmer=rep(0,9),nlme=rep(0,9),
                           INLA=rep(0,9),  MCMCglmm=rep(0,9))
rownames(coefficients_2a)<-c('beta_0^1', 'beta_x^1', 'beta_t^1', 'beta_0^2', 'beta_x^2', 'beta_t^2', 'u_1', 'u_2', 'u_12')

################################################################################
#Model 2: Only random intercept
################################################################################
set.seed(2022)
### Total of N individuals
N <- 100
### n correlated measurements per individual
n <- 2

data_2a<- data.frame(
  id=rep(1:N, each=n),
  x=rnorm(N*n), 
  Period=rep(0:(n-1), times=N), 
  Period2=factor(rep(0:(n-1), times=N)))


errors1 <- matrix(rnorm(N*n), N)
errors2 <- matrix(rnorm(N*n), N)
u<-matrix(rnorm(N*2), N) %*% chol(U)
#matrix2latex(round(cov(u), 2))
coefficients_2a$true<-c(true_betas, diag(cov(u)), cov(u)[1,2])
data_2a$y1 <- true_betas[1] + true_betas[2]*data_2a$x + true_betas[3]*data_2a$Period+rep(u[,1], each=n)+as.vector(t(errors1))
data_2a$y2<- true_betas[4]+true_betas[5]*data_2a$x+true_betas[6]*data_2a$Period+rep(u[,2], each=n)+as.vector(t(errors2))

################################################################################
#Model 2a using nlme
################################################################################
#General information on LMER, NLME and MCMCglmm
#http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#singular-models-random-effect-variances-estimated-as-zero-or-correlations-estimated-as---1
data_melted_2a<-melt(data_2a, measure.vars=c('y1', 'y2'))
data_melted_2a<-data_melted_2a[with(data_melted_2a,order(id)),]

model_lme_2a<-lme(value~variable+variable:x+variable:Period-1,
               random = ~0+variable | id, 
               data=data_melted_2a, na.action=na.omit)


coefficients_2a$nlme<-c(model_lme_2a$coefficients$fixed[c(1,3,5, 2,4,6)],
                             diag(getVarCov(model_lme_2a, individual=1, type='random.effects')),
                             getVarCov(model_lme_2a, individual=1, type='random.effects')[1,2] )

################################################################################
#Model 2a using lmer
################################################################################
model_lmer_2a<-lmer(value~ variable+variable:x+variable:Period-1+(variable-1|id), na.action=na.omit, 
                 data=data_melted_2a)
sum_lmer_2a<-summary(model_lmer_2a)

coefficients_2a$lmer<-c(sum_lmer_2a$coefficients[c(1,3,5, 2,4, 6),1], 
                       sum_lmer_2a$varcor$id[1, 1], sum_lmer_2a$varcor$id[2, 2], 
                       sum_lmer_2a$varcor$id[1,2])


################################################################################
#Model 2a using MCMCglmm
################################################################################

#using an inverse Wishart distribution with parameter R11, R12, R21, R22, v
#Check out the Inverse Wishart priors!!!! Where does the variance come from?
prior = list(R = list(V = diag(2), n = 4),
             G = list(G1=list(V = matrix(c(2,-1,-1, 3), nrow=2), n = 4)))
MCMC_2a<-MCMCglmm(cbind(y1, y2) ~ trait+trait:x+trait:Period - 1,
                 random = ~ us(trait):units,  rcov = ~ idh(trait):units, 
                 family = rep("gaussian", 2),prior=prior, nitt = 20000, burnin = 2000,
                 thin=25, data = data_2a)
coefficients_2a$MCMCglmm<-c(MCMCsummary(MCMC_2a$Sol)[c(1,3,5, 2,4,6),1],MCMCsummary(MCMC_2a$VCV)[c(1,4,3),1])

################################################################################
#Model 2a using INLA 
################################################################################
INLA_data_2a<-data_2a
#Modelling as 2 separate likelihoods
fixed.effects_2a<-list(Intercept1=c(rep(1, n*N), rep(NA, n*N)),
                    Intercept2=c(rep(NA, n*N), rep(1, n*N)),
                    x1=c(INLA_data_2a$x, rep(NA, n*N)),
                    x2=c(rep(NA, n*N), INLA_data_2a$x), 
                    Period1=c(data_2a$Period, rep(NA, n*N)),
                    Period2=c(rep(NA, n*N), data_2a$Period))

random.effects_2a<-list(ID_all=c(INLA_data_2a$id, INLA_data_2a$id+N))
y_final_2a<-list(c(INLA_data_2a$y1, rep(NA, n*N)),
              c(rep(NA, n*N),INLA_data_2a$y2))
final_data_2a<-c(fixed.effects_2a, random.effects_2a)
final_data_2a$Y<-y_final_2a

formula.model_2a=Y~-1+f(Intercept1, model='linear')+  f(Intercept2, model='linear')+
  f(x1, model='linear')+  f(x2, model='linear')+
  f(Period1, model='linear')+f(Period2, model='linear')+
  f(ID_all, model="iid2d", n=2*N,hyper = list(theta1 = list(prior = "wishart2d", param =c(4,1,1, 0))))
final_model_2a<-inla(formula.model_2a, family = c("gaussian","gaussian"),
                    data = final_data_2a,verbose=TRUE)


coefficients_2a$INLA[1:8]<-c(final_model_2a$summary.fixed[c(1,3,5,2,4,6),1],
                       inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                     final_model_2a$internal.marginals.hyperpar[[3]]))$mean,
                       inla.zmarginal(inla.tmarginal(function(x) 1/exp(x),
                                                     final_model_2a$internal.marginals.hyperpar[[4]]))$mean)

coefficients_2a$INLA[9]<-summary(final_model_2a)$hyperpar[5, 1]*sqrt(coefficients_2a$INLA[5]*coefficients_2a$INLA[6])
coefficients_2a

#xtable(coefficients_2a)

