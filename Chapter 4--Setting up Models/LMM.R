source('create_data_single_models.R')

################################################################################
################################################################################
#Fitting the LMM (see model 4.1.1, pages 19-21)
################################################################################
################################################################################

################################################################################
#Parameters
################################################################################
#Setting the variance-covariance matrices D1 and D2.
#Variance-covarince matrix D1 is needed to simulate the exogenous covariate v.
#Our interest when modelling lies only in correctly modelling variance-covariance 
#matrix D2 (relating to outcome y).
parameters<-list(D1=matrix(c(2,1,
                             1,3),2, byrow=TRUE),
                 D2=matrix(c(3,1.5,
                             1.5,4),2, byrow=TRUE),
                 beta_v=1.2)
#Setting model coefficients. For interpretation see further
true_betas<-c(2,4,2.5,3, 1.5, 3.5)
#Combining everything
LMM_coefficients<-list(coefficients=data.frame(real=c(parameters$beta_v, true_betas[1:3]), INLA=rep(NA, 4)),
                      parameters=parameters,
                      INLA_parameters=list(D2=diag(2)))
#The interpretation of the model coefficients true_betas
rownames(LMM_coefficients$coefficients)<-c('beta_v^y', 'beta_0^y', 'beta_w^y','beta_t^y')

################################################################################
#Simulating data from LMM
################################################################################
### Total of N individuals
N <- 250
### Max # of measurements per individual
n <- 10
### Probabilities p1 & p2 that outcome y and covariate v respectively are observed
p1<-1; p2<-1
#Simulate data from LMM
LMM_data<-Create_LMM_data(true_betas,
                LMM_coefficients$parameters,
                N,n,p=c(p1,p2))$data


################################################################################
#LMM using INLA 
################################################################################
#As our covariate x is an exogenous covariate we shall be calling it v instead of x
colnames(LMM_data)[colnames(LMM_data) %in% c("observed_x", "x")] <- c("observed_v", "v")

N_obs<-dim(LMM_data)[1]
#We create 2 independent models for both outcome y and exogenous covariate v
fixed.effects<-list(Intercept=rep(1, N_obs),
                      t=LMM_data$t,
                      w=LMM_data$w,
                      v=LMM_data$v,
                      t=LMM_data$t)
random.effects<-list(Intercept_random=c(LMM_data$id),
                       t_random=c(LMM_data$id+N))

INLA_data<-c(fixed.effects, random.effects)
INLA_data$Y<-LMM_data$y

INLA_formula=Y~-1+
  f(v, model='linear')+f(Intercept, model='linear')+f(w, model='linear')+f(t, model='linear')+ #Model  fixed effects
  f(Intercept_random, model="iid2d", n=2*N)+f(t_random, t, copy="Intercept_random") #Model 1 random effects

LMM_INLA<-inla(INLA_formula, family = "gaussian",
                    data = INLA_data, verbose=TRUE,
                    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))

#Extracting model coefficients
LMM_coefficients$coefficients$INLA<-LMM_INLA$summary.fixed[,1]
#Extracting random effects variance-covariance matrices D1 and D2
LMM_coefficients$INLA_parameters$D2<-obtainVarCov(LMM_INLA, c(2:4), 2)$Cov


################################################################################
#Results
################################################################################
#One can see that INLA correctly computes the model coeffcients
LMM_coefficients$coefficients
#The covariance matrixes D1 and D2 (see equation 5.1 in thesis) are also well approximated
LMM_coefficients$parameters$D2
LMM_coefficients$INLA_parameters$D2

