source('create_data_single_models.R')

################################################################################
################################################################################
#Fitting the Multivariate Joint Model (see definition 4.4.1, pages 21-23)
################################################################################
################################################################################

################################################################################
#Parameters
################################################################################
#Model coefficients
true_betas<-c(2,4,2.5,3, 1.5, 3.5)
#Setting up the residual errors variance-covariance matrix
Sigma1 <- matrix(c(5,2,2,4), 2)
Sigma2<-matrix(c(4,1,1,4),2)
Sigma3<-matrix(c(-1, -0.5, -0.5, -1.5), 2)
Sigma<-rbind(cbind(Sigma1, Sigma3), cbind(t(Sigma3), Sigma2))

Multivariate_coefficients<-list(coefficients=data.frame(true=true_betas, gls=rep(0,6), 
                                             MCMCglmm=rep(0,6), INLA=rep(0,6)),
                     VarCov=list(true=Sigma, gls=NA, MCMCglmm=NA, INLA=NA))
rownames(Multivariate_coefficients$coefficients)<-c( 'beta_0^x', 'beta_w^x', 'beta_t^x', 'beta_0^y', 'beta_w^y', 'beta_t^y')

################################################################################
#Model 1: Simulate data
################################################################################
set.seed(2022)
### Total of N individuals
N <- 1250
### n correlated measurements per individual
n <- 2
### Probabilities p1 & p2 that outcome y and covariate x respectively are observed
p1<-1; p2<-1

Multivariate_data<-Create_multivariate_data(true_betas, 
                                            Multivariate_coefficients$VarCov,
                                            N,n,p=c(p1,p2))$data

################################################################################
#Multivariate using gls
################################################################################
#First melt data
data_melted<-melt(Multivariate_data, measure.vars=c('x', 'y'))
data_melted<-data_melted[with(data_melted,order(id)),]
#Fit gls model
Multivariate_gls<-gls(value~variable+variable:w+variable:t-1,
                      correlation=corSymm(form=~1|id), 
                      weights= varIdent(form=~1|variable*t),
                      data=data_melted, na.action=na.omit, method='REML')

Multivariate_coefficients$coefficients$gls<-summary(Multivariate_gls)$tTable[c(1,3,5, 2,4,6),1]
Multivariate_coefficients$VarCov$gls<-getVarCov(Multivariate_gls)

################################################################################
#Multivariate using MCMCglmm
################################################################################
prior = list(B = list(mu = rep(0,6),V = 100*diag(6)), 
             R = list(V = diag(4), n = 1))

#MCMCglmm needs a factor in the residual covariance definition
Multivariate_data$t2<-factor(Multivariate_data$t)
Multivariate_MCMCglmm<-MCMCglmm(cbind(x, y) ~ trait+trait:w+trait:t - 1,
                 rcov = ~ us(t2:trait):units, 
                 family = rep("gaussian", 2), nitt = 10000, burnin = 1000,
                 prior=prior, thin=25, data = Multivariate_data)
Multivariate_coefficients$coefficients$MCMCglmm<-MCMCsummary(Multivariate_MCMCglmm$Sol)[c(1,3,5, 2,4,6),1]
Multivariate_coefficients$VarCov$MCMCglmm<-matrix(MCMCsummary(Multivariate_MCMCglmm$VCV)[,1], nrow=4, byrow=FALSE)

################################################################################
#Multivariate using INLA
################################################################################
fixed.effects<-list(Intercept_x=c(rep(1, n*N), rep(NA, n*N)),
                    Intercept_y=c(rep(NA, n*N), rep(1, n*N)),
                    w_x=c(Multivariate_data$w, rep(NA, n*N)),
                    w_y=c(rep(NA, n*N), Multivariate_data$w),
                    t_x=c(Multivariate_data$t , rep(NA, n*N)),
                    t_y=c(rep(NA, n*N), Multivariate_data$t))
INLA_data<-c(fixed.effects)
INLA_data$Y<-list(c(Multivariate_data$x, rep(NA, n*N)),
                  c(rep(NA, n*N),Multivariate_data$y))
#setting indixes
INLA_data$i<-c(1:(2*n*N))
#One can get rid of the Gaussian noise term by setting fixed=TRUE in the control family parameter
#https://www.r-inla.org/faq question 5
set_residuals = list(list(initial=15, fixed=T), list(initial=15, fixed=T))
INLA_formula=Y~-1+f(Intercept_x, model='linear')+f(Intercept_y, model='linear')+
  f(w_x, model='linear')+f(w_y, model='linear')+
  f(t_x, model='linear')+f(t_y, model='linear')+
  f(i, model="iid4d", n=2*n*N, hyper = list(theta1 = list(prior = "wishart4d")))
#As we need to write the error structure as random effects and since INLA can handle a 
#maximum of 5 iid random effects ("iid4d"), we can not include more than 2 measurements per subject
Multivariate_INLA<-inla(INLA_formula, family = c("gaussian","gaussian"),
                  data = INLA_data, verbose=TRUE, control.family=set_residuals,
                  control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))


#Extracting model coefficients
Multivariate_coefficients$coefficients$INLA<-Multivariate_INLA$summary.fixed[c(1,3,5, 2,4, 6),1]
#Extracting residual errors variance-covariance matrix
Multivariate_coefficients$VarCov$INLA<-obtainVarCov(Multivariate_INLA, c(1:10), 4)$Cov

################################################################################
#Results
################################################################################
#All methods give very similar model coefficients
Multivariate_coefficients$coefficients
#The difference in residual variance-covariance matrices between the models is quite large
Multivariate_coefficients$VarCov
 