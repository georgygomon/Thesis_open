source('create_data_single_models.R')

################################################################################
################################################################################
#Fitting the JMM (see model 4.6, pages 23-25)
################################################################################
################################################################################

################################################################################
#Parameters
################################################################################
true_betas<-c(3,2,1.5, 4,3, 2.5)
#D: This is the Variance Covariance matrix D, see equation 4.7 on page 24 
D1<-matrix(c(2,-1,
             -1,3),2, byrow=TRUE)
D2<-matrix(c(3,1.5,
             1.5,4),2, byrow=TRUE)
D3<-matrix(c(1.25,1,-0.5,0.3),2, byrow=TRUE)
D<-rbind(cbind(D1, D3), cbind(t(D3), D2))

JMM_coefficients<-list(coefficients=data.frame(true= c(true_betas), 
                                              lmer=rep(0,6),nlme=rep(0,6),
                                              INLA=rep(0,6),  MCMCglmm=rep(0,6)),
                      VarCov=list(D=D, nlme=diag(4), lmer=diag(4), MCMCglmm=diag(4), INLA=diag(4)))
rownames(JMM_coefficients$coefficients)<-c('beta_0^x', 'beta_w^x','beta_t^x', 'beta_0^y', 'beta_w^y','beta_t^y')
################################################################################
#Simulating data from JMM
################################################################################
### Total of N individuals
N <- 100
### n correlated measurements per individual
n <- 10
### Probabilities p1 & p2 that outcome y and covariate x respectively are observed
p1<-1; p2<-1

#Simulate data from JMM
JMM_data<-Create_JMM_data(true_betas,
                JMM_coefficients$VarCov,
                N,n,p=c(p1,p2))$data

################################################################################
#JMM using nlme
################################################################################
#Putting data in long format
JMM_data_melted<-melt(JMM_data, measure.vars=c('x', 'y'))
JMM_data_melted<-JMM_data_melted[with(JMM_data_melted,order(id)),]

#https://stats.stackexchange.com/questions/58669/specifying-multiple-separate-random-effects-in-lme
JMM_lme<-lme(value~variable+variable:w+variable:t-1,
                  random =list(id=~variable+variable:t+0), 
                  data=JMM_data_melted, na.action=na.omit)

#Extracting model coefficients
JMM_coefficients$coefficients$nlme<-JMM_lme$coefficients$fixed[c(1,3,5,2,4,6)]
#Extracting variance covariance matrix D
JMM_coefficients$VarCov$nlme<-matrix(getVarCov(JMM_lme), nrow=4)
################################################################################
#JMM using lmer
################################################################################
JMM_lmer<-lmer(value~ variable+variable:w+variable:t-1+(variable+0+variable:t|id), 
               na.action=na.omit, data=JMM_data_melted, 
                    control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

#Extracting model coefficients
JMM_coefficients$coefficients$lmer<-summary(JMM_lmer)$coefficients[c(1,3,5,2,4,6),1]
#Extracting variance covariance matrix D
JMM_coefficients$VarCov$lmer<-matrix(summary(JMM_lmer)$varcor$id[c(1:4),c(1:4)], nrow=4)
################################################################################
#JMM using MCMCglmm
################################################################################
#Using an inverse Wishart distribution with parameter R11, R12, R21, R22, v
prior = list(G = list(G1=list(V = JMM_coefficients$VarCov$D, n = 4)),
             R = list(V = diag(2*n), n = 4))

#MCMCglmm needs a factor in the residual covariance definition
JMM_data$t2<-factor(JMM_data$t)
JMM_MCMCglmm<-MCMCglmm(cbind(x, y) ~ trait+trait:w +trait:t- 1,
                   random = ~ us(trait:t+trait):units,  rcov = ~ idh(trait:t2):units, 
                   family = rep("gaussian", 2), prior=prior, nitt = 10000, burnin = 2000,
                   thin=25, data = JMM_data)

#Extracting model coefficients
JMM_coefficients$coefficients$MCMCglmm<-c(MCMCsummary(JMM_MCMCglmm$Sol)[c(1,3,5,2,4,6),1])
#Extracting variance covariance matrix D
JMM_coefficients$VarCov$MCMCglmm<-matrix(MCMCsummary(JMM_MCMCglmm$VCV)[c(1:16),1], nrow=4)

################################################################################
#JMM using INLA
################################################################################
#Use copy feature to create dependence between Intercept and Slope
#https://arxiv.org/pdf/1210.0333.pdf
fixed.effects<-list(Intercept_x=c(rep(1, n*N), rep(NA, n*N)),
                        Intercept_y=c(rep(NA, n*N), rep(1, n*N)),
                        w_x=c(JMM_data$w, rep(NA, n*N)),
                        w_y=c(rep(NA, n*N), JMM_data$w), 
                        t_x=c(JMM_data$t, rep(NA, n*N)),
                        t_y=c(rep(NA, n*N), JMM_data$t),
                        t=c(JMM_data$t, JMM_data$t))
#The 4 random effects should be jointly distributed, so they are indexed from 1 till 4N
#With N being the total number of subjects (as both intercept and slope are unique to a subject)
random.effects<-list(Random_Intercept=c(JMM_data$id, JMM_data$id+N),
                         Random_Slope=c(JMM_data$id+2*N, JMM_data$id+3*N))

INLA_data<-c(fixed.effects, random.effects)
INLA_data$Y<-list(c(JMM_data$x, rep(NA, n*N)),
                  c(rep(NA, n*N),JMM_data$y))

INLA_formula=Y~-1+f(Intercept_x, model='linear')+  f(Intercept_y, model='linear')+
  f(w_x, model='linear')+  f(w_y, model='linear')+f(t_x, model='linear')+f(t_y, model='linear')+
  f(Random_Intercept, model="iid4d", n=4*N)+  f(Random_Slope, t,copy="Random_Intercept")

JMM_INLA<-inla(INLA_formula, family = c("gaussian","gaussian"),
                      data = INLA_data,verbose=TRUE, control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))


#Extracting model coefficients
JMM_coefficients$coefficients$INLA<-JMM_INLA$summary.fixed[c(1,3,5,2,4,6),1]
#Extracting variance-covariance matrix D
JMM_coefficients$VarCov$INLA<-obtainVarCov(JMM_INLA, c(3:12), 4)$Cov


################################################################################
#Results
################################################################################
#All methods give very similar model coefficients
JMM_coefficients$coefficients
#Also, all models give similar variance covariance matrices
JMM_coefficients$VarCov


