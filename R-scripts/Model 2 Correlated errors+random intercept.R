source('Load_data.R')
coefficients_21<-list(coefficients=data.frame(gls=rep(0,12), gls_sd=rep(0,12),MCMCglmm=rep(0,12), MCMCglmm_sd=rep(0,12),
                                             INLA=rep(0,12), INLA_sd=rep(0,12)))
rownames(coefficients_21$coefficients)<-c( 'beta_0^t', 'beta_sf^t','beta_ses^t', 
                                          'beta_0^r', 'beta_sf^r', 'beta_ses^r',
                                          'u_t^2','u_r^2', 'rho^u_rt',
                                          'sigma^2_t', 'sigma^2_r', 'rho^e_tr')



################################################################################
#Model 2.1 using nlme
################################################################################

control <- lmeControl(maxIter=300, msMaxIter=300)
model_lme<-lme(value~variable+variable:sex+variable:ses-1,
               random = ~0+variable:Period | id, 
               data=data_long_melted, na.action=na.omit, method='REML')

sum_lme<-summary(model_lme)
sum_lme$tTable[,2]

?gls
gls_unstructured<-gls(value~variable+variable:ses+variable:sex-1,
                      correlation=corSymm(form=~1|Period/id), 
                      data=data_long_melted, na.action=na.omit, method='REML')

cov2cor(getVarCov(gls_unstructured, individual=3))

coefficients_2$nlme<-c(model_lme$coefficients$fixed[c(1,3,5,2,4,6)], 
                       VarCorr(model_lme)[1:2,1], VarCorr(model_lme)[2,3])


?gls
coefficients_2$nlme_sd<-c(sum_lme$tTable[c(1,3,5,2,4,6),2], rep(NA, 3))


################################################################################
#Model 2.1 using MCMCglmm
################################################################################

#Gaussian prior for fixed effects B with mean mu and covariance matrix V
#Inverse Wishart prior for R
R_prior<-diag(8)
R_prior[c(1:4), c(1:4)]<-0.5
R_prior[c(5:8), c(5:8)]<-0.5
prior = list(B = list(mu = rep(0,6),V = 100*diag(6)), 
             R = list(V = R_prior, n = 8))


m_pain<-MCMCglmm(cbind(pain_tolerance, pain_rating) ~ trait+trait:sex+trait:ses - 1,
                 rcov = ~ us(Period:trait):units, 
                 family = rep("gaussian", 2),prior=prior, nitt = 10000, burnin = 1000,
                 thin=25, data = data_long)

coefficients_1$coefficients[,3:4]<-MCMCsummary(m_pain$Sol)[c(1,3,5,2,4,6),c(1,2)]
MCMCsummary(m_pain$VCV)
matrix(MCMCsummary(m_pain$VCV)[,1], nrow=8, byrow=FALSE)
coefficients_1$CovCor$MCMCglmm<-cov2cor(matrix(MCMCsummary(m_pain$VCV)[,1], nrow=8, byrow=TRUE))