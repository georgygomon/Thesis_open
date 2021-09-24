source('Load_data.R')
coefficients_1<-list(coefficients=data.frame(gls=rep(0,6), gls_sd=rep(0,6),MCMCglmm=rep(0,6), MCMCglmm_sd=rep(0,6),
                           INLA=rep(0,6), INLA_sd=rep(0,6)),
                     CovCor=list(gls=NA, MCMCglmm=NA, INLA=NA))
rownames(coefficients_1$coefficients)<-c( 'beta_0^t', 'beta_sf^t','beta_ses^t', 'beta_0^r', 'beta_sf^r', 'beta_ses^r')

xtable(head(data_long[,c(1,2,5,6,13,14)]))

################################################################################
#Model 1 using gls
################################################################################

#Need to find out how to fit mor elaborate correlation structures!!!
#Now fitted a simple unstructured correlation matrix between the (4+4) measurements
#on every individual
gls_unstructured<-gls(value~variable+variable:ses+variable:sex-1,
               correlation=corSymm(form=~1|id), 
               weights= varIdent(form=~1|Period),
               data=data_long_melted, na.action=na.omit, method='REML')

coefficients_1$coefficients[,1:2]<-summary(gls_unstructured)$tTable[c(1,5,3,2,6,4), c(1,4)]
coefficients_1$CovCor$gls<-cov2cor(getVarCov(gls_unstructured, individual=63))

gls_structured<-gls(value~variable+variable:ses+variable:sex-1,
                      correlation=corAR1(form=~1|id),
                      data=data_long_melted, na.action=na.omit, method='REML')

cov2cor(getVarCov(gls_structured, individual=63))

gls_structured$
getVarCov(gls_structured, individual=3)

################################################################################
#Model 1 using MCMCglmm
################################################################################

#Gaussian prior for fixed effects B with mean mu and covariance matrix V
#Inverse Wishart prior for R
R_prior<-diag(8)
R_prior[c(1:4), c(1:4)]<-0.5
R_prior[c(5:8), c(5:8)]<-0.5
prior = list(B = list(mu = rep(0,6),V = 100*diag(6)), 
             R = list(V = R_prior, n = 1))


m_pain<-MCMCglmm(cbind(pain_tolerance, pain_rating) ~ trait+trait:sex+trait:ses - 1,
                 rcov = ~ us(Period:trait):units, 
                 family = rep("gaussian", 2),prior=prior, nitt = 100000, burnin = 10000,
                 thin=25, data = data_long)

coefficients_1$coefficients[,3:4]<-MCMCsummary(m_pain$Sol)[c(1,3,5,2,4,6),c(1,2)]
MCMCsummary(m_pain$VCV)
matrix(MCMCsummary(m_pain$VCV)[,1], nrow=8, byrow=FALSE)
coefficients_1$CovCor$MCMCglmm<-cov2cor(matrix(MCMCsummary(m_pain$VCV)[,1], nrow=8, byrow=TRUE))

xtable(coefficients_1$coefficients[,c(1:4)])
coefficients_1$CovCor$gls
coefficients_1$CovCor$MCMCglmm


################################################################################
#Model 1 using INLA
################################################################################
# fixed.effects<-list(Intercept1=c(rep(1, nLong), rep(NA, nLong)),
#                     Intercept2=c(rep(NA, nLong), rep(1, nLong)),
#                     ses1=c(data_INLA$ses, rep(0, nLong)),
#                     ses2=c(rep(0, nLong), data_INLA$ses),
#                     sex1=c(data_INLA$sex, rep(0, nLong)),
#                     sex2=c(rep(0, nLong), data_INLA$sex))
# 
# random.effects<-list(ID1=c(data_INLA$id, data_INLA$id),
#                      Period=c(data_INLA$Period, data_INLA$Period))
# 
# y_final<-list(c(y$pain_tolerance, rep(NA, nLong)),
#               c(rep(NA, nLong),y$pain_rating))
# final_data<-c(fixed.effects, random.effects)
# final_data$Y<-y_final
# 
# #Fixed effects with Gaussian prior with mean 0 and precision 0.001
# #Random effect with Inverse-Wishart prior. Values set are (r, R11, R22, R12)
# formula.model_easy=Y~-1+f(Intercept1, model='linear', mean.linear = 0, prec.linear = 0.001)+
#   f(Intercept2, model='linear', mean.linear = 0, prec.linear = 0.001)+
#   f(ses1, model='linear', mean.linear = 0, prec.linear = 0.001)+
#   f(ses2, model='linear', mean.linear = 0, prec.linear = 0.001)+
#   f(sex1, model='linear', mean.linear = 0, prec.linear = 0.001)+
#   f(sex2, model='linear', mean.linear = 0, prec.linear = 0.001)+
#   f(ID1, model="iid2d", n=2*63, hyper = list(theta = list(prior = "wishart2d",   
#                                                           param =c(4,10,0.45, -0.50))))+
#   f(Period, model = "ar1")
# final_model_easy<-inla(formula.model_easy, family = c("gaussian","gaussian"),
#                        data = final_data, verbose=TRUE)
# 
# summary(final_model_easy)
# final_model_easy$summary.hyperpar
# final_model_easy$summary.random
