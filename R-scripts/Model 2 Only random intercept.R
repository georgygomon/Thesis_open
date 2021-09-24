source('Load_data.R')
coefficients_2<-data.frame(lmer=rep(0,9), lmer_sd=rep(0,9), nlme=rep(0,9), nlme_sd=rep(0,9),
                           INLA=rep(0,9), INLA_sd=rep(0,9),  MCMCglmm=rep(0,9), MCMCglmm_sd=rep(0,9))
rownames(coefficients_2)<-c( 'beta_0^t', 'beta_sf^t','beta_ses^t', 'beta_0^r', 'beta_sf^r', 'beta_ses^r',
                             'u_t^2','u_r^2', 'rho_rt')

################################################################################
#Model 2 using lmer without correlation between terms
################################################################################
model_lmer<-lmer(value~ variable+variable:sex+variable:ses-1+(variable-1|id), na.action=na.omit, 
                 data=data_long_melted)
sum_lmer<-summary(model_lmer)

coefficients_2$lmer<-c(sum_lmer$coefficients[c(1,3,5,2,4,6),1], 
                       sum_lmer$varcor$id[1, 1], sum_lmer$varcor$id[2, 2], 
                       attributes(sum_lmer$varcor$id)$correlation[1,2])
coefficients_2$lmer_sd<-c(sum_lmer$coefficients[c(1,3,5,2,4,6),2], rep(NA, 3))


################################################################################
#Model 2 using nlme without correlation between terms
################################################################################


model_lme<-lme(value~variable+variable:sex+variable:ses-1,
               random = ~0+variable | id, data=data_long_melted, na.action=na.omit)

sum_lme<-summary(model_lme)
sum_lme$tTable[,2]
coefficients_2$nlme<-as.numeric(c(model_lme$coefficients$fixed[c(1,3,5,2,4,6)], 
                       VarCorr(model_lme)[1:2,1], VarCorr(model_lme)[2,3]))


coefficients_2$nlme_sd<-c(sum_lme$tTable[c(1,3,5,2,4,6),2], rep(NA, 3))


################################################################################
#Model 2 using INLA without correlation between terms
################################################################################
fixed.effects<-list(Intercept1=c(rep(1, nLong), rep(NA, nLong)),
                    Intercept2=c(rep(NA, nLong), rep(1, nLong)),
                    ses1=c(data_INLA$ses, rep(0, nLong)),
                    ses2=c(rep(NA, nLong), data_INLA$ses),
                    sex1=c(data_INLA$sex, rep(NA, nLong)),
                    sex2=c(rep(NA, nLong), data_INLA$sex))

random.effects<-list(ID1=c(data_INLA$id, data_INLA$id),
                     Period=c(data_INLA$Period, data_INLA$Period))

y_final<-list(c(y$pain_tolerance, rep(NA, nLong)),
              c(rep(NA, nLong),y$pain_rating))
final_data<-c(fixed.effects, random.effects)
final_data$Y<-y_final

#Fixed effects with Gaussian prior with mean 0 and precision 0.001
#Random effect with Inverse-Wishart prior. Values set are (r, R11, R22, R12)
formula.model_2=Y~-1+f(Intercept1, model='linear', mean.linear = 0, prec.linear = 0.001)+
  f(Intercept2, model='linear', mean.linear = 0, prec.linear = 0.001)+
  f(ses1, model='linear', mean.linear = 0, prec.linear = 0.001)+
  f(ses2, model='linear', mean.linear = 0, prec.linear = 0.001)+
  f(sex1, model='linear', mean.linear = 0, prec.linear = 0.001)+
  f(sex2, model='linear', mean.linear = 0, prec.linear = 0.001)+
  f(ID1, model="iid2d", n=2*252, hyper = list(theta = list(prior = "wishart2d",   
                                                          param =c(1,10,0.45, -0.50))))
final_model_2<-inla(formula.model_2, family = c("gaussian","gaussian"),
                       data = final_data, verbose=TRUE)
summary(final_model_2)
final_model_2$summary.hyperpar
final_model_2$summary.random

#Check why INLA random effects are so different!!!!! Focus on priors for the random effects!!
#https://inla.r-inla-download.org/r-inla.org/doc/latent/iid.pdf
coefficients_2$INLA<-c(final_model_2$summary.fixed[c(1,5,3,2,6,4),1], 
                       1/final_model_2$summary.hyperpar[c(3,4),6], 
                       final_model_2$summary.hyperpar[5,6])

coefficients_2$INLA_sd<-c(final_model_2$summary.fixed[c(1,5,3,2,6,4),2], rep(NA, 3))


################################################################################
#Model 2 using MCMCglmma without correlation between terms
################################################################################
head(data_long)

#using an inverse Wishart distribution with parameter R11, R12, R21, R22, v
#Check out the Inverse Wishart priors!!!! Where does the variance come from?
prior = list(R = list(V = diag(2), n = 4),
             G = list(G1=list(V = matrix(c(10,-0.50,-0.50, 0.45), nrow=2), n = 1)))

m_pain<-MCMCglmm(cbind(pain_tolerance, pain_rating) ~ trait+trait:sex+trait:ses - 1,
                 random = ~ us(trait):units, rcov = ~ idh(trait):units, 
                 family = rep("gaussian", 2),prior=prior, nitt = 10000, burnin = 1000,
                 thin=25, data = data_long)



coefficients_2$MCMCglmm<-c(MCMCsummary(m_pain$Sol)[c(1,3,5,2,4,6),1],MCMCsummary(m_pain$VCV)[c(1,4,3),1])
coefficients_2$MCMCglmm_sd<-c(MCMCsummary(m_pain$Sol)[c(1,3,5,2,4,6),2], rep(NA, 3))



coefficients_2$MCMCglmm
coefficients_2$INLA_sd<-c(final_model_easy$summary.fixed[c(1,5,3,2,6,4),2], rep(NA, 3))


#MCMCtrace(m_pain$Sol, pdf=FALSE)
#MCMCplot(m_pain$Sol)
coefficients_2[1,1]
colnames(coefficients_2)
round(as.data.frame(coefficients_2))
typeof(coefficients_2)

xtable(coefficients_2)
?xtable
round(coefficients_2)
coefficients_2$nlme
