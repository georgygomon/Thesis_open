source('Load_data.R')

coefficients_3<-data.frame(INLA=rep(0,6), INLA_sd=rep(0,6))
rownames(coefficients_3)<-c( 'beta_0^t', 'b_ses^t', 'b_0^t', 'beta_0^r', 'gamma_1', 'gamma_2')
################################################################################
#Model 3 specified wrongly
################################################################################

#https://groups.google.com/g/r-inla-discussion-group/c/ClNVlx1lgwY
#https://arxiv.org/pdf/1210.0333.pdf
#One can copy only random effects. If one wants to copy fixed effect simply write it as random effect
fixed.effects<-list(Intercept1=c(rep(1, nLong), rep(NA, nLong)),
                    Intercept2=c(rep(NA, nLong), rep(1, nLong)),
                    ses1=c(data_INLA$ses, rep(NA, nLong)),
                    ses2=c(rep(NA, nLong), data_INLA$ses),
                    sex1=c(data_INLA$sex, rep(NA, nLong)),
                    sex2=c(rep(NA, nLong), data_INLA$sex), 
                    ses=c(data_INLA$ses, data_INLA$ses), 
                    sex=as.factor(c(data_INLA$sex,data_INLA$sex)))

sestol = c(rep(1, nLong), rep(NA, nLong))
sesrat = c(rep(NA, nLong), rep(1, nLong))
idtol = c(rep(1, nLong), rep(NA, nLong))
idrat = c(rep(NA, nLong), rep(1, nLong))

random.effects<-list(ID1=c(data_INLA$id, rep(NA, nLong)),
                     ID2=c(rep(NA, nLong), data_INLA$id),
                     ID=c(data_INLA$id, data_INLA$id),
                     Period=c(data_INLA$Period, data_INLA$Period))

y_final<-list(c(y$pain_tolerance, rep(NA, nLong)),
              c(rep(NA, nLong),y$pain_rating))
final_data<-c(fixed.effects, random.effects)
final_data$Y<-y_final

#Fixed effects with Gaussian prior with mean 0 and precision 0.001
#Random effect with Inverse-Wishart prior. Values set are (r, R11, R22, R12)
formula.model_3=Y~-1+f(Intercept1, model='linear', mean.linear = 0, prec.linear = 0.001)+
  f(Intercept2, model='linear', mean.linear = 0, prec.linear = 0.001)+
  f(sestol) +
  f(sesrat, copy="sestol",hyper = list(beta = list(fixed=TRUE)))+
  f(idtol, ID) +
  f(idrat, ID, copy="idtol",hyper = list(beta = list(fixed=TRUE)))
final_model_3<-inla(formula.model_3, family = c("gaussian","gaussian"),
                       data = final_data, verbose=TRUE)


summary(final_model_3)
final_model_3$summary.random
coefficients_3$INLA[c(1,4)]<-final_model_3$summary.fixed[,1]
coefficients_3$INLA[c(5,6)]<-final_model_3$summary.hyperpar[c(5,6),1]
coefficients_3$INLA[c(2,3)]<-1/final_model_3$summary.hyperpar[c(3,4),1]

coefficients_3$INLA_sd[c(1,4)]<-final_model_3$summary.fixed[,2]
coefficients_3$INLA_sd[c(5,6)]<-final_model_3$summary.hyperpar[c(5,6),2]

xtable(coefficients_3, digits=5)
?xtable


# fixed.effects<-list(Intercept1=c(rep(1, nLong), rep(NA, nLong)),
#                     Intercept2=c(rep(NA, nLong), rep(1, nLong)),
#                     ses1=c(data_INLA$ses, rep(NA, nLong)),
#                     ses2=c(rep(NA, nLong), data_INLA$ses),
#                     sex1=c(data_INLA$sex, rep(NA, nLong)),
#                     sex2=c(rep(NA, nLong), data_INLA$sex),
#                     pain_tolerance=c(rep(NA, nLong), y$pain_tolerance))
# 
# random.effects<-list(ID1=c(data_INLA$id, rep(NA, nLong)),
#                      ID2=c(rep(NA, nLong), data_INLA$id),
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
#   f(pain_tolerance, model='linear', mean.linear = 0, prec.linear = 0.001)+
#   f(ID1, model="iid")+  f(ID2, model="iid")
# final_model_easy<-inla(formula.model_easy, family = c("gaussian","gaussian"),
#                        data = final_data, verbose=TRUE)
# 
# summary(final_model_easy)





