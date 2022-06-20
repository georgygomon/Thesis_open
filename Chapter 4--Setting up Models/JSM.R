source('create_data_single_models.R')

################################################################################
################################################################################
#Fitting the JSM (see definition 4.6.1, pages 25-27)
################################################################################
################################################################################

################################################################################
#Parameters
################################################################################
#Setting the variance covariance matrices D1 and D2 (see page 25 of thesis)
D1<-matrix(c(2,1.5,
             1.5,3),2, byrow=TRUE)
D2<-matrix(c(3,1.5,
             1.5,4),2, byrow=TRUE)
#Setting scaling factor gamma
gamma<-1.3
#Setting model coefficients. For interpretation see further
true_betas<-c(2,4,2.5,3, 1.5, 3.5)
#Combining everything
JSM_coefficients<-list(coefficients=data.frame(true= c(true_betas, gamma),  INLA=rep(0,7), 
                                               true_combi=c(true_betas[1:3], true_betas[1:3]*gamma+true_betas[4:6], NA), INLA_combi=c(rep(0,6), NA)),
                       parameters=list(D1=D1, D2=D2, gamma=gamma),
                       INLA=list(D1=diag(2), D2=diag(2)))
#The interpretation of the model coefficients true_betas
rownames(JSM_coefficients$coefficients)<-c( 'beta_0^x', 'beta_w^x','beta_t^x', 'beta_0^y', 'beta_w^y','beta_t^y', 'gamma')
################################################################################
#Simulating data from JSM
################################################################################ 
### Total of N individuals
N <- 200
### Max # of measurements per individual
n <- 10
### Probabilities p1 & p2 that outcome y and covariate x respectively are observed
p1<-1; p2<-1

#Creating data
JSM_data<-Create_JSM_data(true_betas,
                          JSM_coefficients$parameters,
                          N,n,p=c(p1,p2))$data

################################################################################
#JSM using INLA
################################################################################
#Total number of observations
N_obs<-dim(JSM_data)[1]
#Fixed effects for outcome y
fixed.effects<-list(w=c(JSM_data$w, JSM_data$w), 
                    t=c(JSM_data$t, JSM_data$t),
                    Intercept_y=c(rep(NA, N_obs), rep(1, N_obs)),
                    t_y=c(rep(NA, N_obs),JSM_data$t),
                    w_y=c(rep(NA, N_obs),JSM_data$w))
#Random effects 
#Remember that to fit the JSM the fixed effects of the endogenous covariate x
#are written down as random effects
random.effects<-list(Intercept_x=c(rep(1, N_obs), rep(NA, N_obs)),
                     Intercept_x_scaled=c(rep(NA, N_obs), rep(1, N_obs)),
                     idw_x=c(rep(1, N_obs),rep(NA, N_obs)),
                     idw_x_scaled=c(rep(NA, N_obs),rep(1, N_obs)),
                     t_x=c(rep(1, N_obs),rep(NA, N_obs)),
                     t_x_scaled=c(rep(NA, N_obs),rep(1, N_obs)),
                     Random_intercept_x=c(JSM_data$id, rep(NA, N_obs)),
                     Random_intercept_x_scaled=c(rep(NA, N_obs), JSM_data$id),
                     Random_intercept_y=c(rep(NA, N_obs), JSM_data$id),
                     Random_Slope_x=c(JSM_data$id+N, rep(NA, N_obs)),
                     Random_Slope_x_scaled=c(rep(NA, N_obs), JSM_data$id+N),
                     Random_Slope_y=c(rep(NA, N_obs), JSM_data$id+N))

INLA_data<-c(fixed.effects, random.effects)
INLA_data$Y<-list(c(JSM_data$x, rep(NA, n*N)),
                  c(rep(NA, n*N),JSM_data$y))

INLA_formula=Y~-1+
  f(Intercept_x)+ #Fixed intercept in endogenous covariate function
  f(Intercept_x_scaled, copy="Intercept_x",
    hyper = list(beta = list(fixed=FALSE)))+ #Scaling that fixed intercept into outcome y
  f(idw_x, w)+
  f(idw_x_scaled, w, copy="idw_x", same.as = 'Intercept_x_scaled',
    hyper = list(beta = list(fixed=FALSE)))+
  f(t_x, t)+
  f(t_x_scaled, t, copy="t_x", same.as = 'Intercept_x_scaled',
    hyper = list(beta = list(fixed=FALSE)))+
  f(Intercept_y, model='linear')+f(w_y, model='linear')+f(t_y, model='linear')+ #Outcome fixed effects
  f(Random_intercept_x, model="iid2d", n=2*N)+f(Random_Slope_x, t, copy="Random_intercept_x")+ #Endogenous covariate random effects
  f(Random_intercept_x_scaled, copy='Random_intercept_x', same.as='Intercept_x_scaled', fixed=FALSE)+
  f(Random_Slope_x_scaled, t, copy='Random_intercept_x', same.as='Intercept_x_scaled', fixed=FALSE)+ #Scaled endogenous random effects
  f(Random_intercept_y, model="iid2d", n=2*N)+f(Random_Slope_y, t, copy="Random_intercept_y") #Outcome y random effects

JSM_INLA<-inla(INLA_formula, family = c("gaussian","gaussian"),
               data = INLA_data,verbose=TRUE,
               control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))

#Extract model coefficients
#Note that the fixed model coefficients in the linear predictor of the endogenous covariate x 
#are modelled as random effects and thus can be found in summary.random
x_coefficients=rbind(JSM_INLA$summary.random$Intercept_x[c(2,3,4,6)],
                     JSM_INLA$summary.random$idw_x[c(2,3,4,6)],
                     JSM_INLA$summary.random$t_x[c(2,3,4,6)])
rownames(x_coefficients)<-c('Intercept_x', 'w_x', 't_x')
gamma<-summary(JSM_INLA)$hyperpar[12,1]
JSM_coefficients$coefficients$INLA<-c(x_coefficients[,1],JSM_INLA$summary.fixed[,1], gamma)
JSM_coefficients$coefficients$INLA_combi[1:6]<-c(x_coefficients[,1], x_coefficients[,1]*gamma+JSM_INLA$summary.fixed[,1])

#Extracting variance covariance matrixes D1 and D2
JSM_coefficients$INLA$D1<-obtainVarCov(JSM_INLA, c(6:8), 2)$Cov
JSM_coefficients$INLA$D2<-obtainVarCov(JSM_INLA, c(9:11), 2)$Cov

################################################################################
#Results
################################################################################
#All methods give very similar model coefficients
JSM_coefficients$coefficients
#Also, all models give similar variance covariance matrices
JSM_coefficients[c(2,3)]

