source('Make_models.R')
source('Fit_models.R')


#Data of already performed analysis is saved!
#For the large simulation: N=750, n=2, 8 models
load('Model_comparison_big.RData')
Model_comparison_big
#For the small simulation: N=500, n=4, 6 models
load('Model_comparison_smaller.RData')
Model_comparison_big


#Overall parameters
true_betas<-c(2,4,2.5,3, 1.5, 3.5)
N<-1000
n<-2

df_standard=data.frame(true=c(true_betas, rep(NA, 4)),
                       Model_0=c(rep(NA, 10)),
                       Model_1A=c(rep(NA, 10)),
                       Model_2A=c(rep(NA, 10)),
                       Model_2B=c(rep(NA, 10)),
                       Model_2C1=c(rep(NA, 10)),
                       Model_2C2=c(rep(NA, 10)),
                       Model_3A1=c(rep(NA, 10)),
                       Model_3B1=c(rep(NA, 10)))
rownames(df_standard)<-c( 'beta_0^1', 'beta_x^1','beta_t^1', 'beta_0^2', 'beta_x^2','beta_t^2', 'MLIK', 'DIC', 'WAIC', 'CPO')


Model_comparison_big<-list(Model_0=list(Model_0_df=df_standard,
                                        parameters=list(U1=matrix(c(2,1,
                                                                    1,3),2, byrow=TRUE),
                                                        U2=matrix(c(3,1.5,
                                                                    1.5,4),2, byrow=TRUE))),
                           Model_1A=list(Model_1A_df=df_standard,
                                         parameters=list(Sigma=rbind(cbind(matrix(c(5,2,2,4), 2), matrix(c(-1, -0.5, -0.5, -1.5), 2)), 
                                                                     cbind(t(matrix(c(-1, -0.5, -0.5, -1.5), 2)), matrix(c(4,1,1,4),2))))),
                           Model_2A=list(Model_2A_df=df_standard,
                                         parameters=list(U_0=matrix(c(2,-1,
                                                                      -1,3),2, byrow=TRUE))),
                           Model_2B=list(Model_2B_df=df_standard,
                                         parameters=list(U_0=matrix(c(2,-1,
                                                                      -1,3),2, byrow=TRUE),
                                                         Sigma=rbind(cbind(matrix(c(5,2,2,4), 2), matrix(c(-1, -0.5, -0.5, -1.5), 2)), 
                                                                     cbind(t(matrix(c(-1, -0.5, -0.5, -1.5), 2)), matrix(c(4,1,1,4),2))))),
                           Model_2C1=list(Model_2C1_df=df_standard,
                                          parameters=list(U=rbind(cbind(matrix(c(2,-1,-1,3),2, byrow=TRUE), matrix(c(0,0,0,0),2, byrow=TRUE)), 
                                                                  cbind(t(matrix(c(0,0,0,0),2, byrow=TRUE)), matrix(c(3,1.5,1.5,4),2, byrow=TRUE))))),
                           Model_2C2=list(Model_2C2_df=df_standard,
                                          parameters=list(U=rbind(cbind(matrix(c(2,-1,-1,3),2, byrow=TRUE), matrix(c(1.25,1,-0.5,0.3),2, byrow=TRUE)), 
                                                                  cbind(t(matrix(c(1.25,1,-0.5,0.3),2, byrow=TRUE)), matrix(c(3,1.5,1.5,4),2, byrow=TRUE))))),
                           Model_3A1=list(Model_3A1_df=df_standard,
                                          parameters=list(U1=matrix(c(2,0,
                                                                       0,3),2, byrow=TRUE),
                                                          U2=matrix(c(3,1.5,
                                                                       1.5,4),2, byrow=TRUE),
                                                          gamma=0.45)),
                           Model_3B1=list(Model_3B1_df=df_standard,
                                          parameters=list(U1=matrix(c(2,0,
                                                                       0,3),2, byrow=TRUE),
                                                          U2=matrix(c(3,1.5,
                                                                       1.5,4),2, byrow=TRUE),
                                                          gamma1=0.45,
                                                          gamma2=1.2))
)

#Model 0
data<-make_model_0(true_betas, Model_comparison_big$Model_0$parameters$U1, Model_comparison_big$Model_0$parameters$U2, N, n)
print('Data-set 1/8')
Model_comparison_big$Model_0$Model_0_df<-fit_all_models(data$data, N, n, Model_comparison_big$Model_0$Model_0_df)
#Model 1A
data<-make_model_1A(true_betas, Model_comparison_big$Model_1A$parameters$Sigma, N, n)
print('Data-set 2/8')
Model_comparison_big$Model_1A$Model_1A_df<-fit_all_models(data$data, N, n, Model_comparison_big$Model_1A$Model_1A_df)
#Model 2A
data<-make_model_2A(true_betas, Model_comparison_big$Model_2A$parameters$U_0, N, n)
print('Data-set 3/8')
Model_comparison_big$Model_2A$Model_2A_df<-fit_all_models(data$data, N, n, Model_comparison_big$Model_2A$Model_2A_df)
#Model 2B
data<-make_model_2B(true_betas, Model_comparison_big$Model_2B$parameters$U_0,Model_comparison_big$Model_2B$parameters$Sigma, N, n)
print('Data-set 4/8')
Model_comparison_big$Model_2B$Model_2B_df<-fit_all_models(data$data, N, n, Model_comparison_big$Model_2B$Model_2B_df)
#Model 2C1
data<-make_model_2C(true_betas, Model_comparison_big$Model_2C1$parameters$U, N, n)
print('Data-set 5/8')
Model_comparison_big$Model_2C1$Model_2C1_df<-fit_all_models(data$data, N, n, Model_comparison_big$Model_2C1$Model_2C1_df)
#Model 2C2
data<-make_model_2C(true_betas, Model_comparison_big$Model_2C2$parameters$U, N, n)
print('Data-set 6/8')
Model_comparison_big$Model_2C2$Model_2C2_df<-fit_all_models(data$data, N, n, Model_comparison_big$Model_2C2$Model_2C2_df)
#Model 3A1
data<-make_model_3A(true_betas, Model_comparison_big$Model_3A1$parameters$U1,
                    Model_comparison_big$Model_3A1$parameters$U2, Model_comparison_big$Model_3A1$parameters$gamma, N, n)
print('Data-set 7/8')
Model_comparison_big$Model_3A1$Model_3A1_df<-fit_all_models(data$data, N, n, Model_comparison_big$Model_3A1$Model_3A1_df)
#Model 3B1
data<-make_model_3B(true_betas, Model_comparison_big$Model_3B1$parameters$U1,
                    Model_comparison_big$Model_3B1$parameters$U2, Model_comparison_big$Model_3B1$parameters$gamma1,
                    Model_comparison_big$Model_3B1$parameters$gamma2, N, n)
print('Data-set 8/8')
Model_comparison_big$Model_3B1$Model_3B1_df<-fit_all_models(data$data, N, n, Model_comparison_big$Model_3B1$Model_3B1_df)


load('Model_comparison_big.RData')
Model_comparison_big
load('Model_comparison_smaller.RData')
Model_comparison_big
