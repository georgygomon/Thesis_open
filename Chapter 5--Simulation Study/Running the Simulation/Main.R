source('Make_models.R') #Simulates data
source('Fit_models.R')  #Fits models

################################################################################
#Model Parameters
################################################################################
#D are the variance-covariance matrices
#LMM parameters
LMM_parameters=list(D1=matrix(c(2,1.5,1.5,3),2, byrow=TRUE),
                    D2=matrix(c(3,2.5,2.5,4),2, byrow=TRUE),
                    beta_x=1.2)
rownames(LMM_parameters$D1)<-colnames(LMM_parameters$D1)<-c('d_0x', 'd_tx')
rownames(LMM_parameters$D2)<-colnames(LMM_parameters$D2)<-c('d_0y', 'd_ty')
#JMM parameters
JMM_parameters=list(D=rbind(cbind(matrix(c(2,1.5,1.5,3),2, byrow=TRUE), 
                                         matrix(c(1.75,1.6,2,2.5),2, byrow=TRUE)),
                                   cbind(t(matrix(c(1.75,1.6,2,2.5),2, byrow=TRUE)), 
                                         matrix(c(3,2.6,2.6,4),2, byrow=TRUE))))
rownames(JMM_parameters$D)<-colnames(JMM_parameters$D)<-c('d_x0', 'd_y0', 'd_xt', 'd_yt')
#JSM parameters
JSM_parameters=list(D1=matrix(c(2,1.5,1.5,3),2, byrow=TRUE),
                    D2=matrix(c(3,2.5,2.5,4),2, byrow=TRUE),
                    gamma=1.2)
rownames(JSM_parameters$D1)<-colnames(JSM_parameters$D1)<-c('d_x0', 'd_xt')
rownames(JSM_parameters$D2)<-colnames(JSM_parameters$D2)<-c('d_y0', 'd_yt')
#The same model coefficients are used in all models
true_betas<-c(5,3,1,6, 4, 2.2)
names(true_betas)<-c('b_x0', 'b_xw', 'b_xt', 'b_y0', 'b_yw', 'b_yt')

################################################################################
#Other parameters
################################################################################
### Probabilities that measurements of outcome y and endogenous covariate x respectively are observed
p<-c(1, 1)
### Number of measurements per individual: Training set
n<-10
#Ratio of extra measurements used for subsequent measurements on subjects in training set
n_share<-0.5
#Number of subjects
N<-100
#Number of extra subjects for test set
N_extra<-max(floor(0.5*N), 20)
#Obtain indexes for training set, subsequent measurements on training set and test set
indexes<-get_indexes(n, n_share, N, N_extra)



################################################################################
#Simulation study
################################################################################
#Short code to see whether everything is working
failed<-0
temp<-Create_LMM_data(true_betas, LMM_parameters, p, indexes)
temp<-Create_JMM_data(true_betas, JMM_parameters, p, indexes)
temp<-Create_JSM_data(true_betas, JSM_parameters, p, indexes)
LMM_results<-Fit_LMM(temp$data, indexes, failed)
JMM_results<-Fit_JMM(temp$data, indexes, failed)
JSM_results<-Fit_JSM(temp$data, indexes, failed)

# #Actual simulation takes days!!! Thus the code is commented out
# #The results have been saved and are attached. The analyses of the results can be found in the folder Simulation Results
# #The simulation study can best be run from the terminal, supplying extra arguments (which are read by commandArgs)
# #The simulation study is divided into 3 segments: Per segment data is simulated according to just 1 model (LMM, JMM and JSM). This is set via which_model
# args <- commandArgs(trailingOnly = TRUE)
# which_model<- as.integer(args[1])
# CV <- as.integer(args[2])
# n1 <- as.integer(args[3])
# n2 <- as.integer(args[4])
# CV<-250
# N_s<-c(10,25,50, 75, 100, 150, 250, 500, 750)
# increasing_N_comparison<-list()
# 
# for (i in 1:length(N_s)){
#   N<-N_s[i]
#   cat('\n At N', N, '\n')
#   cv<-1
#   while (cv<=CV){
#     failed<-0
#     cat('At N:', N, ', CV:', cv, '\n')
#     print(Sys.time())
#     N_extra<-floor(0.5*N)
#     indexes<-get_indexes(n, n_share, N, N_extra)
#     if (which_model==0){
#       temp<-Create_LMM_data(true_betas,model_0_parameters, p, indexes, cv)
#     }
#     else if (which_model==2){
#     temp<-Create_JMM_data(true_betas,model_2_parameters, p, indexes, cv)
#     }
#     else if (which_model==3){
#     temp<-Create_JSM_data(true_betas,model_3_parameters, p, indexes)
#     }
#     else {print('Can only simulate from models 2 and 3')}
#     LMM_fitted<-Fit_LMM(temp$data, indexes, failed)
#     cat('Model 0 done \n')
#     if (failed==1){
#       cat('Model failed: redoing cv \n')
#       next
#     }
#     JMM_fitted<-Fit_JMM(temp$data, indexes, failed)
#     cat('Model 2 done \n')
#     if (failed==1){
#       next
#     }
#     JSM_fitted<-Fit_JSM(temp$data, indexes, failed)
#     cat('Model 3 done \n')
#     if (failed==1){
#       next
#     }
#     increasing_N_comparison[[paste0('N', N)]][[paste0('CV',cv)]][['GOF']]<-
#       cbind(unlist(LMM_fitted[2:8]), unlist(JMM_fitted[2:8]), unlist(JSM_fitted[2:8]))
#     increasing_N_comparison[[paste0('N', N)]][[paste0('CV',cv)]][['Coefficients']]<-
#       rbind(LMM_fitted[[1]], JMM_fitted[[1]], JSM_fitted[[1]])
#     increasing_N_comparison[[paste0('N', N)]][[paste0('CV',cv)]][['Model 0']]<-LMM_fitted[-c(1:8)]
#     increasing_N_comparison[[paste0('N', N)]][[paste0('CV',cv)]][['Model 2']]<-JMM_fitted[-c(1:8)]
#     increasing_N_comparison[[paste0('N', N)]][[paste0('CV',cv)]][['Model 3']]<-JSM_fitted[-c(1:8)]
#     cv<-cv+1
#   }
# }
# # #Save results!!!!!!!!!
# # save(increasing_N_comparison, file=paste0('Simulation_2_final_', n1, '_', n2,'.RData'))


