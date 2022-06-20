source('Read_data.R')
source('Fit_models_LUMC.R')

################################################################################
################################################################################
#Fitting the LUMC data
################################################################################
################################################################################

################################################################################
#Parameters
################################################################################
comparison_all<-list()
data<-create_data()
data$y<-data$severity_score_pseudo

#Numer of subjects to include into the test set.
#Needed only when computing the MSE on the test set
#Default value=0: Uses entire data-set to train models
Test_set_subjects<-0

#Number of cross-validations to perform
#Default value=1 (in combination with Test_set_subjects=0): Entire dataset is used for training,
#No need for cross-validation
CV<-1


################################################################################
#Fitting the models
################################################################################
for (cv in c(1:CV)){
  cat('\n At loop:', cv, '\n')
  #In total there are 3 cytokines
  for (i in c(1:3)){
    #Determine indexes of train and test sets
    indexes<-test_indexes(data, Test_set_subjects)
    #Select appropriate cytokine
    data$exo<-data[[i+6]]
    #Fit models
    LMM<-Fit_LMM(data, indexes)
    JMM<-Fit_JMM(data, indexes)
    JSM<-Fit_JSM(data, 0, indexes)
    # #The JSM can be fit with lagg up to order 5. See below JSM's with lagg of order 2 & 5
    # JSM_lagg_2<-Fit_JSM(data, 2, indexes)
    # JSM_lagg_5<-Fit_JSM(data, 5, indexes)
    
    #Save results
    comparison_all[[paste0('Cytokine',i)]][[paste0('CV', cv)]][['LMM']]<-LMM
    comparison_all[[paste0('Cytokine',i)]][[paste0('CV', cv)]][['JMM']]<-JMM
    comparison_all[[paste0('Cytokine',i)]][[paste0('CV', cv)]][['JSM']]<-JSM
    # comparison_all[[paste0('Cytokine',i)]][[paste0('CV', cv)]][['JSM_lagg_2']]<-JSM_lagg_2
    # comparison_all[[paste0('Cytokine',i)]][[paste0('CV', cv)]][['JSM_lagg_5']]<-JSM_lagg_5
    comparison_all[[i]][[paste0('CV', cv)]]$results<-rbind(unlist(LMM[c(2:8)]), unlist(JMM[c(2:8)]),
                                       unlist(JSM[c(2:8)]))
  }
}

# save(comparison_all, file="comparison_all.RData")

