


library(tidyverse)


svr_wrapper_parse_Jennifer <- function(svm_warpper_object){
  cors<-svm_warpper_object$wrapperout %>% dplyr::group_by(y) %>% dplyr::summarise(cvcor=cor(cvpred,testyval),cvrsq=rsquared_function(testyval,error))
  print(cors)
  return(cors)
}


SVRfunction <- function(baselineclean){
  
  source("DIRECTORY WHERE YOU SAVE THE 0X_SVRFUNC.R CODE")
  
  
  #xcols<-grep("Var",names(baselineclean),value=TRUE)
  
  xcols<-names(baselineclean)
  # so here I am grabbing all variables in my large dataframe called baselineclean that contain alpha, theta, beta, gamma, in the names
  # these are the variables you want to feed into the SVR 
  xcols <- xcols[grep("Alpha|Theta|Beta|Gamma",xcols)] 
  
  
  # here i am saying remove all variables from the large dataframe (baselineclean), that include inverse, age, group
  # you want to remove the variabel you are trying to predict 
  xcols <- xcols[!grepl("inverse|age|Group",xcols)]
  
  #foldavar ### column name that you want to split data for cross-validation (probably use "subject"_
  #df = dataframe with all data
  #ys = column name or names that you want to predict
  #xcols = column names that you want to build the model from (variables doign the predicting) #this will likely include all of your measures
  
  #tune (hyperparameter tuning) for SVR set to FALSE for now
  #tunefolds number of folds for tuning (no action if tune is set to FALSE) leave at 0
  #outerfold cores how many folds in parallel (depends on your computer), if running locally probably set to 1
  #tune cores internal cores used for various feature selection steps (typically set to the same as outerfold cores
  
  ###unifeatselect (should univariate feature selection be used? ##if so how many features (nfeatures)? can set to FALSE right now
  ###PCA should PCA be used ###if so what perecent of variance (numbers less than one), all of the components = 1, a specific number of components (numbers greater than 1) ### I would recommend running it with PCA=TRUE and propvarretain=1 for now
  ####
  
  
  ###run model#####
  
  SVRage <- svm_wrapper(foldvar="Subject", df=baselineclean, ys="age", xcols=xcols, tune=FALSE, tunefolds=0, outerfoldcores=1, tunecores=1, unifeatselect=FALSE, nfeatures=0, PCA=TRUE, propvarretain=1)
  
  #####print out cross-validation info
  
  foldparse <- svr_wrapper_parse_Shane(SVRage)
  print(foldparse)
  
  plot(SVRage$wrapperout$testyval, SVRage$wrapperout$cvpred, main = "Correlation between Predicted Age and Actual Age", ylab = "Predicted Age", xlab = "Actual Age")
  
  varWeights <- SVRage$weightmean                    
  absVarWeights <- abs(varWeights[1:60])
  
  absValueDF <- gather(absVarWeights, key = "measure", value = "absWeights")
  
  
  ggplot(data = absValueDF, aes(x = reorder(measures, -absWeights), absWeights)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
}











