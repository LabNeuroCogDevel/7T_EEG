



library(tidyverse)


svr_wrapper_parse_Shane <- function(svm_warpper_object){
  cors<-svm_warpper_object$wrapperout %>% dplyr::group_by(y) %>% dplyr::summarise(cvcor=cor(cvpred,testyval),cvrsq=rsquared_function(testyval,error))
  print(cors)
  return(cors)
}


SVRfunctionAge <- function(baselineclean){
  
  
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/0x_svrfuncs.R")
  
  
  #xcols<-grep("Var",names(baselineclean),value=TRUE)
  
  xcols<-names(baselineclean)
  xcols <- xcols[grep("Alpha|Theta|Beta|Gamma|mgsLatency|absPositionError",xcols)]
  xcols <- xcols[!grepl("inverse|age|Group|Trial",xcols)]
  
  
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
  knitr::kable(foldparse, format = "markdown")

  
  varWeights <- SVRage$weightmean   
  cols <- ncol(varWeights)
  absVarWeights <- abs(varWeights[1:cols-1])
  
  absValueDF <- gather(absVarWeights, key = "measure", value = "absWeights")
  
  
  print(ggplot(data = absValueDF, aes(x = reorder(measure, -absWeights), absWeights)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
 
 

}

SVRfunctionInverseAge <- function(baselineclean){
  
  
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/0x_svrfuncs.R")
  
  
  #xcols<-grep("Var",names(baselineclean),value=TRUE)
  
  xcols<-names(baselineclean)
  xcols <- xcols[grep("Alpha|Theta|Beta|Gamma|mgsLatency|absPositionError",xcols)]
  xcols <- xcols[!grepl("inverse|age|Group|Trial",xcols)]
  
  
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
  
  SVRage <- svm_wrapper(foldvar="Subject", df=baselineclean, ys="inverseAge", xcols=xcols, tune=FALSE, tunefolds=0, outerfoldcores=1, tunecores=1, unifeatselect=FALSE, nfeatures=0, PCA=TRUE, propvarretain=1)
  
  #####print out cross-validation info
  
  foldparse <- svr_wrapper_parse_Shane(SVRage)
  print(foldparse)
  
  #plot(SVRage$wrapperout$testyval, SVRage$wrapperout$cvpred, main = "Correlation between Predicted Inverse Age and Actual Inverse Age", ylab = "Predicted Inverse Age", xlab = "Actual Inverse Age")
  
  varWeights <- SVRage$weightmean   
  cols <- ncol(varWeights)
  absVarWeights <- abs(varWeights[1:cols-1])
  
  absValueDF <- gather(absVarWeights, key = "measure", value = "absWeights")
  
  
  print(ggplot(data = absValueDF, aes(x = reorder(measure, -absWeights), absWeights)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  
  
}

SVRfunctionBehaviorPosition <- function(baselineclean){
  
  
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/0x_svrfuncs.R")
  
  
  #xcols<-grep("Var",names(baselineclean),value=TRUE)
  
  xcols<-names(baselineclean)
  xcols <- xcols[grep("Alpha|Theta|Beta|Gamma",xcols)]
  xcols <- xcols[!grepl("inverse|Group|Trial|log_mgsLatency|age|.1",xcols)]
  
  
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
  
  SVRage <- svm_wrapper(foldvar="Subject", df=baselineclean, ys="log_absPositionError", xcols=xcols, tune=FALSE, tunefolds=0, outerfoldcores=1, tunecores=1, unifeatselect=FALSE, nfeatures=0, PCA=TRUE, propvarretain=1)
  
  #####print out cross-validation info
  
  foldparse <- svr_wrapper_parse_Shane(SVRage)
  print(foldparse)
  
  plot(SVRage$wrapperout$testyval, SVRage$wrapperout$cvpred, main = "Predicted Position Error vs Actual Position Error", ylab = "Predicted Position Error", xlab = "Actual Position Error")
  
  varWeights <- SVRage$weightmean   
  cols <- ncol(varWeights)
  absVarWeights <- abs(varWeights[1:cols-1])
  
  absValueDF <- gather(absVarWeights, key = "measure", value = "absWeights")
  
  
  print(ggplot(data = absValueDF, aes(x = reorder(measure, -absWeights), absWeights)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  
  
}

SVRfunctionBehaviorLatency <- function(baselineclean){
  
  
  source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/0x_svrfuncs.R")
  
  
  #xcols<-grep("Var",names(baselineclean),value=TRUE)
  
  xcols<-names(baselineclean)
  xcols <- xcols[grep("Alpha|Theta|Beta|Gamma",xcols)]
  xcols <- xcols[!grepl("inverse|Group|Trial|age|log_absPOsitionError|.1",xcols)]
  
  
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
  
  SVRage <- svm_wrapper(foldvar="Subject", df=baselineclean, ys="log_mgsLatency", xcols=xcols, tune=FALSE, tunefolds=0, outerfoldcores=1, tunecores=1, unifeatselect=FALSE, nfeatures=0, PCA=TRUE, propvarretain=1)
  
  #####print out cross-validation info
  
  foldparse <- svr_wrapper_parse_Shane(SVRage)
  print(foldparse)
  
  plot(SVRage$wrapperout$testyval, SVRage$wrapperout$cvpred, main = "Predicted mgsLatency vs Actual mgsLatency", ylab = "Predicted mgsLatency", xlab = "Actual mgsLatency")
  
  varWeights <- SVRage$weightmean   
  cols <- ncol(varWeights)
  absVarWeights <- abs(varWeights[1:cols-1])
  
  absValueDF <- gather(absVarWeights, key = "measure", value = "absWeights")
  
  
  print(ggplot(data = absValueDF, aes(x = reorder(measure, -absWeights), absWeights)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  
  
}


combineSVRs <- function() {
  model <- lm(GammaSVRage$plot_env$SVRage$wrapperout$testyval ~ GammaSVRage$plot_env$SVRage$wrapperout$cvpred)
  modsum <- print(summary(model))
  r2 = modsum$adj.r.squared
  my.p = modsum$coefficients[2,4]
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  
  
 betamodel <- lm(BetaSVRage$plot_env$SVRage$wrapperout$testyval ~ BetaSVRage$plot_env$SVRage$wrapperout$cvpred)
  betamodsum <- print(summary(betamodel))
  betar2 = betamodsum$adj.r.squared
  betamy.p = betamodsum$coefficients[2,4]
  betarp = vector('expression',2)
  betarp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                     list(MYVALUE = format(betar2,dig=3)))[2]
  betarp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(betamy.p, digits = 2)))[2]
  
  

plot(GammaSVRage$plot_env$SVRage$wrapperout$testyval, GammaSVRage$plot_env$SVRage$wrapperout$cvpred, main = "Correlation between Predicted Age and Actual Age", ylab = "Predicted Age", xlab = "Actual Age", xlim = c(10, 35), ylim = c(5, 40), pch = 1)
  abline(lm(GammaSVRage$plot_env$SVRage$wrapperout$testyval ~ GammaSVRage$plot_env$SVRage$wrapperout$cvpred), col = "red")
  legend('topright', legend = rp, bty = 'n')
    par(new=TRUE)
 plot(BetaSVRage$plot_env$SVRage$wrapperout$testyval, BetaSVRage$plot_env$SVRage$wrapperout$cvpred, main = "Correlation between Predicted Age and Actual Age", ylab = "Predicted Age", xlab = "Actual Age", xlim = c(10, 35), ylim = c(5, 40), pch = 8)
  abline(lm(BetaSVRage$plot_env$SVRage$wrapperout$testyval ~ BetaSVRage$plot_env$SVRage$wrapperout$cvpred), col = "blue")
    legend('topleft', legend = betarp, bty = 'n')
  
  
  
  
}




