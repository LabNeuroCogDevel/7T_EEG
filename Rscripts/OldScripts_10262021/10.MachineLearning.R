
## Going to try SVR from Caret Package

# GammaVars <- grep("Gamma|age|visitno",names(alldata),value=TRUE) 
# GammaVars <- GammaVars[!grepl("Trial.x|Trial.y|Trial_Power", GammaVars)]
# 
# allGamma <- alldata[,GammaVars]
# allGamma <- allGamma[complete.cases(allGamma),]
# allGamma <- allGamma[!duplicated(allGamma),]
# 
# Gamma_baselineclean <- allGamma[allGamma$visitno < 2,]
# Gamma_baselineclean <- Gamma_baselineclean[1:20]
# 
# train_control <- trainControl(method = "repeatedcv", number = 15, repeats = 5)
# 
# # linear Kernel
# svm1 <- train(age~., data = Gamma_baselineclean, method = "svmLinear", trControl = train_control, preProcess = c("center", "scale"), fitBest = FALSE, returnData = TRUE)
# svm1
# 
# res1 <- as_tibble(svm1$results[which.min(svm1$results[,2]),])
# finalMod <- svm1$finalModel
# 
# 
# bothModels <- list(linear = svm1)
# agePredicted <- extractPrediction(bothModels, testX = Gamma_baselineclean[20])
# 
# ggplot(data = agePredicted, aes(x = obs, y = pred)) + geom_point() + xlim(10,32) + ylim(10, 32) + geom_smooth(method = "lm") + ggtitle("Obs. vs Pred. Ages: Gamma Band")
# plot(varImp(svm1, useModel = FALSE, nonpara = FALSE))
# 
# 
# # non-linear kernel: radial
# svm2 <- train(age~., data = Gamma_baselineclean, method = "svmRadial", trControl = train_control, preProcess = c("center", "scale"), tuneLength = 10)
# # Print the best tuning parameter sigma and C that maximizes model accuracy 
# svm2$bestTune
# svm2
# 
# res2 <- as_tibble(svm2$results[1,])
# res2
# 
# linearPred <- predict(svm2)
# 
# plot(Gamma_baselineclean$age, linearPred) + abline(lm(linearPred ~ Gamma_baselineclean$age), col= "red")
# plot(varImp(svm2))
# 
# # non-linear kernel: polynomial
# svm3 <- train(age~., data = Gamma_baselineclean, method = "svmPoly", trControl = train_control, preProcess = c("center", "scale"), tuneLength = 4)
# # Print the best tuning parameter sigma and C that maximizes model accuracy 
# svm3$bestTune
# svm3
# 
# res3 <- as_tibble(svm3$results[18,])
# res3
# 
# linearPred <- predict(svm3)
# 
# plot(Gamma_baselineclean$age, linearPred) + abline(lm(linearPred ~ Gamma_baselineclean$age), col= "red")
# plot(varImp(svm3))
# 
# #display all results
# 
# df <- tibble(Model = c('SVM Linear', 'SVM Radial','SVM Poly'), Rsquared = c(res1$Rsquared, res2$Rsquared, res3$Rsquared))
# df %>% arrange(Rsquared)


linearSVRmodel <- function(baselineclean) {
  
  baselineclean_impute <- kNN(baselineclean, variable = colnames(baselineclean), k = 10, impNA = T, imp_var = F, imp_suffix = F)
  
  train_control <- trainControl(method = "repeatedcv", number = 15, repeats = 5)
  
  # linear Kernel
  # preprocess = c("center", "scale") normalized the variables to make their scales comparable 
  svm1 <- train(age~., data = baselineclean_impute, method = "svmLinear", trControl = train_control, preProcess = c("center", "scale"), na.action = na.pass)
  
  
  return(svm1)

  
  
}

linearSVRmodel_PE <- function(baselineclean) {
  
  train_control <- trainControl(method = "repeatedcv", number = 15, repeats = 5)
  
  # linear Kernel
  svm1 <- train(absPositionError~., data = baselineclean, method = "svmLinear", trControl = train_control, preProcess = c("center", "scale"))
  
  
  return(svm1)
  
  
}


linearSVRmodel_Lat <- function(baselineclean) {
  
  train_control <- trainControl(method = "repeatedcv", number = 15, repeats = 5)
  
  # linear Kernel
  svm1 <- train(mgsLatency~., data = baselineclean, method = "svmLinear", trControl = train_control, preProcess = c("center", "scale"))
  
  
  return(svm1)
  
  
}

