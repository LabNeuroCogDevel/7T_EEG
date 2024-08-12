
mediationAnalysis <- function(DF, other_vars) {

  always_want <- c("Subject","age.x")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- DF[,c(always_want, other_vars)]
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  
#the effect of age on behavior (c)
model0Formula <- as.formula(sprintf("%s ~ age.x", other_vars[2]))
model.0 <- lm(model0Formula, data = mediationMatrix)
print(car::Anova(model.0))
print(summary(model.0))
#anova(model.00, model.0) # significant improvement in fit when you include age 


#the effect of age on eeg activity  (a)
modelMformula <- as.formula(sprintf("%s ~ age.x", other_vars[1]))
model.M <- lm(modelMformula, data = mediationMatrix)
print(car::Anova(model.M))
print(summary(model.M))

#the effect of eeg activity on behavior (b)
modelYformula <- as.formula(sprintf("%s ~ age.x +  %s", other_vars[2],  other_vars[1]))
model.Y <- lm(modelYformula, data = mediationMatrix)
print(car::Anova(model.Y))
print(summary(model.Y))


results <- mediation::mediate(model.M, model.Y, treat = 'age.x', mediator = other_vars[1], boot = FALSE, sims = 1000)
#pvalue <- results$d0.p

#plot(results, group.plots = T)


print(summary(results))
}

mediationAnalysis_MRS_crossSec <- function(DF, other_vars) {
  
  always_want <- c("luna","age")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- DF[,c(always_want, other_vars)]
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  
  #the effect of age on behavior (c)

  model0Formula <- as.formula(sprintf("%s ~ age + %s + %s + age*%s + age*%s + (1|luna)", other_vars[2], other_vars[3], other_vars[4], other_vars[3], other_vars[4]))
  model.0 <- lmer(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ age + %s+ %s+ (1|luna)", other_vars[1], other_vars[4], other_vars[3]))
  model.M <- lmer(modelMformula, data = mediationMatrix )
  print(car::Anova(model.M))
  print(summary(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~ age + %s+ %s + %s + %s*%s + %s*%s + (1|luna)", other_vars[1], other_vars[2], other_vars[3], other_vars[4], other_vars[3], other_vars[2],other_vars[4], other_vars[2]))
  model.Y <- lmer(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  
  
  results <- mediation::mediate(model.M, model.Y, treat = 'age', mediator = other_vars[1], boot = FALSE, sims = 1000)
  #pvalue <- results$d0.p
  
  #plot(results, group.plots = T)
  
  
  print(summary(results))
}


mediationAnalysis_longitud <- function(DF, other_vars) {
  
  always_want <- c("Subject","age.x")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- DF[,c(always_want, other_vars)]
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ age.x + (1|Subject)", other_vars[2]))
  model.0 <- lme4::lmer(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ age.x + (1|Subject)", other_vars[1]))
  model.M <- lmer(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  print(summary(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~ age.x + %s + (1|Subject)", other_vars[2],  other_vars[1]))
  model.Y <- lmer(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  
  
  results <- mediation::mediate(model.M, model.Y, treat = 'age.x', mediator = other_vars[1], boot = FALSE, sims = 1000)
  #pvalue <- results$d0.p
  
  #plot(results, group.plots = T)
  
  
  print(summary(results))
}

mediationAnalysis_invage <- function(DF, other_vars) {
  
  always_want <- c("Subject","invage")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- DF[,c(always_want, other_vars)]
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ invage", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ invage", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  print(summary(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~ invage +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  
  
  results <- mediation::mediate(model.M, model.Y, treat = 'invage', mediator = other_vars[1], boot = FALSE, sims = 1000)
  #pvalue <- results$d0.p
  
  #plot(results, group.plots = T)
  
  
  print(summary(results))
}

mediationAnalysis_invage_behaviorvsFreq <- function(DF, other_vars) {
  
  always_want <- c("Subject","inverseAge")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- DF[,c(always_want, other_vars)]
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ 1 + inverseAge", other_vars[2]))
  model.0 <- lm(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ 1 + inverseAge", other_vars[1]))
  model.M <- lm(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  print(summary(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~ 1 + inverseAge +  %s", other_vars[2],  other_vars[1]))
  model.Y <- lm(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  
  
  results <- mediation::mediate(model.M, model.Y, treat = 'inverseAge', mediator = other_vars[1], boot = FALSE, sims = 1000)
  #pvalue <- results$d0.p
  
  #plot(results, group.plots = T)
  
  
  print(summary(results))
}

mediationAnalysis_Epoch <- function(DF, other_vars) {
  
  always_want <- c("Subject","age")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- DF[,c(always_want, other_vars)]
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ age + %s + (1|Subject)", other_vars[2], other_vars[3]))
  model.0 <- lmer(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ age + %s+ (1|Subject)", other_vars[1], other_vars[3]))
  model.M <- lmer(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  print(summary(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~ age + %s + %s+ (1|Subject)", other_vars[2],  other_vars[1], other_vars[3]))
  model.Y <- lmer(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  
  
  results <- mediation::mediate(model.M, model.Y, treat = 'age', mediator = other_vars[1], boot = FALSE, sims = 1000)
  #pvalue <- results$d0.p
  
  #plot(results, group.plots = T)
  
  
  print(summary(results))
}
mediationAnalysis_Epoch_invage <- function(DF, other_vars) {
  
  always_want <- c("Subject","invage")
  
  #other_vars is defined in rmarkdown 
  mediationMatrix <- DF[,c(always_want, other_vars)]
  mediationMatrix <- mediationMatrix[complete.cases(mediationMatrix), ]
  
  
  #the effect of age on behavior (c)
  model0Formula <- as.formula(sprintf("%s ~ invage + %s + (1|Subject)", other_vars[2], other_vars[3]))
  model.0 <- lmer(model0Formula, data = mediationMatrix)
  print(car::Anova(model.0))
  print(summary(model.0))
  #anova(model.00, model.0) # significant improvement in fit when you include age 
  
  
  #the effect of age on eeg activity  (a)
  modelMformula <- as.formula(sprintf("%s ~ invage + %s+ (1|Subject)", other_vars[1], other_vars[3]))
  model.M <- lmer(modelMformula, data = mediationMatrix)
  print(car::Anova(model.M))
  print(summary(model.M))
  
  #the effect of eeg activity on behavior (b)
  modelYformula <- as.formula(sprintf("%s ~ invage + %s + %s+ (1|Subject)", other_vars[2],  other_vars[1], other_vars[3]))
  model.Y <- lmer(modelYformula, data = mediationMatrix)
  print(car::Anova(model.Y))
  print(summary(model.Y))
  
  
  results <- mediation::mediate(model.M, model.Y, treat = 'invage', mediator = other_vars[1], boot = FALSE, sims = 1000)
  #pvalue <- results$d0.p
  
  #plot(results, group.plots = T)
  
  
  print(summary(results))
}

