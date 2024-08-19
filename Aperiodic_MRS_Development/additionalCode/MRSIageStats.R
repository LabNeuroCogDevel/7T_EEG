
#!/usr/bin/env Rscript

MSRIageStats <- function () {
  
# Libraries ----

library(LNCDR)
library(data.table)
library(dplyr)
library(factoextra)
library(ggplot2)
library(e1071)
library(caret)
attach(mtcars)
library(grid)
library(gridExtra)
library(plotrix)
library(mgcv)
library(readxl)
library(lme4)
library(lubridate)
library(checkmate)
library(lmerTest)
library(tidyr)
library(jtools)

# Load Dataframe ----

MRSlong <- read.csv("/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/allSubjectsDLPFCMRSMeasures_20230613.csv")

# MRSI stats ----

## GLU VS AGE ----

lunaize(ggplot(data = MRSlong, 
               aes(x = age, y = Glu_gamadj, by = luna, shape = Region, linetype = Region)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
          geom_smooth(aes(group = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.6,size=1) + 
          scale_color_manual(values=c("gold3", "blue4"))) + xlab("Age") +ylab("Glutamate")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

#ggsave("dlpfc_glu_age.png", height=6, width = 6, dpi = 300)
gam.model <-  gamm(Glu_gamadj ~ s(age, k = 3) + Region, data = MRSlong, random=list(luna=~1))
summary(gam.model$gam)

#test <-  lmer(Glu_gamadj ~ inverseAge * Region + (1|luna), data = MRSlong)
#summary(test)


MRSlong$oRegion <- ordered(MRSlong$Region, levels = c("LDLPFC","RDLPFC")) # LDLFPC will be the reference group
model_formula <- as.formula("Glu_gamadj ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = MRSlong )
summary(model$gam)

#see which hemisphere has a stronger age association
gam.model <-  gamm(Glu_gamadj ~ s(age, k = 3), data = MRSlong %>% filter(Region == 'RDLPFC'), random=list(luna=~1))
summary(gam.model$gam)


## GABA VS AGE ----
lunaize(ggplot(data = MRSlong,
               aes(x = age, y = GABA_gamadj, by =luna, linetype = Region)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(alpha=.5) + 
          geom_smooth(aes(group = Region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.6,size=1)) + 
  xlab("Age") +ylab("GABA") + theme(text = element_text(size = 30)) + theme(legend.position='none')


ggsave("dlpfc_gaba_age.png", height=6, width = 6, dpi = 300)

gam.model <-  gamm(GABA_gamadj ~ s(age, k = 3) + Region, data = MRSlong, random=list(luna=~1))
summary(gam.model$gam)

MRSlong$oRegion <- ordered(MRSlong$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC will be the reference group
model_formula <- as.formula("gabaRes.residuals ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = MRSlong )
summary(model$gam)


## RATIO VS AGE ----
lunaize(ggplot(data = MRSlong, 
               aes(x = age, y = Ratio_gamadj, by =luna, linetype = Region)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(aes(shape=Region),alpha=.5) + 
          geom_smooth(aes(group = Region, alpha = 0.1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.6,size=1)) + 
  xlab("Age") +ylab("Glu/GABA Ratio")+ theme(text = element_text(size = 30))+ theme(legend.position='none')


gam.model <-  gamm(Ratio_gamadj ~ s(age, k = 3) + Region, data = MRSlong, random=list(luna=~1))
summary(gam.model$gam)

MRSlong$oRegion <- ordered(MRSlong$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC will be the reference group
model_formula <- as.formula("Ratio_gamadj ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = MRSlong%>% filter(age < 26) )
summary(model$gam)



## GLU GABA Imbalance VS AGE ----
lunaize(ggplot(data = MRSlong, 
               aes(x = age, y = GluGABAimbalanceABS, by = luna, shape = Region, linetype = Region)) + 
          geom_line(aes(group=interaction(luna,Region)), alpha = 0.2) + geom_point(aes(shape= Region),alpha=.5) + 
          geom_smooth(aes(group = Region, alpha = 0.1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T),alpha=.6,size=1)) + 
  xlab("Age") +ylab("Glu GABA Imbalance")+ theme(text = element_text(size = 30))+ theme(legend.position='none')

gam.model <-  gamm(GluGABAimbalanceABS ~ s(age, k = 3) +Region , data = MRSlong, random=list(luna=~1))
summary(gam.model$gam)

MRSlong$oRegion <- ordered(MRSlong$Region, levels = c("LDLPFC","RDLPFC")) # eyes open will be the reference group
model_formula <- as.formula("GluGABAimbalanceABS ~ oRegion + s(age, k = 4, fx = T) + s(age, by = oRegion, k = 4, fx = T)")
# Note we keep fx = T for reliable p-values.
model <- gamm(model_formula,
              random = list(luna=~1),
              data = MRSlong )
summary(model$gam)



## MRSI Loop ----

### MRSI vs Age ----
mrsiVars <- c('Ratio_gamadj','GluGABAimbalanceABS','Glu_gamadj','GABA_gamadj')

dim(MRSlong$Ratio_gamadj) <- NULL
dim(MRSlong$Glu_gamadj) <- NULL
dim(MRSlong$GABA_gamadj) <- NULL

output <- c()  
for (mrsiVar in mrsiVars) {
  
  model <- paste0(mrsiVar, ' ~ ', 's(age, k = 3) + Region')
  model.out <- (mgcv::gamm(as.formula(model), random=list(luna=~1), data = MRSlong ))
  F <- summary(model.out$gam)$s.table[1,3]
  p <- summary(model.out$gam)$s.table[1,4]
  RegionF <- summary(model.out$gam)$pTerms.table[1,2]
  Regionp <- summary(model.out$gam)$pTerms.table[1,3]
  
  output <- rbind(output, data.frame(mrsiVar, F, p, RegionF, Regionp))
}
output  

# region interaction
MRSlong$oRegion <- ordered(MRSlong$Region, levels = c("LDLPFC","RDLPFC")) # LDLPFC reference group

output <- c()  
for (mrsiVar in mrsiVars) {
  
  model_formula <- as.formula(paste0(mrsiVar, ' ~ ', ' oRegion + s(age, k = 3, fx = T) + s(age, by = oRegion, k = 3, fx = T)')) 
  
  model.out <- (mgcv::gamm(model_formula,
                           random = list(luna=~1),
                           data = MRSlong))
  
  RegionF <- summary(model.out$gam)$s.table[2,3]
  Regionp <- summary(model.out$gam)$s.table[2,4]
  
  output <- rbind(output, data.frame(mrsiVar, RegionF, Regionp))
}
output  


}
