library(tvem)

foof <- read.csv('H:/Projects/7TBrainMech/scripts/eeg/Shane/fooof/Results/fooofMRSbehavior_20230313.csv')

foof$luna <- as.factor(foof$luna)

# Run separately for trial type given three-way interaction (above)
foof2 <- dplyr::select(foof[foof$Region!="MPFC",], age.x, age, luna, absBestError, visitno,Exponent,Condition,Region)
# foof2 <- foof2[complete.cases(foof2),]
foof2$age_new <- round(foof2$age.x)

# Finn's if you want to average across conditions (keeping region separate). 
#The one I'm running below is putting both conditions (so 2 datapts per sub)
# foof2 <- read.csv('/Users/ashley/Downloads/fooofMRSbehavior_20230313.csv') %>% 
#   filter(Region %in% c('LDLPFC','RDLPFC'))%>% 
#   group_by(luna, visitno, age,Region) %>% 
#   summarize(absBestError = mean(absBestError, na.rm=T),
#             Exponent = mean(Exponent, na.rm=T)) %>%
#   ungroup 

#Model 1 for R DLPFC 
model1 <- tvem(data=foof2[foof2$Region=="RDLPFC",], 
               formula=Exponent~absBestError,
               id=luna, 
               time=age_new,
               num_knots=15)

ages <- model1$time_grid
coefs <- model1$grid_fitted_coefficients$absBestError$estimate
upper <- model1$grid_fitted_coefficients$absBestError$upper
lower <- model1$grid_fitted_coefficients$absBestError$lower

test_df <- cbind(ages,coefs,upper,lower)
test_df <- as.data.frame(test_df)
test_df$Region <- "RDLPFC"

#Model 2 for L DLPFC 
model2 <- tvem(data=foof2[foof2$Region=="LDLPFC",], 
               formula=Exponent~absBestError,
               id=luna, 
               time=age_new,
               num_knots=15)

ages <- model2$time_grid
coefs <- model2$grid_fitted_coefficients$absBestError$estimate
upper <- model2$grid_fitted_coefficients$absBestError$upper
lower <- model2$grid_fitted_coefficients$absBestError$lower

test_df2 <- cbind(ages,coefs,upper,lower)
test_df2 <- as.data.frame(test_df2)
test_df2$Region <- "LDLPFC"

mega_df <- rbind(test_df, test_df2)
mega_df$Region= as.factor(mega_df$Region)

#plot them both together 
age_perf2 <- ggplot(mega_df, aes(x=ages, y=coefs, color=Region))+
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=Region, color=NULL), alpha=.15)+
  geom_line()+ geom_hline(yintercept=0)+
  theme(aspect.ratio=1, plot.title = element_text(size=20,color="black"), 
        axis.text=element_text(size=18,color="black"), axis.title=element_text(size=18,color="black"),
        panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  labs(y="Exponent by Abs Error", x="Age")
print(age_perf2)
