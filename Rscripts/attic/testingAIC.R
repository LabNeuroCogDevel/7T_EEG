

# Testing AIC

delay_alldata_gamma <- delay_alldata %>% filter(delay_alldata$Band == "Gamma")
behaviorAge <- merge(Behavior, agefile, by = c("Subject", "age"))


#generate quadratic variables
delay_alldata_gamma$quadAge <- delay_alldata_gamma$age^2

delay_alldata_gamma$quadAge <- (delay_alldata_gamma$age - mean(delay_alldata_gamma$age))^2

behaviorAge$quadAge <- (behaviorAge$age - mean(behaviorAge$age))^2


linearModel <- lm(data = behaviorAge[], absBestError ~ age) 

inverseModel <- lm(data = behaviorAge[], absBestError ~ inverseAge) 

quadModel <- lm(data = behaviorAge, absBestError ~ age + quadAge)


AIC(linearModel,inverseModel, quadModel) #check to see if age vs inverse age vs quad age, which is the better model


#plot all these types of graphs
#linear
lunaize(ggplot(data = behaviorAge[], aes(x = age, y = absBestError)) + geom_point() + stat_smooth(method = "lm") + ggtitle("Linear") + xlab("Age") + ylab("Best Saccade")) +theme(plot.title = element_text(hjust = 0.5))

#inverse
lunaize(ggplot(data = behaviorAge[], aes(x = age, y = absBestError)) + geom_point() + stat_smooth(method = "lm",formula = 'y~I(1/x)') + ggtitle("Inverse") + xlab("Age") + ylab("Best Saccade")) +theme(plot.title = element_text(hjust = 0.5))

#quadratic 
lunaize(ggplot(data = behaviorAge[], aes(x = age, y = absBestError)) + geom_point() + stat_smooth(method = "lm",formula = 'y ~ x + I(x^2)') + ggtitle("Quadratic") + xlab("Age") + ylab("Best Saccade")) +theme(plot.title = element_text(hjust = 0.5))
