# animate interactions

library(viridis)
library(gganimate)
library(magick)

m <- fooofMRSbehavior %>% 
  filter(Region %in% c('LDLPFC','RDLPFC'))%>% 
  group_by(luna, visitno, age) %>% 
  summarize(mgsLatency = mean(mgsLatency, na.rm=T),
            Offset = mean(Offset, na.rm=T)) %>%
  ungroup 

#### set parameters
ageRange <- seq(10,30,length.out=100) # 100 steps from min to max age
sigma <- 2 # SD of sliding window


#### make moving interaction plot

# set up data frame
movdata <- merge(m, 
                 data.frame(thisage = ageRange)) %>%
  mutate(w = exp(-(age - thisage)^2/(2*(sigma^2)))) %>% 
  group_by(thisage) %>% mutate(w = w/max(w)) %>% mutate(w = pmax(w, .2*max(w))) %>% ungroup()


intPlot <- lunaize(ggplot(data = movdata, aes(x = mgsLatency, y = Offset, color = age)) +
                     geom_smooth(method='lm', aes(weight=w, color=thisage)) +
                     scale_color_viridis(option='H') +
                     geom_point(aes(alpha = w, size=w)) +
                     transition_time(thisage) +
                     ease_aes('linear')) + theme(legend.position='none')



int_gif <- animate(intPlot, width = 480, height = 480)
anim_save("interaction.gif", animation = int_gif)

# make age vs var1 plot
agePlot <- lunaize(ggplot(data = movdata, aes(x = age, y = Offset, color = age)) +
                     stat_smooth(method='lm', formula = y ~ I(1/x)) +
                     geom_vline(aes(xintercept = thisage)) +
                     scale_color_viridis(option='H') +
                     geom_point(aes(alpha = w, size=w)) +
                     transition_time(thisage) +
                     labs(x = 'Age', y = 'Offset') +
                     ease_aes('linear')) + theme(legend.position='none')

age_gif <- animate(agePlot, width = 480, height = 480)
anim_save("age.gif", animation = age_gif)

# stitch together
a_mgif <- image_read('interaction.gif')
b_mgif <- image_read('age.gif')
new_gif <- image_append(c(a_mgif[1], b_mgif[1]))

for(i in 2:100){
  combined <- image_append(c(a_mgif[i], b_mgif[i]))
  new_gif <- c(new_gif, combined)
}
new_gif
anim_save("combined.gif", animation = new_gif)