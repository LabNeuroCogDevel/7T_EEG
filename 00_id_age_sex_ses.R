#!/usr/bin/env Rscript
# use LNCDDB to get eeg 7T id_vist, age sex and session
# see Makefile (rule for id_age_sex_ses.csv)

# dbsession is what we recorded from the scan sheet
# but it's sometimes missing. so calcutate visit number by ranking visitdate (sort within id)
#
# 20230112 - WF init
pacman::p_load(dplyr)
pacman::p_load(tidyr)
dbcmd <- "selld8 l |sed 's/\t/,/g'|grep 'eeg.*BrainM'|cut -d, -f1-3,7"
ages <- read.csv(text=system(intern=T, dbcmd), header=F, col.names=c("ld8","age","sex","dbses"))
# calc visit by ordering
ages %>% separate(ld8,c("lunaid","vdate")) %>%
   arrange(lunaid,vdate) %>%
   group_by(lunaid) %>%
   mutate(visitno=rank(vdate)) %>%
   unite("ld8",lunaid,vdate,sep="_") %>%
   write.csv(stdout(), row.names=F, quote=F)
