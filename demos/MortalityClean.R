rm(list=ls())
library(dplyr)

# do some age assignments with my bare hands
DFAge <- data.frame(age_group_id=c(2:20, 30:33, 44:45)) %>%
    mutate(age_mean=c(7/356/2, 21/356/2, 328/356/2, 3, seq(7.5, 107.5, 5))) %>%
    mutate(age_time=c(7/356, 21/356, 328/356, 4, rep(5, 21)))

# read and filter the adult data
DFAdult <- "~/Downloads/ProbabilityOfDeath_estimates.csv" %>%
    read.csv(stringsAsFactors=FALSE) %>%
    filter(sex=="Both" & year==2016) %>% 
    filter(age_group_id %in% c(6:20, 30:33, 44:45)) %>%
    unique %>% rename(px=mean) %>% 
    select(location, year, age_group, age_group_id, sex, px)

# load in the child data and ombine with adult and the age DFs
DFDeath <-"~/Downloads/5q0Results_estimates.csv" %>%
    read.csv(stringsAsFactors=FALSE) %>%
    filter(sex=="Both" & year==2016 & age_group_id %in% 2:5) %>%
    unique %>% rename(px=mean) %>% 
    select(location, year, age_group, age_group_id, sex, px) %>%
    bind_rows(DFAdult) %>% left_join(DFAge, by="age_group_id") %>%
    arrange(location, age_group_id) %>%
    mutate(qx=1-px, hx=1-((qx)^(1/age_time))) %>%
    group_by(location) %>% mutate(Sx=cumprod(qx), Fx=1-Sx) %>% as.data.frame

save(DFDeath, file="../data/DFDeath.Rda")
