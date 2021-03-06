rm(list=ls())
library(dplyr)

# do some age assignments with my bare hands
DFAge <- data.frame(age_group_id=c(2:20, 30:33, 44:45)) %>%
    mutate(age_mean=c(7/356/2, 21/356/2, 328/356/2, 3, seq(7.5, 107.5, 5))) %>%
    mutate(age_time=c(7/356, 21/356, 328/356, 4, rep(5, 21))) %>%
    mutate(age_end=c(7/356, 28/356, 1, seq(5, 110, 5))) %>%
    mutate(age_start=c(0, 7/356, 28/356, 1, seq(5, 105, 5)))

# read and filter the adult data
DFAdult <- "~/Documents/DemographicSimulation/data/ProbabilityOfDeath_estimates.csv" %>%
    read.csv(stringsAsFactors=FALSE) %>%
    filter(sex=="Both") %>% 
    filter(age_group_id %in% c(6:20, 30:33, 44:45)) %>%
    unique %>% rename(qx=mean) %>% 
    select(location, year, age_group, age_group_id, sex, qx)

# load in the child data and ombine with adult and the age DFs
DFDeath <-"~/Documents/DemographicSimulation/data/5q0Results_estimates.csv" %>%
    read.csv(stringsAsFactors=FALSE) %>%
    filter(sex=="Both")  %>%
    filter(age_group_id %in% 2:5) %>%
    unique %>% rename(qx=mean) %>% 
    select(location, year, age_group, age_group_id, sex, qx) %>%
    bind_rows(DFAdult) %>% left_join(DFAge, by="age_group_id") %>%
    arrange(location, year, sex, age_group_id) %>%
    mutate(px=1-qx, hx=1-((px)^(1/age_time))) %>%
    group_by(location, year, sex) %>% mutate(Sx=cumprod(px), Fx=1-Sx) %>% 
    as.data.frame

save(DFDeath, file="~/Documents/DemographicSimulation/data/DFDeath.RData")
