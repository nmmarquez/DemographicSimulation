# Recording Mortality

# we want to record each individual death that occurs with an upper and
# lower bound at age of death but we also want to move populations forward as
# often as possible. It seems like having an object with a mortality curve,
# a population that has an indicator for age, and a template to record the
# information when someone dies.
rm(list=ls())
library(dplyr)
library(parallel)
library(ggplot2)

load("~/Documents/DemographicSimulation/data/USProbDeath.Rda")

ages <- c(7/365/2, mean(c(28/365, 7/365)), mean(c(1 ,28/365)), 3, 
          seq(7.5, 107.5, 5))

DFAge <- data.frame(age=ages, age_group_id=c(2:20, 30:33, 44, 45))

DFAdultDeath <- "~/Downloads/ProbabilityOfDeath_estimates.csv" %>%
    read.csv(stringsAsFactors=FALSE) %>%
    filter(year == 2016 & sex == "Both") %>%
    filter(age_group_id %in% 6:20 | age_group_id %in% c(30:33, 44, 45)) %>%
    unique %>% select(location, age_group, age_group_id, year, sex, mean)

ProbDeath <- "~/Downloads/5q0Results_estimates.csv" %>%
    read.csv(stringsAsFactors=FALSE) %>%
    filter(year==2016 & age_group_id %in% 2:5) %>%
    unique %>% 
    select(location, age_group, age_group_id, year, sex, mean) %>%
    bind_rows(DFAdultDeath) %>% left_join(DFAge) %>%
    mutate(mean=ifelse(age_group_id == 2, mean * 365/7, mean)) %>%
    mutate(mean=ifelse(age_group_id == 3, mean * 365/21, mean)) %>%
    mutate(mean=ifelse(age_group_id == 4, mean * 365/337, mean)) %>%
    rename(px=mean)

save(ProbDeath, file="~/Documents/DemographicSimulation/data/ProbDeath.Rda")

ggplot(ProbDeath, aes(x=age, y=px, color=location, group=location)) + 
    geom_line() + 
    coord_trans(y="log")


mort_ledger <- function(start_pop, mort_func, time=0){
    ledgerDF <- data.frame(
        age_start=vector("numeric"),
        age_end=vector("numeric"),
        cycle_death=vector("numeric"))
    ledger <- list(pop=start_pop, mortality=mort_func, ledger=ledgerDF, time=0)
    class(ledger) <- "MortalityLedger"
    ledger
}

mort_unif_func <- function(age, x=.19){
    # Given an age return the yearly rate of mortality
    rep(x, length(age))
}

mortUnif <- function(x=.19){
    proto=list(
        mortFunc = function(x) rep(.19, length(x)),
        logSurvInt = function(a, b) log(.81) * (b-a),
        periodMort = function(a, b) 1 - exp(log(.81) * (b-a))
    )
    class(proto) <- "Mortality"
    proto
}

progress_pop <- function(ledger, stepsize=1){
    if(class(ledger$mortality) == "Mortality"){
        pdeath <- ledger$mortality$periodMort(ledger$pop, ledger$pop+stepsize)
    }
    else{
        logsurv <- function(y){
           log(1 - ledger$mortality(y))
        }
        pdeath <- sapply(ledger$pop, function(a)
           integrate(logsurv, a, a+stepsize)$value %>% exp %>% `-`(1, .))
    }
    deathVec <- as.logical(rbinom(length(ledger$pop), 1, pdeath))
    deathN <- sum(deathVec)
    newDeaths <- data.frame(
        age_start=ledger$pop[deathVec],
        age_end=ledger$pop[deathVec] + stepsize,
        cycle_death=rep(ledger$time, deathN))
    ledger$ledger <- rbind(ledger$ledger, newDeaths)
    ledger$pop <- ledger$pop[!deathVec] + stepsize
    ledger$time <- ledger$time + stepsize
    return(ledger)
}


sims <- 1000
pSize <- 1000

## comapare analytical integeration vs numerical

## Integration through R's integrate function
system.time(lapply(1:sims, function(x) 
    progress_pop(mort_ledger(rep(0, pSize), mort_unif_func), 1)$pop %>% 
        length) %>% 
        unlist) # about 92 seconds

## Integration done by hand
system.time(lapply(1:sims, function(x) 
    progress_pop(mort_ledger(rep(0, pSize), mortUnif()), 1)$pop %>% length) %>% 
    unlist) # about half a second

## comapre probabilities of death when taking many little steps vs one big step
## to ensure consistancy

mledger <- mort_ledger(rep(0, pSize), mortUnif())

lilsteps <- lapply(1:sims, function(x){
    small_jump <- mledger
    for(i in 1:100){
        small_jump <- progress_pop(small_jump, stepsize=1/100)
    }
    length(small_jump$pop)
}) %>% unlist

bigsteps <- sapply(1:sims, function(x) length(progress_pop(mledger)$pop))

data.frame(pop=c(bigsteps, lilsteps), model=rep(c("big","small"), each=sims)) %>%
    ggplot(aes(x=pop, group=model, fill=model)) + geom_density(alpha=.3)

siler_func <- function(x, a=.01, b=.01, c=.001, d=.002, f=.040){
    a * exp(-b * x) + c + d * exp(f * x)
}

ggplot(data.frame(x=c(.001, 100)), aes(x=x)) + 
    stat_function(fun=siler_func) +
    coord_trans(y="log")

ggplot(data.frame(x=c(-3,3)), aes(x=x)) + 
    stat_function(fun=mledger$mortality$mortFunc)
