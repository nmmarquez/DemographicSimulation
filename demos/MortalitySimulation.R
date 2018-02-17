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

load("../data/ProbDeath.Rda")

ages <- c(7/365/2, mean(c(28/365, 7/365)), mean(c(1 ,28/365)), 3, 
          seq(7.5, 107.5, 5))
yearsPassed <- c(7/365, 21/365, (365-28)/365, 4, rep(5, 21))

ggplot(ProbDeath, aes(x=age, y=px, color=location, group=location)) + 
    geom_line() + 
    coord_trans(y="log")

USDeath <- filter(ProbDeath, location=="United States")

siler_func <- function(x, params){
    a <- params[1]
    b <- params[2]
    c <- params[3]
    d <- params[4]
    f <- params[5]
    a * exp(-b * x) + c + d * exp(f * x)
}

lik_func <- function(params, ages=USDeath$age, y=USDeath$px){
    modelp <- exp(params)
    silerfit <- siler_func(ages, modelp)
    return(sqrt(mean((log(y)-log(silerfit))^2)))
}

sparams <- log(c(11.3915130, 1.944036, .000000000003, 3.835808e-04, 7.967152e-02))

optSiler <- optim(sparams, lik_func)

opt_siler_func <- function(x){
    siler_func(x, exp(optSiler$par))
}


ggplot(data.frame(x=c(.001, 100)), aes(x=x)) + 
    stat_function(fun=function(y) opt_siler_func(y)) +
    coord_trans(y="log") + 
    geom_point(data=USDeath, mapping=aes(x=age, y=px)) +
    xlim(c(0.001, 100)) +
    labs(x="Age", y="px", title="US Probability of Death by Age")


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


silerLedger <- mort_ledger(rep(0, 10000), opt_siler_func)

silerLedgerList <- list(mort_ledger(rep(0, 10000), opt_siler_func))
steps <- c(7/356, 21/356, (356-21-7)/356, 4, rep(5, 14))

for(i in 1:length(steps)){
    print(i)
    silerLedgerList[[i + 1]] <- progress_pop(silerLedgerList[[i]], steps[i])
}

estPx <- silerLedgerList[[19]]$ledger %>% 
    mutate(age=(age_start + age_end)/2) %>%
    group_by(age) %>% summarise(deaths=n()) %>%
    mutate(pop=10000) %>%
    mutate(pop=pop- cumsum(deaths) + deaths) %>%
    mutate(rawp=deaths/pop) %>%
    mutate(time=yearsPassed[1:length(steps)]) %>%
    mutate(px=1-((1-rawp)^((time)^-1)))


ggplot(data.frame(x=c(.001, 100)), aes(x=x)) + 
    stat_function(fun=function(y) opt_siler_func(y)) +
    coord_trans(y="log") + 
    geom_point(data=USDeath, mapping=aes(x=age, y=px)) +
    xlim(c(0.001, 100)) +
    labs(x="Age", y="px", title="US Probability of Death by Age") +
    geom_point(data=estPx, mapping=aes(x=age, y=px), color=2)

ggplot(data.frame(x=c(.001, 100)), aes(x=x)) + 
    stat_function(fun=function(y) opt_siler_func(y)) +
    geom_point(data=USDeath, mapping=aes(x=age, y=px)) +
    xlim(c(0.001, 100)) +
    labs(x="Age", y="px", title="US Probability of Death by Age") +
    geom_point(data=estPx, mapping=aes(x=age, y=px), color=2)

# logSurvSiler <- function(y){
#     log(1 - opt_siler_func(y))
# }
# 
# my.uniroot <- function(x) uniroot(logSurvSiler, c(0, 100), tol = 0.0001)$root
# x <- runif(100)
# system.time({
#     r <- vapply(x, my.uniroot, numeric(1))
# })

#msteps <-