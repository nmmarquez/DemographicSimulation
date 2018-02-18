# Mortality fits
rm(list=ls())
library(rgr)
library(dplyr)
library(parallel)
library(ggplot2)
library(splines)
library(GoFKernel)

load("../data/DFDeath.Rda")
USDeath <- DFDeath %>%
    filter(location == "United States" & sex == "Both" & year == 2016)

ggplot(USDeath, aes(x=age_end, y=Fx)) + 
    geom_point() +
    labs(x="Age", y="F(x)", title="Observed Period Failure Rate: 2016")

ggplot(USDeath, aes(x=age_end, y=Fx)) +
    geom_point() +
    coord_trans(y="log") +
    labs(x="Age", y="F(x)", title="Observed Period Failure Rate: 2016")

ggplot(USDeath, aes(x=age_end, y=hx)) + 
    geom_point() + 
    coord_trans(y="log") +
    labs(x="Age", y="h(x)", title="Instantaneous Period Hazard: 2016")

fitfunc <- function(x, params){
    expp <- exp(params)
    pgamma(x, expp[1], expp[2])# * p
}

likfunc <- function(params, datur=USDeath$Fx, ages=USDeath$age_end){
    fit <- fitfunc(ages, params)
    return(sqrt(mean(((datur)-(fit))^2)))
}

Opt <- optim(c(0,0), likfunc)
Wfit <- exp(Opt$par)

ggplot(USDeath, aes(x=age_end, y=Fx)) + 
    geom_point() +
    labs(x="Age", y="F(x)", title="Observed Period Failure Rate: 2016") + 
    stat_function(fun=function(y) fitfunc(y, Opt$par))

ggplot(USDeath, aes(x=age_end, y=Fx)) + 
    geom_point() +
    labs(x="Age", y="F(x)", title="Observed Period Failure Rate: 2016") + 
    stat_function(fun=function(y) fitfunc(y, Opt$par)) +
    xlim(c(0,5)) + ylim(c(0,.05))

linInterp <- splinefun(
    c(-1, 0, USDeath$age_end, 115, 120), 
    c(0, 0, USDeath$Fx, 1, 1), 
    method="monoH.FC")

ggplot(USDeath, aes(x=age_end, y=Fx)) + 
    geom_point() +
    labs(x="Age", y="F(x)", title="Observed Period Failure Rate: 2016") + 
    stat_function(fun=function(x) linInterp(x))

ggplot(USDeath %>% filter(age_end < 1.1), aes(x=age_end, y=Fx)) + 
    geom_point() +
    labs(x="Age", y="F(x)", title="Observed Period Failure Rate: 2016") + 
    stat_function(fun=function(x) linInterp(x))

hxfunc <- function(x){
    linInterp(x, deriv=1)/(1-linInterp(x))
}
    

data.frame(x=USDeath$age_end) %>%#seq(0.01, 105, .01)) %>% 
    mutate(y=linInterp(x)) %>%
    ggplot(aes(x, y)) + geom_line() + 
    labs(title="Interpolated CDF of Failure Rate")


data.frame(x=seq(0.025, 115, .025)) %>% 
    mutate(y=linInterp(x, deriv=1)) %>%
    ggplot(aes(x, y)) + geom_line() + 
    labs(title="First Derivative of Interpolated CDF(Interpolated PDF)")

data.frame(x=USDeath$age_end) %>%#seq(0.01, 105, .01)) %>% 
    mutate(y=hxfunc(x)) %>%
    ggplot(aes(x, y)) + geom_line() + 
    coord_trans(y="log") +
    labs(title="Interpolated Hazard Function")

rSx <- function(y, max_age=140){
    f_ <- function(x) sqrt((linInterp(x)-y)^2)
    optim(0, f_, method="Brent", lower=0, upper=max_age)$par
}

simUSAMort <- function(n, mc.cores=1, max_age=140){
    if(mc.cores == 1){
        return(sapply(runif(n), rSx, max_age=max_age))
    }
    else{
        return(unlist(mclapply(runif(n), rSx, max_age=max_age,
                               mc.cores=mc.cores)))
    }
}

set.seed(123)
m <- 10000
simDeaths <- lapply(1:100, function(y) simUSAMort(m, 6))

aggData <- function(sims, sim_num, m_=m){
    USDeath %>% select(age_group_id, age_time, age_end) %>%
        mutate(ldeaths=sapply(age_end, function(a) sum(sims < a))) %>%
        mutate(deaths=ldeaths-lag(ldeaths)) %>%
        mutate(deaths=ifelse(is.na(deaths), ldeaths, deaths)) %>%
        mutate(pop_size=m_-lag(ldeaths)) %>%
        mutate(pop_size=ifelse(is.na(pop_size), m_, pop_size)) %>%
        mutate(px=deaths/pop_size, qx=1-px) %>%
        mutate(hx=1-(qx^(1/age_time))) %>% 
        mutate(Sx=cumprod(qx), Fx=1-Sx) %>%
        mutate(simulation=sim_num)
}

simDF <- bind_rows(lapply(1:100, function(i) aggData(simDeaths[[i]], i))) 

simDF %>% filter(age_end < 115 & px < 1 & px > 0) %>%
    ggplot(aes(x=age_end, y=px, color=simulation, group=simulation)) + 
    geom_line(alpha=.3) + 
    geom_line(aes(x=age_end, y=px, group=1), data=USDeath, color="red") + 
    coord_trans(y="logit") + 
    labs("Simulated Probability of Death")

simDF %>% filter(age_end < 115 & hx != 0) %>%
    ggplot(aes(x=age_end, y=hx, color=simulation, group=simulation)) + 
    geom_line(alpha=.3) + 
    geom_line(aes(x=age_end, y=hx, group=1), data=USDeath, color="red") + 
    coord_trans(y="log") + 
    labs("Simulated Instantaneous Hazard")

simDF %>% filter(age_end < 115 & Sx < 1 & Sx > 0) %>%
    ggplot(aes(x=age_end, y=Sx, color=simulation, group=simulation)) + 
    geom_line(alpha=.3) + 
    geom_line(aes(x=age_end, y=Sx, group=1), data=USDeath, color="red") +
    coord_trans(y="logit") +
    labs(title="Simulated Survival Curves")

