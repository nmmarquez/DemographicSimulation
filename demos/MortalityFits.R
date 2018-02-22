# Mortality fits
rm(list=ls())
library(rgr)
library(dplyr)
library(parallel)
library(ggplot2)
library(splines)
library(GoFKernel)
library(sn)

load("../data/DFDeath.RData")
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

data.frame(x=seq(0.025, 115, .025)) %>% 
    mutate(y=linInterp(x, deriv=1)) %>%
    ggplot(aes(x, y)) + geom_line() + 
    labs(title="First Derivative of Interpolated CDF(Interpolated PDF)") +
    xlim(c(0,3))

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
system.time(simDeaths <- lapply(1:100, function(y) simUSAMort(m, 6)))

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

data.frame(death=unlist(simDeaths[1:9]), simulation=rep(1:9, each=m)) %>%
    ggplot(aes(x=death)) + geom_density() + 
    facet_wrap(~simulation)


transparams <- function(params){
    exp(params)
}

doublikfunc <- function(params, datur=simDeaths[[1]]){
    shape1 <- transparams(params)[1]
    scale1 <- transparams(params)[2]
    #p <- transparams(params)[3]
    alpha2 <- transparams(params)[3]
    beta2 <- transparams(params)[4]
    dlo <- datur[datur < 5]
    dhi <- datur[datur >= 5]
    nll <- sum(dgamma(dlo, shape=shape1, scale=scale1, log=TRUE) * -1)
    nll <- sum(nll + dbeta(dhi/110, alpha2, beta2, log=TRUE) * -1)
    nll
}

Opt <- optim(c(log(c(1,1)), log(c(3, 1.8))), doublikfunc)

simDist <- function(n, params, p=.007){
    ydeath <- as.logical(rbinom(n, 1, p))
    shape1 <- transparams(params)[1]
    scale1 <- transparams(params)[2]
    alpha2 <- transparams(params)[3]
    beta2 <- transparams(params)[4]
    c(rgamma(n, shape=shape1, scale=scale1)[ydeath], 
      rbeta(n, alpha2, beta2)[!ydeath] * 110)
}

data.frame(deaths=c(simDist(m, Opt$par), simDeaths[[1]])) %>%
    mutate(type=rep(c("Parametric", "Non-Parametric"), each=m)) %>%
    ggplot(aes(x=deaths, group=type, fill=type)) + 
    geom_density(alpha=.2)

data.frame(deaths=c(simDist(m, Opt$par), simDeaths[[1]])) %>%
    mutate(type=rep(c("Parametric", "Non-Parametric"), each=m)) %>%
    filter(deaths<=5) %>%
    ggplot(aes(x=deaths, group=type, fill=type)) + 
    geom_density(alpha=.2)


doublikfunc2 <- function(params, datur=simDeaths[[1]]){
    shape1 <- transparams(params)[1]
    scale1 <- transparams(params)[2]
    omega2 <- transparams(params)[3]
    alpha2 <- -transparams(params)[4]
    xi2 <- params[5]
    dlo <- datur[datur < 5]
    dhi <- datur[datur >= 5]
    nll <- sum(dgamma(dlo, shape=shape1, scale=scale1, log=TRUE) * -1)
    nll <- sum(nll + dsn(dhi, xi2, omega2, alpha2, log=TRUE) * -1)
    nll
}

Opt2 <- optim(c(log(c(1,1)), log(c(3, 1.8)), 80), doublikfunc2)

# simDist2 <- function(n, params, p=.007){
#     ydeath <- as.logical(rbinom(n, 1, p))
#     shape1 <- transparams(params)[1]
#     scale1 <- transparams(params)[2]
#     omega2 <- transparams(params)[3]
#     alpha2 <- -transparams(params)[4]
#     xi2 <- params[5]
#     sim_ <- c(rgamma(n, shape=shape1, scale=scale1)[ydeath], 
#               rsn(n, xi2, omega2, alpha2)[!ydeath])
#     sim_[sim_ <=0] <- sample(sim_[sim_ >0], size=sum(sim_<=0), replace=T)
#     sim_
# }
# 
# data.frame(deaths=c(simDist2(m, Opt2$par), simDeaths[[1]])) %>%
#     mutate(type=rep(c("Parametric", "Non-Parametric"), each=m)) %>%
#     ggplot(aes(x=deaths, group=type, fill=type)) + 
#     geom_density(alpha=.2)
# 
# data.frame(deaths=c(simDist2(m, Opt2$par), simDeaths[[1]])) %>%
#     mutate(type=rep(c("Parametric", "Non-Parametric"), each=m)) %>%
#     filter(deaths<=5) %>%
#     ggplot(aes(x=deaths, group=type, fill=type)) + 
#     geom_density(alpha=.2)
# 
# simParamDeaths <- lapply(1:100, function(x) simDist2(m, Opt2$par))
# 
# simParamDF <- bind_rows(lapply(1:100, function(i) 
#     aggData(simParamDeaths[[i]], i)))
# 
# simParamDF %>% filter(age_end < 115 & hx != 0) %>%
#     ggplot(aes(x=age_end, y=hx, color=simulation, group=simulation)) + 
#     geom_line(alpha=.3) + 
#     geom_line(aes(x=age_end, y=hx, group=1), data=USDeath, color="red") + 
#     coord_trans(y="log") + 
#     labs("Simulated Instantaneous Hazard")
# 
# simParamDF %>% filter(age_end < 115 & Sx < 1 & Sx > 0) %>%
#     ggplot(aes(x=age_end, y=Sx, color=simulation, group=simulation)) + 
#     geom_line(alpha=.3) + 
#     geom_line(aes(x=age_end, y=Sx, group=1), data=USDeath, color="red") +
#     coord_trans(y="logit") +
#     labs(title="Simulated Survival Curves")

paramFitFunc <- function(
    start_params=c(0, 0, 1.1, .57, 80),
    location_="United States", 
    year_=2016, 
    sex_="Both",
    max_age=140,
    returnParams=FALSE
    ){
    DFsub <- subset(DFDeath,location==location_ & year == year_ & sex == sex_)
    childDF <- DFsub %>% filter(age_group_id <= 5)
    adultDF <- DFsub %>% filter(age_group_id > 5)
    childCDF <- childDF$Fx / childDF$Fx[nrow(childDF)]
    adultCDF <- c(adultDF$Fx - childDF$Fx[nrow(childDF)], 1)
    
    doublikfunc2 <- function(params){
        shape1 <- exp(params)[1]
        scale1 <- exp(params)[2]
        omega2 <- exp(params)[3]
        alpha2 <- -exp(params)[4]
        xi2 <- params[5]
        childCDFHat <- pgamma(childDF$age_end, shape=shape1, scale=scale1)
        adultCDFHat <- psn(c(adultDF$age_end, max_age), xi2, omega2, alpha2)
        adultWeight <- diff(c(childDF$Fx[nrow(childDF)], adultDF$Fx, 1))
        childLoss <- (childCDF - childCDFHat)^2
        adultLoss <- (adultCDF - adultCDFHat)^2 * adultWeight
        sqrt(mean(c(childLoss, adultLoss)))
    }
    
    Opt2 <- optim(start_params, doublikfunc2)
    
    p <- 1 - cumprod(childDF$qx)[nrow(childDF)]
    shape1 <- exp(Opt2$par)[1]
    scale1 <- exp(Opt2$par)[2]
    omega2 <- exp(Opt2$par)[3]
    alpha2 <- -exp(Opt2$par)[4]
    xi2 <- Opt2$par[5]
    
    #plot(c(adultDF$age_end, max_age), adultCDF)
    #lines(seq(5, 140, .01), psn(seq(5, 140, .01), xi2, omega2, alpha2))
    
    if(returnParams){
        return(c(shape=shape1, scale=scale1, xi=xi2, 
                 omega=omega2, alpha=alpha2, p=p))
    }
    
    simDist2 <- function(n){
        ydeath <- as.logical(rbinom(n, 1, p))
        sim_ <- c(rgamma(n, shape=shape1, scale=scale1)[ydeath], 
                  rsn(n, xi2, omega2, alpha2)[!ydeath])
        sim_[sim_ <=0] <- sample(sim_[sim_ >0], size=sum(sim_<=0), replace=T)
        sim_
    }
    
    return(simDist2)
}

funcTest <- paramFitFunc()

simDF <- data.frame(Deaths=c(
    unlist(lapply(1:9, function(x) funcTest(m))),
    unlist(lapply(1:9, function(x) simUSAMort(m, 6))))) %>%
    mutate(simulation=rep(rep(1:9, 2), each=m)) %>%
    mutate(Type=rep(c("Parametric", "Non-Parametric"), each=m*9))


simDF %>% 
    ggplot(aes(x=Deaths, group=Type, fill=Type)) + 
    geom_density(alpha=.2) +
    facet_wrap(~simulation)
    

simDF %>% 
    filter(Deaths <= 5) %>%
    ggplot(aes(x=Deaths, group=Type, fill=Type)) + 
    geom_density(alpha=.2) +
    facet_wrap(~simulation)

simParamDeaths <- lapply(1:100, function(x) funcTest(m))

simParamDF <- bind_rows(lapply(1:100, function(i) 
    aggData(simParamDeaths[[i]], i)))

simParamDF %>% filter(age_end < 115 & hx != 0) %>%
    ggplot(aes(x=age_end, y=hx, color=simulation, group=simulation)) + 
    geom_line(alpha=.3) + 
    geom_line(aes(x=age_end, y=hx, group=1), data=USDeath, color="red") + 
    coord_trans(y="log") + 
    labs("Simulated Instantaneous Hazard")

simParamDF %>% filter(age_end < 12.5 & hx != 0) %>%
    ggplot(aes(x=age_end, y=hx, color=simulation, group=simulation)) + 
    geom_line(alpha=.3) + 
    geom_line(aes(x=age_end, y=hx, group=1), 
              data=USDeath %>% filter(age_end < 12.5), color="red") + 
    coord_trans(y="log") + 
    labs("Simulated Instantaneous Hazard")


simParamDF %>% filter(age_end < 115 & Sx < 1 & Sx > 0) %>%
    ggplot(aes(x=age_end, y=Sx, color=simulation, group=simulation)) + 
    geom_line(alpha=.3) + 
    geom_line(aes(x=age_end, y=Sx, group=1), data=USDeath, color="red") +
    labs(title="Simulated Survival Curves")
