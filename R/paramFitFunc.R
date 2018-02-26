#' Generate Parametric Mortality Smulator 
#'
#' @description Generates a parametric mortality simulator using a gamma 
#' distribution for under 5 mortality and a skewed normal for above five mortality
#' where the probability of puuling from the gamma distribution is generated from the 
#' observed data. 
#' 
#' @param start_params vector of length five of starting values for parameters 
#' @param location_ character, country name
#' @param year_  1970 <= int <= 2016, year of mortality function to generate 
#' @param sex_ character, either 'Female', 'Male', or 'Both'
#' @param max_age float, float for the max age allowed for simulatuins generally a value >= 120 
#' @param return_params logical, whether to return the fitted parameter values rather than the simulation function
#'
#' @return List of Demographic Mortality and Simulation Functions
#' 
#' Generates a parametric mortality simulator using a gamma 
#' distribution for under 5 mortality and a skewed normal for above five mortality
#' where the probability of puuling from the gamma distribution is generated from the 
#' observed data. 
#'
#' @examples
#' # Show list of countries we can simulate from 
#' unique(DFDeath$location)
#' 
#' # run plotting and simulation code
#' require(dplyr)
#' require(ggplot2)
#' # Generate the demographic functions for Mexico 1980
#' MX2016sim <- paramFitFunc(location_="Mexico", year_=1980, sex_="Both")
#' 
#' 
#' m <- 10000
#' system.time(simDeaths <- lapply(1:10, function(y) MX2016sim(m)))
#' 
#' MXDeath <- DFDeath %>%
#'     filter(location=="Mexico" & year==1980 & sex=="Both")
#' 
#' aggData <- function(sims, sim_num, m_=m){
#'     MXDeath %>% select(age_group_id, age_time, age_end) %>%
#'         mutate(ldeaths=sapply(age_end, function(a) sum(sims < a))) %>%
#'         mutate(deaths=ldeaths-lag(ldeaths)) %>%
#'         mutate(deaths=ifelse(is.na(deaths), ldeaths, deaths)) %>%
#'         mutate(pop_size=m_-lag(ldeaths)) %>%
#'         mutate(pop_size=ifelse(is.na(pop_size), m_, pop_size)) %>%
#'         mutate(px=deaths/pop_size, qx=1-px) %>%
#'         mutate(hx=1-(qx^(1/age_time))) %>% 
#'         mutate(Sx=cumprod(qx), Fx=1-Sx) %>%
#'         mutate(simulation=sim_num)
#' }
#' 
#' simDF <- bind_rows(lapply(1:10, function(i) aggData(simDeaths[[i]], i))) 
#' 
#' simDF %>% filter(age_end < 115 & hx != 0) %>%
#'     ggplot(aes(x=age_end, y=hx, color=simulation, group=simulation)) + 
#'     geom_line(alpha=.3) + 
#'     geom_line(aes(x=age_end, y=hx, group=1), data=MXDeath, color="red") + 
#'     coord_trans(y="log") + 
#'     labs(title="Parametric Simulated Instantaneous Hazard", x="Age", y="Hazard")
#' 
#' 
#' @export
#' 

paramFitFunc <- function(
    start_params=c(0, 0, 1.1, .57, 80),
    location_="United States", 
    year_=2016, 
    sex_="Both",
    max_age=140,
    returnParams=FALSE
){
    library(sn)
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
