#' Generate Mortality Functions For a Location Year
#'
#' @description Generates Survial, Hazard, Failure(CDF), 
#' Inverse Failure(invCDF), derivative of failure, and simulation functions for
#' a specified country, year, sex.
#' 
#' @param location_ character, country name
#' @param year_  1970 <= int <= 2016, year of mortality function to generate 
#' @param sex_ character, either 'Female', 'Male', or 'Both'
#' @param max_age float, float for the max age allowed for simulatuins generally a value >= 120 
#'
#' @return List of Demographic Mortality and Simulation Functions
#' 
#' Generates Typical Demographic Mortality Functions Based on Global Burden of 
#' Disease Estimate Data for Survival, Failure(CDF), Hazard, Inverse Failure, 
#' and the Age of Death Distribution(Derivative of Failure) using a non 
#' Parametric Spline interpolation of the the (CDF) with strictly >= 0 for
#' derivatives. Also inlcudes a convience function to simulate data. A list of 
#' names of possible countries to use for the simulation are shown in the code
#' example.
#'
#' @examples
#' # Show list of countries we can simulate from 
#' unique(DFDeath$location)
#' 
#' # run plotting and simulation code
#' require(dplyr)
#' require(ggplot2)
#' # Generate the demographic functions for Mexico 1980
#' MX2016F <- GenPopFuncs(location_="Mexico", year_=1980, sex_="Both")
#' 
#' 
#' # plot some of the functions
#' data.frame(Age=seq(.01, 120, .01)) %>%
#'     mutate(CDF=MX2016F$CDF(Age)) %>%
#'     ggplot(aes(x=Age, y=CDF)) + 
#'     geom_line() + 
#'     coord_trans(y="log") +
#'     labs(title="Failure Function of Mexico Mortality: 1980", x="Age", y="Failure")
#' 
#' data.frame(Age=seq(.1, 100, .005)) %>%
#'     mutate(Hazard=MX2016F$Hxfunc(Age)) %>%
#'     ggplot(aes(x=Age, y=Hazard)) + 
#'     geom_line() + 
#'     coord_trans(y="log") +
#'     labs(title="Hazard Function of Mexico Mortality: 1980", x="Age", y="Hazard")
#' 
#' m <- 10000
#' system.time(simDeaths <- lapply(1:10, function(y) MX2016F$simPop(m, 1)))
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
#'     labs(title="Non-Parametric Simulated Instantaneous Hazard", x="Age", y="Hazard")
#' 
#' 
#' @export
#' 

GenPopFuncs <- 
    function(location_="United States", year_=2016, sex_="Both", max_age=140){
    library(splines)
    library(parallel)
    DFsub <- subset(DFDeath,location==location_ & year == year_ & sex == sex_)
    CDF <- splinefun(
        c(-1, 0, DFsub$age_end, max_age, max_age+1), 
        c(0, 0, DFsub$Fx, 1, 1), 
        method="monoH.FC")
    Sxfunc <- function(x){
        1 - CDF(x)
    }
    
    PDF <- function(x){
        CDF(x, deriv=1)
    }
    
    Hxfunc <- function(x){
        PDF(x) / Sxfunc(x)
    }
    
    invCDF <- function(x){
        f_ <- function(y) sqrt((CDF(y)-x)^2)
        optim(0, f_, method="Brent", lower=0, upper=max_age+2)$par
    }
    
    simPop <- function(n, mc.cores=1){
        if(mc.cores == 1){
            return(sapply(runif(n), invCDF))
        }
        else{
            return(unlist(mclapply(runif(n), invCDF,mc.cores=mc.cores)))
        }
    }
    
    list(CDF=CDF, PDF=PDF, Sxfunc=Sxfunc, Hxfunc=Hxfunc, 
         invCDF=invCDF, simPop=simPop)
}
