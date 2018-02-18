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

# MX1982 <- GenPopFuncs(location="Mexico", year=1982)
# MX2016 <- GenPopFuncs(location="Mexico", year=2016)
# MsimsOld <- MX1982$simPop(10000, 1)
# MsimsNew <- MX2016$simPop(10000, 1)
# data.frame(age=c(MsimsOld, MsimsNew)) %>%
#     mutate(year=rep(c("1982", "2016"), each=10000)) %>% 
#     ggplot(aes(x=age, group=year, fill=year)) + geom_density(alpha=.2)
