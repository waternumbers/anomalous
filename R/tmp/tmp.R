rm(list=ls())
devtools::load_all()




## this replicates Example 2 in the test_anomaly.R script of the anomlay package

set.seed(0)
x <- rnorm(5000)
x[401:500] <- rnorm(100, 4, 1)
x[1601:1800] <- rnorm(200, 0, 0.01)
x[3201:3500] <- rnorm(300, 0, 10)
x[c(1000, 2000, 3000, 4000)] <- rnorm(4, 0, 100)
x <- (x - median(x))/mad(x)


fCost <- gaussMean(x)


system.time({
    tmp <- rep(NA,length(x))
    for(ii in 1:length(x)){
        tmp[ii] <- collectiveCost(fCost,ii,length(x),0)
    }
})


tmp <- partition(3*log(length(x)),3*log(length(x)),10)

addCollective(tmp,1,400,fCost)
addBase(tmp,1,400,fCost)
addPoint(tmp,401,fCost)

gc()


opt <- capa(tmp,x=fCost)

##htmlwidgets::saveWidget(p, "profile.html")
     
## Can open in browser from R
##browseURL("profile.html",browser="firefox")

## opt is the obect to look at....

## ############################################
## compare
orig <- readRDS("test_output.rds")
##orig$pa$location
tinytest::expect_equal(orig$ca, collective_anomalies(opt))
tinytest::expect_equal(orig$pa, point_anomalies(opt))
       
