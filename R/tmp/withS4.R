rm(list=ls())
##load("./inst/tinytest/anomaly_test_objects.RData")


## this replicates Example 2 in the test_anomaly.R script of the anomlay package

set.seed(0)
x <- rnorm(5000)
x[401:500] <- rnorm(100, 4, 1)
x[1601:1800] <- rnorm(200, 0, 0.01)
x[3201:3500] <- rnorm(300, 0, 10)
x[c(1000, 2000, 3000, 4000)] <- rnorm(4, 0, 100)
x <- (x - median(x))/mad(x)


## set up the cost class and methods
setGeneric("collectiveCost",function(obj,...){standardGeneric("collectiveCost")})
setGeneric("baseCost",function(obj,...){standardGeneric("baseCost")})
setGeneric("pointCost",function(obj,...){standardGeneric("pointCost")})


setClass("gaussCost", representation(summaryStats = "matrix", min_length="numeric", maxT="numeric"))

gaussCost <- function(x,min_length){
    new("gaussCost",summaryStats = cbind(cumsum(x),cumsum(x^2),x),
        min_length = min_length,maxT = length(x))
}

setMethod("collectiveCost","gaussCost",function(obj,a,b,pen){
    a <- a-1
    n <- b-a
    if(n < obj@min_length){ return(Inf) }
    if(a<1){
        Syy <- obj@summaryStats[b,2]
        Sy <- obj@summaryStats[b,1]
    }else{
        Syy <- obj@summaryStats[b,2] - obj@summaryStats[a,2]
        Sy <- obj@summaryStats[b,1] - obj@summaryStats[a,1]
    }
    ( n*log(2*pi) + Syy - (Sy^2)/n ) + pen
})

setMethod("baseCost","gaussCost",function(obj,a,b,pen=0){
    a <- a-1
    n <- b-a
    if(a<1){
        Syy <- obj@summaryStats[b,2]
    }else{
        Syy <- obj@summaryStats[b,2] - obj@summaryStats[a,2]
    }
    ( n*log(2*pi) + Syy ) + pen
})

setMethod("pointCost","gaussCost",function(obj,a,pen){
    if(a<2){
        Syy <- obj@summaryStats[a,2]
    }else{
        Syy <- obj@summaryStats[a,2] - obj@summaryStats[a-1,2]
    }
    gamma <- exp(-(1+pen))
    log(2*pi) + log(gamma + Syy) + 1 +pen
    ###-2 * dnorm(obj@summaryStats[a,3],0,abs(obj@summaryStats[a,3]),log=TRUE) + pen
})

fCost <- gaussCost(x,10)


system.time({
    tmp <- rep(NA,length(x))
    for(ii in 1:length(x)){
        tmp[ii] <- collectiveCost(fCost,ii,length(x),0)
    }
})


## ########################################
## set up the partition class and methods in S3

## set up generic methods

addCollective <- function(obj,...){ UseMethod("addCollective",obj) }
addBase <- function(obj,...){ UseMethod("addBase",obj) }
addPoint <- function(obj,...){ UseMethod("addPoint",obj) }
collective_anomalies <- function(obj,...){ UseMethod("collective_anomalies",obj) }
point_anomalies <- function(obj,...){ UseMethod("point_anomalies",obj) }

partition <- function(beta,betaP){ structure(list(ca=list(),pa=list(),beta=beta,betaP=betaP,cost=0,last_time=0),class="partition") }

addCollective.partition <- function(p,x,s,e){
    cst <- collectiveCost(x,s,e,p$beta)
    p$ca[[length(p$ca)+1]] <- c(start=s,end=e,cost=cst)
    p$cost <- p$cost + cst
    p$last_time <- e
    return(p)
}

addBase.partition <- function(p,x,s,e){
    p$cost <- p$cost + baseCost(x,s,e,0)
    p$last_time <- e
    return(p)
}

addPoint.partition <- function(p,x,s){
    cst <- pointCost(x,s,p$betaP)
    p$pa[[length(p$pa)+1]] <- c(location=s,cost=cst) 
    p$cost <- p$cost + cst
    p$last_time <- s
    return(p)
}

collective_anomalies <- function(p){ as.data.frame( do.call(rbind,p$ca) ) }
point_anomalies <- function(p){ as.data.frame( do.call(rbind,p$pa) ) }


tmp <- partition(3*log(length(x)),3*log(length(x)))

addCollective(tmp,fCost,1,400)
addBase(tmp,fCost,1,400)
addPoint(tmp,fCost,401)

gc()

prune=TRUE
##system.time({
##p <- profvis::profvis({
ctlg <- list()
ctlg[[1]] <- partition(3*log(length(x)),3*log(length(x))) ###new("partition") ## offset by 1 versus time!!
##METHOD <- selectMethod(addCollective, class(ctlg[[1]]))

cnst <- max(ctlg[[1]]$beta,ctlg[[1]]$betaP)

for(tt in 1:length(x)){
    if(tt %% 100==0) {
        ## Print on the screen some message
        cat(paste0("time step: ", tt, "\n"))
    }
    ##browser()
    n <- length(ctlg)

    ## compute C2 from paper
    opt <- addBase(ctlg[[n]],fCost,tt,tt)
    
    ## compute C3 from paper
    tmp <- addPoint(ctlg[[n]],fCost,tt)
    if(tmp$cost < opt$cost){ opt <- tmp }

    ## loop ctlg
    ctlgCost <- rep(-Inf,length(ctlg))
    for(ii in 1:length(ctlg)){
        tmp <- addCollective(ctlg[[ii]], fCost,ctlg[[ii]]$last_time + 1,tt)
        if(tmp$cost < opt$cost){ opt <- tmp }
        if(ctlg[[ii]]$last_time < tt-10){
            ctlgCost[ii] <- tmp$cost
        }
    }
    if(prune){ ctlg <- ctlg[ ctlgCost <= opt$cost+cnst ] }
    
    ctlg[[length(ctlg)+1]] <- opt  
}
##})

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
       
