rm(list=ls())
load("./inst/tinytest/anomaly_test_objects.RData")

## this replicates Example 2 in the test_anomaly.R script of the anomlay package

set.seed(0)
x <- rnorm(5000)
x[401:500] <- rnorm(100, 4, 1)
x[1601:1800] <- rnorm(200, 0, 0.01)
x[3201:3500] <- rnorm(300, 0, 10)
x[c(1000, 2000, 3000, 4000)] <- rnorm(4, 0, 100)
x <- (x - median(x))/mad(x)

## #################################################
## run the anomlay package and save output
## The anomaly package called is from the master branch of https://github.com/grosed/anomaly
## taken after commit b24aad5
anomOut <- anomaly::capa(x,type="mean") ##part$beta,part$betaP,"mean")
orig <- list(ca = anomaly::collective_anomalies(anomOut),
             pa = anomaly::point_anomalies(anomOut))

## #################################################
## Very direct coding of the capa algorithm with no pruning....
## feel the slow.....

gaussMean <- function(x,a,b,pen){
    -2 * sum(dnorm(x[a:b],mean(x[a:b]),1,log=TRUE)) + pen
}

gaussKnown <- function(x,a,b,pen){
    -2 * sum(dnorm(x[a:b],0,1,log=TRUE)) + pen
}

gaussPoint <- function(x,a,pen){
    ## TODO apply correction to stop point anomalies at 0 by inflating the cost
    -2 * dnorm(x[a],0,abs(x[a]),log=TRUE) + pen
}

## create a blank record of partitions with the parameters taken from the anomaly output
part <- list(cost = 0,
             beta = anomOut@beta[1],
             betaP = anomOut@beta_tilde,
             min_length = anomOut@min_seg_len,
             ca = data.frame(start=numeric(0),end=numeric(0),cost=numeric(0)),
             pa = data.frame(location=numeric(0),cost=numeric(0))
             )

addCollective <- function(p,x,s,e){
    tmp <- data.frame(start = s,end = e,cost = gaussMean(x,s,e,p$beta))
    if(tmp$end-tmp$start+1 < p$min_length){ tmp$cost <- Inf } ## ensure segments are at least 10 long
    p$ca <- rbind(p$ca, tmp)
    p$cost <- p$cost + tmp$cost
    return(p)
}

addBase <- function(p,x,s,e){
    p$cost <- p$cost + gaussKnown(x,s,e,0)
    return(p)
}

addPoint <- function(p,x,s){
    tmp <- data.frame(location=s, cost = gaussPoint(x,s,p$betaP))
    p$pa <- rbind(p$pa,tmp)
    p$cost <- p$cost + tmp$cost
    return(p)
}

ctlg <- list()
ctlg[[1]] <- part ## offset by 1 versus time!!
for(tt in 1:length(x)){
    if(tt %% 100==0) {
      # Print on the screen some message
      cat(paste0("time step: ", tt, "\n"))
    }
    
    ## compute C2 from paper
    opt <- addBase(ctlg[[tt]],x,tt,tt)
    
    ## compute C3 from paper
    tmp <- addPoint(ctlg[[tt]],x,tt)
    if(tmp$cost < opt$cost){ opt <- tmp }
    
    ## loop ctlg
    for(ii in 1:tt){
        tmp <- addCollective(ctlg[[ii]], x,ii,tt)
        if(tmp$cost < opt$cost){ opt <- tmp }
    }

    ctlg[[tt+1]] <- opt  
}
## opt is the obect to look at....

## ############################################
## compare
orig$ca[,c("start","end")] == opt$ca[,c("start","end")]

setequal(orig$pa$location,opt$pa$location) ## not equal
all( opt$pa$location %in% orig$pa$location ) ## the points in the basic algorithm are a subset of those found in anomaly
idx <- setdiff( orig$pa$location, opt$pa$location ) # the extra points returned by anomaly

## to test the extra point anomalies
## From the paper the cost must be less then if they were from N(0,1)
## else they would not be selected as optimal
## Assuming the cost functions above are correct
print(idx)
## [1] 3249 3277 3313 3353 3354 3408 3416 3433 3451 3483 3496
bc <- sapply(idx,function(ii){gaussKnown(x,ii,ii,0)}) ## cost of each point if background i.e. N(0,1)
pc <- sapply(idx,function(ii){gaussPoint(x,ii,part$betaP)}) ## cost of each point if a point anomaly
all(bc<pc)
       
saveRDS(opt,"test_output.rds")
