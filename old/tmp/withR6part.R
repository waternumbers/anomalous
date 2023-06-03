rm(list=ls())
library(R6)
##load("./inst/tinytest/anomaly_test_objects.RData")


## this replicates Example 2 in the test_anomaly.R script of the anomlay package

set.seed(0)
x <- rnorm(5000)
x[401:500] <- rnorm(100, 4, 1)
x[1601:1800] <- rnorm(200, 0, 0.01)
x[3201:3500] <- rnorm(300, 0, 10)
x[c(1000, 2000, 3000, 4000)] <- rnorm(4, 0, 100)
x <- (x - median(x))/mad(x)



gaussCost <- R6Class("gaussCost",
                 public = list(
                     summaryStats = NULL,
                     beta = NULL,
                     betaP = NULL,
                     min_length = NULL,
                     initialize = function(x,beta,betaP,min_length){
                         self$summaryStats <- cbind(cumsum(x),cumsum(x^2),x)
                         self$beta <- beta
                         self$betaP <- betaP
                         self$min_length <- min_length
                     },
                     collective = function(a,b){
                         a <- a-1
                         n <- b-a
                         if(n < self$min_length){ return(Inf) }
                         ##browser()
                         if(a<1){
                             Syy <- self$summaryStats[b,2]
                             Sy <- self$summaryStats[b,1]
                         }else{
                             Syy <- self$summaryStats[b,2] - self$summaryStats[a,2]
                             Sy <- self$summaryStats[b,1] - self$summaryStats[a,1]
                         }
                         ##print(Syy)
                         ##print( Sy )
                         ##print(n)
                         ##print( n*log(2*pi) )
                         ( n*log(2*pi) + Syy - (Sy^2)/n ) + self$beta
                     },
                     collective2 = function(a,b){
                         if(b-a+1 < self$min_length){ return(Inf) }
                         -2 * sum(dnorm(self$summaryStats[a:b,3],mean(self$summaryStats[a:b,3]),1,log=TRUE)) + self$beta
                     },
                     base = function(a,b){
                         a <- a-1
                         n <- b-a
                         if(a<1){
                             Syy <- self$summaryStats[b,2]
                         }else{
                         ##browser()
                             Syy <- self$summaryStats[b,2] - self$summaryStats[a,2]
                         }
                         ( n*log(2*pi) + Syy )
                     },
                     base2 = function(a,b){
                         -2 * sum(dnorm(self$summaryStats[a:b,3],0,1,log=TRUE))
                     },
                     point = function(a){
                         -2 * dnorm(self$summaryStats[a],0,abs(self$summaryStats[a]),log=TRUE) + self$betaP
                     }
                 ))

fCost <- gaussCost$new(x,3*log(length(x)),3*log(length(x)),10)


partition <- R6Class("partition",
                 public = list(
                     cost = 0,
                     ca = list(),
                     pa = list(),
                     initialize = function(){
                         invisible(self)
                     },
                     addCollective = function(x,s,e){
                         cst <- x$collective(s,e)
                         self$ca[[length(self$ca)+1]] <- c(start=s,end=e,cost=cst)
                         self$cost <- self$cost + cst
                         invisible(self)
                     },
                     addBase = function(x,s,e){
                         cst <- x$base(s,e)
                         self$cost <- self$cost + cst
                         invisible(self)
                     },
                     addPoint = function(x,a){
                         cst <- x$point(a)
                         self$pa[[length(self$pa)+1]] <- c(location=a, cost=cst)
                         self$cost <- self$cost + cst
                         invisible(self)                             
                     }
                 ))



p <- profvis::profvis({
ctlg <- list()
ctlg[[1]] <- partition$new() ## offset by 1 versus time!!
for(tt in 1:1000){##length(x)){
    if(tt %% 100==0) {
      # Print on the screen some message
      cat(paste0("time step: ", tt, "\n"))
    }
    
    ## compute C2 from paper
    opt <- ctlg[[tt]]$clone()$addBase(fCost,tt,tt)
    
    ## compute C3 from paper
    tmp <- ctlg[[tt]]$clone()$addPoint(fCost,tt)
    if(tmp$cost < opt$cost){ opt <- tmp }
    
    ## loop ctlg
    for(ii in 1:tt){
        tmp <- ctlg[[ii]]$clone() ## This is very slow.....
        tmp$addCollective(fCost,ii,tt)
        if(tmp$cost < opt$cost){ opt <- tmp }
    }
    
    ctlg[[tt+1]] <- opt  
}
})

htmlwidgets::saveWidget(p, "profile.html")
     
     # Can open in browser from R
browseURL("profile.html",browser="firefox")

## opt is the obect to look at....

## ############################################
## compare
##orig <- readRDS("test_output.rds")
##all(orig$ca==opt$ca)
##all(opt$pa==orig$pa)
       
