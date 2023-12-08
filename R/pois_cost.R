#' @export
poisCost <- R6Class("poisCost",
                     public=list(
                         summaryStats = NULL,
                         maxT = 0,
                         initialize = function(x,rate=1){
                             S <- cbind(x,lfactorial(x),rate,x*log(rate),1)
                             S[is.na(x),] <- NA
                             self$summaryStats <- apply(S,2,cumsumNA)
                             self$maxT <- length(x)
                             invisible(self)
                         },
                         isValid = function(a,b){
                             ## start and end must be valid time steps
                             !( any(is.na(self$summaryStats[b,])) | any(is.na(self$summaryStats[a,])) )
                         },
                         fixa = function(a){                             
                             ##If the period starts at a we want the last finite value of sumStats before it
                             a <- a-1
                             if(a>1){ while( any(is.na(self$summaryStats[a,])) & a>0 ){ a <- a-1 } }
                             a
                         },
                         baseCost = function(a,b,pen=0){
                             if( !self$isValid(a,b) ){ return(Inf) }
                             a <- self$fixa(a)
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             2*sumStat[3] - 2*sumStat[4] + 2*sumStat[2] + pen
                         },
                         pointCost = function(b,pen){
                             if( !self$isValid(b,b) ){ return(Inf) }
                             a <- self$fixa(b)
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             if( sumStat[1] == 0){
                                 xlgx <- 0
                             }else{
                                 xlgx <- sumStat[1]*log(sumStat[1])
                             }
                             2*(sumStat[1] - xlgx  + sumStat[2]) + pen
                         },
                         collectiveCost = function(a,b,pen){
                             if( !self$isValid(a,b) ){ return(Inf) }
                             a <- self$fixa(a)
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             rhat <- sumStat[1] / sumStat[3]
                             ##if(rhat < .Machine$double.eps){ return(0) }
                             rhat <- max(rhat, .Machine$double.eps)
                             
                             2*rhat*sumStat[3] - 2*sumStat[1]*log(rhat) - 2*sumStat[4] + 2*sumStat[2] + pen
                         }
                     )
                    )
