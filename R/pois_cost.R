#' @export
poisCost <- R6Class("poisCost",
                     public=list(
                         summaryStats = NULL,
                         maxT = 0,
                         rate=NA,
                         initialize = function(x,rate=1){
                             self$summaryStats <- apply( cbind(x,lfactorial(x)), 2, cumsum )
                             self$maxT <- length(x)
                             self$rate <- rate 
                             invisible(self)
                         },
                         baseCost = function(a,b,pen=0){
                             a <- a-1
                             n <- b-a
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             2*n*self$rate - 2*sumStat[1]*log(self$rate) + 2*sumStat[2] + pen
                         },
                         pointCost = function(a,pen){
                             stop("Not implimented")
                         },
                         collectiveCost = function(a,b,pen){
                             a <- a-1
                             n <- b-a
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             rhat <- sumStat[1] / n
                             if(rhat < .Machine$double.eps){ return(0) }
                             2*n*rhat - 2*sumStat[1]*log(rhat) + 2*sumStat[2] + pen
                         }
                     )
                    )
