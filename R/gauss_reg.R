#' @export
gaussCost <- R6Class("gaussCost",
                     public=list(
                         summaryStats = NULL,
                         maxT = 0,
                         initialize = function(x,m=0,s=1){
                             self$summaryStats <- apply(cbind(1/s,log(s),(x-m)/s,((x-m)^2)/s),2,cumsum)
                             self$maxT <- length(x)
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
                             n*log(2*pi) + sumStat[2] + sumStat[4] + pen
                         },
                         pointCost = function(a,pen){
                             if(a<2){
                                 sumStat <- self$summaryStats[a,]
                             }else{
                                 sumStat <- self$summaryStats[a,] - self$summaryStats[a-1,]
                             }
                             gamma <- max(.Machine$double.xmin, exp(-(1+pen)))
                             log(2*pi) + sumStat[2] + log(gamma + sumStat[4]) + 1 + pen
                         }
                     ),
                     private=list(
                         meanChange = function(a,b,pen){
                             a <- a-1
                             n <- b-a
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             mhat <- sumStat[3] / sumStat[1]
                             n*log(2*pi) + sumStat[2] + sumStat[4] - (mhat^2)*sumStat[1] + pen
                         },
                         varChange = function(a,b,pen){
                             a <- a-1
                             n <- b-a
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             shat <- sumStat[4] / n ## TODO add catch for close to zero
                             shat <- max(shat,.Machine$double.xmin)
                             n*log(2*pi*shat) + sumStat[2] + n + pen
                         },
                         meanVarChange = function(a,b,pen){
                             a <- a-1
                             n <- b-a
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             mhat <- sumStat[3] / sumStat[1]
                             shat <- (sumStat[4] - (mhat^2)*sumStat[1])/n
                             shat <- max(shat,.Machine$double.xmin)
                             n*log(2*pi*shat) + sumStat[2] + n + pen
                         }
                     )
                     )
                     
#' @export
gaussMean <- R6Class("gaussMean",
                     inherit = gaussCost,
                     public = list(
                         collectiveCost = function(a,b,pen){ private$meanChange(a,b,pen) }
                     )
                     )

#' @export
gaussVar <- R6Class("gaussVar",
                    inherit = gaussCost,
                    public = list(
                        collectiveCost = function(a,b,pen){ private$varChange(a,b,pen) }
                    )
                    )

#' @export
gaussMeanVar <- R6Class("gaussMeanVar",
                        inherit = gaussCost,
                        public = list(
                            collectiveCost = function(a,b,pen){ private$meanVarChange(a,b,pen) }
                        )
                        )



