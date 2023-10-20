#' @export
gaussRepCost <- R6Class("gaussRepCost",
                     public=list(
                         summaryStats = NULL,
                         maxT = 0,
                         initialize = function(x,m=0,s=1){
                             if(is.list(x)){
                                 tmp <- matrix(NA,length(x),4)
                                 for(ii in 1:length(x)){
                                     tmp[ii,] <- c(length(x[[ii]]),
                                                   length(x[[ii]])/s[ii],
                                                   length(x[[ii]])*log(s[ii]),
                                                   sum(x[[ii]]-m[ii])/s[ii],
                                                   sum( (x[[ii]]-m[ii])^2 )/s[ii]
                                                   )
                                 }
                                 self$summaryStats <- apply(tmp,2,cumsum)
                             }else{
                                 self$summaryStats <- apply(cbind(1,1/s,log(s),(x-m)/s,((x-m)^2)/s),2,cumsum)
                             }
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
                             
                             sumStat[1]*log(2*pi) + sumStat[3] + sumStat[5] + pen
                         },
                         pointCost = function(a,pen){
                             private$meanVarChange(a,a,pen)
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
                             mhat <- sumStat[4] / sumStat[2]
                             sumStat[1]*log(2*pi) + sumStat[3] + sumStat[5] - (mhat^2)*sumStat[2] + pen
                         },
                         varChange = function(a,b,pen){
                             a <- a-1
                             n <- b-a
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             shat <- sumStat[5] / n ## TODO add catch for close to zero
                             shat <- max(shat,.Machine$double.xmin)
                             sumStat[1]*log(2*pi*shat) + sumStat[3] + sumStat[1] + pen
                         },
                         meanVarChange = function(a,b,pen){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             mhat <- sumStat[4] / sumStat[2]
                             shat <- (sumStat[5] - (mhat^2)*sumStat[2])/sumStat[1]
                             shat <- max(shat,.Machine$double.xmin)
                             
                             sumStat[1]*log(2*pi*shat) + sumStat[3] + sumStat[1] + pen
                         }
                     )
                     )
                     
#' @export
gaussRepMean <- R6Class("gaussMean",
                     inherit = gaussRepCost,
                     public = list(
                         collectiveCost = function(a,b,pen){ private$meanChange(a,b,pen) }
                     )
                     )

#' @export
gaussRepVar <- R6Class("gaussVar",
                    inherit = gaussRepCost,
                    public = list(
                        collectiveCost = function(a,b,pen){ private$varChange(a,b,pen) }
                    )
                    )

#' @export
gaussRepMeanVar <- R6Class("gaussMeanVar",
                        inherit = gaussRepCost,
                        public = list(
                            collectiveCost = function(a,b,pen){ private$meanVarChange(a,b,pen) }
                        )
                        )



