#' @export
gaussCost <- R6Class("gaussCost",
                     public=list(
                         summaryStats = NULL,
                         maxT = 0,
                         point_type = "var",
                         initialize = function(x,m=0,s=1,point_type=c("var","mean")){
                             self$point_type <- match.arg(point_type)
                             self$maxT <- length(x)
                             S <- cbind( 1/s,log(s),(x-m)/s,((x-m)^2)/s,1 )
                             S[is.na(x),] <- NA
                             self$summaryStats <- apply(S,2,cumsumNA)
                             invisible(self)
                         },
                         baseCost = function(a,b,pen=0){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             sumStat[5]*log(2*pi) + sumStat[2] + sumStat[4] + pen
                         },
                         pointCost = function(b,pen){
                             a <- b-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             gamma <- max(.Machine$double.xmin, exp(-(1+pen)))
                             mhat <- sumStat[3] / sumStat[1]
                             n <- 1
                             cst <- switch(
                                 self$point_type,
                                 "var" = log(2*pi) + sumStat[2] + log(gamma + sumStat[4]) + 1,
                                 "mean" = n*log(2*pi) + sumStat[2] + sumStat[4] - (mhat^2)*sumStat[1],
                                 stop("Unknown point anomaly type")
                             )
                             cst + pen
                         }
                     ),
                     private=list(
                         meanChange = function(a,b,pen,len){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             if( is.na(sumStat[5]) | sumStat[5]<len ){ return(NA) } ## check length and if NA
                             
                             mhat <- sumStat[3] / sumStat[1]
                             sumStat[5]*log(2*pi) + sumStat[2] + sumStat[4] - (mhat^2)*sumStat[1] + pen
                         },
                         varChange = function(a,b,pen,len){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             if( is.na(sumStat[5]) | sumStat[5]<len ){ return(NA) } ## check length and if NA
                             shat <- sumStat[4] / sumStat[5]
                             shat <- max(shat,.Machine$double.xmin)
                             sumStat[5]*log(2*pi*shat) + sumStat[2] + sumStat[5] + pen
                         },
                         meanVarChange = function(a,b,pen,len){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             if( is.na(sumStat[5]) | sumStat[5]<len ){ return(NA) } ## check length and if NA
                             mhat <- sumStat[3] / sumStat[1]
                             shat <- (sumStat[4] - (mhat^2)*sumStat[1])/sumStat[5]
                             shat <- max(shat,.Machine$double.xmin)
                             sumStat[5]*log(2*pi*shat) + sumStat[2] + sumStat[5] + pen
                         },
                         paramMean = function(a,b){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             sumStat[3] / sumStat[1]
                         },
                         paramVar = function(a,b){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             shat <- sumStat[4] / sumStat[5]
                             max(shat,.Machine$double.xmin)
                         },
                         paramMeanVar = function(a,b){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             mhat <- sumStat[3] / sumStat[1]
                             shat <- (sumStat[4] - (mhat^2)*sumStat[1])/sumStat[5]
                             shat <- max(shat,.Machine$double.xmin)
                             c(mhat,shat)
                         }

             

                     )
                     )
                     
#' @export
gaussMean <- R6Class("gaussMean",
                     inherit = gaussCost,
                     public = list(
                         collectiveCost = function(a,b,pen,len){ private$meanChange(a,b,pen,len) },
                         param = function(a,b){ private$paramMean(a,b) }
                     )
                     )

#' @export
gaussVar <- R6Class("gaussVar",
                    inherit = gaussCost,
                    public = list(
                        collectiveCost = function(a,b,pen,len){ private$varChange(a,b,pen,len) },
                        param = function(a,b){ private$paramVar(a,b) }
                    )
                    )

#' @export
gaussMeanVar <- R6Class("gaussMeanVar",
                        inherit = gaussCost,
                        public = list(
                            collectiveCost = function(a,b,pen,len){ private$meanVarChange(a,b,pen,len) },
                            param = function(a,b){ private$paramMeanVar(a,b) }
                        )
                        )



