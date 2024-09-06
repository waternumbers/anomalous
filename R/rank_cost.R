#' @export
rankCost <- R6Class("rankCost",
                    public=list(
                        summaryStats = NULL,
                        nStep = NULL,
                        Sr = NULL,
                        m=NULL,
                        maxT = 0,
                        initialize = function(x,m=0){
                            x <- as.matrix(x)
                            x <- x-m
                            self$maxT <- nrow(x)
                            S <- apply(x,2,rank,na.last="keep") - (nrow(x)+1)/2
                            self$summaryStats <- apply(S,2,cumsumNA)
                            self$nStep <- apply(is.finite(S),2,cumsumNA)
                            self$Sr <- var(S,na.rm=T)
                            invisible(self)
                        },
                        baseCost = function(a,b,pen=0){
                            self$collectiveCost(a,b,pen,1)
                        },
                        pointCost = function(b,pen){
                            self$collectiveCost(b,b,pen,1)
                        },
                        collectiveCost = function(a,b,pen,len){
                            a <- a-1
                            if(a<1){
                                sumStat <- self$summaryStats[b,]
                                nS <- self$nStep[b,]
                            }else{
                                sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                                nS <- self$nStep[b,] - self$nStep[a,]
                            }
                            ##browser()
                            if( any(is.na(nS) | nS<len) ){ return(NA) } ## check length and if NA
                            
                            rbar <- sumStat/nS
                            idx <- is.finite(rbar)
                            if(!any(idx)){ return(NA) }

                            -(b-a-1)* rbar %*% solve(self$Sr,rbar) + pen
                        },
                        param = function(a,b,type){ return(NULL) }
                    ))

                     
