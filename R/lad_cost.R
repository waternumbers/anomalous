#' @export
ladCost <- R6Class("ladCost",
                   public=list(
                       summaryStats = NULL,
                       maxT = 0,
                       tau = NA,
                       chck = function(x){x*(self$tau - (x<0))},
                       initialize = function(x,m=0,tau=0.5){
                           self$summaryStats <- x-m
                           self$maxT <- length(x)
                           self$tau <- tau
                           invisible(self)
                       },
                       baseCost = function(a,b,pen=0){
                           sumStat <- self$summaryStats[a:b]
                           sumStat <- sumStat[ is.finite(sumStat) ]
                           if( length(sumStat) == 0 ){ return(Inf) }
                           2*( sum(self$chck(sumStat)) + length(sumStat)*log(self$tau*(1-self$tau)) ) + pen
                       },
                       pointCost = function(a,pen){ 
                           sumStat <- self$summaryStats[a]
                           sumStat <- sumStat[ is.finite(sumStat) ]
                           if( length(sumStat) == 0 ){ return(Inf) }
                           theta <- sumStat
                           2*( sum(self$chck(sumStat-theta)) + length(sumStat)*log(self$tau*(1-self$tau)) ) + pen
                           
                       },
                       collectiveCost = function(a,b,pen,len){
                           if( (b-1+1)<len ){ return(NA) } ## check length and if NA
                           sumStat <- self$summaryStats[a:b]
                           sumStat <- sumStat[ is.finite(sumStat) ]
                           if( length(sumStat) == 0 ){ return(Inf) }
                           if( self$tau == 0.5 ){ theta <- fast_med(sumStat) } ##median(sumStat) }
                           else{ theta <- quantile(sumStat,self$tau) }
                           2*( sum(self$chck(sumStat-theta)) + length(sumStat)*log(self$tau*(1-self$tau)) ) + pen
                       },
                       param = function(a,b){
                           quantile( self$summaryStats[a:b],tau )
                       }
                       
                   )
                   )

