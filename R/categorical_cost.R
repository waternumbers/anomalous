#' @export
categoricalCost <- R6Class("categoricalCost",
                     public=list(
                         summaryStats = NULL,
                         m=NULL,
                         maxT = 0,
                         initialize = function(x,m=rep(1/ncol(x),ncol(x))){
                             self$maxT <- nrow(x)
                             S <- x
                             S[is.na(rowSums(S)),] <- NA
                             self$summaryStats <- apply(S,2,cumsumNA)
                             nS <- rep(1,nrow(S))
                             nS[is.na(S[,1])] <- NA
                             self$nStep <- cumsumNA(nS)
                             self$m <- log(m) #matrix(m,nrow(S),ncol(S),byrow=T)
                             invisible(self)
                         },
                         baseCost = function(a,b,pen=0){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             n <- sum(sumStat)
                             -2*sum( sumStat*self$m ) + pen
                         },
                         pointCost = function(b,pen){
                             self$collectiveCost(b,b,pen)
                         },
                         collectiveCost = function(a,b,pen,len){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                                 nS <- self$nStep[b]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                                 nS <- self$nStep[b] - self$nStep[a] + 1
                             }
                             if( is.na(nS) | nS<len ){ return(NA) } ## check length and if NA
                             
                             n <- sum(sumStat)
                             lambda <- sumStat/n
                             lambda[lambda==0] <- 1e-10 ## this is to get finite logs which are multiplied by 0
                             -2*sum( sumStat*log(lambda) ) + pen
                         },
                         param = function(a,b,type){
                             if( type=="background" ){ return(self$m) }
                            
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             n <- sum(sumStat)
                             sumStat/n
                         }
                     ))

                     
