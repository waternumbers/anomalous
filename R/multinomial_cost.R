#' @export
multinomialCost <- R6Class("multinomialCost",
                     public=list(
                         summaryStats = NULL,
                         m=NULL,
                         maxT = 0,
                         initialize = function(x,m=rep(1/ncol(x),ncol(x))){
                             self$maxT <- nrow(x)
                             S <- x
                             S[is.na(rowSums(S)),] <- NA
                             self$summaryStats <- apply(S,2,cumsumNA)
                             self$m <- matrix(m,nrow(S),ncol(S),byrow=T)
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
                             -2*( lfactorial(n) + sum( sumStat*log(m) - lfactorial(sumStat) ) ) + pen
                         },
                         pointCost = function(b,pen){
                             a <- b-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             n <- sum(sumStat)
                             m <- sumStat/n
                             m[m==0] <- 1e-10 ## this is to get finite logs which are multiplied by 0
                             -2*( lfactorial(n) + sum( sumStat*log(m) - lfactorial(sumStat) ) ) + pen
                         },
                         collectiveCost = function(a,b,pen=0){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             n <- sum(sumStat)
                             m <- sumStat/n
                             m[m==0] <- 1e-10 ## this is to get finite logs which are multiplied by 0
                             -2*( lfactorial(n) + sum( sumStat*log(m) - lfactorial(sumStat) ) ) + pen
                         },
                         param = function(a,b){
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

                     
