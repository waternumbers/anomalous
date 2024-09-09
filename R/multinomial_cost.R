#' @export
multinomialCost <- R6Class("multinomialCost",
                     public=list(
                         summaryStats = NULL,
                         summaryStatsSlow = NULL,
                         nStep = NULL,
                         m=NULL,
                         maxT = 0,
                         initialize = function(x,m=rep(1/ncol(x),ncol(x))){
                             self$maxT <- nrow(x)
                             S <- cbind(x,
                                        lfactorial(rowSums(x)),
                                        rowSums(lfactorial(x)),
                                        1 )
                                        #                             S <- x
                             S[is.na(rowSums(S)),] <- NA
                             self$summaryStats <- apply(S,2,cumsumNA)
                             self$summaryStatsSlow <- x
                             nS <- rep(1,nrow(S))
                             nS[is.na(S[,1])] <- NA
                             self$nStep <- cumsumNA(nS)
                             self$m <- m #matrix(m,nrow(S),ncol(S),byrow=T)
                             invisible(self)
                         },
                         baseCost = function(a,b,pen=0){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             nm <- length(self$m)
                             -2*( sumStat[nm+1] - sumStat[nm+2] + sum(sumStat[1:nm]*log(self$m)) ) + pen
                         },
                         baseCostSlow = function(a,b,pen=0){
                             nm <- length(self$m)
                             lp <- sapply(a:b,function(i){dmultinom(self$summaryStatsSlow[i,1:nm], prob=self$m,log=T)})
                             -2*sum(lp) + pen
                         },
                         pointCost = function(b,pen){
                             ## a <- b-1
                             ## if(a<1){
                             ##     sumStat <- self$summaryStats[b,]
                             ## }else{
                             ##     sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             ## }
                             ## n <- sum(sumStat)
                             ## m <- sumStat/n
                             ## m[m==0] <- 1e-10 ## this is to get finite logs which are multiplied by 0
                             ## -2*( lfactorial(n) + sum( sumStat*log(m) - lfactorial(sumStat) ) ) + pen
                             self$collectiveCost(b,b,pen,1)
                         },
                         collectiveCost = function(a,b,pen,len){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             nm <- length(self$m)
                             if( is.na(sumStat[nm+3]) | sumStat[nm+3]<len ){ return(NA) } ## check length and if NA
                             m <- sumStat[1:nm]/sum(sumStat[1:nm])
                             idx <- (1:nm)[m>0]
                             -2*( sumStat[nm+1] - sumStat[nm+2] + sum(sumStat[idx]*log(m[idx])) ) + pen
                             
                         },
                         collectiveCostSlow = function(a,b,pen,len){
                             m <- colSums(self$summaryStatsSlow[a:b,,drop=F],na.rm=TRUE)
                             m <- m/sum(m)
                             lp <- sapply(a:b,function(i){dmultinom(self$summaryStatsSlow[i,], prob=m,log=T)})
                             -2*sum(lp) + pen
                         },
                         param = function(a,b){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             nm <- length(self$m)
                             sumStat[1:nm]/sum(sumStat[1:nm])
                             
                             ## a <- a-1
                             ## if(a<1){
                             ##     sumStat <- self$summaryStats[b,]
                             ## }else{
                             ##     sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             ## }
                             ## n <- sum(sumStat)
                             ## sumStat/n
                         }
                     ))

                     
