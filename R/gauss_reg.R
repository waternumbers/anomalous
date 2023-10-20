#' @export
gaussRegCost <- R6Class("gaussRegCost",
                     public=list(
                         summaryStats = NULL,
                         maxT = 0,
                         theta0 = NA,
                         sigma0 = NA,
                         initialize = function(x,theta0=NULL,sigma0=1){
                             self$summaryStats <- list(XtX=list(),Xty=list(),
                                                       yty=rep(NA,length(x)),
                                                       n=rep(NA,length(x)))
                             XtX <- matrix(0,ncol(x[[1]]$X),ncol(x[[1]]$X))
                             Xty <- matrix(0,ncol(x[[1]]$X),1)
                             yty <- 0
                             n <- 0
                             for(ii in 1:length(x)){
                                 XtX <- XtX + crossprod(x[[ii]]$X)
                                 Xty <- Xty + crossprod(x[[ii]]$X,x[[ii]]$y)
                                 n <- n + length(x[[ii]]$y)
                                 yty <- yty + sum(x[[ii]]$y^2)
                                 self$summaryStats$XtX[[ii]] <- XtX
                                 self$summaryStats$Xty[[ii]] <- Xty
                                 self$summaryStats$n[ii] <- n
                                 self$summaryStats$yty[ii] <- yty
                             }
                             self$maxT <- length(x)
                             if(is.null(theta0)){ theta0 <- rep(0,ncol(XtX)) }
                             if(length(theta0) != ncol(XtX)){ stop("theta0 incorrect size") }
                             self$theta0 <- theta0
                             self$sigma0 <- sigma0
                             invisible(self)
                         },
                         baseCost = function(a,b,pen=0){
                             a <- a-1
                             if(a<1){
                                 XtX <- self$summaryStats$XtX[[b]]
                                 Xty <- self$summaryStats$Xty[[b]]
                                 n <- self$summaryStats$n[b]
                                 yty <- self$summaryStats$yty[b]
                             }else{
                                 XtX <- self$summaryStats$XtX[[b]] - self$summaryStats$XtX[[a]]
                                 Xty <- self$summaryStats$Xty[[b]] - self$summaryStats$Xty[[a]]
                                 n <- self$summaryStats$n[b] - self$summaryStats$n[a]
                                 yty <- self$summaryStats$yty[b] - self$summaryStats$yty[a]
                             }
                             sigma <- self$sigma0
                             theta <- self$theta0
                             as.numeric( n*log(2*pi*sigma) + (1/sigma)*( yty - 2*(t(theta)%*%Xty) + t(theta)%*%XtX%*%theta ) +pen )
                         },
                         pointCost = function(a,pen){
                             self$collectiveCost(a,a,pen)
                         },
                         collectiveCost = function(a,b,pen){
                             p <- self$param(a,b)
                             as.numeric( p$n*log(2*pi*p$sigma) + p$n + pen)
                         },
                         param = function(a,b){
                             a <- a-1
                             if(a<1){
                                 XtX <- self$summaryStats$XtX[[b]]
                                 Xty <- self$summaryStats$Xty[[b]]
                                 n <- self$summaryStats$n[b]
                                 yty <- self$summaryStats$yty[b]
                             }else{
                                 XtX <- self$summaryStats$XtX[[b]] - self$summaryStats$XtX[[a]]
                                 Xty <- self$summaryStats$Xty[[b]] - self$summaryStats$Xty[[a]]
                                 n <- self$summaryStats$n[b] - self$summaryStats$n[a]
                                 yty <- self$summaryStats$yty[b] - self$summaryStats$yty[a]
                             }

                             ## Check for non singular regressio
                             QR <- qr(XtX)                 # Get the QR decomposition
                             vars <- QR$pivot[seq_len(QR$rank)]     # Variable numbers that are OK to include
                             ## solve- set other values to 0
                             theta <- rep(0,ncol(XtX))
                             theta[vars] <- solve(XtX[vars,vars],Xty[vars])
                             sigma <- ( yty - 2*(t(theta)%*%Xty) + t(theta)%*%XtX%*%theta ) / n
                             list(theta=theta, sigma=sigma, n=n)
                         }
                     )
                     )

