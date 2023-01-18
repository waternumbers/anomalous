#' An R implimentation of the segmented search algorithm
#' @param y univariate data series
#' @param segType type of segment
#' @param beta Penalisation term for a new segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#'@export
op <- function(y,segType,beta,min_length=2,max_length=.Machine$integer.max){
    
    n <- length(y)
    opt <- rep(list(NA),n) ## optimal partions
    cst <- rep(NA,n)
    ctlg <- rep(list(NA),n) ## catalog of partitions to try

    
    fC <- function(p){
        sum( sapply(p,cost) )
    }

    ## handle tt == 1
    ctlg[[1]] <- list( createSegment(segType,0,min_length,max_length) ) ## no penalty for the first segment
    ctlg[[1]][[1]] <- update(ctlg[[1]][[1]],y[1])
    opt[[1]] <- ctlg[[1]]
    cst[1] <- fC(ctlg[[1]])
    
    for(tt in 2:length(y)){
        cst[tt] <- Inf

        ## create new segment starting with current obs
        p <- opt[[tt-1]]
        p <- c(p, createSegment(segType,beta,min_length,max_length) )
        ctlg[[tt]] <- p

        ## loop to see which is optimal
        for(tau in 1:tt){
            p <- ctlg[[tau]]
            p[[length(p)]] <- update( p[[length(p)]], y[tt] )
            ctlg[[tau]] <- p
            
            tmp <- fC(p)
            if( tmp < cst[tt] ){
                opt[[tt]] <- ctlg[[tau]]
                cst[tt] <- tmp
            }
        }
    }
    return(opt[[tt]])
}

            
##         ## first case is to add a new segment
##         p <- ctlg[[tt-1]]
##         p <- c(p, createSegment(segType,beta,min_length,max_length) )
##         p[[length(p)]] <- update( p[[length(p)]], y[tt] )
##         F[tt] <- fC(p)
##         ctlg[[tt]] <- p

##         ## search impact of extending all previous groups
##         for(tau in 1:(tt-1)){
##             p <- ctlg[[tau]]
##             p[[length(p)]] <- update( p[[length(p)]], y[tt] )
##             tmp <- fC(p)
            
##             ctlg[[tau]] <- p
##             if( tmp < F[tt] ){
##                 #browser()
##                 F[tt] <- tmp
##                 ctlg[[tt]] <- p
##             }
##         }
##     }
##     browser()
##     ctlg <- unlist(tail(ctlg,1)) ## final classification
##     return(ctlg)
##     ##return( which(diff(ctlg)>0) )
## }



## op <- function(y,fC,Beta){
##     n <- length(y)
##     F <- rep(NA,n+1)
##     ctlg <- rep(list(NULL),n+1)
##     F[1] <- -Beta
##     for(tt in 1:n){
##         Fvec <- rep(NA,tt)
##         for(tau in 0:(tt-1)){
##             jj <- tau+1
##             Fvec[jj] <- F[jj] + fC(y[jj:tt]) + Beta
##         }
##         tauhat <- which.min(Fvec)-1
##         ctlg[[tt+1]] <- c(ctlg[[tauhat+1]],tauhat)
##         F[tt+1] <- Fvec[tauhat+1]
##     }
##     ctlg <- unlist(tail(ctlg,1))[-1] ## final value and trim intial 0
##     return(ctlg)
## }


