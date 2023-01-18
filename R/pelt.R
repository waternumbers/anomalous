#' An R implimentation of the pelt agorithm
#' @param y univariate data series
#' @param segType type of segment
#' @param beta Penalisation term
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation of the PELT algorihm
#'@export
pelt <- function(y,segType,beta,min_length=2,max_length=.Machine$integer.max){
    
    n <- length(y)
    opt <- rep(list(NA),n) ## optimal partions
    cst <- rep(NA,n)
    ctlg <- list(NULL) ## catalog of partitions to try
    
    fC <- function(p){
        sum( sapply(p,cost) )
    }

    fN <- function(p){
        p[[length(p)]]@n < p[[length(p)]]@min_length
    }
    
    
    ## handle tt ==1
    ctlg[[1]] <- list( createSegment(segType,0,min_length,max_length) ) ## no penalty for the first segment
    ctlg[[1]][[1]] <- update(ctlg[[1]][[1]],y[1])
    opt[[1]] <- ctlg[[1]]
    cst[1] <- fC(ctlg[[1]])
    
    
    for(tt in 2:n){
        
        ## append to current catalog
        p <- opt[[tt-1]]
        p <- c(p, createSegment(segType,beta,min_length,max_length) )
        ctlg[[ length(ctlg)+1 ]] <- p

        print(paste(tt,length(ctlg)))
        ## update catalog
        for(tau in 1:length(ctlg)){
            p <- ctlg[[tau]]
            p[[length(p)]] <- update( p[[length(p)]], y[tt] )
            ctlg[[tau]] <- p
        }
        

        ## evaluate the catalog
        ##browser()
        Cvec <- sapply(ctlg,fC)
        Nvec <- sapply(ctlg,fN)

        ##browser()
        ## ## update catalog
        ## ctlg <- lapply(ctlg, function(x,ii){x[ii] <- x[ii-1]; return(x)},ii=tt) ## extend current periods
        ## ## new change point starting at current value
        ## m <- length(ctlg)+1
        ## ctlg[[ m ]] <- ctlg[[ 1 ]]
        ## ctlg[[m]][tt] <- ctlg[[m]][tt]+1

        
        ## evaluate catalog
        idx <- order(Cvec) ## index of lowest value first
        
        ## copy min value over
        opt[[tt]] <- ctlg[[ idx[1] ]]
        cst[tt] <- Cvec[ idx[1] ]

        ## trim catalog
        idx <- idx[ (Cvec[idx] <= cst[tt]+beta) | Nvec ]
        ctlg <- ctlg[ idx ]
        
        ##print(paste(tt,length(ctlg)))
    }

    return(opt[[tt]])
    ##return( which(diff(cp)>0) )
}

## pelt <- function(y,fC,beta){
##     n <- length(y)
##     F <- rep(NA,n+1)
##     cp <- rep(list(rep(NA,n)),n)
    
##     ## handle tt ==1
##     cp[[1]][1] <- 1
##     F[1] <- fC(y,cp[[1]],beta)
    
##     Fvec <- rep(Inf,n)
##     R <- 1
    
##     for(tt in 2:n){
##         Fvec[] <- Inf
        
##         ## first case is to keep the same
##         cl <- cp[[tt-1]]
##         cl[tt] <- cl[tt-1]
##         Fvec[tt] <- fC(y,cl,beta)
##         ##F[tt] <- fC(y,cl,beta)
##         ##cp[[tt]] <- cl
        
##         ## search impact of extending all previous groups
##         for(tau in R){
##             cl <- cp[[tau]]
##             cl[(tau+1):tt] <- cl[tau]+1
##             Fvec[tau] <- fC(y,cl,beta)
##             ##if( Fvec[tau] < F[tt] ){
##                 #browser()
##             ##    F[tt] <- Fvec[tau]
##             ##    cp[[tt]] <- cl
##             ##}
##         }
##         ## work out best estimate
##         tau <- which.min(Fvec) ## [1:(tt-min_length+1)])-1
        
##         F[tt] <- Fvec[tau]
##         if(tau==tt){
##             cl <- cp[[tt-1]]
##             cl[tt] <- cl[tt-1]
##         }else{
##             cl <- cp[[tau]]
##             cl[(tau+1):tt] <- cl[tau]+1
##         }
##         cp[[tt]] <- cl

##         R <- which( Fvec <= F[tt]+beta )
##         print(paste(tt,length(R)))
##     }
##     cp <- unlist(tail(cp,1)) ## final classification
##     return(cp)
##     ##return( which(diff(cp)>0) )
## }

## pelt <- function(y,fC,Beta,min_length=2){
##     n <- length(y)
##     n <- length(y)
##     F <- rep(NA,n+1)
##     cp <- rep(list(NULL),n+1)
##     F[1] <- -Beta
##     R <- 0
##     Fvec <- F
##     for(tt in min_length:n){
##         ##print(tt)
##         Fvec[] <- Inf
##         for(tau in R){#0:(tt-1)){
##             jj <- tau+1
##             Fvec[jj] <- F[jj] + fC(y[jj:tt]) + Beta
##         }
##         tauhat <- which.min(Fvec[1:(tt-min_length+1)])-1
##         cp[[tt+1]] <- c(cp[[tauhat+1]],tauhat)
##         F[tt+1] <- Fvec[tauhat+1]
        
##         idx <- which( Fvec <= F[tt+1]+Beta )
##         R <- c(idx-1,tt)
        
##     }
##     #browser()
##     cp <- unlist(tail(cp,1))[-1] ## final value and trim intial 0
##     return(cp)
## }
