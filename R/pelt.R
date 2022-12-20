#' An R implimentation of the pelt agorithm
#' @param y univariate data series
#' @param fC cost function
#' @param Beta Penalisation term
#' @param min_length minimum length of a segment
#' 
#' @details Basic r implimentation of the PELT algorihm
#'@export
pelt <- function(y,fC,beta){
    n <- length(y)
    F <- rep(NA,n+1)
    cp <- rep(list(rep(NA,n)),n)
    
    ## handle tt ==1
    
    cp[[1]] <- fC$new(y)
    F[1] <- cp[[1]]$cost(beta) ##cost(cp[[1]],beta)  ##fC(y,cp[[1]],beta)

    ## form first catalog
    ctlg <- list(cp[[1]])
    
    
    for(tt in 2:n){
        ## update partitions
        
        ctlg <- c( lapply(ctlg, update,y=y) , list(update(ctlg[[1]],y,new=TRUE) ))
        
        ## evaluate the catalog
        Fvec <- sapply(ctlg,cost,beta=beta)
        
        ##browser()
        ## ## update catalog
        ## ctlg <- lapply(ctlg, function(x,ii){x[ii] <- x[ii-1]; return(x)},ii=tt) ## extend current periods
        ## ## new change point starting at current value
        ## m <- length(ctlg)+1
        ## ctlg[[ m ]] <- ctlg[[ 1 ]]
        ## ctlg[[m]][tt] <- ctlg[[m]][tt]+1

        
        ## evaluate catalog
        ## Fvec <- sapply(ctlg,fC, y=y,beta=beta) ##unction(x,y,beta){ return( fC(y,x,beta) ) }, y=y,beta=beta)
        idx <- order(Fvec) ## index of lowest value first
        
        ## copy min value over
        F[tt] <- Fvec[ idx[1] ]
        cp[[tt]] <- ctlg[[ idx[1] ]]

        ## trim catalog
        idx <- idx[ (Fvec[idx] <= F[tt]+beta) | is.na(Fvec[idx]) ]
        ctlg <- ctlg[ idx ]
        
        ##print(paste(tt,length(ctlg)))
    }
    
    cp <- tail(cp,1)[[1]] ## final classification
    return(cp)
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
