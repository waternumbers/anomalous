#' An R implimentation of the pelt agorithm
#' @param y univariate data series
#' @param fC cost function
#' @param Beta Penalisation term
#' @param min_length minimum length of a segment
#' 
#' @details Basic r implimentation of the PELT algorihm
#'@export
peltR <- function(y,fC,Beta,min_length=2){
    n <- length(y)
    n <- length(y)
    F <- rep(NA,n+1)
    cp <- rep(list(NULL),n+1)
    F[1] <- -Beta
    R <- 0
    Fvec <- F
    for(tt in min_length:n){
        ##print(tt)
        Fvec[] <- Inf
        for(tau in R){#0:(tt-1)){
            jj <- tau+1
            Fvec[jj] <- F[jj] + fC(y[jj:tt]) + Beta
        }
        tauhat <- which.min(Fvec[1:(tt-min_length+1)])-1
        cp[[tt+1]] <- c(cp[[tauhat+1]],tauhat)
        F[tt+1] <- Fvec[tauhat+1]
        
        idx <- which( Fvec <= F[tt+1]+Beta )
        R <- c(idx-1,tt)
        
    }
    #browser()
    cp <- unlist(tail(cp,1))[-1] ## final value and trim intial 0
    return(cp)
}
