#' An R implimentation of the pelt agorithm
#' @param y univariate data series
#' @param segType function to generate an object of class Segment
#' @param beta Penalisation term
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation of the PELT algorihm
#'@export
pelt <- function(y,mu,sigma,segType,beta,min_length=2,max_length=.Machine$integer.max){
    
    n <- length(y)
    opt <- rep(list(NULL),n) ## optimal partions
    cst <- rep(NA,n)
    ctlg <- list(NULL) ## catalog of partitions to try
    
    ## handle tt ==1
    ctlg[[1]] <- partition(min_length,max_length)
    ctlg[[1]] <- update(ctlg[[1]],y[1],mu[1],sigma[1],segType(beta,1))
    opt[[1]] <- ctlg[[1]]
    cst[1] <- ctlg[[1]]@cost
    
    ## loop time
    for(tt in 2:n){

        ## update catalog
        for(tau in 1:length(ctlg)){
            ctlg[[tau]] <- update( ctlg[[tau]], y[tt], mu[tt], sigma[tt] )
        }

        ## append new break to catalog
        ctlg[[ length(ctlg)+1 ]] <- update( opt[[tt-1]],y[tt], mu[tt], sigma[tt], segType(beta,tt) )

        ## compute costs and validitiy
        cstVec <- isValid <- rep(NA,length(ctlg))
        for(ii in 1:length(ctlg)){
            cstVec[ii] <- ctlg[[ii]]@cost
            isValid[ii] <- ctlg[[ii]]@is_valid
        }

        ## find minimum
        idx <- which.min(cstVec)

        ## copy min value over
        opt[[tt]] <- ctlg[[ idx ]]
        cst[tt] <- cstVec[ idx ]

        ## trim catalog
        idx <- (cstVec <= cst[tt]+beta) | !isValid 
        ctlg <- ctlg[ idx ]
    }
    
    return(opt[[tt]])
}
