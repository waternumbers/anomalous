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

    ctlg <- list(NULL) ## catalog of partitions to try

    ## initialise
    p <- partition() ##min_length,max_length)
    for(tt in 1:min_length){
        if(tt==1){ seg <- segType(beta,1) }else{seg <- NULL}
       
        p <- update(p,y[1],mu[1],sigma[1],seg)
    }

    ctlg[[1]] <- p
    opt[[min_length]] <- p
    
    ## loop time
    for(tt in (min_length+1):n){
        last_n <- NULL
        cst <- NULL
        
        ## update catalog
        for(tau in 1:length(ctlg)){
            ctlg[[tau]] <- update( ctlg[[tau]], y[tt], mu[tt], sigma[tt] )
            last_n <- c(last_n ,ctlg[[tau]]@last_n )
            cst <- c(cst, ctlg[[tau]]@cost)
        }

        ## append new break to catalog
        nc <- length(ctlg)+1
        ctlg[[ nc ]] <- update( opt[[tt-1]],y[tt], mu[tt], sigma[tt], segType(beta,tt) )
        last_n <- c(last_n ,ctlg[[nc]]@last_n )
        cst <- c(cst, ctlg[[nc]]@cost)

        ## find minimum
        not_yet_valid <- min_length > last_n
        cst[ not_yet_valid | max_length < last_n ] <- Inf
               
        ii <- which.min(cst) ## the optimal choice
        

        ## copy min value over
        opt[[tt]] <- ctlg[[ii]]

        ## trim catalog
        idx <- (cst <= ctlg[[ii]]@cost + beta) | not_yet_valid
        ctlg <- ctlg[ idx ]
    }
    
    return(opt[[tt]])
}
