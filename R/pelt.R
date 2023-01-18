#' An R implimentation of the pelt agorithm
#' @param y univariate data series
#' @param segType function to generate an object of class Segment
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
        ## for loop appears quicker then sum(sapply(...))
        out <- 0
        for(ii in 1:length(p)){out <- out + p[[ii]]@cost}
        return(out)
    }

    fN <- function(p){
        p[[length(p)]]@n < p[[length(p)]]@min_length
    }
    
    
    ## handle tt ==1
    ctlg[[1]] <- list( segType(y[1],0,min_length,max_length) )
    opt[[1]] <- ctlg[[1]]
    cst[1] <- fC(ctlg[[1]])
    
    ## loop time
    for(tt in 2:n){

        ## update catalog
        for(tau in 1:length(ctlg)){
            p <- ctlg[[tau]]
            p[[length(p)]] <- update( p[[length(p)]], y[tt] )
            ctlg[[tau]] <- p
        }

        ## append new break to catalog
        p <- opt[[tt-1]]
        p <- c(p, segType(y[tt],beta,min_length,max_length) )
        ctlg[[ length(ctlg)+1 ]] <- p

        ## evaluate the catalog
        Cvec <- sapply(ctlg,fC)
        Nvec <- sapply(ctlg,fN)
        
        ## evaluate catalog
        idx <- order(Cvec) ## index of lowest value first
        
        ## copy min value over
        opt[[tt]] <- ctlg[[ idx[1] ]]
        cst[tt] <- Cvec[ idx[1] ]

        ## trim catalog
        idx <- idx[ (Cvec[idx] <= cst[tt]+beta) | Nvec ]
        ctlg <- ctlg[ idx ]
        

    }

    return(opt[[tt]])
}
