#' An R implimentation of the pelt agorithm
#' @param y univariate data series
#' @param baseType function to generate an object of class Segment for base distribution
#' @param anomType function to generate an object of class Segment for anomalies
#' @param beta penalisation term for collective anomaly segments
#' @param beta_point penalisation term for point anomaly segments
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation of the PELT algorihm
#'@export
pelt <- function(y,baseType,anomType,beta,beta_point,min_length=2,max_length=.Machine$integer.max){
    
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

    ## work out if should be kept due to short segment length
    fN <- function(p){
        ii <- length(p)
        while( inherits(p[[ii]],"pointSegment") & ii > 0){
            ii <- ii-1
        }
        if(ii==0){ stop("Whoa shouldn't get here...") } ## all segments are point anonalies Ahhh TODO: tidy up
        else{ out <- (p[[ii]]@n < p[[ii]]@min_length) |
                  (p[[1]]@n < p[[1]]@min_length) }## since p[[1]] is base case 

        return(out)
    }
    
    
    ## handle tt ==1 - assume base case TODO - change so not the case
    ctlg[[1]] <- list( baseType(y[1],0,min_length,max_length) )
    opt[[1]] <- ctlg[[1]]
    cst[1] <- fC(ctlg[[1]])
    
    ## loop time
    for(tt in 2:n){

        ## update catalog adding to last segment
        updatedBase <- FALSE
        for(tau in 1:length(ctlg)){
            p <- ctlg[[tau]]
            if( length(p) == 1){ updatedBase <- TRUE }
            p[[length(p)]] <- update( p[[length(p)]], y[tt] )
            ctlg[[tau]] <- p
        }

        ## append new possible steps to catalog
        p <- opt[[tt-1]]
        p <- c(p, anomType(y[tt],beta,min_length,max_length) )
        ctlg[[ length(ctlg)+1 ]] <- p
        
        p <- opt[[tt-1]]
        p <- c(p, pointSegment(y[tt],beta_point) )
        ctlg[[ length(ctlg)+1 ]] <- p

        if(!updateBase){
            p <- opt[[tt-1]]
            p[[1]] <- update( p[[1]], y[tt] )
            ctlg[[ length(ctlg)+1 ]] <- p
        }
        

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
