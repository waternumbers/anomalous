#' An R implimentation of the pelt agorithm for partition structures
#' @param y univariate data series
#' @param mu background mean
#' @param sigma background variance
#' @param segType defined types of segment allowed
#' @param beta Penalisation term for each type of segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation of the PELT algorihm
#'@export
cpPelt <- function(y,mu,sigma,cpType,beta,min_length=2,max_length=.Machine$integer.max){

    
    n <- length(y)
    opt <- rep(list(NA),n) ## optimal partions
    cst <- rep(NA,n)
    ctlg <- list() ## catalog of partitions to try
    state <- NULL

    ## initial timestep
    ctlg[[1]] <- update(gaussPartition(beta,min_length,max_length),y[1],mu[1],sigma[1],1,cpType)
    opt[[1]] <- ctlg[[1]]
    cst[1] <- ctlg[[1]]@cost
            
    for(tt in 2:n){

        ## update catalog
        for(tau in 1:length(ctlg)){
            ctlg[[tau]] <- update(ctlg[[tau]],y[tt],mu[tt],sigma[tt],tt,"current")
        }

        ## append new break to catalog
        p <- opt[[tt-1]]
        p <- update(p,y[tt],mu[tt],sigma[tt],tt,cpType)
        ctlg[[ length(ctlg)+1 ]] <- p

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
        idx <- (cstVec <= cst[tt]+beta[cpType]) | !isValid 
        ctlg <- ctlg[ idx ]
    }
    
    return(opt[[tt]])
}
