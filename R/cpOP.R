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
cpOP <- function(y,mu,sigma,cpType,beta,min_length=2,max_length=.Machine$integer.max){

    
    n <- length(y)
    opt <- rep(list(NA),n) ## optimal partions
    cst <- rep(NA,n)
    ctlg <- list() ## catalog of partitions to try
    state <- NULL

    ctlg[[1]] <- update(gaussPartition(beta,min_length,max_length),y[1],mu[1],sigma[1],1,cpType)
    opt[[1]] <- ctlg[[1]]
    cst[1] <- ctlg[[1]]@cost
            
    for(tt in 2:n){

        ctlg[[tt]] <- update(opt[[tt-1]],y[tt],mu[tt],sigma[tt],tt,cpType)
        opt[[tt]] <- ctlg[[tt]]
        cst[tt] <- ctlg[[tt]]@cost

        for(tau in 1:(tt-1)){
            
            p <- update(ctlg[[tau]], y[tt],mu[tt],sigma[tt],tt,"current")
            ctlg[[tau]] <- p

            
            if( p@cost < cst[tt] ){
                opt[[tt]] <- p
                cst[tt] <- p@cost
            }
        }
    }

    return(opt[[tt]])
}
