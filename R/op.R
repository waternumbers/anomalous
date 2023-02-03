#' An R implimentation of the segmented search algorithm
#' @param y univariate data series
#' @param segType function to generate an object of class Segment
#' @param beta Penalisation term for a new segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#'@export
op <- function(y,mu,sigma,segType,beta,min_length=2,max_length=.Machine$integer.max){
    ##browser()
    n <- length(y)
    opt <- rep(list(NULL),n) ## optimal partions
    cst <- rep(NA,n)
    ctlg <- rep(list(NULL),n) ## catalog of partitions to try

    ## handle tt == 1
    ctlg[[1]] <- partition(min_length,max_length)
    ctlg[[1]] <- update(ctlg[[1]],y[1],mu[1],sigma[1],segType(beta,1))
    opt[[1]] <- ctlg[[1]]
    cst[1] <- ctlg[[1]]@cost

    ## loop time
    for(tt in 2:n){

        ctlg[[tt]] <- update(opt[[tt-1]], y[tt],mu[tt],sigma[tt], segType(beta,tt))
        opt[[tt]] <- ctlg[[tt]]
        cst[tt] <- ctlg[[tt]]@cost

        for(tau in 1:(tt-1)){

            ctlg[[tau]] <- update( ctlg[[tau]], y[tt], mu[tt],sigma[tt] )

            tmp <- ctlg[[tau]]@cost
            if( tmp < cst[tt] ){
                opt[[tt]] <- ctlg[[tau]]
                cst[tt] <- tmp
            }
        }
    }
    
    return(opt[[tt]])
}
