#' An R implimentation of the segmented search algorithm
#' @param y univariate data series
#' @param segType function to generate an object of class Segment
#' @param beta Penalisation term for a new segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#'@export
op <- function(y,segType,beta,min_length=2,max_length=.Machine$integer.max){
    ##browser()
    n <- length(y)
    opt <- rep(list(NA),n) ## optimal partions
    cst <- rep(NA,n)
    ctlg <- rep(list(NA),n) ## catalog of partitions to try

    
    fC <- function(p){
        out <- 0
        ## for loop appears quicker then sum(sapply(...))
        for(ii in 1:length(p)){out <- out + p[[ii]]@cost}
        return(out)
    }

    ## handle tt == 1
    ctlg[[1]] <- list( segType(y[1],0,min_length,max_length) )
    opt[[1]] <- ctlg[[1]]
    cst[1] <- fC(ctlg[[1]])

    ## loop time
    for(tt in 2:length(y)){

        ctlg[[tt]] <- c(opt[[tt-1]], segType(y[tt],beta,min_length,max_length))
        opt[[tt]] <- ctlg[[tt]]
        cst[tt] <- fC(ctlg[[tt]])

        for(tau in 1:(tt-1)){
            
            p <- ctlg[[tau]]
            p[[length(p)]] <- update( p[[length(p)]], y[tt] )
            ctlg[[tau]] <- p
            
            tmp <- fC(p)
            if( tmp < cst[tt] ){
                opt[[tt]] <- ctlg[[tau]]
                cst[tt] <- tmp
            }
        }
    }
    
    return(opt[[tt]])
}
