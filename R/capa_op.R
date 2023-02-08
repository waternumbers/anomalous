#' An R implimentation of the segmented search algorithm
#' @param y univariate data series
#' @param segType function to generate an object of class Segment
#' @param beta Penalisation term for a new segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#'@export
capa_op <- function(y,mu,sigma,segType,beta,betaP, min_length=2,max_length=.Machine$integer.max){
    ##browser()
    n <- length(y)
    optRec <- rep(list(NULL),n) ## optimal partions
    optAnomRec <- rep(NA,n)

    
    ctlg <- rep(list(NULL),n) ## catalog of partitions to try
    cst <- rep(NA,n)
    last_n <- rep(NA,n)

    opt <- partition()
    optAnom <- TRUE
    for(tt in 1:n){
        
        ## update the catalog
        last_n <- NULL
        cst <- NULL
        if(tt>1){
            for(tau in 1:(tt-1)){
                ctlg[[tau]] <- update(ctlg[[tau]],y[tt],mu[tt],sigma[tt])
                last_n[tau] <- ctlg[[tau]]@last_n
                cst[tau] <-  ctlg[[tau]]@cost
            }
        }
        ## append to ctlg with new anomaly point
        ctlg[[tt]] <- update(opt,y[tt],mu[tt],sigma[tt],segType(beta,tt))
        last_n[tt] <- ctlg[[tt]]@last_n
        cst[tt] <-  ctlg[[tt]]@cost

        ## find minimal which is C1 in the paper
        idx <- min_length > last_n | max_length < last_n
        cst[idx] <- Inf
        ii <- which.min(cst)

        nopt <- ctlg[[ii]] ## initial new optimal value
        ncst <- cst[ii]
        nAnom <- TRUE
        
        ## update optimal wih base mode
        if( optAnom ){
            p <- update(opt, y[tt],mu[tt],sigma[tt], gaussFixed(0,tt))
        }else{
            p <- update(opt, y[tt],mu[tt],sigma[tt])
        }
        if(p@cost < ncst){
            nopt <- p
            nAnom <- FALSE
            ncst <- p@cost
        }

        p <- update(opt, y[tt],mu[tt],sigma[tt], gaussPoint(betaP,tt),isPoint=TRUE)
        if(p@cost < ncst){
            nopt <- p
            nAnom <- optAnom
            ncst <- p@cost
        }
        
        opt <- nopt
        optAnom <- nAnom
        optRec[[tt]] <- opt
        optAnomRec[tt] <- optAnom
    }

    return(opt)
}
