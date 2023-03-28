#' An R implimentation of the segmented search algorithm
#' @param y univariate data series
#' @param mu mean of univariate series y
#' @param sigma variance of univariate series y
#' @param segType function to generate an object of class Segment
#' @param beta Penalisation term for a new segment
#' @param betaP Cost of new point anomaly
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#'@export
capa <- function(y,mu,sigma,segType,beta,betaP, min_length=2,max_length=.Machine$integer.max){
    ##browser()
    n <- length(y)
    cnst <- max(beta,betaP)
    
    ## optimal partions
    optRec <- rep(list(NULL),n) ## optimal partions
    optAnomRec <- rep(NA,n)

    opt <- partition()
    optAnom <- TRUE

    ctlg <- list()
    ## loop time
    for(tt in 1:n){
        ##print(tt)
        ##if(tt==351){browser(); pjs <- 99}
        ## update the catalog
        last_n <- NULL
        cst <- NULL
        nc <- length(ctlg)
        if(nc>0){
            for(tau in 1:nc){
                if(is.null(ctlg[[tau]])){browser(); pjs <- 99}
                ctlg[[tau]] <- update(ctlg[[tau]],y[tt],mu[tt],sigma[tt])
                last_n[tau] <- ctlg[[tau]]@last_n
                cst[tau] <-  ctlg[[tau]]@cost
            }
        }
        ## append to ctlg with new anomaly point
        if(is.null(opt)){browser(); pjs <- 99}
        ctlg[[nc+1]] <- update(opt,y[tt],mu[tt],sigma[tt],segType(beta,tt))
        last_n[nc+1] <- ctlg[[nc+1]]@last_n
        cst[nc+1] <-  ctlg[[nc+1]]@cost
        
        ## find minimal which is C1 in the paper
        idx <- min_length > last_n | max_length < last_n
        cst[idx] <- Inf
        ii <- which.min(cst)

        nopt <- ctlg[[ii]] ## initial new optimal value
        ncst <- cst[ii]
        nAnom <- TRUE
        
        ## update optimal with base mode
        if(is.null(opt)){browser(); pjs <- 99}
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
        if(is.null(opt)){browser(); pjs <- 99}
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
        
        ## trim catalog
        ##browser()
        idx <- cst <= (opt@cost + cnst) | last_n < min_length
        if( !all(is.finite(idx)) ){ browser() }
        ##if(length(idx)>length(ctlg)){ browser();pjs <- 98}
        ctlg <- ctlg[idx]
        ##print(paste(tt,length(ctlg)))
        
    }
    
    return(opt)
}
