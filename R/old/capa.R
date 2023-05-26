#' An R implimentation of the segmented search algorithm
#' @param y univariate data series
#' @param mu mean of univariate series y
#' @param sigma variance of univariate series y
#' @param baseSeg function to generate an instance of the base segment
#' @param collectiveSeg function to generate an instance of the collectove anomlay segment
#' @param pointSeg function to generate an instance of the point anomlay segment
#' @param beta Penalisation term for a new segment
#' @param betaP Cost of new point anomaly
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#'@export
capa <- function(y,mu,sigma,baseSeg, collectiveSeg, pointSeg,beta,betaP, min_length=2,max_length=.Machine$integer.max){
    ##browser()
    n <- length(y)
    cnst <- max(beta,betaP)


    opt <- partition()

    ctlg <- list()
    ##browser()
    ## loop time
    for(tt in 1:n){
        ## update ctlg
        ctlg[[ length(ctlg)+1 ]] <- addCollective(opt,collectiveSeg,beta,tt)
        
        nc <- length(ctlg)
        last_n <- rep(NA,nc)
        cst <- rep(NA,nc)

        if(nc>0){
            for(tau in 1:nc){
                ctlg[[tau]] <- update(ctlg[[tau]],y[tt],mu[tt],sigma[tt])
                last_n[tau] <- ctlg[[tau]]@last_n
                cst[tau] <-  ctlg[[tau]]@cost
            }
        }

        ## find minimal which is C1 in the paper
        idx <- (min_length > last_n) | (max_length < last_n)
        cst[idx] <- Inf
        ii <- which.min(cst)
        nopt <- ctlg[[ii]]

        ## update optimal with base model C2 in paper
        p <- opt
        if(opt@isAnom){ p <- addBase(baseSeg,0,tt) }
        p <- update(p,x[tt],mu[tt],sigma[tt])

        if( p@cost < nopt@cost ){ nopt <- p }

        ## update opt by adding a point anomaly
        p <- addPoint(opt,y[tt],mu[tt],sigma[tt],pointSeg,betaP,tt)
        if( p@cost < nopt@cost ){ nopt <- p }

        opt <- nopt

        ## trim catalog
        idx <- (cst <= (opt@cost + cnst)) | (last_n < min_length)
        ctlg <- ctlg[idx]
        ##print(paste(tt,length(ctlg)))
    }
    
    return(opt)
}
