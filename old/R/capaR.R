#' An R implimentation of the segmented search algorithm
#' @param y univariate data series
#' @param mu mean of univariate series y
#' @param sigma variance of univariate series y
#' @param segType the names of the segment class to use
#' @param beta Penalisation term for a new segment
#' @param betaP Cost of new point anomaly
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#'@export
capaR <- function(y,mu,sigma,segType,beta,betaP, min_length=2,max_length=.Machine$integer.max){
    ##browser()
    n <- length(y)
    cnst <- max(beta,betaP)

    ## generator functions
    fb <- function(beta,start){ new(gaussKnown,beta,start) }
    fc <- function(beta,start){ new(gaussMean,beta,start) }
    fp <- function(beta,start){ new(gaussKnown,beta,start) } ##TODO change to gaussPoint
   
    ## optimal partions
    optRec <- rep(list(NULL),n) ## optimal partions

    opt <- addBase(partition(),fb,0,1)
    

    ctlg <- list()
    ctlg[[1]] <- opt
    ## loop time
    for(tt in 1:n){
        print(tt)
        ##if(tt==351){browser(); pjs <- 99}
        ## update the catalog
        nc <- length(ctlg)
        last_n <- rep(NA,nc+2)
        cst <- rep(NA,nc+2)

        if(nc>0){
            for(tau in 1:nc){
                if(is.null(ctlg[[tau]])){browser(); pjs <- 99}
                ctlg[[tau]] <- update(ctlg[[tau]],y[tt],mu[tt],sigma[tt])
                last_n[tau] <- ctlg[[tau]]@last_n
                cst[tau] <-  ctlg[[tau]]@cost
            }
        }
        
        ## update opt by changing collective to base or visa versa
        if(opt@isAnom){
            ## in a collective anomaly period so try changing to baseline
            ctlg[[nc+1]] <- addBase(opt,fb,0,tt)
        }else{
            ## in base so start a collective anomaly
            ctlg[[nc+1]] <- addCollective(opt,fc,0,tt)
        }
        ctlg[[nc+1]] <- update(ctlg[[nc+1]],y[tt],mu[tt],sigma[tt])
        last_n[nc+1] <- ctlg[[nc+1]]@last_n
        cst[nc+1] <-  ctlg[[nc+1]]@cost

        ## update opt by adding a point anomaly
        ctlg[[nc+2]] <- addPoint(opt,y[tt],mu[tt],sigma[tt],fc,betaP,tt)
        last_n[nc+2] <- ctlg[[nc+2]]@last_n
        cst[nc+2] <-  ctlg[[nc+2]]@cost

        ## work out new optimum
        browser()
        tmp <- cst
        tmp[(last_n < min_length) | (last_n > max_length)] <- Inf
        odx <- which.min(tmp)

        opt <- ctlg[[odx]]
        ctlg <- ctlg[ (cst<=(tmp[odx]+cnst)) | (last_n < min_length) ]
        
    }
    
    return(opt)
}
