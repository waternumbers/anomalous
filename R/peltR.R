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
peltP <- function(y,mu,sigma,moveType,beta,min_length=2,max_length=.Machine$integer.max){
    
    n <- length(y)
    opt <- rep(list(NA),n) ## optimal partions
    cst <- rep(NA,n)
    ctlg <- list(NULL) ## catalog of partitions to try

    ## handle tt = 1 - can do any move except current
    ctlg[[1]] <- gaussPartition(beta,min_length,max_length)

    
    for(tt in 1:n){
        print(tt)
        ## update the ctlg
        octlg <- ctlg
        nc <- length(octlg) ## append to this so define length explicitly
        ctlg <- list()
        for(ii in 1:nc){
            for(mv in setdiff(moveType,"current")){
                ctlg[[length(ctlg)+1]] <- update(octlg[[ii]],y[tt],mu[tt],sigma[tt],mv)
            }
            if(is.finite(octlg[[ii]]@current_grp)){
                ctlg[[ii]] <- update(octlg[[ii]],y[tt],mu[tt],sigma[tt],"current")
            }
        }

        ## extract cst and validity
        ctlgCst <- isValid <- rep(NA,length(ctlg))
        for(ii in 1:length(ctlg)){
            ctlgCst[ii] <- ctlg[[ii]]@cost
            isValid[ii] <- ctlg[[ii]]@is_valid
        }

        if(tt == 10){ browser() }
        #browser()
##        print(range(ctlgCst))
        idx <- which.min(ctlgCst)
        opt[[tt]] <- ctlg[[idx]]
        cst[tt] <- ctlgCst[idx]
        idx <- (ctlgCst <= cst[tt]+max(beta)) | !isValid
        ctlg <- ctlg[idx]
        print(length(ctlg))
        
    }
    browser()
    return(opt[[tt]])
}
