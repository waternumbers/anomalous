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
capaP <- function(y,mu,sigma,anomType,beta,min_length=2,max_length=.Machine$integer.max){
    
    n <- length(y)
    opt <- rep(list(NA),n) ## optimal partions
    cst <- rep(NA,n)
    ctlg <- list() ## catalog of partitions to try
    state <- NULL

    ctlg[[1]] <- gaussPartition(beta,min_length,max_length)
    
    ## handle tt = 1
    ctlg <- list() ## catalog of partitions to try
    state <- NULL
    tmp <- gaussPartition(beta,min_length,max_length)
    for(ii in  c("background",anomType)){
        ctlg[[length(ctlg)+1]] <- update(tmp,y[1],mu[1],sigma[1],1,ii)
        state[length(state)+1] <- ii
    }
    
    for(tt in 2:n){
        ##print(tt)
        ##if(tt==100){ browser() }
        ## update the ctlg
        octlg <- ctlg
        ostate <- state
        nc <- length(octlg) ## append to this so define length explicitly
        ctlg <- list()
        state <- NULL
        for(ii in 1:nc){
            moveTo <- switch(ostate[ii],
                             "background" = c("current",anomType),
                             c("current","background"))
            if( ostate[ii] %in% moveTo){
                print(c(ostate[ii],moveTo))
            }
            
            for(mv in moveTo){
                ctlg[[length(ctlg)+1]] <- update(octlg[[ii]],y[tt],mu[tt],sigma[tt],tt,mv)
                state[length(state)+1] <- ifelse(mv=="current",ostate[ii],mv)
            }
        }
        #browser()
        ## extract cst and validity
        ctlgCst <- isValid <- rep(NA,length(ctlg))
        for(ii in 1:length(ctlg)){
            ctlgCst[ii] <- ctlg[[ii]]@cost
            isValid[ii] <- ctlg[[ii]]@is_valid
        }

        if(tt == 30){ browser() }
        #browser()
##        print(range(ctlgCst))
        idx <- which.min(ctlgCst)
        opt[[tt]] <- ctlg[[idx]]
        cst[tt] <- ctlgCst[idx]
        idx <- (ctlgCst <= cst[tt]+max(beta)) | !isValid
        ctlg <- ctlg[idx]
        state <- state[idx]
        print(paste(tt, length(ctlg)))
        
    }

    return(opt[[tt]])
}
