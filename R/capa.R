#' An R implimentation of the segmented search algorithm
#' @param y univariate data series
#' @param segType function to generate an object of class Segment
#' @param beta Penalisation term for a new segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#'@export
capa <- function(y,mu,sigma,segType,beta,betaP, min_length=2,max_length=.Machine$integer.max){
    ##browser()
    n <- length(y)
    cnst <- max(beta,betaP)

    opt <- rep(list(NULL),n) ## optimal partions
    isAnom <- rep(NA,n)
    
    ctlg <- rep(list(NULL),n) ## catalog of partitions to try

    ## initialise the catalog
    ctlg <- list(partition(),partition())
    for(tt in 1:min_length){
        if(tt==1){
            seg1 <- gaussFixed(0,1)
            seg2 <- segType(beta,1)
        }else{
            seg1 <- seg2 <- NULL
        }
        ctlg[[1]] <- update(ctlg[[1]],y[tt],mu[tt],sigma[tt],seg1)
        ctlg[[2]] <- update(ctlg[[2]],y[tt],mu[tt],sigma[tt],seg2)
    }

    cAnom <- c(FALSE,TRUE)
    last_n <- c(min_length,min_length)
    
    cst <- sapply(ctlg,function(x){x@cost})
    ii <- which.min(cst)
    
    opt[[min_length]] <- ctlg[[ ii ]]
    isAnom <- cAnom[ii]

    idx <- cst <= (cst[ii] + cnst)
    ctlg <- ctlg[idx]
    cAnom <- cAnom[idx]

    browser()
    ## loop time
    for(tt in (min_length+1):n){

        last_n <- NULL
        cst <- NULL
        nc <- length(ctlg)

        
        ## update catalog
        for(tau in 1:nc){
            ## ## try a point
            
            ## nadd <- length(ctlg)+1
            ## ctlg[[ nadd ]] <- update( ctlg[[tau]], y[tt], mu[tt], sigma[tt], gaussFixed(0,tt) )
            ##         cAnom[nadd] <- FALSE
            ##         last_n[nadd] <- ctlg[[ nadd ]]@last_n
            ##         cst[nadd] <- ctlg[[ nadd ]]@cost
            ## ## we can update of we can try a point or we can updata
            
            ## if(ctlg[[tau]]@last_n >= min_length){
            ##     ## try move to new collective anomaly
            ##     nadd <- length(ctlg)+1
            ##     ctlg[[ nadd ]] <- update( ctlg[[tau]], y[tt], mu[tt], sigma[tt], segType(beta,tt) )
            ##     cAnom[nadd] <- TRUE
            ##     last_n[nadd] <- ctlg[[ nadd ]]@last_n
            ##     cst[nadd] <- ctlg[[ nadd ]]@cost
                
            ##     if(cAnom[tau]){ ## if anomlay alrady try move to background type
            ##         nadd <- length(ctlg)+1
            ##         ctlg[[ nadd ]] <- update( ctlg[[tau]], y[tt], mu[tt], sigma[tt], gaussFixed(0,tt) )
            ##         cAnom[nadd] <- FALSE
            ##         last_n[nadd] <- ctlg[[ nadd ]]@last_n
            ##         cst[nadd] <- ctlg[[ nadd ]]@cost
            ##     }
            ## }
            ## update current structure
            ctlg[[ tau ]] <- update( ctlg[[tau]], y[tt], mu[tt], sigma[tt])
            last_n[tau] <- ctlg[[ tau ]]@last_n
            cst[tau] <- ctlg[[ tau ]]@cost
            
            ## try a point anomlay
            nadd <- length(ctlg)+1
            ctlg[[ nadd ]] <- update( ctlg[[tau]], y[tt], mu[tt], sigma[tt], gaussVar(betaP,tt),isPoint=TRUE )
            cAnom[nadd] <- cAnom[tau]
            last_n[nadd] <- ctlg[[ nadd ]]@last_n
            cst[nadd] <- ctlg[[ nadd ]]@cost
            
        }

        ## try a change for the optimum
        if( isAnom ){
            seg <-  gaussFixed(0,tt)
            segA <- FALSE
        }else{
            seg <-  segType(beta,tt)
            segA <- TRUE
        }
        
        nadd <- length(ctlg)+1
        ctlg[[ nadd ]] <- update( ctlg[[tau]], y[tt], mu[tt], sigma[tt], seg)
        cAnom[nadd] <- segA
        last_n[nadd] <- ctlg[[ nadd ]]@last_n
        cst[nadd] <- ctlg[[ nadd ]]@cost

        ## find minimum
        if(tt==100){browser()}
        not_yet_valid <- min_length > last_n
        cst[ not_yet_valid | max_length < last_n ] <- Inf
               
        ii <- which.min(cst) ## the optimal choice

        ## copy min value over
        opt[[tt]] <- ctlg[[ii]]
        isAnom <- cAnom[ii]

        
        ## trim catalog
        idx <- cst <= (cst[ii] + cnst) | not_yet_valid
        ctlg <- ctlg[ idx ]
        cAnom <- cAnom[idx]

        print(paste(tt,length(ctlg)))
        
    }
    
    return(opt[[tt]])
}
