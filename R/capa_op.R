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
    opt <- rep(list(NULL),n) ## optimal partions
    isAnom <- rep(NA,n)
    
    ctlg <- rep(list(NULL),n) ## catalog of partitions to try
    cst <- rep(NA,n)
    last_n <- rep(NA,n)
    cAnom <- rep(NA,n)

    ## initialise
    p <- list(partition(),partition())
    for(tt in 1:min_length){
        if(tt==1){
            seg1 <- gaussFixed(0,1)
            seg2 <- segType(beta,1)
        }else{
            seg1 <- seg2 <- NULL
        }
        p[[1]] <- update(p[[1]],y[tt],mu[tt],sigma[tt],seg1)
        p[[2]] <- update(p[[2]],y[tt],mu[tt],sigma[tt],seg2)
    }

    if( p[[1]]@cost < p[[2]]@cost ){
        ctlg[[min_length]] <- p[[1]]
        cAnom[min_length] <- FALSE
    }else{
        ctlg[[min_length]] <- p[[2]]
        cAnom[min_length] <- TRUE
    }
    
    opt[[min_length]] <- ctlg[[min_length]]
    isAnom[min_length] <- cAnom[min_length]
    
    ## loop time
    for(tt in (min_length+1):n){

        ## get the optimal choce from the catalog
        for(tau in min_length:(tt-1)){
            ctlg[[tau]] <- update( ctlg[[tau]], y[tt], mu[tt],sigma[tt] )
            last_n[tau] <- ctlg[[tau]]@last_n
            cst[tau] <- ctlg[[tau]]@cost
        }
        if(isAnom[tt-1]){
            ctlg[[tt]] <- update(opt[[tt-1]], y[tt],mu[tt],sigma[tt], gaussFixed(0,tt))
        }else{
            ctlg[[tt]] <- update(opt[[tt-1]], y[tt],mu[tt],sigma[tt], segType(beta,tt))
        }
        last_n[tt] <- ctlg[[tt]]@last_n
        cst[tt] <- ctlg[[tt]]@cost
        cAnom[tt] <- !isAnom[tt-1]

        idx <- min_length > last_n | max_length < last_n
        cst[idx] <- Inf
        ii <- which.min(cst)

        #browser()
        opt[[tt]] <- ctlg[[ii]]
        isAnom[[tt]] <- cAnom[ii]

        ## try an extension to the last optimal
        p <- update(opt[[tt-1]], y[tt],mu[tt],sigma[tt])
        if( p@cost < opt[[tt]]@cost ){
            opt[[tt]] <- p
            isAnom[[tt]] <- isAnom[tt-1]
        }
        
        ## try adding a point change to the last optimal
        p <- update(opt[[tt-1]], y[tt],mu[tt],sigma[tt],gaussVar(betaP,tt),isPoint=TRUE)
        if( p@cost < opt[[tt]]@cost ){
            opt[[tt]] <- p
            isAnom[[tt]] <- isAnom[tt-1]
        }
        
    }
    
    return(opt[[tt]])
}
