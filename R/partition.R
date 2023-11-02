## set up generic methods
#' @export
addCollective <- function(p,...){ UseMethod("addCollective",p) }
#' @export
addBase <- function(p,...){ UseMethod("addBase",p) }
#' @export
addPoint <- function(p,...){ UseMethod("addPoint",p) }
#' @export
collective_anomalies <- function(p){ UseMethod("collective_anomalies",p) }
#' @export
point_anomalies <- function(p){ UseMethod("point_anomalies",p) }

#' @export
partition <- function(beta,betaP,min_length){
    structure(list(
        ca=list(),
        pa=list(),
        beta=beta,
        betaP=betaP,
        min_length=min_length,
        baseCost=0,
        cost=0,last_time=0),
        class="partition")
}

#' @export
addCollective.partition <- function(p,x,s,e,...){
    ##cst <- collectiveCost(x,s,e,p$beta)
    cst <- x$collectiveCost(s,e,p$beta)
    p$ca[[length(p$ca)+1]] <- c(start=s,end=e,cost=cst)
    p$cost <- p$cost + cst
    p$last_time <- e
    return(p)
}

#' @export
addBase.partition <- function(p,x,s,e,...){
    cst <- x$baseCost(s,e,0)
    p$baseCost <- p$baseCost + cst
    p$cost <- p$cost + cst
    p$last_time <- e
    return(p)
}

#' @export
addPoint.partition <- function(p,x,s,...){
    cst <- x$pointCost(s,p$betaP)
    ## cst <- pointCost(x,s,p$betaP)
    p$pa[[length(p$pa)+1]] <- c(location=s,cost=cst) 
    p$cost <- p$cost + cst
    p$last_time <- s
    return(p)
}

#' @export
collective_anomalies.partition <- function(p){ as.data.frame( do.call(rbind,p$ca) ) }

#' @export
point_anomalies.partition <- function(p){ as.data.frame( do.call(rbind,p$pa) ) }


#' @export
plot.partition <- function(p,xx,yy,...){
    
    ca <- collective_anomalies(p)
    pa <- point_anomalies(p)
    showRegions <- TRUE
    lenx <- max(c(ca$end,pa$location))

    if( missing(xx) ){ xx <- 1:lenx }
    if( length(xx) < lenx ){ stop("xx is to short") }
    
    if( missing(yy) ){
        ## set y to be the score
        yy <- rep(NA,lenx)
        for(ii in 1:nrow(ca)){
            yy[ ca$start[ii]:ca$end[ii] ] <- ca$cost[ii]
        }
        yy[ pa$location ] <- pa$cost
        showRegions <- FALSE
    }

    plot(xx,yy,...)
    if(showRegions){
        for(ii in 1:nrow(ca)){
            rect(xleft = xx[ ca$start[ii] ], xright = xx[ ca$end[ii] ],
                 ybottom = par("usr")[3], ytop = par("usr")[4], 
                 border = NA, col = adjustcolor("blue", alpha = 0.3))
        }
        points( xx[pa$location], yy[pa$location], pch=23, col = "blue" )
    }
}
        
        ## then we will plot the
    
