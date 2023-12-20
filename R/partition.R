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
plot.partition <- function(x,...){
    
    ca <- collective_anomalies(x)
    pa <- point_anomalies(x)
    showRegions <- TRUE
    lenx <- max(c(ca$end,pa$location))

    z <- list(...)

    ## if( "xx" %in% names(z) ){ xx <- z$xx; z$xx <- NULL}else{ xx <- 1:lenx }
    ## if( length(xx) < lenx ){ stop("xx is to short") }


    
    ## if( "yy" %in% names(z) ){
    ##     yy <- z$yy
    ##     z$yy <- NULL
    ## }else{
    ##     ## set y to be the score
    ##     yy <- rep(NA,lenx)
    ##     for(ii in 1:nrow(ca)){
    ##         yy[ ca$start[ii]:ca$end[ii] ] <- ca$cost[ii]
    ##     }
    ##     yy[ pa$location ] <- pa$cost
    ##     showRegions <- FALSE
    ## }

    if( "xx" %in% names(z) ){ z$x <- z$xx; z$xx <- NULL}else{ z$x <- 1:lenx }
    if( "yy" %in% names(z) ){
        z$y <- z$yy
        z$yy <- NULL
    }else{
        ## set y to be the score
        yy <- rep(NA,lenx)
        for(ii in 1:nrow(ca)){
            yy[ ca$start[ii]:ca$end[ii] ] <- ca$cost[ii]
        }
        yy[ pa$location ] <- pa$cost
        showRegions <- FALSE
        z$y <- yy
    }
    if(!("xlab" %in% names(z))){z$xlab=""}
    if(!("ylab" %in% names(z))){z$ylab=""}
    
    do.call(plot,z)
    ##plot(x=xx,y=yy,z) ##...)
    if(showRegions){
        if(nrow(ca)>0){
            for(ii in 1:nrow(ca)){
                graphics::rect(xleft = z$x[ ca$start[ii] ], xright = z$x[ ca$end[ii] ],
                               ybottom = graphics::par("usr")[3],
                               ytop = graphics::par("usr")[4], 
                               border = NA, col = grDevices::adjustcolor("blue", alpha = 0.3))
            }
        }
        graphics::points( z$x[pa$location], z$y[pa$location], pch=23, col = "red" )
    }
}
    
