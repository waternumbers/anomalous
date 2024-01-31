## set up generic methods
#' @export
collective_anomalies <- function(p){ UseMethod("collective_anomalies",p) }
#' @export
point_anomalies <- function(p){ UseMethod("point_anomalies",p) }

#' @export
partition <- function(beta,betaP,min_length){
    structure(list(
        endPoint=NULL,
        cost=NULL,
        type=NULL,
        beta=beta,
        betaP=betaP,
        min_length=min_length),
        class="partition")
}

#' @export
collective_anomalies.partition <- function(p,t){
    tmp <- summary(p,t)
    tmp[tmp$type=="collective",]
}

#' @export
point_anomalies.partition <- function(p,t){
    tmp <- summary(p,t)
    tmp$end <- NULL
    names(tmp) <- gsub("start","location",names(tmp))
    tmp[tmp$type=="point",]
}


#' @export
summary.partition <- function(object,t){
    if(missing(t)){ t <- length(object$endPoint) }
    
    tmp <- t
    while(tmp[1] > 0){ tmp <- c( object$endPoint[tmp[1]], tmp) }
    
    data.frame(start = head(tmp+1, -1),
               end = tmp[-1],
               type = object$type[ tmp[-1] ],
               cost = diff( c(0,object$cost[ tmp[-1] ]) )
               )
}

#' @export
plot.partition <- function(x,...){
    
    ca <- collective_anomalies(x)
    pa <- point_anomalies(x)
    showRegions <- TRUE
    lenx <- max(c(ca$end,pa$location))

    z <- list(...)

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
    
