## set up generic methods
#' @export
collective_anomalies <- function(p,t){ UseMethod("collective_anomalies",p) }
#' @export
point_anomalies <- function(p,t){ UseMethod("point_anomalies",p) }

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
collective_anomalies.partition <- function(p,t=NULL){
    tmp <- summary(p,t)
    tmp[tmp$type=="collective",]
}

#' @export
point_anomalies.partition <- function(p,t=NULL){
    tmp <- summary(p,t)
    tmp$end <- NULL
    names(tmp) <- gsub("start","location",names(tmp))
    tmp[tmp$type=="point",]
}


#' @export
summary.partition <- function(object,...){

##    if(missing(t)){ t <- length(object$endPoint) }
    t <- list(...)$t
    if( is.null(t) ){ t <- tail(which(!is.na(object$cost)), 1) } ## length(object$endPoint) }
    if( length(t) == 0 ){
        return( data.frame(start = integer(0),
                           end = integer(0),
                           type = character(0),
                           cost = numeric(0) ) )
    }
    
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

    eli <- list(...)
    t <- eli$t
    showRegions <- ifelse(is.null(eli$showRegions[1]),TRUE,as.logical(eli$showRegions[1]))
    
    sm <- summary(x,t)
    if(length(sm)==0){ return() }
    lenx <- max(sm$end)

    z <- list(...)

    if( "xx" %in% names(z) ){ z$x <- z$xx; z$xx <- NULL }else{ z$x <- 1:lenx }
    if( "yy" %in% names(z) ){ z$y <- z$yy; z$yy <- NULL }
    else{        
        ## set y to be the score
        z$y <- rep(NA,lenx)
        for(ii in 1:nrow(sm)){
            z$y[ sm$start[ii]:sm$end[ii] ] <- sm$cost[ii] / (sm$end[ii]-sm$start[ii]+1)
        }
        z$ylab = "Cost per observation"
    }
    if(!("xlab" %in% names(z))){z$xlab=""}
    if(!("ylab" %in% names(z))){z$ylab=""}
    
    do.call(plot,z)
    ##plot(x=xx,y=yy,z) ##...)
    if(showRegions){
        for(ii in which(sm$type=="collective")){
            graphics::rect(xleft = z$x[ sm$start[ii] ], xright = z$x[ sm$end[ii] ],
                           ybottom = graphics::par("usr")[3],
                           ytop = graphics::par("usr")[4], 
                           border = NA, col = grDevices::adjustcolor("blue", alpha = 0.3))
        }
        idx <- sm$start[ which(sm$type=="point") ]
        graphics::points( z$x[idx], z$y[idx], pch=23, col = "red" )
    }
}

#' @export
coef.partition <- function(object,...){

    eli <- list(...)
    t <- eli$t
    fCost <- eli$cost
    if(is.null(fCost)){ stop("A cost function is required") }
    if( !("param" %in%  names(fCost)) ){ stop("Parameter method not available for cost function") }

    sm <- summary(object,t)

    tmp <- fCost$param(sm$start[1],sm$end[1],sm$type[1])
    out <- matrix(NA,nrow(sm),length(tmp),dimnames=list(NULL,names(tmp)))
    
    for(ii in 1:nrow(sm)){
        out[ii,] <- fCost$param(sm$start[ii],sm$end[ii],sm$type[ii])
    }
    return(out)
}
