## set up generic methods

addCollective <- function(p,...){ UseMethod("addCollective",obj) }
addBase <- function(p,...){ UseMethod("addBase",obj) }
addPoint <- function(p,...){ UseMethod("addPoint",obj) }
collective_anomalies <- function(p,...){ UseMethod("collective_anomalies",obj) }
point_anomalies <- function(p,...){ UseMethod("point_anomalies",obj) }

#' @export
partition <- function(beta,betaP,min_length){
    structure(list(
        ca=list(),
        pa=list(),
        beta=beta,
        betaP=betaP,
        min_length=min_length,
        cost=0,last_time=0),
        class="partition")
}

#' @export
addCollective.partition <- function(p,s,e,x){
    ##cst <- collectiveCost(x,s,e,p$beta)
    cst <- x$collectiveCost(s,e,p$beta)
    p$ca[[length(p$ca)+1]] <- c(start=s,end=e,cost=cst)
    p$cost <- p$cost + cst
    p$last_time <- e
    return(p)
}

#' @export
addBase.partition <- function(p,s,e,x){
    ##p$cost <- p$cost + baseCost(x,s,e,0)
    p$cost <- p$cost + x$baseCost(s,e,0)
    p$last_time <- e
    return(p)
}

#' @export
addPoint.partition <- function(p,s,x){
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
