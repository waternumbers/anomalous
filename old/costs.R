
#' Univariate cost functions
#'
#' @description Basic R implimentation of various univariate cost function
#' 
#' @param y univariate data series
#'
#' @name UnivariateCost
#' 
#' @details `fCmn` Gaussian mean, common known variance
#' @rdname UnivariateCost
#'@export
fCmn <- function(y,cl,beta){
    y <- split(y,cl)
    out <- sapply(y,function(x){ sum( (x-mean(x))^2 ) })
    out <- sum(out) + (length(out)-1)*beta
    return( out )
}

## fCmn <- function(y){ return(sum( (y-mean(y))^2 )) } # Gaussian mean, known variance

#' @details `fCv` Gaussian variance, common known mean
#' @rdname UnivariateCost
#'@export
fCv <- function(y,cl,beta){
    y <- split(y,cl)
    out <- sapply(y,function(x){ length(x)*(log(2*pi) + log( mean(x^2) ) +1 ) })
    out <- sum(out) + (length(out)-1)*beta
    return( out )
}

##fCv <- function(y){  return( length(y)*(log(2*pi) + log( mean(y^2) ) +1 ) ) } # Gaussian known mean, unknown variance

#' @details `fCmnv` Gaussian mean and variance
#' @rdname UnivariateCost
#'@export
fCmnv <- function(y,cl,beta){
    y <- split(y,cl)
    out <- sapply(y,function(x){ length(x)*(log(2*pi) + log( mean( (x-mean(x))^2 ) ) +1 ) })
    out <- sum(out) + (length(out)-1)*beta
    return( out )
}
##fCmnv <- function(y){ return( length(y)*(log(2*pi) + log( mean( (y-mean(y))^2 ) ) +1 ) ) }

#' @details `fCpois` Poisson
#' @rdname UnivariateCost
#' @export
fCpois <- function(y,cl,beta){
    y <- split(y,cl)
    out <- sapply(y,function(x){
        s <- sum(x)
        if(s==0){return(0)}
        else{ return( 2*s*(log(length(x))-log(s))) }
    })
    out <- sum(out) + (length(out)-1)*beta
    return( out )
}

## fCpois <- function(y){
##     x <- sum(y)
##     if(x==0){return(0)}
##     return( 2*x*(log(length(y))-log(x)) )
## }
