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
fCmn <- function(y){ return(sum( (y-mean(y))^2 )) } # Gaussian mean, known variance

#' @details `fCv` Gaussian variance, common known mean
#' @rdname UnivariateCost
#'@export
fCv <- function(y){  return( length(y)*(log(2*pi) + log( mean(y^2) ) +1 ) ) } # Gaussian known mean, unknown variance

#' @details `fCmnv` Gaussian mean and variance
#' @rdname UnivariateCost
#'@export
fCmnv <- function(y){ return( length(y)*(log(2*pi) + log( mean( (y-mean(y))^2 ) ) +1 ) ) }

#' @details `fCpois` Poisson
#' @rdname UnivariateCost
#' @export
fCpois <- function(y){
    x <- sum(y)
    if(x==0){return(0)}
    return( 2*x*(log(length(y))-log(x)) )
}
