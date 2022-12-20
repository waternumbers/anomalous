partition <-
    R6::R6Class("Partition",
                public = list(
                    initialize = function(x,n){
                        self$cl <- rep(NA,n)
                        self$update(x,TRUE)
                    },
                    update = function(x,new=FALSE){
                        browser()
                        ## TODO add error checking
                        n <- which.max(is.na(self$cl))
                        if(new){
                            self$nseg <- self$nseg + 1
                            self$cl[n] <- self$nseg
                        }else{
                            self$cl[n] <- self$cl[n-1]
                        }
                        self$updateSummary(self$cl[n],x)
                        invisible(self)
                    },
                    cost = function(beta, min_seg_length=2){
                        ## TODO write checks
                        return( self$computeCost(beta, min_seg_length) )
                    },
                    cl = NULL,
                    ss = NULL,
                    nseg = 0,
                    updateSummary = function(seg,x){ invisible(self) },
                    computeCost = function(beta){ NA }
                )
                )

meanCP <-
    R6::R6Class("meanCP",
                inherit = partition,
                public = list(
                    updateSummary = function(seg,x){
                        if( length(self$ss$n) < seg ){
                            m <- seg - length(self$ss$n)
                            self$ss <- rbind(self$ss,
                                        data.frame(n=rep(0,m),
                                                   sx=rep(0,m),
                                                   sxx=rep(0,m)))
                        }                        
                        self$ss$n[seg] <- self$ss$n[seg] + 1
                        self$ss$sx[seg] <- self$ss$sx[seg] + x
                        self$ss$sxx[seg] <- self$ss$sxx[seg] + x*x
                        invisible(self)
                    },
                    computeCost = function(beta,min_seg_length){
                        out <- self$ss$sxx - ( (self$ss$sx^2)/self$ss$n )
                        jj <- length(out)
                        out <- sum(out) + (jj-1)*beta
                        ## apply the segment length rules
                        if( self$ss$n[jj] < min_seg_length ){ out <- NA }
                        if( any( self$ss$n[-jj] < min_seg_length) ){ out <- Inf }
                        return(out)
                    }
                )
                )
                        
varCP <-
    R6::R6Class("varCP",
                inherit = partition,
                public = list(
                    updateSummary = function(seg,x){
                        if( length(self$ss$n) < seg ){
                            m <- seg - length(self$ss$n)
                            self$ss <- rbind(self$ss,
                                        data.frame(n=rep(0,m),
                                                   sxx=rep(0,m)))
                        }                        
                        self$ss$n[seg] <- self$ss$n[seg] + 1
                        self$ss$sxx[seg] <- self$ss$sxx[seg] + x*x
                        invisible(self)
                    },
                    computeCost = function(beta,min_seg_length){
                        out <- self$ss$n*(log(2*pi) + log( self$ss$sxx/self$ss$n ) +1 )
                        jj <- length(out)
                        out <- sum(out) + (jj-1)*beta
                        ## apply the segment length rules
                        if( self$ss$n[jj] < min_seg_length ){ out <- NA }
                        if( any( self$ss$n[-jj] < min_seg_length) ){ out <- Inf }
                        return(out)
                    }
                )
                )


meanvarCP <-
    R6::R6Class("meanvarCP",
                inherit = partition,
                public = list(
                    updateSummary = function(seg,x){
                        if( length(self$ss$n) < seg ){
                            m <- seg - length(self$ss$n)
                            self$ss <- rbind(self$ss,
                                        data.frame(n=rep(0,m),
                                                   sxx=rep(0,m)))
                        }                        
                        self$ss$n[seg] <- self$ss$n[seg] + 1
                        self$ss$sxx[seg] <- self$ss$sxx[seg] + x*x
                        invisible(self)
                    },
                    computeCost = function(beta,min_seg_length){
                        out <- self$ss$sxx - ( (self$ss$sx^2)/self$ss$n )
                        out <- self$ss$n*(log(2*pi) + log( out / self$ss$n ) +1 )
                        jj <- length(out)
                        out <- sum(out) + (jj-1)*beta
                        ## apply the segment length rules
                        if( self$ss$n[jj] < min_seg_length ){ out <- NA }
                        if( any( self$ss$n[-jj] < min_seg_length) ){ out <- Inf }
                        return(out)
                    }
                )
                )




prt_mn <- function(y){
    cl <- rep(NA,length(y))
    cl[1] <- 1
    ss <- data.frame(cl = 1,
                     n = 1,
                     sx = y[1],
                     sxx = y[1]^2)
    return(list(cl=cl,ss=ss))
}

update <- function(prt,y,new=FALSE){
    n <- which.max(is.na(prt$cl))
    prt$cl[n] <- prt$cl[n-1]
    if(new){
        prt$cl[n] <- prt$cl[n]+1
        
        prt$ss <- rbind(prt$ss,
                        data.frame(cl = prt$cl[n],
                                   n = 1,
                                   sx = y[n],
                                   sxx = y[n]^2))
    }else{
        ii <- prt$cl[n]
        prt$ss$n[ii] <- prt$ss$n[ii] + 1
        prt$ss$sx[ii] <- prt$ss$sx[ii] + y[n]
        prt$ss$sxx[ii] <- prt$ss$sxx[ii] + y[n]^2
    }
    return(prt)
}

cost <- function(prt,beta,min_seg_length=10){
    out <- prt$ss$sxx - ( (prt$ss$sx^2)/prt$ss$n )
    jj <- length(out)
    out <- sum(out) + (jj-1)*beta
    ## apply the segment length rules
    if( prt$ss$n[jj] < min_seg_length ){ out <- NA }
    if( any( prt$ss$n[-jj] < min_seg_length) ){ out <- Inf }
    return(out)
}


    
    






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
