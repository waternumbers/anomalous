#' @export
segment <-
    R6::R6Class("segment",
                public = list(
                    initialize = function(x,beta){
                        self$beta <- beta
                        self$update(x)
                    },
                    n = 0,
                    ss = NA,
                    beta = NA,
                    update = function(x){
                        self$n <- self$n + 1
                        self$updateSS(x)
                        invisible(self)
                    },
                    cost = function(min_length=2, max_length=Inf){
                        if( (self$n < min_length) | (self$n > max_length) ){ return(Inf) }
                        out <- self$computeCost() + self$beta
                        return(out)
                    },
                    computeCost = function(){Inf},
                    updateSS = function(){}
                 )
                )

#' @export
seg_mean <-
    R6::R6Class("seg_mean",
                inherit = segment,
                public = list(
                    ss = rep(0,2),
                    updateSS = function(x){
                        self$ss <- self$ss + c(x, x*x)
                    },
                    computeCost = function(){
                        return( self$ss[2] - ( (self$ss[1]^2)/self$n ) )
                    }
                )
                )

#' @export
seg_var <-
    R6::R6Class("seg_var",
                inherit = segment,
                public = list(
                    ss = 0,
                    updateSS = function(x){
                        self$ss <- self$ss + x*x
                    },
                    computeCost = function(){
                        return( self$n*(log(2*pi) + log( self$ss/self$n ) +1 ) )
                    }
                )
                )
#' @export
seg_mean_var <- 
    R6::R6Class("seg_mean_var",
                inherit = segment,
                public = list(
                    ss = rep(0,2),
                    updateSS = function(x){
                        self$ss <- self$ss + c(x,x*x)
                    },
                    computeCost = function(){
                        out <- self$ss[2] - ( (self$ss[1]^2)/self$n )
                        out <- self$n*(log(2*pi) + log( out / self$n ) +1 )
                        return(out)
                    }
                )
                )
