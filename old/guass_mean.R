#' Univariate Gaussian Mean estimation functions
#'
#' @description Basic R implimentation of functions for estimating the guassian mean
#' 
#' @param y new data point
#' @param s data frame of summary statistics
#' @param pen penality value
#'
#' @name GuassianMean
#' 
#' @details `guass_mean_update` Add new observation to summary statistics
#' @rdname GuassianMean
#'@export
guass_mean <- function(obj,y,grp,pen){
    s <- obj$summaryStat
    if(grp<1){
        s <- rbind(s,data.frame(n=1,sum_x=y,sum_x2=y^2))
    }else{
        s$n[grp] <- s$n[grp]+1
        s$sum_x[grp] <- s$sum_x[grp] + y
        s$sum_x2[grp] <- s$sum_x2[grp] + y^2
    }
    obj$summaryStat <- s
    obj$cost <- sum( s$sum_x2 - (s$sum_x^2)/s$n ) + nrow(s)*pen
    
    return(s)
}


