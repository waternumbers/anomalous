#' An R implimentation of the segmented search algorithmpelt algorithm
#' 
#' @param part the starting partition
#' @param fCost the cost function
#' @param prune logical, should pruning be used
#' @param verbose logical, print out progress
#'
#' @return the optimal partition
#' 
#' @details Basic R implimentation of pelt - not efficent
#' @export
capa <- function(part,fCost,prune = TRUE,verbose=FALSE,...){
    
    cnst <- max(part$beta, part$betaP)
    
    endPoints <- c(0) ## end points to search over
    
    for(tt in 1:fCost$maxT){ ##fCost$validTimes){ ##maxT){
        if(verbose && (tt %% 100==0)) {
            ## Print on the screen some message
            cat(paste0("time step: ", tt, "\n"))
        }
        
        endPointCosts <- rep(NA,length(endPoints)) ## costs at those end Points
        endPointType <- rep(NA,length(endPoints)) ## what type is it?

        
        for(ii in seq_along(endPoints)){
            jj <- endPoints[ii]
            if(jj == 0){ jjCost <- 0 }else{ jjCost <- part$cost[jj] }


            endPointCosts[ii] <- jjCost + fCost$baseCost(jj+1,tt,0)
            endPointType[ii] <- "background"
            
            if( jj == tt-1 ){ ## then could be a point anomaly
                tmp <- jjCost + fCost$pointCost(tt,part$betaP)
                if(tmp < endPointCosts[ii]){
                    endPointCosts[ii] <- tmp
                    endPointType[ii] <- "point"
                }
            }else if(jj > tt - part$min_length){
                ## could be collective but to short a period
                endPointCosts[ii] <- NA
            }else{
                endPointCosts[ii] <- jjCost + fCost$collectiveCost(jj+1,tt,part$beta) 
                endPointType[ii] <- "collective"
            }
        }
        

        
        ## for(ii in seq_along(endPoints)){
        ##     jj <- endPoints[ii]
        ##     if(jj == 0){ jjCost <- 0 }else{ jjCost <- part$cost[jj] }
            
        ##     if( jj == tt-1 ){
        ##         ## test adding as base or point anomaly
        ##         ## compute C2 from paper
        ##         endPointCosts[ii] <- jjCost + fCost$baseCost(tt,tt,0)
        ##         endPointType[ii] <- "background"
        ##         ## compute C3 from paper
        ##         tmp <- jjCost + fCost$pointCost(tt,part$betaP)
        ##         if(tmp < endPointCosts[ii]){
        ##             endPointCosts[ii] <- tmp
        ##             endPointType[ii] <- "point"
        ##         }
        ##     }else if(jj > tt - part$min_length){
        ##         ## could be collective but to short a period
        ##         endPointCosts[ii] <- NA
        ##     }else{
        ##         endPointCosts[ii] <- jjCost + fCost$collectiveCost(jj+1,tt,part$beta) 
        ##         endPointType[ii] <- "collective"
        ##     }
        ## }
        

        if( all(is.na(endPointCosts)) ){ next } ## can't evaluate at tt

        ## find minimum
        idx <- which.min(endPointCosts) ## ignores NA
        part$endPoint[tt] <- endPoints[idx]
        part$cost[tt] <- endPointCosts[idx]
        part$type[tt] <- endPointType[idx]

        if(prune){
            idx <- is.na(endPointCosts) | (endPointCosts <= part$cost[tt] + cnst)
            endPoints <- endPoints[ idx ]
        }

        endPoints <- c(endPoints, tt)
        
    }
    
    return(part)
}

