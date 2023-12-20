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
pelt <- function(part,fCost,prune = TRUE,verbose=FALSE,...){
    ctlg <- list()
    ctlg[[1]] <- part ## offset by 1 versus time!!

    cnst <- max(ctlg[[1]]$beta)
    
    for(tt in fCost$validTimes){ ##maxT){
        if(verbose && (tt %% 100==0)) {
            ## Print on the screen some message
            cat(paste0("time step: ", tt, "\n"))
        }

        opt <- addCollective(ctlg[[1]], fCost, ctlg[[1]]$last_time + 1, tt,...)
        ## loop ctlg
        ctlgCost <- rep(-Inf,length(ctlg))
        for(ii in 1:length(ctlg)){
            if(ctlg[[ii]]$last_time > tt - ctlg[[ii]]$min_length){ next }
            
            tmp <- addCollective(ctlg[[ii]], fCost, ctlg[[ii]]$last_time + 1, tt,...)
            ctlgCost[ii] <- tmp$cost
            if(tmp$cost < opt$cost){ opt <- tmp }

        }
        
        if(prune){ ctlg <- ctlg[ ctlgCost <= opt$cost+cnst ] }
        
        ctlg[[length(ctlg)+1]] <- opt  
    }
    return(opt)
}

