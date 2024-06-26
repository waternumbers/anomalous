#' An R implimentation of the segmented search algorithm
#' @param part the starting partition
#' @param fCost the cost function
#' @param prune logical, should pruning be used
#' @param verbose logical, print out progress
#'
#' @return the optimal partition
#' 
#' @details Basic R implimentation of capa - not efficent
#' @export
capa <- function(part,fCost,prune = TRUE,verbose=FALSE,...){
    
    ctlg <- list()
    part$last_time <- 0 #fCost$validTimes[1]-1
    ctlg[[1]] <- part ## offset by 1 versus time!!

    cnst <- max(ctlg[[1]]$beta,ctlg[[1]]$betaP)

    for(tt in 1:fCost$maxT){
        if(verbose && (tt %% 100==0)) {
            ## Print on the screen some message
            cat(paste0("time step: ", tt, "\n"))
        }
        
        ## compute C2 from paper
        opt <- addBase(ctlg[[length(ctlg)]],fCost,tt,tt,...)

        if( is.na(opt$cost) ){next} ## since this is returned when tt is missing data
        
        ## compute C3 from paper
        tmp <- addPoint(ctlg[[length(ctlg)]],fCost,tt,...)
        if(tmp$cost < opt$cost){ opt <- tmp }
        
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

