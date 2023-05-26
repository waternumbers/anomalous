#' An R implimentation of the segmented search algorithm
#' @param cost a cost function
#'
#' @return the optimal partition
#' 
#' @details Basic R implimentation - not efficent
#'@export
capa <- function(fCost,prune = FALSE){
    ##p <- profvis::profvis({
    ctlg <- list()
    ctlg[[1]] <- new("partition") ## offset by 1 versus time!!
    for(tt in 1:length(x)){
        if(tt %% 100==0) {
            ## Print on the screen some message
            cat(paste0("time step: ", tt, "\n"))
        }
        
        ## compute C2 from paper
        opt <- addBase(ctlg[[tt]],fCost,tt,tt)
        
        ## compute C3 from paper
        tmp <- addPoint(ctlg[[tt]],fCost,tt)
        if(tmp@cost < opt@cost){ opt <- tmp }
        
        ## loop ctlg
        for(ii in 1:tt){
            tmp <- addCollective(ctlg[[ii]], fCost,ii,tt)
            if(tmp@cost < opt@cost){ opt <- tmp }
        }
        
        ctlg[[tt+1]] <- opt  
    }
    return(opt)
}

