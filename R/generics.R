#' Return the cost of a segment
#' @param obj a segment class object
#' @rdname cost-methods
#' @export
setGeneric("cost", 
           function(obj) standardGeneric("cost"),
           signature = "obj"
           )

#' Update a segment with new observations
#' @param obj a segment class object
#' @param x new observations to add
#' @rdname update-methods
#' @export
setGeneric("update",
           function(obj,x,...) standardGeneric("update"),
           signature = "obj"
           )
