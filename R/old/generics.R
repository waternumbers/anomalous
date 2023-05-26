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
#' @param ... other inputs
#' @rdname update-methods
#' @export
setGeneric("update",
           function(obj,x,...) standardGeneric("update"),
           signature = "obj"
           )

#' Add a point segment
#' @param obj a segment class object
#' @param x new observations to add
#' @param ... other inputs
#' @rdname update-methods
#' @export
setGeneric("addPoint",
           function(obj,...) standardGeneric("addPoint"),
           signature = "obj"
           )

#' Add a Collective segment
#' @param obj a segment class object
#' @param ... other inputs
#' @rdname update-methods
#' @export
setGeneric("addCollective",
           function(obj,...) standardGeneric("addCollective"),
           signature = "obj"
           )

#' Add a Baseline period
#' @param obj a segment class object
#' @param ... other inputs
#' @rdname update-methods
#' @export
setGeneric("addBase",
           function(obj,...) standardGeneric("addBase"),
           signature = "obj"
           )
