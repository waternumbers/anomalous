## basic segment class
setClass("Segment", 
         slots = c(start = "integer",
                   n = "integer",
                   summaryStats = "numeric",
                   param = "numeric",
                   penalty = "numeric",
                   cost = "numeric"),
         prototype = list(
             start = NA_integer_,
             n = integer(1),
             summaryStats = NA_real_,
             param = NA_real_,
             penalty = NA_real_,
             cost = NA_real_)
         )

##setMethod("cost","Segment",function(obj){return(obj@cost)})
