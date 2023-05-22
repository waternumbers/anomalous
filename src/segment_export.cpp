// segment_export.cpp

// Include Rcpp system header file (e.g. <>)
#include <Rcpp.h>

// Include our definition of the student file (e.g. "")
#include "segment.h"

// Expose gaussMean class
RCPP_MODULE(segments) {             // Name used to "loadModule" in R script
  Rcpp::class_<gaussKnown>("gaussKnown")
    .field("cost", &gaussKnown::cost)
    .field("n", &gaussKnown::n)
    .field("start", &gaussKnown::start)
    .field("beta", &gaussKnown::beta)
    .field("param", &gaussKnown::param)
    .constructor<double, int>()
    .method("update", &gaussKnown::update);
  Rcpp::class_<gaussMean>("gaussMean")       // This must be the C++ class name.
    .field("cost", &gaussMean::cost)
    .field("n", &gaussMean::n)
    .field("start", &gaussMean::start)
    .field("beta", &gaussMean::beta)
    .field("param", &gaussMean::param)
    .constructor<double, int>()
    .method("update", &gaussMean::update);
}
