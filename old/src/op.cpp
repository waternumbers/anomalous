#include <Rcpp.h>
#include <vector>
#include <math.h>
#include "segment.h"
#include "partition.h"

// [[Rcpp::export]]
Rcpp::List opc(std::vector<double> &x,std::vector<double> &mu,std::vector<double> &sigma,
			double &beta, double &minLen, double &maxLen){

  long unsigned int n = x.size();
  std::vector<partition> ctlg(n);
  partition opt;

  for(long unsigned int tt = 0; tt < n; tt++){
    ctlg[tt] = opt;
    ctlg[tt].addNew(beta,tt);

    double minCost = std::numeric_limits<double>::max();
    for(long unsigned int ii = 0; ii <= tt; ii++){
      ctlg[ii].update(x[tt],mu[tt],sigma[tt], minLen, maxLen);
	
      if( (ii==0) or (ctlg[ii].cost < minCost) ){
	opt = ctlg[ii];
	minCost = ctlg[ii].cost;
      }
    }
  }

  Rcpp::List out;
  for(long unsigned int ii = 0; ii < opt.part.size(); ii++){
    Rcpp::List L = Rcpp::List::create(Rcpp::Named("start") = opt.part[ii].start + 1.0,
				      Rcpp::Named("n") = opt.part[ii].n,
				      Rcpp::Named("cost") = opt.part[ii].cost,
				      Rcpp::Named("beta") = opt.part[ii].beta,
				      Rcpp::Named("param") = opt.part[ii].param,
				      Rcpp::Named("summaryStats") = opt.part[ii].summaryStats);
    out.push_back( L );
  }

  
  
  return(out);

};
