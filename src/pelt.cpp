#include <Rcpp.h>
#include <vector>
#include <math.h>
#include <map>
#include "segment.h"
#include "partition.h"

// [[Rcpp::export]]
Rcpp::List peltc(std::vector<double> &x,std::vector<double> &mu,std::vector<double> &sigma,
		 double &beta, double &minLen, double &maxLen){
  
  long unsigned int n = x.size();
  std::vector<partition> ctlg;
  partition opt;

  for(long unsigned int tt = 0; tt < n; tt++){
    // add new break point
    ctlg.push_back( opt );
    long unsigned int ii = ctlg.size()-1;
    ctlg[ii].addNew(beta,tt);

    // update costs
    double minCost = std::numeric_limits<double>::max();
    for(long unsigned int ii = 0; ii < ctlg.size(); ii++){
      ctlg[ii].update(x[tt],mu[tt],sigma[tt], minLen, maxLen);
	
      if( (ii==0) or (ctlg[ii].cost < minCost) ){
	opt = ctlg[ii];
	minCost = ctlg[ii].cost;
      }
    }
    
    // trim
    double thold = minCost + beta;
    auto it = std::remove_if(ctlg.begin(), ctlg.end(), [thold](partition x){return x.cost > thold;});
    ctlg.erase(it, ctlg.end());
    
  }

  // //std::vector<double> c = {1.0,2,3,4,5,6,7};
  // double thold = 395.0;
  
  // auto it = std::remove_if(ctlg.begin(), ctlg.end(), [thold](partition x){return x.cost < thold;});
  // ctlg.erase(it, ctlg.end());

  // Rcpp::Rcout << thold << std::endl;
  // for(uint ii = 0; ii < ctlg.size(); ii++){
  //   Rcpp::Rcout << ctlg[ii].cost << std::endl;
  // }
  
  // // could use erase_if in C++20
  // //std::erase_if(c,[](double x){x<6;});
  
  // int out = ctlg.size(); //ctlg.size();
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
