#ifndef SEG
#define SEG

#include <vector>
#include <utility>
#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>

// generic class
class segment {
public:
   // initialisation
  segment();
  int start;
  double n = 0.0, cost = 0.0;
  double beta;
  std::vector<double> summaryStats;
  std::vector<double> param;
  virtual void update(double&, double&, double&); // update with new observation
};

// gaussMean segment
class gaussMean: public segment {
public:
  gaussMean(double, int);
  std::vector<double> param = {-999.9,1.0};
  std::vector<double> summaryStats = {0.0,0.0,0.0,0.0};
  void update(double&, double&, double&);
};


#endif
