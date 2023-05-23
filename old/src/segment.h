#ifndef SEG
#define SEG

#include <vector>
#include <utility>
#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>

// ideally we'd use a polymorphic class here
// containing start,n,cost,beta,summaryStats, beta
// but I can't get it to play nicely with Rcpp
// so everything is declared within the classes....

// // generic class
// class segment {
// public:
//    // initialisation
//   segment();
//   int start;
//   double n = 0.0, cost = 0.0;
//   double beta;
//   std::vector<double> summaryStats;
//   std::vector<double> param;
//   virtual void update(double&, double&, double&); // update with new observation
// };

// gaussKnown segment
class gaussKnown { //: public segment {
public:
  int start;
  double n = 0.0, cost = 0.0;
  double beta;
  gaussKnown(double, int);
  std::vector<double> param = {0.0,1.0};
  std::vector<double> summaryStats = {0.0,0.0,0.0,0.0};
  void update(double&, double&, double&);
};

// gaussMean segment
class gaussMean { //: public segment {
public:
  int start;
  double n = 0.0, cost = 0.0;
  double beta;
  gaussMean(double, int);
  std::vector<double> param = {-999.9,1.0};
  std::vector<double> summaryStats = {0.0,0.0,0.0,0.0};
  void update(double&, double&, double&);
};


#endif
