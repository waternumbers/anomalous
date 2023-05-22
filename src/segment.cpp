#include "segment.h"

// template class functions
//segment::segment(){ }
//void segment::update(double& x, double& mu, double &sigma ){};

// Gaussian with known mean and variance
gaussKnown::gaussKnown(double beta_, int t_){
  beta = beta_;
  start = t_;
};

void gaussKnown::update( double& x, double& mu, double &sigma ){
  // update summary statistics
  n += 1.0;
  summaryStats[0] += 1.0/sigma;
  summaryStats[1] += std::log(sigma);
  summaryStats[2] += (x-mu)/sigma;
  summaryStats[3] += std::pow(x-mu,2.0)/sigma;
    
  double kappa = (summaryStats[3] - 2.0*param[0]*summaryStats[2] +
		  std::pow(param[0],2.0)*summaryStats[0] ) / n;

  cost = n *std::log( 2.0*M_PI*param[1]) + summaryStats[1] +
    (n*kappa/param[1]) + beta;

};


// Gaussian with unknown mean shift but known variance
gaussMean::gaussMean( double beta_, int t_){
  beta = beta_;
  start = t_;
};

void gaussMean::update( double& x, double& mu, double &sigma ){
  // update summary statistics
  n += 1.0;
  summaryStats[0] += 1.0/sigma;
  summaryStats[1] += std::log(sigma);
  summaryStats[2] += (x-mu)/sigma;
  summaryStats[3] += std::pow(x-mu,2.0)/sigma;
  param[0] = summaryStats[2] / summaryStats[0];
    
  double kappa = (summaryStats[3] - 2.0*param[0]*summaryStats[2] +
		  std::pow(param[0],2.0)*summaryStats[0] ) / n;
  cost = n *std::log( 2.0*M_PI*param[1]) + summaryStats[1] +
    (n*kappa/param[1]) + beta;

};









// class gaussMeanVarCpp {
// public:
//   int start;
//   double n, penalty, cost;
//   std::vector<double> summaryStats = {0.0,0.0,0.0,0.0};
//   std::vector<double> param = {-999.9,-999.9};
//   gaussMeanVarCpp(int s_, double p_): start(s_), penalty(p_){}
//   void update(double &x, double &mu, double &sigma){
//     this->n += 1.0;
//     this->summaryStats[0] += 1.0/sigma;
//     this->summaryStats[1] += std::log(sigma);
//     this->summaryStats[2] += (x-mu)/sigma;
//     this->summaryStats[3] += std::pow(x-mu,2.0)/sigma;

//     this->param[0] = this->summaryStats[2] / this->summaryStats[0];
	
//     double kappa = (this->summaryStats[3] - 2.0*this->param[0]*this->summaryStats[2] +
// 		    std::pow(this->param[0],2.0)*this->summaryStats[0] ) / this->n;
//     kappa = std::max(kappa, std::numeric_limits<double>::min());
//     this->param[1] = kappa;

//     this->cost = this->n *std::log( 2.0*M_PI*this->param[1]) + this->summaryStats[1] +
//       (this->n*kappa/this->param[1]) + this->penalty;
//   }
// };


// class gaussMeanCpp {
// public:
//   int start;
//   double n, penalty, cost;
//   std::vector<double> summaryStats = {0.0,0.0,0.0,0.0};
//   std::vector<double> param = {-999.9,1.0};
//   gaussMeanCpp(int s_, double p_): start(s_), penalty(p_){}
//   void update(double &x, double &mu, double &sigma){
//     this->n += 1.0;
//     this->summaryStats[0] += 1.0/sigma;
//     this->summaryStats[1] += std::log(sigma);
//     this->summaryStats[2] += (x-mu)/sigma;
//     this->summaryStats[3] += std::pow(x-mu,2.0)/sigma;

//     this->param[0] = this->summaryStats[2] / this->summaryStats[0];
    
//     double kappa = (this->summaryStats[3] - 2.0*this->param[0]*this->summaryStats[2] +
// 		    std::pow(this->param[0],2.0)*this->summaryStats[0] ) / this->n;
    
    
//     this->cost = this->n *std::log( 2.0*M_PI*this->param[1]) + this->summaryStats[1] +
//       (this->n*kappa/this->param[1]) + this->penalty;
//   }
// };

// class gaussVarCpp {
// public:
//   int start;
//   double n, penalty, cost;
//   std::vector<double> summaryStats = {0.0,0.0,0.0,0.0};
//   std::vector<double> param = {0.0,-999.9};
//   gaussVarCpp(int s_, double p_): start(s_), penalty(p_){}
//   void update(double &x, double &mu, double &sigma){
//     this->n += 1.0;
//     this->summaryStats[0] += 1.0/sigma;
//     this->summaryStats[1] += std::log(sigma);
//     this->summaryStats[2] += (x-mu)/sigma;
//     this->summaryStats[3] += std::pow(x-mu,2.0)/sigma;
//     double kappa = (this->summaryStats[3] - 2.0*this->param[0]*this->summaryStats[2] +
// 		    std::pow(this->param[0],2.0)*this->summaryStats[0] ) / this->n;
//     kappa = std::max(kappa, std::numeric_limits<double>::min());
//     this->param[1] = kappa;
    
//     this->cost = this->n *std::log( 2.0*M_PI*this->param[1]) + this->summaryStats[1] +
//       (this->n*kappa/this->param[1]) + this->penalty;
//   }
// };


// class gaussFixedCpp {
// public:
//   int start;
//   double n, penalty, cost;
//   std::vector<double> summaryStats = {0.0,0.0,0.0,0.0};
//   std::vector<double> param = {0.0,1.0};
//   gaussFixedCpp(int s_, double p_): start(s_), penalty(p_){}
//   void update(double &x, double &mu, double &sigma){
//     this->n += 1.0;
//     this->summaryStats[0] += 1.0/sigma;
//     this->summaryStats[1] += std::log(sigma);
//     this->summaryStats[2] += (x-mu)/sigma;
//     this->summaryStats[3] += std::pow(x-mu,2.0)/sigma;
    
//     double kappa = (this->summaryStats[3] - 2.0*this->param[0]*this->summaryStats[2] +
// 		    std::pow(this->param[0],2.0)*this->summaryStats[0] ) / this->n;
    
//     this->cost = this->n *std::log( 2.0*M_PI*this->param[1]) + this->summaryStats[1] +
//       (this->n*kappa/this->param[1]) + this->penalty;
//   }
// };


// class gaussPointCpp {
// public:
//   int start;
//   double n, penalty, cost;
//   std::vector<double> summaryStats = {0.0,0.0,0.0,0.0};
//   std::vector<double> param = {0.0,-999.9};
//   gaussPointCpp(int s_, double p_): start(s_), penalty(p_){}
//   void update(double &x, double &mu, double &sigma){
//     this->n += 1.0;
//     this->summaryStats[0] += 1.0/sigma;
//     this->summaryStats[1] += std::log(sigma);
//     this->summaryStats[2] += (x-mu)/sigma;
//     this->summaryStats[3] += std::pow(x-mu,2.0)/sigma;
//     double kappa = (this->summaryStats[3] - 2.0*this->param[0]*this->summaryStats[2] +
// 		    std::pow(this->param[0],2.0)*this->summaryStats[0] ) / this->n;
//     kappa = std::max(kappa, std::numeric_limits<double>::min());
//     this->param[1] = kappa;
    
//     this->cost = std::log( 2.0*M_PI) + std::log(sigma) + std::log( std::exp(-this->penalty) + kappa ) + 1 + this->penalty;
//   }
// };


// //RCPP_MODULE(gaussSeg){
// //using namespace Rcpp;
  
//   // class_<gaussMeanVarCpp>("gaussMeanVarCpp")
//   //   .constructor<int,double>()
//   //   .method("update", &gaussMeanVarCpp::update)
//   //   .field("start", &gaussMeanVarCpp::start)
//   //   .field("param", &gaussMeanVarCpp::param)
//   //   .field("cost", &gaussMeanVarCpp::cost)
//   //   .field("n", &gaussMeanVarCpp::n)
//   //   .field("summaryStats", &gaussMeanVarCpp::summaryStats)
//   //   ;
//   //class_<gaussMean>("gaussMean")
//     .constructor<double,int>()
//     .method("update", &gaussMean::update)
//     .field("start", &gaussMean::start)
//     .field("param", &gaussMean::param)
//     .field("cost", &gaussMean::cost)
//     .field("n", &gaussMean::n)
//     .field("summaryStats", &gaussMean::summaryStats)
//     ;
//   // class_<gaussVarCpp>("gaussVarCpp")
//   //   .constructor<int,double>()
//   //   .method("update", &gaussVarCpp::update)
//   //   .field("start", &gaussVarCpp::start)
//   //   .field("param", &gaussVarCpp::param)
//   //   .field("cost", &gaussVarCpp::cost)
//   //   .field("n", &gaussVarCpp::n)
//   //   .field("summaryStats", &gaussVarCpp::summaryStats)
//   //   ;
//   // class_<gaussFixedCpp>("gaussFixedCpp")
//   //   .constructor<int,double>()
//   //   .method("update", &gaussFixedCpp::update)
//   //   .field("start", &gaussFixedCpp::start)
//   //   .field("param", &gaussFixedCpp::param)
//   //   .field("cost", &gaussFixedCpp::cost)
//   //   .field("n", &gaussFixedCpp::n)
//   //   .field("summaryStats", &gaussFixedCpp::summaryStats)
//   //   ;
//   // class_<gaussPointCpp>("gaussPointCpp")
//   //   .constructor<int,double>()
//   //   .method("update", &gaussPointCpp::update)
//   //   .field("start", &gaussPointCpp::start)
//   //   .field("param", &gaussPointCpp::param)
//   //   .field("cost", &gaussPointCpp::cost)
//   //   .field("n", &gaussPointCpp::n)
//   //   .field("summaryStats", &gaussPointCpp::summaryStats)
//   //   ;
// }
