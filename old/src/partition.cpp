#include "partition.h"


partition::partition(){
  cost = 0.0;
};
void partition::addNew(double beta, int t){
  part.push_back( gaussMean(beta,t) );
  int ii = part.size() - 1;
  cost += part[ii].cost;
}

void partition::update(double& x, double& mu, double &sigma, double &minLen, double &maxLen){
  uint ii = part.size() - 1;
  part[ii].update(x,mu,sigma);

  // this is ugly ( and maybe also slow .....)
  cost = 0.0;
  for(long unsigned int ii=0; ii < part.size(); ii++){
    cost += part[ii].cost;
    if( (part[ii].n < minLen) or (part[ii].n>maxLen) ){
      cost = std::numeric_limits<double>::max();
      ii = part.size();
    }
  }
  
};

