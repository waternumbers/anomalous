#ifndef PART
#define PART

#include <vector>
#include "segment.h"

// generic class
class partition {
public:
   // initialisation
  partition();
  std::vector<gaussMean> part;
  double cost;
  void addNew(double, int);
  void update(double&, double&, double&, double&, double&); // update with new observation
};


#endif
