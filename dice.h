#ifndef _DICE_H_
#define _DICE_M_

#include <vector>
#include "myrand.h"

class Dice{
  private:
    int space;
    std::vector<double> distVector;
    MyRand random;
  public:
    Dice();
    Dice(std::vector<double> distVector);
    Dice(double* array, int length);
	Dice(int* array, int length);
      
    void setDistribution(std::vector<double>& distVector);

    double sample();
    double sample(double offset);
};

#endif
