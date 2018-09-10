#include <cstdio>
#include "dice.h"

Dice::Dice()
{
    space = 0;
}

Dice::Dice(std::vector<double> _distVector)
{
    this->setDistribution(_distVector);
}

Dice::Dice(double* array, int length)
{
    std::vector<double> vec(array, array + length);
    this->setDistribution(vec);
}

Dice::Dice(int* array, int length)
{
	std::vector<double> vec;
	for(int i = 0; i < length; i++)
		vec.push_back((double)array[i]);

    this->setDistribution(vec);
}

void Dice::setDistribution(std::vector<double>& _distVector)
{
    distVector = _distVector;
    space = (int)distVector.size();

    double total = 0.0;

    for(int i = 0; i < space; i++)
    {
        total += distVector[i];
    }
    if(total <= 0.0)
    {
      //error
    }else{
        for(int i = 0; i < space; i++)
          distVector[i] = distVector[i] / total;
    }
}

double Dice::sample()
{
    return this->sample(0.0);
}

double Dice::sample(double offset)
{
    double partialSum = 0.0;
    double pointer = random.uniform();
    int i = 0;
    for(i = 0; i < space; i++)
    {
        partialSum += distVector[i] + offset;
        if(partialSum >= pointer)
        {  
          break;
        }
    }

    return i;
}

