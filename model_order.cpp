#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "model_order.h"
#include "global.h"
#include "dice.h"
#include "sort.h"

Model_Order::Model_Order() :
histogram(*(new Array2D<double>()))
{
	init(0, NULL, 0);
}

Model_Order::Model_Order(int n_ell, Chromosome* pop, int size):
histogram(*(new Array2D<double>(n_ell, n_ell+1)))
{
	init(n_ell, pop, size);
}

Model_Order::~Model_Order()
{
	delete &histogram;
	
	delete[] firstPosition;
}

void Model_Order::init(int n_ell, Chromosome* pop, int size)
{
    MODEL_NAME = "ORDER";
    problemSize = n_ell;
    population = pop;
    n_size = size;
	firstPosition = new int[problemSize];
}

void Model_Order::clear()
{
	for(int i = 0; i < problemSize; i++)
	{
		for(int j=0; j < problemSize+1; j++)
		{
			histogram[i][j] = 0.0;  	    
		}
		
		firstPosition[i] = 0;
	}	
}

void Model_Order::build()
{
	if(problemSize == 0)
	{
		printf("ERROR: Problem size of Model can't be 0.");
		exit(1);
	}
    epsilon = (n_size / 2) * Bratio;
	clear();
	
	for(int i = 0; i < n_size; i++)
	{
		const Chromosome* t = &population[i];
		for(int k = 0; k < problemSize; k++)
		{
			if(k == 0) continue;
			
			for(int l = 0; l < k; l++)
			{
				int a = (int) t->getVal(l);
				int b = (int) t->getVal(k);		
				
				histogram[a][b]++;
			}
		}	
		
		int first = (int) t->getVal(0);
		firstPosition[first] += 1;
	}

  for(int i = 0; i < problemSize; i++)
  {
    double margin = 0.0;
    for(int j = 0; j < problemSize; j++)
    {
      margin += histogram[i][j];
    }
    histogram[i][problemSize] = margin;
  }
}

Chromosome Model_Order::resample()
{
	
	int templateIndex = random.uniformInt(0, n_size-1);
	return this->resampleWithTemplate(&population[templateIndex]);
}

Chromosome Model_Order::resampleWithTemplate(const Chromosome* temp)
{
	Chromosome t;
	t.clone(*temp);
	
	int *cutPoints = new int[NUM_CUTPOINTS];
	random.uniformArray(cutPoints, NUM_CUTPOINTS, 1, problemSize-1);
	Sort::bubble(cutPoints, NUM_CUTPOINTS);

	int resampledSegIndex = random.uniformInt(0, NUM_CUTPOINTS-1);
	
	int startIndex = cutPoints[resampledSegIndex];
	int endIndex = (resampledSegIndex == NUM_CUTPOINTS-1) ? 
			problemSize - 1 : cutPoints[resampledSegIndex+1];
	int numResampled = endIndex - startIndex + 1;
	
	//generate window (if -1: used, ignored by following code)	
	int *window = new int[numResampled];
	bool *used = new bool [problemSize];
	
	for(int g = 0; g < problemSize; g++)
		used[g] = true;
	for(int g = 0; g < numResampled; g++)
	{
		int index = (startIndex+g) % problemSize;
		window[g] = t.getVal(index);
		t.setVal(index, -1);
		used[window[g]] = false;
	}
	
	for(int i = startIndex; i <= endIndex; i++)
	{
		int index = -1;
		if(i == 0)
		{
      std::vector<double> v(firstPosition, firstPosition + n_size);

      for(int j = 0; j < problemSize; j++)
      {
        if(used[j] == true)
          v[j] = 0.0;
        else
          v[j] += epsilon;
      }

            Dice d(v);
            index = d.sample();
		}else
		{
      std::vector<double> v;
      for(int j = 0; j < problemSize; j++)
      {
        double p = 1.0;
        if(used[j] == true)
          p = 0.0;
        else
        {
          for(int k = 0; k < i; k++)
          {
            int a = t.getVal(k);
            double margin = histogram[a][problemSize];
            double s = histogram[a][j]; 
            p *= s/margin; 

          }
        }

        v.push_back(p + epsilon);
      }

        
      Dice d(v);
      index = d.sample();
		}
			
		t.setVal(i, index);
		used[index] = true;
	}	
	
	delete[] used;
	delete[] window;
	delete[] cutPoints;


	
	
	return t;
}


void Model_Order::update()
{

}




void Model_Order::dump()
{
    std::cout << "Model_Order:\n";
    for(int i = 0; i < this->problemSize; i++)
    {
		for(int j = 0; j < this->problemSize; j++)
		{
			std::cout<< std::setw(5)<< histogram[i][j];
		}

		std::cout << std::endl;
    }
}
