#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include "model_Hybrid.h"
#include "global.h"
#include "statistics.h"
#include "sort.h"

Model_Hybrid::Model_Hybrid()
{
    init(0, NULL, 0);
}

Model_Hybrid::Model_Hybrid(double _r1, double _r2, int n_ell, Chromosome* _pop, int p_size) 
: Model_EHBSA(n_ell, _pop, p_size), Model_NHBSA(n_ell, _pop, p_size)
{
	problemSize = n_ell;
    init(n_ell, _pop, p_size);
	r1 = _r1;
	r2 = _r2;
}

void Model_Hybrid::init(int n_ell, Chromosome* _pop, int p_size)
{
    MODEL_NAME = "Hybrid";
	problemSize = n_ell;
    population = _pop;
    n_size = p_size;
}

void Model_Hybrid::build()
{
	Model_EHBSA::build();
	Model_NHBSA::build();
}

void Model_Hybrid::update()
{
	Model_EHBSA::update();
	Model_NHBSA::update();
}

Chromosome Model_Hybrid::resample()
{
	int templateIndex = random.uniformInt(0, n_size-1);
	return this->resampleWithTemplate(&population[templateIndex]);
}

Chromosome Model_Hybrid::resampleWithTemplate(const Chromosome* temp)
{
	return weightResample(temp);
}

Chromosome Model_Hybrid::weightResample(const Chromosome* temp)
{
    MyRand random;

    double token = random.uniform();
    if(token <= r1)
    {
		return Model_EHBSA::resampleWithTemplate(temp);
    }
    else{
		return Model_NHBSA::resampleWithTemplate(temp);
    }
}

Chromosome Model_Hybrid::bitwiseResample(const Chromosome* temp)
{
	Chromosome templateX;
	
	templateX.clone(*temp);

	int *cutPoints = new int[NUM_CUTPOINTS];
	random.uniformArray(cutPoints, NUM_CUTPOINTS, 1, problemSize-1);
	Sort::bubble(cutPoints, NUM_CUTPOINTS);	

	int resampledSegIndex = random.uniformInt(0, NUM_CUTPOINTS-1);
	
	int startIndex = cutPoints[resampledSegIndex];
	int endIndex = (resampledSegIndex == NUM_CUTPOINTS-1) ? 
			problemSize + cutPoints[0] - 1 : cutPoints[resampledSegIndex+1];
	int numResampled = endIndex - startIndex + 1;  //Lwt
	int *r = new int[problemSize];  // random permutation of position index
	random.uniformArray(r, problemSize, 0, problemSize-1);

	
	//generate candidate (if -1: used, ignored by following code)	
	int *candidate = new int[numResampled];
	bool *candidateUsed = new bool [problemSize];
	
	for(int g = 0; g < problemSize; g++)
		candidateUsed[g] = true;
	for(int g = 0; g < numResampled; g++)
	{
		int index = r[(startIndex+g) % problemSize];
		candidate[g] = templateX.getVal(index);
		templateX.setVal(index, -1);
		candidateUsed[candidate[g]] = false;
	}
	

	
	//- roulette wheel for each empty position-//
	for(int i = startIndex; i <= endIndex; i++)
	{
		int position = r[i % problemSize];
		double totalHistogram  = 0.0;
		double* pArray = new double[numResampled];
		for(int j = 0; j < numResampled; j++)
		{
			if(candidateUsed[candidate[j]] == true) continue;
			else{
				
				double nodeP = Model_NHBSA::NodeHistogramMatrix[position][candidate[j]] + Model_NHBSA::epsilon;
				
				double leftP = 1.0, rightP = 1.0;
				int left = (position-1 + problemSize) % problemSize;
				int right = (position+1) % problemSize;
				
				int leftIndex = templateX.getVal(left);
				int rightIndex = templateX.getVal(right);
				
				int centreIndex = candidate[j];

				if(leftIndex != -1)
					leftP *= (EdgeHistogramMatrix[ leftIndex ][centreIndex] + EdgeHistogramMatrix[ centreIndex ][ leftIndex ] + Model_EHBSA::epsilon);
					
				if(rightIndex != -1)
					rightP *= (EdgeHistogramMatrix[ rightIndex ][centreIndex] + EdgeHistogramMatrix[ centreIndex ][ rightIndex ] + Model_EHBSA::epsilon);
				
				double egdeP = leftP * rightP;
				
				//Linear Combination
				double P = r1 * egdeP + r2 * nodeP;
				
				totalHistogram += P;
				pArray[j] = P;
			}
		}
		
		double pointer = totalHistogram * random.uniform();
		double partialSum = 0.0;
		int index = -1;
		
		for(int j = 0; j < numResampled; j++)
		{
			if(candidateUsed[candidate[j]] == true) continue;
			partialSum += pArray[j];
			
			if(partialSum >= pointer)
			{
				index = candidate[j];
				break;
			}
		}   
		delete[] pArray;		
		
		templateX.setVal(position, index);
		candidateUsed[index] = true;	
	}	
	
	delete[] candidateUsed;
	delete[] candidate;
	delete[] cutPoints;
	delete[] r;
	

	return templateX;	
}
void Model_Hybrid::dump()
{
	Model_EHBSA::dump();
	Model_NHBSA::dump();
}


void Model_Hybrid::clearHM()
{
	Model_EHBSA::clearEHM();
	Model_NHBSA::clearNHM();
}



