#include <stdio.h>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "model_NHBSA.h"
#include "global.h"
#include "statistics.h"
#include "sort.h"


Model_NHBSA::Model_NHBSA(): 
NodeHistogramMatrix(*(new Array2D<double>()))
{
    init(0, NULL, 0);
}

Model_NHBSA::Model_NHBSA(int n_ell, Chromosome* pop, int size): 
NodeHistogramMatrix(*(new Array2D<double>(n_ell, n_ell)))
{
    init(n_ell, pop, size);
}

void Model_NHBSA::init(int n_ell, Chromosome* pop, int size)
{
	MODEL_NAME = "NHBSA";
	
	problemSize = n_ell;
	totalHistogram = 0.0;
    population = pop;
    n_size = size;
	initNHM();
}

Model_NHBSA::~Model_NHBSA()
{

	delete &NodeHistogramMatrix;
}

void Model_NHBSA::initNHM()
{

}

void Model_NHBSA::clearNHM()
{
	totalHistogram = 0.0;
	epsilon = ( n_size/ problemSize) * Bratio; 
	
	for(int i = 0; i < problemSize; i++)
		for(int j=0; j < problemSize; j++)
			NodeHistogramMatrix[i][j] = 0.0;	
}

void Model_NHBSA::build()
{
	if(problemSize == 0)
	{
		printf("ERROR: Problem size cann't be 0.");
		exit(1);
	}
	clearNHM();
	
	for(int i = 0; i < n_size; i++)
	{

		const Chromosome* templateX = &population[i];

		for(int j = 0; j < problemSize; j++)
		{

			int index = (int) templateX->getVal(j);
			double& cell = NodeHistogramMatrix[j][index];
			cell = cell+1;
			totalHistogram++;
		}

	}
}

Chromosome Model_NHBSA::resample()
{
	int templateIndex = random.uniformInt(0, n_size-1);
	return this->resampleWithTemplate(&population[templateIndex]);

}

Chromosome Model_NHBSA::resampleWithTemplate(const Chromosome* temp)
{
	Chromosome templateX;
	templateX.clone(*temp);


	int *cutPoints = new int[NUM_CUTPOINTS];
	random.uniformArray(cutPoints, NUM_CUTPOINTS, 1, problemSize-1);
	Sort::bubble(cutPoints, NUM_CUTPOINTS);	

	int resampledSegIndex = random.uniformInt(0, NUM_CUTPOINTS-1);
	int startIndex = cutPoints[resampledSegIndex];
	int endIndex = (resampledSegIndex == NUM_CUTPOINTS-1) ? 
			problemSize - 1 : cutPoints[resampledSegIndex+1];
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
	

	
	/*- roulette wheel for each empty position-*/
	for(int i = startIndex; i <= endIndex; i++)
	{
		int position = r[i % problemSize];
		double totalHistogram  = 0.0;
		double* pArray = new double[numResampled];
		for(int j = 0; j < numResampled; j++)
		{
			if(candidateUsed[candidate[j]] == true) continue;
			else{
				
				double P = NodeHistogramMatrix[position][candidate[j]] + epsilon;
					
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


void Model_NHBSA::dump()
{
    std::cout << "Model_NHBSA :\n";
    for(int i = 0; i < this->problemSize; i++)
    {
		for(int j = 0; j < this->problemSize; j++)
		{
			std::cout<< std::setw(7)<< NodeHistogramMatrix[i][j];
		}

		std::cout << std::endl;
    }
}


void Model_NHBSA::update()
{
}

/*------------------NHBSAwoTemplate PART-----------------*/
Chromosome Model_NHBSAwoTemplate::resample()
{
	Chromosome templateX;
	templateX.clone(population[0]);

	bool *candidateUsed = new bool [problemSize];
	

	for(int g = 0; g < problemSize; g++)
	{
		templateX.setVal(g, -1);
		candidateUsed[g] = false;
	}
	
	/*- roulette wheel for each empty position-*/
	for(int i = 0; i < problemSize; i++)
	{
		int position = i;
		double totalHistogram  = 0.0;
		double* pArray = new double[problemSize];
		for(int j = 0; j < problemSize; j++)
		{
			if(candidateUsed[j] == true) continue;
			else{
				
				double P = NodeHistogramMatrix[position][j] + epsilon;
					
				totalHistogram += P;
				pArray[j] = P;
			}
		}
		
		double pointer = totalHistogram * random.uniform();
		double partialSum = 0.0;
		int index = -1;
		
		for(int j = 0; j < problemSize; j++)
		{
			if(candidateUsed[j] == true) continue;
			partialSum += pArray[j];
			
			if(partialSum >= pointer)
			{
				index = j;
				break;
			}
		}   
		delete[] pArray;		
		
		templateX.setVal(position, index);
		candidateUsed[index] = true;	
	}	
	
	delete[] candidateUsed;
	
	return templateX;
}
