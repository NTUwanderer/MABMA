#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include "model_EHBSA.h"
#include "global.h"
#include "statistics.h"
#include "sort.h"

Model_EHBSA::Model_EHBSA() :
EdgeHistogramMatrix(*(new Array2D<double>()))
{
    init(0, NULL, 0);
}

Model_EHBSA::Model_EHBSA(int n_ell, Chromosome* _pop, int p_size):
EdgeHistogramMatrix(*(new Array2D<double>(n_ell, n_ell)))
{
	init(n_ell, _pop, p_size);
}


void Model_EHBSA::init(int n_ell, Chromosome* _pop, int p_size)
{
    MODEL_NAME = "EHBSA";
	problemSize = n_ell;
	totalEdges = 0.0; 

    firstPositionMatrix = NULL;
    population = _pop;
    n_size = p_size;
    initEHM();    
}

Model_EHBSA::~Model_EHBSA()
{
	delete &EdgeHistogramMatrix;
	
	delete[] firstPositionMatrix;
}

void Model_EHBSA::initEHM()
{
	firstPositionMatrix = new int[problemSize];
}

void Model_EHBSA::clearEHM()
{
	totalEdges = 0.0;
	epsilon = ( n_size*2/ (problemSize-1)) * Bratio; 
	for(int i = 0; i < problemSize; i++)
	{
		for(int j=0; j < problemSize; j++)
		{
		    if(i == j) EdgeHistogramMatrix[i][j] = 0.0;
		    else{
				EdgeHistogramMatrix[i][j] = 0.0;  
		    }		    
		}
		
		firstPositionMatrix[i] = 0;
	}
}

void Model_EHBSA::build()
{
	if(problemSize == 0)
	{
		printf("ERROR: Problem size of Model can't be 0.");
		exit(1);
	}
	clearEHM();
	
	for(int i = 0; i < n_size; i++)
	{
		const Chromosome* templateX = &population[i];
		for(int k = 0; k < problemSize; k++)
		{
			int j = k % problemSize;
			int index = (int) templateX->getVal(j);
			int l = (j + 1)% problemSize;  
			int indexNext = (int)templateX->getVal(l);
			EdgeHistogramMatrix[index][indexNext] = EdgeHistogramMatrix[index][indexNext] + 1;
		}	
		
		int first = (int) templateX->getVal(0);
		firstPositionMatrix[first] += 1;
	}
}


Chromosome Model_EHBSA::resample()
{
	
	int templateIndex = random.uniformInt(0, n_size-1);
	return this->resampleWithTemplate(&population[templateIndex]);
}

Chromosome Model_EHBSA::resampleWithTemplate(const Chromosome* temp)
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
	int numResampled = endIndex - startIndex + 1;
	
	//generate candidate (if -1: used, ignored by following code)	
	int *candidate = new int[numResampled];
	bool *candidateUsed = new bool [problemSize];
	
	for(int g = 0; g < problemSize; g++)
		candidateUsed[g] = true;
	for(int g = 0; g < numResampled; g++)
	{
		int index = (startIndex+g) % problemSize;
		candidate[g] = templateX.getVal(index);
		templateX.setVal(index, -1);
		candidateUsed[candidate[g]] = false;
	}
	

	
	/*- roulette wheel for each empty position-*/
	for(int i = startIndex; i <= endIndex; i++)
	{
		int position = i % problemSize;
		
		double totalHistogram  = 0.0;
		double* pArray = new double[numResampled];
		for(int j = 0; j < numResampled; j++)
		{
			if(candidateUsed[candidate[j]] == true) continue;
			else{
				double leftP = 1.0, rightP = 1.0;
				int left = (position-1 + problemSize) % problemSize;
				int right = (position+1) % problemSize;
				
				int leftIndex = templateX.getVal(left);
				int rightIndex = templateX.getVal(right);
				
				int centreIndex = candidate[j];
				
				if(leftIndex != -1)
					leftP *= (EdgeHistogramMatrix[ leftIndex ][centreIndex] + EdgeHistogramMatrix[ centreIndex ][ leftIndex ] + epsilon);
					
				if(rightIndex != -1)
					rightP *= (EdgeHistogramMatrix[ rightIndex ][centreIndex] + EdgeHistogramMatrix[ centreIndex ][ rightIndex ] + epsilon);
					
				totalHistogram += leftP * rightP;
				pArray[j] = leftP * rightP;
			}
			
		}
		
		double pointer = totalHistogram * random.uniform();
		double partialSum = 0.0;
		int index = -1;
		
		for(int j = 0; j < numResampled; j++)
		{
			if(candidateUsed[candidate[j]] == true) continue;
			else if(pArray[j] == 0 && totalHistogram != 0) continue;
			
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


	
	return templateX;
}


void Model_EHBSA::dump()
{
    std::cout << "Model_EHBSA:\n";
    for(int i = 0; i < this->problemSize; i++)
    {
		for(int j = 0; j < this->problemSize; j++)
		{
			std::cout<< std::setw(2)<< EdgeHistogramMatrix[i][j];
		}

		std::cout << std::endl;
    }
}


void Model_EHBSA::update()
{
}

/*------------------EHBSAwoTemplate PART-----------------*/
Chromosome Model_EHBSAwoTemplate::resample()
{
	Chromosome templateX;
	templateX.clone(population[0]);
		
	//generate candidate (if -1: used, ignored by following code)	
	bool *candidateUsed = new bool [problemSize];
	
	for(int g = 0; g < problemSize; g++)
	{
		candidateUsed[g] = false;
		templateX.setVal(g, -1);
	}

	
	/*- roulette wheel for each empty position-*/
	for(int i = 0; i < problemSize; i++)
	{
		int position = i;
		
		if(i == 0) //if totally empty, randomly choose an index
		{
			int index = random.uniformInt(0, problemSize-1);
			templateX.setVal(position, index);
			candidateUsed[index] = true;
			continue;
		}
		
		double totalHistogram  = 0.0;
		double* pArray = new double[problemSize];
		for(int j = 0; j < problemSize; j++)
		{
			if(candidateUsed[j] == true) continue;
			else{
				double leftP = 1.0, rightP = 1.0;
				int left = (position-1 + problemSize) % problemSize;
				int right = (position+1) % problemSize;
				
				int leftIndex = templateX.getVal(left);
				int rightIndex = templateX.getVal(right);
				
				if(leftIndex != -1)
					leftP *= EdgeHistogramMatrix[ leftIndex ][j] + EdgeHistogramMatrix[ j ][ leftIndex ] + epsilon;
					
				if(rightIndex != -1)
					rightP *= EdgeHistogramMatrix[ rightIndex ][j] + EdgeHistogramMatrix[ j ][ rightIndex ]  + epsilon;
					
				totalHistogram += leftP * rightP;
				pArray[j] = leftP * rightP;
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


/*------------------Model_EHBSAwCircleResample PART-----------------*/

Chromosome Model_EHBSAwCircleResample::resample()
{
	int templateIndex = 0;
	Chromosome templateX;
	templateX.clone(population[templateIndex]);
	

	int *cutPoints = new int[NUM_CUTPOINTS];
	random.uniformArray(cutPoints, NUM_CUTPOINTS, 0, problemSize-1);
	Sort::bubble(cutPoints, NUM_CUTPOINTS);

	int resampledSegIndex = random.uniformInt(0, NUM_CUTPOINTS-1);

	int startIndex = cutPoints[resampledSegIndex];
	int endIndex = (resampledSegIndex == NUM_CUTPOINTS-1) ? 
			problemSize-1 : cutPoints[resampledSegIndex+1];
	int numResampled = endIndex - startIndex + 1;
	
	//generate candidate (if -1: used, ignored by following code)	
	int *candidate = new int[numResampled];
	bool *candidateUsed = new bool [problemSize];
	
	for(int g = 0; g < problemSize; g++)
		candidateUsed[g] = true;
	for(int g = 0; g < numResampled; g++)
	{
		candidate[g] = templateX.getVal(startIndex+g);
		templateX.setVal(startIndex+g, -1);
		candidateUsed[candidate[g]] = false;
	}
	
	
	/*- roulette wheel for each empty position-*/
	for(int i = startIndex; i <= endIndex; i++)
	{
		int position = i;
		double totalHistogram  = 0.0;
		double* pArray = new double[numResampled];
		for(int j = 0; j < numResampled; j++)
		{
			if(candidateUsed[candidate[j]] == true) continue;
			else{
				double leftP = 1.0, rightP = 1.0;
				int left = (position-1 + problemSize) % problemSize;
				int right = (position+1) % problemSize;
				
				int leftIndex = templateX.getVal(left);
				int rightIndex = templateX.getVal(right);
				
				if(leftIndex != -1)
					leftP *= EdgeHistogramMatrix[ leftIndex ][candidate[j]] + EdgeHistogramMatrix[ candidate[j] ][ leftIndex ] + epsilon;
					
				if(rightIndex != -1)
					rightP *= EdgeHistogramMatrix[ rightIndex ][candidate[j]] + EdgeHistogramMatrix[ candidate[j] ][ rightIndex ] + epsilon;
					
				totalHistogram += leftP * rightP;
				pArray[j] = leftP * rightP;
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

	return templateX;
}
