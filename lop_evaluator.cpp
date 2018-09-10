#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <exception>

#include "global.h"
#include "evaluator.h"



LOP_Evaluator::LOP_Evaluator()
{
	IEvaluator::init();
	
	problemSize = 0;
	best = 0.0;
	isReverse = false;


}


LOP_Evaluator::LOP_Evaluator(const char* fileName)
{
	IEvaluator::init();
	
	best = 0.0;
	problemSize = 0;
	isReverse = false;
	
	readData(fileName);
	readBest(fileName);
}

LOP_Evaluator::~LOP_Evaluator()
{

}


bool LOP_Evaluator::readData(const char* fileName)
{
	std::ifstream file(fileName);
	char str[80];
	int  i,k;
	file.getline(str,10000);	  //Ignore first line (title).	
	file >> str;				  //Read size of problem
	problemSize = atoi(str);	
	
	costMatrix.init(problemSize, problemSize);
	
	i = 0;
	for (k=0; k<problemSize*problemSize; k++)		    
	{
		file >> str;
		costMatrix[i][k%problemSize] = atoi(str);
		if ( (k+1)%problemSize == 0 ) 
			i++;
	}

	file.close();
	return true;
}

void LOP_Evaluator::readBest(const char* fileName)
{
	std::string fileLoc(fileName);
	
	int pos = fileLoc.find("/lolib/") + 6; //at "b"
	fileLoc.insert(pos, "_opt");
	fileLoc.append(".opt");
	
	std::ifstream file(fileLoc.c_str());
	char str[80];
  if (!file) {
    printf("file not found: %s\n", fileName);

    return;
  }

	file >> str;				  //Read size of problem
	
	Chromosome* temp = new Chromosome(problemSize);
	
	for (int k = 0; k < problemSize; k++)		    
	{
		file >> str;
		int t = atoi(str) - 1;
		temp->setVal(k, t);
	}
	
	this->best = this->calc(temp);
	
	delete temp;
	file.close();
}


double LOP_Evaluator::calc(const Chromosome* c)
{
	this->increaseNFE();

	double val = 0;
	int i,j;
	for (i=0; i < problemSize-1; i++)
	{
		for (j=i+1; j < problemSize; j++)
		{
			if (isReverse)
				val += costMatrix[ c->getVal(j) ][ c->getVal(i) ];
			else
				val += costMatrix[ c->getVal(i) ][ c->getVal(j) ];
		}
	}
	return -1 * val;	
}


double LOP_Evaluator::getBest()
{
	return this->best;
}


int LOP_Evaluator::getProblemSize()
{
	return this->problemSize;
}



