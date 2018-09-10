#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <exception>

 #include <unistd.h>

#include "global.h"
#include "evaluator.h"



QAP_Evaluator::QAP_Evaluator()
{
	IEvaluator::init();
	
	problemSize = 0;
	best = 0.0;
	isReverse = false;


}


QAP_Evaluator::QAP_Evaluator(const char* fileName)
{
	IEvaluator::init();
	
	best = 0.0;
	problemSize = 0;
	isReverse = false;
	
	readData(fileName);
	readBest(fileName);
}

QAP_Evaluator::~QAP_Evaluator()
{

}


bool QAP_Evaluator::readData(const char* fileName)
{
	std::ifstream file(fileName);
	char str[1000];
	
	int k = 1;
	int i = 0;
	
	file >> str;				//Read size of problem
	problemSize = atoi(str);	
	
	distanceMatrix.init(problemSize, problemSize);
	flowMatrix.init(problemSize, problemSize);
	
	
	i = 0;
	for (k=0; k<problemSize*problemSize; k++)		    
	{
	    file >> str;
	    flowMatrix[i][k%problemSize] = atoi(str);
	    if ( (k+1)%problemSize == 0 )
	    {
		i++;
	    }	
	}
	i = 0;
	for (k=0; k<problemSize*problemSize; k++)		    
	{
	    file >> str;
	    distanceMatrix[i][k%problemSize] = atoi(str);
	    if ( (k+1)%problemSize == 0 )
	    {
		i++;
	    }	
	}
	
	file.close();
	return true;
}


void QAP_Evaluator::readBest(const char* fileName)
{
	std::string fileLoc(fileName);
	
	int pos = fileLoc.find("/qap/") + 5; //at "/"
	fileLoc.insert(pos, "sol/");
	pos = fileLoc.find(".dat");
	fileLoc.replace(pos, 4, ".sln");
	
	std::ifstream file(fileLoc.c_str());
	char str[80];
	
	file >> str;				  //Read size of problem

	
	file >> str;
	
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


double QAP_Evaluator::calc(const Chromosome* c)
{
	this->increaseNFE();

	double val = 0.0;
	for(int i = 0; i < problemSize ; i++){
		for(int j = 0; j < problemSize; j++){
			val += flowMatrix[i][j]*distanceMatrix[ c->getVal(i) ][ c->getVal(j) ];
		}
	}

	return val;	
}


double QAP_Evaluator::getBest()
{
	return this->best;
}


int QAP_Evaluator::getProblemSize()
{
	return this->problemSize;
}

void QAP_Evaluator::dump()
{

	for(int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			printf("%d ", (int)flowMatrix[i][j]);
		}
		
		printf("\n");
	}
	
	printf("\n");printf("\n");
	
	for(int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
			printf("%d ", (int)distanceMatrix[i][j]);

		printf("\n");
	}

}


