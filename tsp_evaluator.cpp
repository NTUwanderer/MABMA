#include <stdio.h>
#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "global.h"
#include "evaluator.h"
#include "tinyxml2.h" 
#include "Array2D.h" 

using namespace tinyxml2;

TSP_Evaluator::TSP_Evaluator(){
	IEvaluator::init();
	
	problemSize = 0;
	best = 0.0;
}

TSP_Evaluator::TSP_Evaluator(const char* fileName)
{
	IEvaluator::init();
	
	readData(fileName);
	readBest(fileName);
}

TSP_Evaluator::~TSP_Evaluator(){
	
}

bool TSP_Evaluator::readData(const char* fileName){
	XMLDocument doc;
	doc.LoadFile(fileName); 
	
	XMLElement *scene=doc.RootElement();  
    XMLElement *graph=scene->FirstChildElement("graph");  
    XMLElement *vertex=graph->FirstChildElement("vertex");  
	
	int count = 0;
	while(vertex)
	{
		count++;
		vertex = (XMLElement*)vertex->NextSibling();
	}
	
	this->problemSize = count;
	
	vertex=graph->FirstChildElement("vertex");
	vDistanceMatrix.init(problemSize, problemSize);
	XMLElement *edge;
	
	int i = 0; 
	int j = 0;
	while(vertex)
	{
		j = 0;
		edge = vertex->FirstChildElement("edge");  
		while(edge)
		{
			if(j == i)
			{
				vDistanceMatrix[i][j] = 0.0;
				j++;
			}else{
				vDistanceMatrix[i][j] = round(edge->DoubleAttribute("cost"));
				edge = (XMLElement*)edge->NextSibling();  
				j++;
			}
		}
		i++;
		vertex = (XMLElement*)vertex->NextSibling();
	}

    return true;
}

void TSP_Evaluator::readBest(const char* fileName)
{
	// ../benchmarks/tsp/eil51.xml
	std::string fileLoc(fileName);
	
	int pos = fileLoc.find("/tsp/") + 5; //at "/"
	fileLoc.insert(pos, "opt/");
	pos = fileLoc.find(".xml");
	fileLoc.replace(pos, 4, ".opt.tour");
	
	std::ifstream file(fileLoc.c_str());
	char str[256];
  if (!file) {
    printf("file not found: %s\n", fileName);
    return;
  }


	file >> str;  //Read size of problem
	

	file.getline(str, 256);
	file.getline(str, 256);

	
	char* pch;
	pch = strtok(str, "( )");
	while(strncmp(pch, "Length", 6) != 0){
	    pch = strtok(NULL, "( )");
	}

	pch = strtok(NULL, "( )");


	
	file.getline(str, 256);
	while(strncmp(str, "TOUR_SECTION", 11) != 0){

		file.getline(str, 256);
	}

	Chromosome* temp = new Chromosome(problemSize);
	file >> str;
	for(int i = 0; i < problemSize ; i++)
	{
		
		temp->setVal(i, atoi(str)-1);
		file >> str;
	}
	
	this->best = this->calc(temp);
	
	delete temp;
	file.close();
}

double TSP_Evaluator::getBest(){
	return this->best;
}

int TSP_Evaluator::getProblemSize(){
	return this->problemSize;
}

double TSP_Evaluator::calc(const Chromosome *c)
{
	this->increaseNFE();
	double tmp = 0.0;
	
	for (int it=0; it<problemSize-1; ++it){
		int city1 = c->getVal(it);
		int city2 = c->getVal(it+1); 
		tmp += vDistanceMatrix[city1][city2];


	}
	int cityFirst = c->getVal(0);
	int cityLast = c->getVal(this->problemSize - 1);
	


	tmp += vDistanceMatrix[cityLast][cityFirst];

	return tmp;
}

//------------------------------------------------------------



TSPFixed_Evaluator::TSPFixed_Evaluator(char* filename)
: TSP_Evaluator(filename)
{
    this->best = this->best - 5.0;
}


double TSPFixed_Evaluator::calc(const Chromosome *c)
{
	double tmp = 0.0;
	tmp = TSP_Evaluator::calc(c);

	if(c->getVal(0) == 8)
	    tmp = tmp - 5;

	return tmp;
}
