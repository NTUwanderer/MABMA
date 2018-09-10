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

FSSP_Evaluator::FSSP_Evaluator()
{
	IEvaluator::init();
	
	problemSize = 0;
	best = 0.0;
	
	job_num = 0;
	machine_num = 0;
	isValid = false;
	job_times = NULL;
	due_times = NULL;
	
	initial_seed = 0;
	current_job = 0;
	
	flowTable = NULL;
}


FSSP_Evaluator::FSSP_Evaluator(const char* fileName)
{
	IEvaluator::init();
	
	best = 0.0;
	readData(fileName);
	
	problemSize = job_num;
	
	flowTable = new int[problemSize];
}

FSSP_Evaluator::~FSSP_Evaluator()
{
	for(int i = 0; i < job_num; i++)
	{
		if(job_times[i] != NULL)
			delete[] job_times[i];
	}
	delete[] job_times;
	delete[] due_times;
	
	delete[] flowTable;
}


bool FSSP_Evaluator::readData(const char* fileName)
{
	try {
		isValid = false;
		std::ifstream sc(fileName);
		if(!sc){
	        printf("ERROR: Failed to Read problem file.");
	        return false;
        }
		
		int state = -1;
		int current_job = 0;
		
		std::string line = "";
		while (getline(sc, line)){
			switch (state) {
			case -1:
				best = atoi(line.c_str());
				state ++;
				break;				
			case 0: // iş sayısını al
				job_num = atoi(line.c_str());
				state ++;
				break;
			case 1: // makine sayısını al. ve değişkenleri başlat
				machine_num = atoi(line.c_str());
				state ++;
				job_times = new int*[job_num];
				for(int i = 0; i < job_num; i++)
					job_times[i] = new int[machine_num];
				
				due_times = new int[job_num];
				break;
			case 2: // başlangıç seed'ini oku (rastgele sayı üreteci için)
				initial_seed = atol(line.c_str());
				state ++;
				break;
			case 3: // işlenecek iş numarasını oku
				current_job = atoi(line.c_str());
				state ++;
				break;
			case 4: // işlenecek iş numarasını oku
				due_times[current_job]= atoi(line.c_str());
				state ++;
				break;
			case 5: // işlenecek iş numarasını oku
				{
					int loc0 = 0;
					int loc1 = (int) line.find(' ');
					// std::cout << "loc1: " << loc1 << std::endl;
					int counter = 0;
					while (loc1 != -1){
						std::string job_time = line.substr(loc0,loc1-loc0);
						// std::cout << "job_time: " << job_time << std::endl;
						job_times[current_job][counter] = atoi(job_time.c_str());
						counter++;
						loc0 = loc1+1;
						loc1 = (int) line.find(' ', loc0);
					}
				}
				state = 3;
				break;
			default:
				break;
			}
		}
		sc.close();
		isValid = true;
	} catch (const std::exception &e) {
		std::cout << "Standard exception: " << e.what() << std::endl;
		isValid = false;
		return false;
	}
	
	return true;
}


double FSSP_Evaluator::calc(const Chromosome* c)
{
	this->increaseNFE();
	for(int j = 0; j < machine_num; j++)
	{
		flowTable[j] = 0;
	}

	for(int i = 0; i < job_num; i++)
	{
		for(int j = 0; j < machine_num; j++)
		{
			int jobIndex = c->getVal(i);
			int thisJobProcessTime = job_times[jobIndex][j];
			if(j == 0)
			{
				flowTable[0] = flowTable[0] + thisJobProcessTime;
			}else{
				//last job end time at this machine: table[j]
				//this job end time at last machine: table[j-1]
				//choose the LARGER one + thisJobProcessTime
				flowTable[j] = ((flowTable[j] > flowTable[j-1]) ? 
						flowTable[j] : flowTable[j-1]) + thisJobProcessTime;				
			}
		}
	}
	
	return flowTable[machine_num-1];
}


double FSSP_Evaluator::getBest()
{
	return this->best;
}


int FSSP_Evaluator::getProblemSize()
{
	return this->problemSize;
}

