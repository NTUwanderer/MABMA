/*
 *  PFSP.cpp
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "PFSP.h"
#include "Tools.h"
#include "EDA.h"
#include <fstream>
#include <iostream>
using std::ifstream;
using std::ofstream;

/*
 *Class constructor.
 */
PFSP::PFSP()
{
	
}

/*
 * Class destructor.
 */
PFSP::~PFSP()
{
	for	(int i=0;i<m_machines;i++)
		delete [] JOBPROCESSINGMATRIX[i];
	delete [] JOBPROCESSINGMATRIX;
}

int PFSP::Evaluate(int * genes)
{
    EVALUATIONS++;
	return -EvaluateFSPTotalFlowtime(genes);
}

						 
int PFSP::Read(string filename)
{
	bool readMatrix=false;
	bool readDimension=false;
	//double **coordinates;
	char line[2048]; // variable for input value
	string data="";
	ifstream indata;
	indata.open(filename.c_str(),ios::in);
	//int num=0;
	while (!indata.eof())
	{
		//LEER LA LINEA DEL FICHERO
		indata.getline(line, 2048);
		stringstream ss;
		string sline;
		ss << line;
		ss >> sline;
		if (strContains(line,"number of jobs")==true && readMatrix==true)
		{
			break;
		}
		else if (strContains(line,"number of jobs")==true)
		{
			readDimension=true;
		}
		else if (readDimension==true)
		{
			//char * pch;
			m_jobs = atoi(strtok (line," "));
			//cout<<"JOB NUM: "<<m_jobs<<endl;
			m_machines = atoi(strtok (NULL, " "));
			//cout<<"MACHINE NUM: "<<m_machines<<endl;
			readDimension=false;
			
		}
		else if (readMatrix)
		{;
			if (data=="")
				data = line;
			else
				data = data+' '+line;
		}
		else if (strContains(line,"processing times :"))
		{
			readMatrix=true;
		}
	}
	indata.close();
	
	//BUILD JOB PROCESSING MATRIX
	//cout << "--> BUILDING JOB PROCESSING MATRIX" << endl;
	JOBPROCESSINGMATRIX = new int*[m_machines];
	for (int i=0;i<m_machines;i++)
	{
		JOBPROCESSINGMATRIX[i]= new int[m_jobs];
	}
	
	//FILL JOB PROCESSING MATRIX
	//cout << "--> FILLING JOB PROCESSING MATRIX: "<<data << endl;
	istringstream iss(data);
	int i=0;
	int j=0;
	do
	{
		string sub;
	    iss >> sub;
	    if (sub!="")
	    {
			//save distance in distances matrix. Save negative distance in order to minimize fitness instead of
			//maximize.
	    	JOBPROCESSINGMATRIX[i][j]= atoi(sub.c_str());
	    	if (j==m_jobs-1)
	    	{
	    		i++;
	    		j=0;
	    	}
	    	else
	    	{
	    		j++;
	    	}
	    }
	    else
	    {
	    	break;
	    }
	} while (iss);
	
	return (m_jobs);
}

double PFSP::EvaluateFSPMakespan(int * genes)
{
	double fitness=0;
	double timeTable[m_machines];
	
	for (int j=0;j<m_machines;j++)
	{
		timeTable[j]=0;
	}
	for (int z=0;z<m_jobs;z++)
	{
		int job=genes[z];
		//cout<<"Job "<<job<<endl;
		for (int machine=0;machine<m_machines;machine++)
		{
			double processingTime=JOBPROCESSINGMATRIX[machine][job];
			if (machine==0)
			{
				timeTable[machine]+=processingTime;
			}
			else
			{
				if (timeTable[machine-1]<timeTable[machine])
				{
					timeTable[machine]=timeTable[machine]+processingTime;
				}
				else
				{
					timeTable[machine]= timeTable[machine-1]+processingTime;
				}
			}
			//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
		}
	}
	
	fitness=timeTable[m_machines-1];
	return fitness;
}

int PFSP::EvaluateFSPTotalFlowtime(int * genes)
{
    
	int timeTable[m_machines];
	int j,z, job;
	int machine;
    int prev_machine=0;
	//int genes[20]={2,16,8,14,13,7,18,12,15,5,6,0,1,3,4,17,19,11,10,9};
	int first_gene=genes[0];
	timeTable[0]=JOBPROCESSINGMATRIX[0][first_gene];
	for (j=1;j<m_machines;j++)
	{
		timeTable[j]=timeTable[j-1]+JOBPROCESSINGMATRIX[j][first_gene];
	}
	
	int fitness=timeTable[m_machines-1];
	for (z=1;z<m_jobs;z++)
	{
		job=genes[z];
		
		//machine 0 is always incremental, so:
		timeTable[0]+=JOBPROCESSINGMATRIX[0][job];
		prev_machine=timeTable[0];
		for (machine=1;machine<m_machines;machine++)
		{			
			timeTable[machine]= MAX(prev_machine,timeTable[machine])+ JOBPROCESSINGMATRIX[machine][job];
			prev_machine=timeTable[machine];
		}
		
		fitness+=timeTable[m_machines-1];
	}
	
	return fitness;
}


