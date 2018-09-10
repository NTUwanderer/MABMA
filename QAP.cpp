/*
 *  QAP.cpp
 *  MixtureMallowsEDA
 *
 *  Created by Josu Ceberio Uribe on 6/26/12.
 *  Copyright 2012 University of the Basque Country. All rights reserved.
 *
 */

#include "QAP.h"
#include "Tools.h"
#include "EDA.h"
#include <fstream>
#include <iostream>
using std::ifstream;
using std::ofstream;
/*
 *Class constructor.
 */
QAP::QAP()
{
	
}

/*
 * Class destructor.
 */
QAP::~QAP()
{
	for (int i=0;i<m_size;i++)
	{
		delete [] m_distance_matrix[i];
		delete [] m_flow_matrix[i];
	}
	delete [] m_flow_matrix;
	delete [] m_distance_matrix;
	
}



int QAP::Read(string filename)
{

	char line[2048]; // variable for input value
	
	ifstream indata;
	indata.open(filename.c_str(),ios::in);
	int num=0;
	
	while (!indata.eof())
	{
		//LEER LA LINEA DEL FICHERO
		indata.getline(line, 2048);
		stringstream ss;
		string sline;
		ss << line;
		ss >> sline;
		
		if (num==0)
		{
			//OBTENER EL TAMAÑO DEL PROBLEMA
			m_size = atoi(sline.c_str());
			//cout << "--> BUILDING DISTANCE MATRIX" << endl;
			//cout << "--> BUILDING FLOW MATRIX" << endl;
			m_distance_matrix = new int*[m_size];
			m_flow_matrix = new int*[m_size];
			for (int i=0;i<m_size;i++)
			{
				m_distance_matrix[i]= new int[m_size];
				m_flow_matrix[i] = new int[m_size];
			}
		}
		else if (1<=num && num<=m_size)
		{
			//LOAD DISTANCE MATRIX
			char * pch;
			pch = strtok (line," ");
			int distance=atoi(pch);
			m_distance_matrix[num-1][0]=distance;
			for (int i=1;i < m_size; i++)
			{
				pch = strtok (NULL, " ,.");
				distance=atoi(pch);
				m_distance_matrix[num-1][i]=distance;
			}
			
		}
		else if (num>m_size && num<=(2*m_size))
		{
			
			//LOAD FLOW MATRIX
			char * pch;
			pch = strtok (line," ");
			int weight=atoi(pch);
			m_flow_matrix[num-m_size-1][0]=weight;
			for (int i=1;i < m_size; i++)
			{
				pch = strtok (NULL, " ,.");
				weight=atoi(pch);
				m_flow_matrix[num-m_size-1][i]=weight;
			}
		}
		else
		{
			break;
		}
		num++;
	}
	indata.close();
	
	return (m_size);
}


int QAP::Evaluate(int * genes)
{
	EVALUATIONS++;

	int fitness=0;
	int FactA, FactB;
	int distAB, flowAB, i ,j;
	for (i=0;i<m_size;i++)
	{
		for (j=0;j<m_size;j++)
		{
			FactA = genes[i];
			FactB = genes[j];
			
			distAB= m_distance_matrix[i][j];
			flowAB= m_flow_matrix[FactA][FactB];
			fitness= fitness+(distAB*flowAB);			
		}
	}
	
	return -fitness;
}
