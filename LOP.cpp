/*
 *  LOP.cpp
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "LOP.h"
#include "Tools.h"
#include "EDA.h"
#include <fstream>
#include <iostream>
using std::ifstream;
using std::ofstream;

/*
 * Class constructor.
 */
LOP::LOP()
{
	
}

/*
 * Class destructor.
 */
LOP::~LOP()
{
	for (int i=0;i<m_problemsize;i++)
		delete [] m_matrix[i];
	delete [] m_matrix;
}

int LOP::Read(string filename)
{
	char line[2048]; // variable for input value
	string data="";
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
		if (sline=="")
		{
			break;
		}
		if (num==0)
		{
			//primera linea
			m_problemsize = atoi(line);
			//cout << "Problem size: "<< line<<endl;
		}
		else
		{
			if (data=="")
				data = line;
			else
				data = data+' '+line;
		}
		num++;
	}
	indata.close();

	//BUILD DISTANCE MATRIX
	m_matrix = new int*[m_problemsize];
	for (int i=0;i<m_problemsize;i++)
	{
		m_matrix[i]= new int[m_problemsize];
	}

	istringstream iss(data);
	int i=0;
	int j=0;
	do
	{
		string sub;
	    iss >> sub;
	    if (sub!=""){
			//save distance in distances matrix.
	    	m_matrix[i][j]= atoi(sub.c_str());
	    	if (j==(m_problemsize-1))
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
    
	return (m_problemsize);
}

int LOP::Evaluate(int * genes)
{
	EVALUATIONS++;
	int fitness=0;
    int i,j;
	for (i=0;i<m_problemsize-1;i++)
		for (j=i+1;j<m_problemsize;j++)
			fitness+= m_matrix[genes[i]][genes[j]];

	return fitness;
}


