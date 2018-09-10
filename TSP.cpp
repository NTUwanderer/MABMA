/*
 *  TSP.cpp
 *  MixtureMallowsEDA
 *
 *  Created by Josu Ceberio Uribe on 6/26/12.
 *  Copyright 2012 University of the Basque Country. All rights reserved.
 *
 */

#include "TSP.h"
#include "Tools.h"
#include "EDA.h"
#include <fstream>
#include <iostream>
using std::ifstream;
using std::ofstream;
/*
 *Class constructor.
 */
TSP::TSP()
{
	
}

/*
 * Class destructor.
 */
TSP::~TSP()
{
	for (int i=0;i<m_size;i++)
		delete [] m_distance_matrix[i];
	delete [] m_distance_matrix;
	
}


double CalculateGEODistance(double latitudeX, double latitudeY, double longitudeX, double longitudeY)
{
	double PI = 3.141592;
	double RRR = 6378.388;
	
	double deg = (double)((int)latitudeX);
	double min = latitudeX - deg;
	double latitudeI = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
	
	deg = (double)((int)latitudeY);
	min = latitudeY - deg;
	double longitudeI = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
	
	deg = (double)((int)longitudeX);
	min = longitudeX - deg;
	double latitudeJ = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
	
	deg = (double)((int)longitudeY);
	min = longitudeY - deg;
	double longitudeJ = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
	
	double q1 = cos(longitudeI - longitudeJ);
	double q2 = cos(latitudeI - latitudeJ);
	double q3 = cos(latitudeI + latitudeJ);
	
	return (int)(RRR * acos(0.5 * ((1.0+q1) * q2 - (1.0 - q1) * q3)) + 1.0);
}

int TSP::Read(string filename)
{
	bool readm_distance_matrix=false;
	bool coordinatesData=false;
	string DISTANCE_TYPE;
	double **coordinates;
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
		//cout<<sline<<endl;
		if (sline=="EOF")
		{
			break;
		}
		if (readm_distance_matrix && coordinatesData==false)
		{
			//cout << "reading distance m_distance_matrix "<<line<< endl;
			if (data=="")
				data = line;
			else
				data = data+' '+line;
		}
		else if (readm_distance_matrix && coordinatesData==true)
		{
			//FILL DISTANCE m_distance_matrix
			char * coordPieces;
			coordPieces = strtok (line," ");
			if (strcmp(coordPieces, " ") == 0)
			{
				coordPieces = strtok (NULL," ");
			}
			int coordNum = atoi(coordPieces);
			coordPieces = strtok (NULL, " ");
			double latitud = atof(coordPieces);
			coordPieces = strtok (NULL, " ");
			double longitud = atof(coordPieces);
			double *coordinate= new double[2];
			coordinate[0]=latitud;
			coordinate[1]=longitud;
			
			coordinates[coordNum-1]= coordinate;
			//cout<<"coordNum "<<coordNum-1<<" latit: "<<latitud<<" long: "<<longitud<<endl;
		}
		
		//COMPROBAR SI SE DEFINE EL TAMAÑO DEL PROBLEMA
		if(strContains(sline,"DIMENSION"))
		{
			char * pch;
			pch = strtok (line," ");
			pch = strtok (NULL, " ");
			if (strcmp(pch,":")==0)
			{
				pch = strtok (NULL, " ");
			}
			m_size = atoi(pch);
			//cout<<"problemsize "<<problemsize<<endl;
		}
		else if (strContains(sline,"EDGE_WEIGHT_TYPE"))
		{
			char * pch;
			pch = strtok (line," ");
			pch = strtok (NULL, " ");
			if (strcmp(pch,":")==0)
			{
				pch = strtok (NULL, " ");
			}
			stringstream s;
			string type;
			s << pch;
			s >> type;
			DISTANCE_TYPE = type;
		}
		//COMPROBAR SI SE DEFINEN LOS DATOS DEL PROBLEMA EN DISTANCIAS
		else if (sline =="EDGE_WEIGHT_SECTION")
		{
			readm_distance_matrix=true;
			coordinatesData=false;
		}
		//COMPROBAR SI SE DEFINEN LOS DATOS DEL PROBLEMA EN COORDENADAS
		else if (sline=="NODE_COORD_SECTION"){
			readm_distance_matrix=true;
			coordinatesData=true;
			coordinates= new double*[m_size];
		}
		
	}
	indata.close();
	
	//BUILD DISTANCE m_distance_matrix
	//cout << "--> BUILDING DISTANCE m_distance_matrix" << endl;
	m_distance_matrix = new int*[m_size];
	for (int i=0;i<m_size;i++)
	{
		m_distance_matrix[i]= new int[m_size];
	}
	
	//FILL DISTANCE m_distance_matrix
	if (coordinatesData==true)
	{
		//CALCULATE EUCLIDEAN DISTANCES
		//cout << "--> CALCULATING EUCLIDEAN DISTANCE m_distance_matrix" << endl;
		for (int i=0;i<m_size;i++)
		{
			//get coordinate A
			double *coordA=coordinates[i];
			double coordAx = coordA[0];
			double coordAy = coordA[1];
			for (int j=i;j<m_size;j++)
			{
				//get coordinate B.
				double *coordB=coordinates[j];
				double coordBx=coordB[0];
				double coordBy=coordB[1];
				double euclidean;
				if (DISTANCE_TYPE=="GEO")
				{
					//calculate geographic distance between A and B.
					euclidean=CalculateGEODistance(coordAx, coordAy, coordBx, coordBy);
				}
				else
				{
					//calculate euclidean distance between A and B.
					double absolute= fabs(pow((coordAx-coordBx),2) + pow((coordAy-coordBy),2));
					euclidean= sqrt(absolute);
				}
				//convert double to string
				std::ostringstream stream;
				stream << euclidean;
				std::string euclideanString = stream.str();
				//save distance in distances m_distance_matrix. Save negative distance in order to minimize fitness instead of
				//maximize.
				m_distance_matrix[i][j]=  (int)euclidean;
				m_distance_matrix[j][i]= (int)euclidean;//<-symmetric m_distance_matrix
			}
		}
	}
	else
	{
		//FILL DISTANCE m_distance_matrix
		//cout << "--> FILLING DISTANCE m_distance_matrix: "<<data << endl;
		istringstream iss(data);
		int i=0;
		int j=0;
		do
		{
			string sub;
		    iss >> sub;
		    if (sub!=""){
				//save distance in distances m_distance_matrix. Save negative distance in order to minimize fitness instead of
				//maximize.
		    	m_distance_matrix[i][j]= atoi(sub.c_str());
		    	m_distance_matrix[j][i]= atoi(sub.c_str());//<-symmetric m_distance_matrix
		    	if (sub=="0")
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
	}

	
		return (m_size);
}


int TSP::Evaluate(int * genes)
{
	EVALUATIONS++;

	//cout<<"here"<<endl;
	int distanceSum=0;
	int distAB=0;
	int IDCityA, IDCityB;
	for (int i=0;i<m_size;i++)
	{
		IDCityA = genes[i];
		IDCityB = genes[0];
		if (i+1<m_size)
		{
			IDCityB = genes[i+1];
		}
		
		distAB = m_distance_matrix[IDCityA][IDCityB];
		distanceSum=distanceSum+distAB;
		
		//security condition
		if (IDCityA==m_size || IDCityB==m_size){
			distanceSum=0;
			break;
		}
	}
	//cout << "evaluation:"<<distanceSum<<endl;
	return -distanceSum;
}
