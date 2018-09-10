/*
 *  TSP.h
 *  MixtureMallowsEDA
 *
 *  Created by Josu Ceberio Uribe on 6/26/12.
 *  Copyright 2012 University of the Basque Country. All rights reserved.
 *
 */

#ifndef _TSP_H__
#define _TSP_H__

#include "PBP.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdio.h>


using std::istream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;

class TSP : public PBP
{
	
public:
	
    /*
     * Matrix of distances between the cities.
     */
	int ** m_distance_matrix;
	
	/*
	 * The number of jobs of the problem.
	 */
	int m_size;
	
	/*
     * The constructor. It initializes a flowshop scheduling problem from a file.
     */
	TSP();
	
    /*
     * The destructor.
     */
    virtual ~TSP();
	
	/*
	 * Read TSP file.
	 */
	int Read(string filename);
	
	/*
	 * This function evaluates the individuals for the TSP problem.
	 */
	int Evaluate(int * genes);
	
private: 
	
};
#endif
