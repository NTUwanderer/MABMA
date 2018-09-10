/*
 *  QAP.h
 *  MixtureMallowsEDA
 *
 *  Created by Josu Ceberio Uribe on 6/26/12.
 *  Copyright 2012 University of the Basque Country. All rights reserved.
 *
 */


#ifndef _QAP_H__
#define _QAP_H__

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

class QAP : public PBP
{
	
public:
	
    /*
     * The matrix of distances between the cities.
     */
	int ** m_distance_matrix;
	
    /*
     * The flow matrix.
     */
	int ** m_flow_matrix;
	
	/*
	 * The number of jobs of the problem.
	 */
	int m_size;
	
	/*
     * The constructor. It initializes a QAP from a file.
     */
	QAP();
	
    /*
     * The destructor.
     */
    virtual ~QAP();
	
	/*
	 * Read QAP file.
	 */
	int Read(string filename);
	
	/*
	 * This function evaluates the individuals for the QAP problem.
	 */
	int Evaluate(int * genes);
	
private: 
	
};
#endif
