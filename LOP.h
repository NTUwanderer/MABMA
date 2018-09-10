/*
 *  LOP.h
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#ifndef _LOP_H__
#define _LOP_H__

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

class LOP : public PBP
{
	
public:
	
    /*
     * Entries matrix of the LOP.
     */
	int ** m_matrix;
    
	/*
	 * The number of jobs of the problem.
	 */
	int m_problemsize;
		
    /*
     * The constructor. It initializes a flowshop scheduling problem from a file.
     */
	LOP();
	
    /*
     * The destructor.
     */
    virtual ~LOP();
	
	/*
	 * Read LOP file.
	 */
	int Read(string filename);
	
	/*
	 * This function evaluates the individuals for the LOP problem.
	 */
	int Evaluate(int * genes);
	
private: 
	
};
#endif





