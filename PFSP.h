/*
 *  PFSP.h
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#ifndef _PFSP_H__
#define _PFSP_H__

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdio.h>
#include "PBP.h"
using std::istream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;

class PFSP : public PBP
{
	
public:
	
	/*
	 * The number of jobs of the problem.
	 */
	int m_jobs;
	
	/*
	 * The number of machines of the problem.
	 */
	int m_machines;
	
	/*
	 * The processing times matrix.
	 */
	int **JOBPROCESSINGMATRIX;


	// The constructor. It initializes a flowshop scheduling problem from a file.
	PFSP();
	
    // The destructor.
    virtual ~PFSP();
	
	/*
	 * Reads and PFSP file of Taillard Data set with the given filename.
	 */
	int Read(string filename);
	
	/*
	 * Evaluates the given solution with the default criterion.
	 */
	int Evaluate(int * genes);
	
	/*
	 * This function evaluates the individuals for the makespan of FSP problem.
	 */
	double EvaluateFSPMakespan(int * genes);

	/*
	 * This function evaluates the individuals for the total flowtime of FSP problem.
	 */
	int EvaluateFSPTotalFlowtime(int * genes);

private: 
	
	
	
};
#endif





