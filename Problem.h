// Problem.h: definition of the functions specific to the problem.
//
// Author:	Ramon Etxeberria
// Date:	1999-09-25
//
// Note:
//		You must implement the following three functions if
//		you want to solve a specific problem.

#ifndef _PROBLEM_
#define _PROBLEM_

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include "EDA.h"
#include "PFSP.h"
#include "LOP.h"
#include "QAP.h"
#include "TSP.h"
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;

// The evaluation function of the individuals.
int Metric(int * genes, int size);

// It gets the information related to the problem being solved.
// At least, this function must initaliaze the size of the
// individuals and the number of states of their genes.
int GetProblemInfo(string problemType, string filename);

// This functions is used to free all the memory required
// for the problem's information.
void RemoveProblemInfo();

#endif



