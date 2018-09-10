#ifndef _EDA_
#define _EDA_
#include "Individual.h"
#include "Variables.h"
#include "PBP.h"
#include <fstream>
#include <iostream>
using std::ifstream;
using std::ofstream;

extern PBP * m_problem;

//The type of the problem to solve.
extern char PROBLEM_TYPE[4];

// Number of generations carried out at the moment.
extern int GEN_NUM;

// The size of the population.
extern int POP_SIZE;

// The number of selected individuals.
extern int SEL_SIZE;

// The offspring size.
extern int OFFSPRING_SIZE;

// Whether elitism is used or not.
extern int ELITISM;

// The seed asigned to the process
extern int SEED;

// The individual size.
extern int IND_SIZE;

// Type of the probabibilistic mode to learn
// 0 ->Mallows
// 1 ->Generalized Mallows
// 2 ->PlackettLuce
extern int MODEL_TYPE;

// Number of evaluations performed.
extern unsigned long long int EVALUATIONS;

// Convergence evaluation of the best fitness.
extern unsigned long long int CONVERGENCE_EVALUATIONS;

// Maximum number of evaluations allowed performed.
extern unsigned long long int MAX_EVALUATIONS;

// Maximum number of generations (used by the stopping criterium).
// When MAX_GENERATIONS=0 -> no limit
extern int MAX_GENERATIONS;

// Best result convergente iteration.
extern int CONVERGENCE_ITERATION;

// Maximum number of generations without improvement.
extern int MAX_GENERATIONS_NO_IMPROVEMENT;

//Maximum number of restarts without improvement.
extern int MAX_RESTARTS_NO_IMPROVEMENT;

// Name of the file where the result will be stored.
extern char OUTPUTFILE[50];

// Name of the file where the log results will be stored (optional).
extern char LOGFILE[50];

// Name of the file where the instances is stored.
extern char DATA_FILE_NAME[50];

// Output file stream of the result.
extern ofstream foutput;

// Output file stream of the log.
extern ofstream flog;

// Best current solution of the EDA.
extern CIndividual * BEST;

#endif
