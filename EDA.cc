/*  ===============================================================================
 *  Program:        EDA (Estimation Distribution Algorithms)
 *  Module:         EDA.cpp : Defines the entry point for the console application.
 *  Date:           23-Dic-2011
 *  Author:         Josu Ceberio
 *  Description:    
 *      Calculates the optimum solution of a problem following an EDA approach. 
 *      All the problem specification is defined in Problem.cpp. The rest of the
 *        files should remain unaltered.
 *      Output is displayed only in the screen. The function WriteSolution can 
 *        be modified in order to generate the desired aoutput format.
 *
 *  History: 
 *
 *		23-Dic-2011 - Mallows probabilistic model for permutations EDA added by Josu Ceberio.
 *					  
 *  =============================================================================== 
 */ 

#include "EDA.h"
#include "Solution.h"
#include "Problem.h"
#include "PermutationTools.h"
#include "Tools.h"
#include "PBP.h"
#include "Seed.h"
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
using std::istream;
using std::ostream;
using namespace std;
#define itoa(a,b,c) sprintf(b, "%d",a)
#include "string.h"

// time calculation variables.
timeval start,stop;

//It is the instance of the problem to optimize.
PBP * m_problem;

//The type of the problem to solve.
char PROBLEM_TYPE[4];

// Number of generations carried out at the moment.
int GEN_NUM;

// The size of the population.
int POP_SIZE;

// The number of selected individuals.
int SEL_SIZE;

// The offspring size.
int OFFSPRING_SIZE;

// Number of restarts carried out during the evolutionary process.
int RESTARTS=0;

// Whether elitism is used or not.
int ELITISM=1;

// The individual size.
int IND_SIZE;

// Type of the probabibilistic mode to learn
// 0 ->Mallows
// 1 ->Generalized Mallows
// 2 ->PlackettLuce
int MODEL_TYPE=0;

// Number of evaluations performed.
unsigned long long int EVALUATIONS = 0;

// Convergence evaluation of the best fitness.
unsigned long long int CONVERGENCE_EVALUATIONS = 0;

// Maximum number of evaluations allowed performed.
unsigned long long int MAX_EVALUATIONS = 0;

// Maximum number of generations (used by the stopping criterium).
// When MAX_GENERATIONS=0 -> no limit
int MAX_GENERATIONS=30;

// Maximum number of generations without improvement.
int MAX_GENERATIONS_NO_IMPROVEMENT=30;

//Maximum number of restarts without improvement.
int MAX_RESTARTS_NO_IMPROVEMENT=30;

// Number of iterations without improvement.
int NO_IMPROVEMENT_GENERATIONS=0;

// Best result convergente iteration.
int CONVERGENCE_ITERATION=0;

// Name of the file where the result will be stored.
char OUTPUTFILE[50];

// Name of the file where the log results will be stored (optional).
char LOGFILE[50];

// Name of the file where the instances is stored.
char DATA_FILE_NAME[50];

// The seed asigned to the process
int SEED;

// Input file stream of the instance.
ifstream fdatafile;

// Output file stream of the result.
ofstream foutput;

// Output file stream of the log.
ofstream flog;

// Best current solution of the EDA
CIndividual * BEST;

/*
 * Get next command line option and parameter
 */
int GetOption (int argc, char** argv, char* pszValidOpts, char** ppszParam)
{
	
    static int iArg = 1;
    char chOpt;
    char* psz = NULL;
    char* pszParam = NULL;
	
    if (iArg < argc)
    {
        psz = &(argv[iArg][0]);
		
        if (*psz == '-' || *psz == '/')
        {
            // we have an option specifier
            chOpt = argv[iArg][1];
			
            if (isalnum(chOpt) || ispunct(chOpt))
            {	
                // we have an option character
                psz = strchr(pszValidOpts, chOpt);
				
                if (psz != NULL)
                {
                    // option is valid, we want to return chOpt
                    if (psz[1] == ':')
                    {
                        // option can have a parameter
                        psz = &(argv[iArg][2]);
                        if (*psz == '\0')
                        {
                            // must look at next argv for param
                            if (iArg+1 < argc)
                            {
                                psz = &(argv[iArg+1][0]);
                                if (*psz == '-' || *psz == '/')
                                {
                                    // next argv is a new option, so param
                                    // not given for current option
                                }
                                else
                                {
                                    // next argv is the param
                                    iArg++;
                                    pszParam = psz;
                                }
                            }
                            else
                            {
                                // reached end of args looking for param
                            }
							
                        }
                        else
                        {
                            // param is attached to option
                            pszParam = psz;
                        }
                    }
                    else
                    {
                        // option is alone, has no parameter
                    }
                }
                else
                {
                    // option specified is not in list of valid options
                    chOpt = -1;
                    pszParam = &(argv[iArg][0]);
                }
            }
            else
            {
                // though option specifier was given, option character
                // is not alpha or was was not specified
                chOpt = -1;
                pszParam = &(argv[iArg][0]);
            }
        }
        else
        {
            // standalone arg given with no option specifier
            chOpt = 1;
            pszParam = &(argv[iArg][0]);
        }
    }
    else
    {
        // end of argument list
        chOpt = 0;
    }
	
    iArg++;
	
    *ppszParam = pszParam;
    return (chOpt);
}

/*
 * Help command output.
 */
void usage(char *progname)
{
          cerr
	      << "Usage: " << progname << " -N pop_size -S sel_size  -a caching -i simulation_type" << endl
	      << "        -l learning_type -c score -f offspring_size -e elitism -g save_structures -o output_file -w data_filename" << endl
	      << endl
	      << "-N pop_size: number of individuals in the population. (def:2000)" << endl
	      << "-S sel_size: number of selected individuals. (def:1000)" << endl
	      << "-c score: type of score (only for EBNA): " << endl
	      << "       BIC (0) or K2 (1). (def:0)" << endl
	      << "-l learning: how the Bayesian network is learned:" << endl
	      << "       UMDA (0), EBNA with B algorithm (1), " << endl
	      << "       EBNA with local search (2), PBIL with alpha = 0.5 (3)," << endl
	      << " 	 BSC (4), TREE(5), MIMIC(6), EBNA K2 + penalization (7)," << endl
	      << "   EBNA PC (8) or EBNA with global search (9)," << endl
	      << "   Permutations Marginal Modelling (11), Mallows (12)." << endl
	      << "-b bias correction: (only for K-order marginals based EDA): [0-1]  (def:1)"<<endl
	      << "-f offspring_size: how many individuals are created in" << endl
	      << "       each generation. (def:1999)" << endl
	      << "-e elitism: how the next population is created: " << endl
	      << "       elitism (1) or no-elitism (0). (def:1)" << endl
	      << "-a caching: whether the values of the evaluated individuals.  (def:1)" << endl	      
          << "-i simulation: type of simulation: PLS (0), " << endl
		  << "       PLS forcing all values LTM (1) and ATM (2)," << endl
		  << "       correction after generation (3) and penalization (4).(def:0)" << endl	
	      << "-g save_structure: whether the structures of each generation " << endl	      
	      << "       have to be saved or not.  (def:0)" << endl
	      << "-P part_time: work in parts, execute only until generation 'part_time'  " << endl	      
	      << "       and store partial population in file 'pop.out'. (def:0,no output)" << endl
          << "-o output_file: name of output_file to store the results" << endl
	      << "       (optional parameter). (def:-)" << endl
		  << "-d log_file: name of log_file to store the log results" << endl
		  << "-t problem type: 'tsp','tsf','qap','fsp' or 'lop'"<< endl;
}

/*
 * Obtaint the execution parameters from the command line.
 */
bool GetParameters(int argc,char * argv[])
{
	char c;
    if(argc==1)
    {
    	usage(argv[0]);
        return false;
    }
	char** optarg;
	optarg = new char*[argc];
    while ((c = GetOption (argc, argv, "F:hN:l:e:o:p:u:x:d:t:w:",optarg)) != '\0')
    {
    	switch (c)
    	{
		
            case 'h' :  usage(argv[0]);
                        return false;
                        break;

            case 'l' :  MODEL_TYPE = atoi(*optarg);
                        break;

            case 'x' :  SEED = atoi(*optarg);
						break;
				
           	case 't':  strcpy(PROBLEM_TYPE, *optarg);
                        break;
                
            case 'o' :  strcpy(OUTPUTFILE, *optarg);
                          //OUTPUTFILE=optarg;
                          // open the output file
                        if (OUTPUTFILE)
                        {
                        	foutput.open(OUTPUTFILE);
                            if (!foutput) {
                                cerr << "Could not open file " << OUTPUTFILE << ". Ignoring this file." << endl;
                                OUTPUTFILE[0]='\0';
                            }
                        }
                        break;
				
			case 'd' :  strcpy(LOGFILE, *optarg);
                        break;
				
			case 'w':  strcpy(DATA_FILE_NAME, *optarg);
                        if (DATA_FILE_NAME){
                        fdatafile.open(DATA_FILE_NAME);
                            if (!fdatafile){
                                cerr << "Could not open file " << DATA_FILE_NAME << ". Ignoring this file." << endl;
                                DATA_FILE_NAME[0]='\0';
                            }
                            fdatafile.close();
                        }
                        break;

		}
     }

     delete [] optarg;

	return true;
}

/*
 * Main function.
 */
int xmain(int argc, char* argv[])
{
    //Start timer.
	gettimeofday(&start,NULL);
    
	//Initialize seed.
	seed();
//    srand(1000);
	//Parse the parameters from the command line.
    cout<<"Parsing parameters..."<<endl;
	if(!GetParameters(argc,argv)) return -1;
	
	//Read the problem instance to optimize.
    cout<<"Reading problem instance..."<<endl;
	if (GetProblemInfo(PROBLEM_TYPE,DATA_FILE_NAME) < 0) return -2;
	
	//Initialize the population.
	cout<<"Initializing population...."<<endl;
	CSolution solution;

	//Print execution parameters previous optimization.
	cout<<"Starting evolutionary algorithm..."<<endl;
	cout<<"Instance to optimize: "<<DATA_FILE_NAME<<endl;
    cout<<"Problem type: "<<PROBLEM_TYPE<<endl;
    cout<<"Model type: "<<MODEL_TYPE<<endl;
	cout<<"Population size: "<<POP_SIZE<<endl;
	cout<<"Offspring size: "<<OFFSPRING_SIZE<<endl;
	cout<<"Selection size: "<<SEL_SIZE<<endl;
	cout<<"Max generations: "<<MAX_GENERATIONS<<endl;
	cout<<"Max generations non improvement: "<<MAX_GENERATIONS_NO_IMPROVEMENT<<endl;
	cout<<"Max restarts non improvement: "<<MAX_RESTARTS_NO_IMPROVEMENT<<endl;
	cout<<"Max evaluations: "<<MAX_EVALUATIONS<<endl;
    cout<<"-----------------------------------------------------"<<endl;

    //Evolutionary process starts.
    for (;EVALUATIONS<MAX_EVALUATIONS ;)
    {
        //Print log if required
        if ((GEN_NUM % 10) ==0)
        {
            float averageSelected= solution.m_population->AverageFitnessPopulation(SEL_SIZE);
            cout<<GEN_NUM<<", "<<EVALUATIONS<<", "<<solution.GetBestPopulationIndividual()->Value()<<", "<<averageSelected<<", "<<BEST->Value()<<", "<<CONVERGENCE_ITERATION<<", "<<CONVERGENCE_EVALUATIONS<<endl;
        }
   
        //Improve solution.
		solution.Improve();
        
        //Increment generation number.
        GEN_NUM++;
	}
    cout<<"-----------------------------------------------------"<<endl;

	//Write solution in the output file.
	solution.WriteSolution();
	int best=BEST->Value();

	RemoveProblemInfo();
	gettimeofday(&stop,NULL);
	return -best;
}


