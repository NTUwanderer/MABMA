// Solution.h: interface for the CSolution class.
//
// Description:
//      The class which solves the problem. It consists of
//      a population which evolves according to EDA.
//
// Author:  Josu Ceberio
// Date:    1999-09-25
//
// Notes: 
//      There will be only one object of this class in the
//      program. All the global variables that define the
//      problem and the optimization parameters must be 
//      initialized before that object is created.

#ifndef _SOLUTION_
#define _SOLUTION_

#include <time.h>

#include "Individual.h"
#include "Population.h"
#include "RankingModel.h"
#include "GeneralizedMallowsModel.h"
#include "MallowsModel.h"
#include "PlackettLuce.h"
#include "Timer.h"

class CSolution 
{
public:

	//This timer is used to calculate the execution time of the whole
    //program.store the time. This is used to 
    //display at the end of the program the time needed by the
    //algorithm to reach the final solution
    CTimer timer;
	
    // It performs a generation of the EDA.
    void Improve();

    // The constructor. It creates the initial population.
    CSolution();

    // The destructor. It frees all the used memory.
    virtual ~CSolution();

	//Creates all the structures used to save the results and populations.
	void CreateStructures();
	
	//Destroyes all the structures used to save the results and populations, freeing memory.
	void DestroyStructures();
	
    // It displays the best solution found so far, together
    // with the information of the optimization process. 
    // Usually used when the algorithm has converged.
    friend ostream & operator<<(ostream & os,CSolution & soluzioa);

    // It returns the best individual for the actual generation, 
    // as well as its fitness value.
    CIndividual  * & GetBestPopulationIndividual();
	
    // It inserts the given individual in the population.
    // The individual is inserted according to its value.
    void AddToPopulation(CIndividual * individual);

    // This function prepares the final solution and writes it
    // in the screen. It can also be modified in order to
    // write the answer also in a file.
    void WriteSolution();

    // Sets the timer to 0 and start it
    void StartTimer();

    // Stops the timer
    void StopTimer();

    // The population of individuals. It will be an ordered
    // list. This way, the best individual will be allways
    // at the head of the list.
    CPopulation * m_population;
    
private:


    
    //This variable defines the model that it is employed in the EDA whatever are its characteristics.
    CRankingModel * m_ranking_model;
	
};

#endif 
