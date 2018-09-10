// Solution.cpp: implementation of the CSolution class.
//

#include "Solution.h"
#include "EDA.h"
#include "PermutationTools.h"
#include "Problem.h"
#include "Tools.h"
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#define itoa(a,b,c) sprintf(b, "%d",a)


CSolution::CSolution()
{
    m_population= new CPopulation(POP_SIZE, OFFSPRING_SIZE,IND_SIZE);
    
    int * permutation= new int[IND_SIZE];
	//Initial population created randomly as usual
	for(int i=0;i<POP_SIZE;i++)
	{
		//Create random individual
		GenerateRandomPermutation(IND_SIZE,permutation);
        m_population->SetToPopulation(permutation, i);
    }
    delete [] permutation;
    
    m_population->SortPopulation(0);
    
    // Save the best individal
	BEST = GetBestPopulationIndividual()->Clone();
    
	CreateStructures();
    
    //starts the timer from 0.
    StartTimer();
}

CSolution::~CSolution()
{
	DestroyStructures();
}

void CSolution::CreateStructures()
{
    if (MODEL_TYPE==0)
    {
        m_ranking_model= new CMallowsModel();
    }
    else if (MODEL_TYPE==1)
    {
        m_ranking_model= new CGeneralizedMallowsModel();
    }
    else
    {
        m_ranking_model = new CPlackettLuceModel();
    }
}

void CSolution::DestroyStructures()
{
    cout<<"Destroying structures..."<<endl;

    
    cout<<"     probabilistic model...";
    delete m_ranking_model;
    cout<<"OK"<<endl;
    
    cout<<"     best solution...";
    delete BEST;
    cout<<"OK"<<endl;

    
    cout<<"     population...";
   delete m_population;
    
    cout<<"OK"<<endl;


}

void CSolution::Improve()
{
    int i;
    
    m_ranking_model->Learn(m_population);
    
    int * genes= new int[IND_SIZE];
    //cout<<"Sample new solutions and update population."<<endl;
	for(i=0;i<OFFSPRING_SIZE && EVALUATIONS<MAX_EVALUATIONS;i++)
	{
        m_ranking_model->Sample(genes);
        m_population->AddToPopulation(genes, i);
    }
    delete [] genes;
    
    m_population->SortPopulation(1);
    
    //This routine executes any adhoc post processes required by each model type.
    m_ranking_model->PostProcesses(m_population);
    
	//checking convergence of the best individual
    double newScore=GetBestPopulationIndividual()->Value();
    if (newScore>BEST->Value())
    {
        BEST->SetGenes(GetBestPopulationIndividual()->Genes());
        BEST->SetValue(GetBestPopulationIndividual()->Value());
        CONVERGENCE_ITERATION=GEN_NUM;
        CONVERGENCE_EVALUATIONS=EVALUATIONS;
    }
}


ostream & operator<<(ostream & os,CSolution & solution)
{

    os
        << "Instance to optimize: "<<DATA_FILE_NAME<<endl
        << "Problem type: "<<PROBLEM_TYPE<<endl
        << "Model type: "<<MODEL_TYPE<<endl
        << "Population size: " << POP_SIZE << endl
        << "Selected individuals: " << SEL_SIZE << endl
        << "Offspring size: " << OFFSPRING_SIZE << endl
        << "Generations: " << GEN_NUM << endl
        << "Convergence generation: "<<CONVERGENCE_ITERATION<<endl
        << "Evaluations: " << EVALUATIONS << endl
		<< "Convergence evaluation: "<<CONVERGENCE_EVALUATIONS<<endl
        << "Total time: " << solution.timer.TimeString() << endl
        << "Population best individual: " << solution.GetBestPopulationIndividual()<<endl
        << "EDA best individual: " << BEST<<endl;
return os;
}


void CSolution::StartTimer()
{
    //set the time at the beginning
    timer.Reset();
}


void CSolution::StopTimer()
{
    //Stop the timer.
    timer.End();
}


CIndividual  * & CSolution::GetBestPopulationIndividual()
{
    return m_population->m_individuals[0];
}


void CSolution::WriteSolution()
{
//  WriteSolutionProblem(*this);

    // Stop the timer
    StopTimer();

    //Write the information of the solution in the screen
    cout << *this << endl;

    //Write the information of the solution in the file (if exists)
    if (foutput)
    {
        foutput << "== final RESULTS =="<< endl
                << *this << endl;
        foutput.close();
    }
}


