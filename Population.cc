//
//  Population.cc
//  GMEDAforTEC
//
//  Created by Josu Ceberio Uribe on 2/7/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "Population.h"

using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;

/*
 * Constructor function.
 */
CPopulation::CPopulation(int pop_size, int offspring_size, int individual_size)
{
    m_pop_size=pop_size;
    m_offspring_size=offspring_size;
    m_size=pop_size+offspring_size;
    m_individuals.resize(m_size);
    
    //Initialize population with empty solutions
    for (int i=0;i<m_size;i++)
    {
        m_individuals[i]= new CIndividual(individual_size);
    }
}

/*
 * Destructor function.
 */
CPopulation::~CPopulation()
{
    for (int i=0;i<m_size;i++)
    {
        delete m_individuals[i];
    }
    m_individuals.clear();
}

/*
 * Function to add an individual to the population
 */
void CPopulation::AddToPopulation(int * genes, int index)
{
    m_individuals[m_pop_size+index]->SetGenes(genes);
}

/*
 * Function to set an individual in a specific position of the population
 */
void CPopulation::SetToPopulation(int * genes, int index)
{
    m_individuals[index]->SetGenes(genes);
}

void CPopulation::SetToPopulation(int* genes, int index, int fitness)
{
    m_individuals[index]->SetGenes(genes, fitness);
}


/*
 * Prints the current population.
 */
void CPopulation::Print()
{
	for(int i=0;i<m_pop_size; i++)
		cout<<m_individuals[i]<<endl;
}

/*
 * Determines if the individuals in the population are equal.
 */
bool CPopulation::Same(int size)
{
	int best=m_individuals[0]->Value();
	for(int i=1; i<size; i++)
	{
		if (m_individuals[i]->Value()!=best)
			return false;
	}
	return true;
}

/*
 * Calculates the average fitness of the first 'size' solutions in the population
 */
float CPopulation::AverageFitnessPopulation(int size)
{
	float result=0;
	for(int i=0;i<size;i++)
    {
        result+=m_individuals[i]->Value();
    }
    return result/size;
    
}


/*
 * Sorts the individuals in the population in decreasing order of the fitness value.
 * mode = 0 means sort only the population.
 * mode = 1 means sort population and offsprings population together.
 */
void CPopulation::SortPopulation(int mode)
{
    if (mode==0)
        sort(m_individuals.begin(), m_individuals.begin()+m_pop_size, Better);
    else
        sort(m_individuals.begin(), m_individuals.end(), Better);
}


