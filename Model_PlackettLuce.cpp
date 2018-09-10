#include "chromosome.h"
#include "myrand.h"
#include "Model_PlackettLuce.h"

Model_PlackettLuce::Model_PlackettLuce()
{
    this->init(0, NULL, 0);
}

Model_PlackettLuce::Model_PlackettLuce(int n_ell, Chromosome* pop, int size)
{
    this->init(n_ell, pop, size);
}

Model_PlackettLuce::~Model_PlackettLuce()
{
	delete plModel;
	delete m_population;
}

void Model_PlackettLuce::init(int n_ell, Chromosome* pop, int size)
{
	MODEL_NAME = "PlackettLuce";
	
	problemSize = n_ell;
    population = pop;
    n_size = size;
    
    plModel = new CPlackettLuceModel(problemSize, size);
    m_population = new CPopulation(n_size, 0, problemSize);
}

void Model_PlackettLuce::build()
{
    int * permutation= new int[problemSize];
	for(int i = 0; i < n_size; i++)
	{
        for(int j = 0; j < problemSize; j++)
            permutation[j] = population[i].getVal(j);
        m_population->SetToPopulation(permutation, i, (int)population[i].getFitness());
    } 
    delete[] permutation;
    m_population->SortPopulation(0);
    plModel->Learn(m_population);
}

Chromosome Model_PlackettLuce::resample()
{
	int templateIndex = random.uniformInt(0, n_size-1);
	return this->resampleWithTemplate(&population[templateIndex]);
}

Chromosome Model_PlackettLuce::resampleWithTemplate(const Chromosome* c)
{
    int * genes= new int[problemSize];
    
    
    plModel->Sample(genes);
    
    Chromosome n;
    n.clone(*c);

    for(int i = 0; i < problemSize; i++)
    {
        n.setVal(i, genes[i]);

    }

    delete[] genes;
    return n;
}

void Model_PlackettLuce::update()
{
    
}

void Model_PlackettLuce::dump()
{

}
