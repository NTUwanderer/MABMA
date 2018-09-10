#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <sys/stat.h> 
#include <ctime>
#include <cstring>

#include "hybrid.h"
#include "global.h"
#include "statistics.h"
#include "myrand.h"
#include "mkpath.h"
#include "model.h"
#include "model_EHBSA.h"
#include "model_NHBSA.h"
#include "Model_PlackettLuce.h"
#include "dice.h"


Hybrid::Hybrid ()
{
    ell = 0;
    nInitial = 0;
    nCurrent = 0;
    fe = 0;
    nNextGeneration = 0;
    maxFe = -1;
    population = NULL;
    evaluator = NULL;
    ModelList.clear();
    ofileList.clear();

}


Hybrid::Hybrid (int n_ell, int n_nInitial, int n_maxFe, IEvaluator *evaluator )
{
    init (n_ell, n_nInitial, n_maxFe, evaluator );
}


Hybrid::~Hybrid ()
{
    for(int i = 0; i < (int)ofileList.size(); i++){
	    ofileList[i]->close();
	    delete ofileList[i];
    }
    
    for(int i = 0; i < (int)ModelList.size(); i++){
		  delete ModelList[i];
    }
    
	
    delete[] population;
}


void
Hybrid::init (int n_ell, int n_nInitial, int n_maxFe, IEvaluator *_evaluator )
{
    int i;
    ell = n_ell;
    nInitial = n_nInitial;
    nCurrent = nInitial;
    maxFe = n_maxFe;
	  evaluator = _evaluator;

    population = new Chromosome[nInitial];


    for (i = 0; i < nInitial; i++) {
        population[i].init (ell, evaluator);
    }

	  initializePopulation ();
	  initializeModelBuilders ();
}

/* Initial the permutations */
void Hybrid::initializePopulation ()
{
    for (int i = 0; i < nInitial; i++)
	{
		int* temp = new int[ell];
		myRand.uniformArray(temp, ell, 0, ell-1);
			
        for (int j = 0; j < ell; j++)
		{
			population[i].setVal(j, temp[j]);
		}
		delete[] temp;
	}
}



void Hybrid::initializeModelBuilders()
{
	ModelList.clear();
	

  ModelList.push_back(new Model_EHBSA(ell, population, nInitial));
  ModelList.push_back(new Model_NHBSA(ell, population, nInitial));
  ModelList.push_back(new Model_PlackettLuce(ell, population, nInitial));
 
	
	int size = ModelList.size();
	ModelWeight.resize(size);
    UCB_mean.resize(size);
    UCB_times.resize(size);
    UCB_total_times = 0;
	for(int i=0; i < size; i++){
		ModelWeight[i] = (double)1 / (double)size;
        UCB_mean.at(i).reset();
        UCB_times.at(i) = 0;
	}
}



bool Hybrid::randomPickAndSampleByIndex(IModel* model, int index)
{
  int t = random.uniformInt(0, nInitial-1);
  int tempt = t;
  Chromosome child = model->resampleWithTemplate(&population[t]);
  if (isRTR) {
      t = FindClosest(population, child, index, tempt);
  } // isRTR
  
  if(child.getFitness() < population[t].getFitness())
  {
        if (tempt != t)
          numOfReplacement++;
        population[t].clone(child);
        return true;
  }
  else{
    return false;
  }
  return false;
}




void Hybrid::ModelBuilding()
{	

	for(int i = 0; i < (int)ModelList.size(); i++)
	{
	    IModel *m = ModelList.at(i);
	    m->build();
	
	}
}



void Hybrid::UCBTunedSample()
{
  int numModels = (int) ModelList.size();
  if(generation == 0)
  {
    for(int i = 0; i < numModels; i++)
    {
      IModel *m = ModelList.at(i);
      UCB_times[i]++;
      UCB_total_times++;
      if (this->randomPickAndSampleByIndex(m, i) == true)
      {
        UCB_mean[i].record(1.0);
      }else
      {
        UCB_mean[i].record(0.0);
      }
    }
  }else
  {
    //find the maximun UCB
    int index = -1;
    double max = -DBL_MAX;
    for(int i = 0; i < numModels; i++)
    {
       double v = UCB_mean[i].getVariance() + sqrt(2*log(UCB_total_times) / (double)UCB_times[i]);
       double V = (v > 0.25)? 0.25 : v;
       double ucb = 
         UCB_mean[i].getMean() + sqrt( 2*log(UCB_total_times) / (double)UCB_times[i] * V);
       if(ucb > max)
       {
         max = ucb;
         index = i;
       }



    }

  
    IModel *m = ModelList.at(index);
    UCB_total_times++;
    UCB_times[index]++;

    if(this->randomPickAndSampleByIndex(m,index) == true)
    {
      UCB_mean[index].record(1.0);
    }
    else
    {
      UCB_mean[index].record(0.0);  
    }

  }
}


void Hybrid::showStatistics ()
{
    printf ("Gen:%d  Fitness:(StD/Max/Mean/Min):%f/%f/%f/%f Chromosome Length:%d\n",
        generation, stFitness.getStdev(),stFitness.getMax (), stFitness.getMean (),
        stFitness.getMin (), population[0].getLength ());
    printf ("best chromosome:");
    population[bestIndex].printf();
}




void Hybrid::oneRun (bool _output)
{
    int i;
    output = _output;

    ModelBuilding();
    UCBTunedSample();
    


    double min = DBL_MAX;
    stFitness.reset ();
    for (i = 0; i < nCurrent; i++) {
        double fitness = population[i].getFitness ();
        if (fitness < min) {
            min = fitness;
            bestIndex = i;
        }
        stFitness.record (fitness);
    }

    if (output)
        showStatistics ();


    generation++;
}


int Hybrid::doIt (bool output)
{
    generation = 0;

    while (!shouldTerminate ()) {
        oneRun (output);
    }
    
    return generation;
}


bool Hybrid::shouldTerminate ()
{

  bool termination = false;

    // Reach maximal # of function evaluations
    if (maxFe != -1) {
       if (evaluator->getNFE() >= (unsigned int)maxFe)
            termination = true;
    }

    // Found a satisfactory solution
    if (stFitness.getMin() - 0.001 <= evaluator->getBest())
      termination = true;

    // The population loses diversity
    if (stFitness.getMax()-1e-6 < stFitness.getMean())
        termination = true;
		
    

    return termination;

}



int Hybrid::FindClosest(Chromosome* p, Chromosome& c, int m, int tempChro)
{
  return FindNearestModelWT(p, c, m, tempChro);
}



int Hybrid::FindNearestModelWT(Chromosome* p, Chromosome& c, int m, int tempChro)
{
  double dis = 9999999;
  int near = 0;
  double temp = 0;



  if (m == 0)
    R1 = 0.5;
  else if (m == 1)
    R1 = 0.2;
  else if (m == 2)
    R1 = 0.1;


  int aint = ell*R1;
  int *RTRarr = new int[aint];
  if (aint > nInitial)
    aint = nInitial;
  myRand.uniformArray(RTRarr, aint, 0, nInitial - 1);


  if (m == 0)
    dis = c.edgeDis(p[tempChro]);
  else if (m == 1)
    dis = c.nodeDis2(p[tempChro]);
  else if (m == 2)
    dis = c.orderDis(p[tempChro]);

  near = tempChro;


  for (int i = 0; i<aint; i++){

    if (m == 0)
      temp = c.edgeDis(p[RTRarr[i]]);
    else if (m == 1)
      temp = c.nodeDis2(p[RTRarr[i]]);
    else if (m == 2)
      temp = c.orderDis(p[RTRarr[i]]);

    if (dis>temp){
      dis = temp;
      near = RTRarr[i];
    } // if
  } // for



  return near;


}




