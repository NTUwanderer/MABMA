#ifndef Hybrid_H
#define Hybrid_H

#include <vector>
#include <fstream>
#include "evaluator.h"
#include "chromosome.h"
#include "statistics.h"
#include "myrand.h"
#include "model.h"


class Hybrid
{

    public:
        Hybrid ();
        Hybrid (int n_ell, int n_nInitial, int n_maxFe, IEvaluator *evaluator );

        ~Hybrid ();

        void init (int n_ell, int n_nInitial, int n_maxFe, IEvaluator *evaluator );

        void initializePopulation ();
        void initializeModelBuilders ();
        void evaluate();

        void selection();
        void tournamentSelection ();
        bool randomPickAndSampleByIndex(IModel*, int index);
        void ModelBuilding ();
        void UCBTunedSample();
	

        void showStatistics ();
        void oneRun (bool output = true);
        int doIt (bool output = true);

        bool shouldTerminate ();
		
        Statistics stFitness;
 
        int FindClosest(Chromosome*, Chromosome&, int m, int tempChros);
        int FindNearestModelWT(Chromosome*, Chromosome&, int m, int tempChro);

    protected:

        int ell;                 // chromosome length
        int nInitial;            // initial population size
        int nCurrent;            // current population size
        int nNextGeneration;     // population size for the next generation
        int selectionPressure;

        Chromosome *population;

        
        int maxFe;
        int repeat;
        int fe;
        int generation;
        int bestIndex;
        IEvaluator *evaluator;
        std::vector<IModel*> ModelList;
        std::vector<double> ModelWeight;
        
        std::vector<Statistics> UCB_mean;
        std::vector<int> UCB_times;
        int UCB_total_times;
        
        std::vector<std::ofstream*> ofileList;
        std::ofstream fin;

    private:
    	  bool output;
		    MyRand random;
};
#endif
