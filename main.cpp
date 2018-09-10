#include <string>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include "statistics.h"
#include "hybrid.h"
#include "chromosome.h"
#include "global.h"

using namespace std;

int main (int argc, char *argv[])
{
    const clock_t begin_time = clock();
    if (argc != 6) {
        printf("Program_name nInitial maxFe repeat Benchmark Problem_type\n");
        return -1;
    }
    
   
    int nInitial = atoi (argv[1]);
    int maxFe = atoi (argv[2]); 
    int repeat = atoi (argv[3]); 
    int bit = atoi(argv[4]);
    char *benchmark = argv[4];
    std::string Problem (argv[5]);
    int i;


    Statistics stGenS, stGenF, stFitness;
    int usedGen;

   
    char benchmarkLoc[100] = "";
    strcpy(benchmarkLoc, benchmark);
    
    IEvaluator *evaluator = NULL;
    if(Problem == "FSSP")   
    	evaluator = new FSSP_Evaluator(benchmarkLoc);
    else if(Problem == "LOP")
	    evaluator = new LOP_Evaluator(benchmarkLoc);
    else if(Problem == "QAP")    
	    evaluator = new QAP_Evaluator(benchmarkLoc);
    else if(Problem == "TSP")
	    evaluator = new TSP_Evaluator(benchmarkLoc);
    else if(Problem == "TSPF")
	    evaluator = new TSPFixed_Evaluator(benchmarkLoc);
    else if(Problem == "FLAT")
	    evaluator = new Flat_Evaluator(bit);
    else if(Problem == "CVRP")
      evaluator = new CVRP_Evaluator(benchmarkLoc);
    else{
	    printf("ERROR: Prblem invalid, %s\n", Problem.c_str());
      exit(1);
    }
        
    int ell = evaluator->getProblemSize();	
    for (i = 0; i < repeat; i++) {
      numOfReplacement = 0;
      evaluator->init();
      Hybrid hybrid ( ell, nInitial, maxFe, evaluator );
      usedGen = hybrid.doIt (false);
      Chromosome ch(ell);
      
      
      
      if (hybrid.stFitness.getMin() - 0.01 > evaluator->getBest()) {
        cout << "-, Fitness: " << hybrid.stFitness.getMin() << ", RPD: " << ((evaluator->getBest() - hybrid.stFitness.getMin()) / evaluator->getBest()) * 100 << ", Generation: " << usedGen << endl;
      }
      else {
        cout << "+, Fitness: " << hybrid.stFitness.getMin() << ", RPD: " << ((evaluator->getBest() - hybrid.stFitness.getMin()) / evaluator->getBest()) * 100 << ", Generation: " << usedGen << endl;
        stGenS.record (usedGen);
      }
      
      stFitness.record( hybrid.stFitness.getMin() );

      fflush (NULL);

    }
    
    printf ("\nAverage Gen of Success: %f\n", stGenS.getMean());
    printf ("Average Gen of Failure: %f\n", stGenF.getMean());
	  printf("Average Fitness: %f\n", stFitness.getMean());
    printf("ARPD: %f\n", 100 * (evaluator->getBest() - stFitness.getMean()) / evaluator->getBest());
    printf ("Execution time: %f\n", float(clock() - begin_time)/ CLOCKS_PER_SEC);


    delete evaluator;

    return EXIT_SUCCESS;
}
