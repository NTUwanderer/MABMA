// Problem.cpp: implementation of the methods specific to the problem.
//


#include "Problem.h"
#include "EDA.h"

/*
 * Reads the problem info of the instance set.
 */
int GetProblemInfo(string problemType, string filename)
{
	if (problemType=="fsp")
	{
		//problem is FSP
		m_problem= new PFSP();
	}
	else if (problemType=="lop")
	{
		//problem is LOP
		m_problem= new LOP();
    }
	else if (problemType=="tsp")
	{
		//problem is TSP
		m_problem= new TSP();
    }
	else if (problemType=="qap")
	{
		//problem is QAP
		m_problem= new QAP();
    }
    else
    {
        cout<<"System error creating problem instance..."<<endl;
        exit(1);
    }
    
    int problemsize= m_problem->Read(filename);
    cout << "Problem type: " << problemType << ", Filename: " << filename << ", Problem size: "<< problemsize<<endl;
    
    POP_SIZE=problemsize*10;
	SEL_SIZE=problemsize;
	OFFSPRING_SIZE=POP_SIZE-1;
	IND_SIZE = problemsize;
	MAX_GENERATIONS=100*IND_SIZE;
	MAX_EVALUATIONS=POP_SIZE*MAX_GENERATIONS;
	MAX_GENERATIONS_NO_IMPROVEMENT=10*IND_SIZE;
	MAX_RESTARTS_NO_IMPROVEMENT=10*IND_SIZE;
    
	return problemsize;
}

/*
 * Frees the resources used by the values of the instance read.
 */
void RemoveProblemInfo()
{
	// The memory allocated in GetProblemInfo() is returned.
	delete m_problem;
}

/*
 * Calculates the objective function of the given genes.
 */
int Metric(int * genes, int size)
{
	return m_problem->Evaluate(genes);
}





