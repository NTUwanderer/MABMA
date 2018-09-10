/*
 *  PlackettLuce.cpp
 *  MixtureMallowsEDA
 *
 *  Created by Josu Ceberio Uribe on 7/16/12.
 *  Copyright 2012 University of the Basque Country. All rights reserved.
 *
 */

#include "PlackettLuce.h"
#include "PlackettLuceInference.h"
#include "Tools.h"
#include "PermutationTools.h"
#include <math.h>
#include "SortingTools.h"
#include <complex>
using namespace std;
ofstream fweights;
/*
 * Class constructor.
 */
CPlackettLuceModel::CPlackettLuceModel()
{
    init(IND_SIZE, SEL_SIZE);
}


CPlackettLuceModel::CPlackettLuceModel(int problem_size, int samples_size)
{
    init(problem_size, samples_size);
}


void CPlackettLuceModel::init(int problem_size, int samples_size)
{
	m_problem_size=problem_size;
    m_samples_size=samples_size;
    
    m_MM_iterations=10000;
    m_twomodels_learning=false;
    m_sampling_by_threshold=false;
    
	m_weights= new float[m_problem_size];
    m_mostProbable= new CIndividual(m_problem_size);

    m_sampling_threshold=0.6;
    
    int aux1[m_problem_size];
    for (int i=0;i<m_problem_size;i++)
        aux1[i]=i;
    m_mostProbable->SetGenes(aux1);
    
    if (m_twomodels_learning)
        m_inverse_mostProbable->SetGenes(aux1);

    
    m_aux1= new int[m_problem_size];
    m_aux2= new int[m_problem_size];
    

    //creating sampling structures.
    distribution= new float[m_problem_size];
    objects= new int[m_problem_size];
    m_last= new int[m_problem_size];


        //create learning structures.
        m_M = m_problem_size;
        m_N = m_samples_size; //number of rankings in the population.
        m_P = m_problem_size;
    
        r2 = new float*[m_M];
        f= new int*[m_P];
        r= new int*[m_M];
        aux= new float*[m_P];
        g= new float*[m_P];
        for (int i=0;i<m_P;i++)
        {
            f[i]= new int[m_N];
            r[i]= new int[m_N];
            r2[i]= new float[m_N];
            aux[i]= new float[m_N];
            g[i]= new float[m_N];
        }
        
        dgamma= new float[m_M];
        gamma= new float[m_M];
        newgamma= new float[m_M];
    
        w= new float[m_M];
        pp= new int[m_N];
        sr2= new float[m_M];
   
    m_p1=new int[m_problem_size];
    m_p2= new int[m_problem_size];
    m_located= new int[m_problem_size];
}

/*
 * Class destructor.
 */
CPlackettLuceModel::~CPlackettLuceModel()
{
    
	delete [] m_weights;
    delete m_mostProbable;
    
    if (m_twomodels_learning)
    {
        delete [] m_inverse_weights;
        delete m_inverse_mostProbable;
    }
    
    delete [] m_last;
    delete [] distribution;
    delete [] objects;
    delete [] m_aux1;
    delete [] m_aux2;
    
        for (int i=0;i<m_P;i++)
        {
            delete [] f[i];
            delete [] r[i];
            delete [] r2[i];
            delete [] aux[i];
            delete [] g[i];
        }
        delete [] f;
        delete [] r;
        delete [] r2;
        delete [] aux;
        delete [] g;
    
        delete [] dgamma;
        delete [] gamma;
        delete [] newgamma;
        delete [] w;
        delete [] pp;
        delete [] sr2;
    
    delete [] m_p1;
    delete [] m_p2;
    delete [] m_located;
}

void cumsum(float * values, int size, float * resul)
{
	resul[0]=values[0];
	for (int i=1;i<size;i++)
	{
		resul[i]=resul[i-1]+values[i];
	}
}
void cumsum(float ** values,int size1, int size2, float ** resul)
{
	for (int j=0;j<size2;j++)
		resul[0][j]=values[0][j];
	
	for(int i=1;i<size1;i++)
		for (int j=0;j<size2;j++)
			resul[i][j]=values[i][j]+resul[i-1][j];
}

float sum (float * values, int size)
{
	float res=0;
	for (int i=0;i<size;i++) { res+=values[i];}
	return res;
}
double norm(float * values, int size)
{
	double result=0;
	for (int i=0;i<size;i++)
		result+=pow(values[i],2);
	return sqrt(result);
}

bool CPlackettLuceModel::Learn(CPopulation * population)
{
    bool result=LearnMM(population,m_samples_size);
   /* if (result==false)
    {
        cout<<"fail"<<endl;
    }
    if (m_twomodels_learning && result)
    {
        result=LearnMMInverse(population, m_samples_size);
    }*/
    return result;
}

/*
 * It estimates the Plackett-Luce model from the given cases. The learning process consists
 * of finding the MLE v parameters. For that task, Minorise-Maximize (MM) algorithm is utilized.
 */
bool CPlackettLuceModel::LearnMM(CPopulation * population, int sel_total)
{
	//This code is a parse of the code of Hunter et al. 2004 of the MM algorithm implemented in matlab.
	
	//Initialization of auxiliary data.
	int i,j, aux_i,aux_j;
//	int M = ProblemSize;
//	int N = sel_total; //number of rankings in the population.
//	int P = ProblemSize;

	int * individua;
	//extract individuals from population
	for (i = 0; i < m_N-1; i++)
	{
		individua = population->m_individuals[i]->Genes();
		for (j = 0; j < m_M; j++)
		{
            //we consider that individuals in the population have an ordering description
			//and these matrix consider ordering description of the solutions.
            f[j][i]=individua[j];
            r[individua[j]][i]= j + m_P*i;
		}
	}
    
    
    GenerateRandomPermutation(m_problem_size, m_last);
    for (j = 0; j < m_M; j++)
    {
        //we consider that individuals in the population have an ordering description
        //and these matrix consider ordering description of the solutions.
        f[j][i]=m_last[j];
        r[m_last[j]][i]= j + m_P*i;
    }
 
    //create a copy of r in r2.
	for (i=0;i<m_M;i++) for(j=0;j<m_N;j++)r2[i][j]=r[i][j];

    
	for (i=0;i<m_M;i++) w[i]=0;
	for (i=0;i<m_N;i++) pp[i]=m_M;

	for (i=0;i<m_N;i++)
		for (j=0;j<m_P-1;j++)
			w[f[j][i]]++;	
		
	
	for (i=0;i<m_N;i++) pp[i]+= i*m_P;
	
	//Starts the MM algorithm

    for (i=0;i<m_M;i++) {gamma[i]=1;dgamma[i]=1;}

	int iterations=0;
	//double ** aux= new double*[P];
	m_norm=1;

	//cout<<"Starting MM algorithm...."<<endl;
	while(m_norm>0.000000001 && iterations<m_MM_iterations)
	{
		//cout<<"norm: "<<norm(dgamma, M)<<"iterations: "<<iterations<<endl;
		//PRIMERA PARTE
		iterations++;

		// for (i=0;i<N;i++) for(j=0;j<N;j++)g[i][j]=f[i][j]; //esta linea en matlab hace falta, pero aqui no.
 
		for(i=0;i<m_P;i++)
			for(j=0;j<m_N;j++)
				g[i][j]=gamma[f[i][j]];
		
		//for (j=0;j<N;j++)
	//		aux[P-1][j]=g[P-1][j];
        
        memcpy(aux[m_P-1], g[m_P-1], sizeof(float)*m_N);
        
		for(i=m_P-2;i>=0;i--)
			for (j=0;j<m_N;j++)
				aux[i][j]=g[i][j]+aux[i+1][j];

		for(i=m_P-1;i>=0;i--)
			memcpy(g[i],aux[i],sizeof(float)*m_N);
	
		
		for(j=0;j<m_N;j++)
			g[m_P-1][j]=0;
	
		
		for(i=0;i<m_M-1;i++)
			for(j=0;j<m_N;j++)
				g[i][j]=(float)1/g[i][j]; //en este for anidado va a dar error porque la ultima fila de g es de zeros
		//como debe, y no estamos seguros d que se controle. Hay que ver si el N-1 esta bien puesto en el primer
		//for
		
		cumsum(g,m_P,m_N,g);

		//SEGUNDA PARTE
		for (j=0;j<m_N;j++)
		{
			for(i=0;i<m_M;i++)
			{
				aux_i= r[i][j] % m_P;
				aux_j= r[i][j] / m_P;
				r2[i][j]=g[aux_i][aux_j];
			}
		}
		
		for (i=0;i<m_M;i++)
		{
			sr2[i] = sum(r2[i],m_N);
			newgamma[i] = w[i] / sr2[i];
			dgamma[i] = newgamma[i] - gamma[i];
		}

		memcpy(gamma, newgamma, sizeof(float)*m_M);

        m_norm=norm(dgamma, m_M);
	}

    
	float sumGamma=sum(gamma,m_M);

    for (int i=0;i<m_M;i++) m_weights[i]=gamma[i]/sumGamma;
//    if (iterations==m_MM_iterations)
//        PrintArray(m_weights,IND_SIZE,"weights: ");

    CalculateMostProbableSolution();
    
   // cout<<"straight iterations: "<<iterations<<endl;
    return true;
}

/*
 * It estimates the Plackett-Luce model from the given cases. The learning process consists
 * of finding the MLE v parameters. For that task, Minorise-Maximize (MM) algorithm is utilized.
 */
bool CPlackettLuceModel::LearnMMInverse(CPopulation * population, int sel_total)
{
	//This code is a parse of the code of Hunter et al. 2004 of the MM algorithm implemented in matlab.
	
	//Initialization of auxiliary data.
	int i,j, aux_i,aux_j;
    //	int M = ProblemSize;
    //	int N = sel_total; //number of rankings in the population.
    //	int P = ProblemSize;
    
	int * individua;
	//extract individuals from population.
	for (i = 0; i < m_N-1; i++)
	{
		individua = population->m_individuals[i]->Genes();
		for (j = 0; j < m_M; j++)
        {
            //we consider that individuals in the population have an ordering description
			//and these matrix consider ordering description of the solutions.
            f[j][i]=individua[m_M-j-1];
            r[individua[m_M-j-1]][i]= j + m_P*i;
		}
	}
    
    for (j = 0; j < m_M; j++)
    {
        //we consider that individuals in the population have an ordering description
        //and these matrix consider ordering description of the solutions.
        f[j][i]=m_last[j];
        r[m_last[j]][i]= j + m_P*i;
    }

    
    //create a copy of r in r2.
	for (i=0;i<m_M;i++) for(j=0;j<m_N;j++)r2[i][j]=r[i][j];
    
    
	for (i=0;i<m_M;i++) w[i]=0;
	for (i=0;i<m_N;i++) pp[i]=m_M;
    
	for (i=0;i<m_N;i++)
		for (j=0;j<m_P-1;j++)
			w[f[j][i]]++;
    
	
	for (i=0;i<m_N;i++) pp[i]+= i*m_P;
	
	//Starts the MM algorithm
    
    for (i=0;i<m_M;i++) {gamma[i]=1;dgamma[i]=1;}
    
	int iterations=0;
	//double ** aux= new double*[P];
	m_norm=1;
    
	//cout<<"Starting MM algorithm...."<<endl;
	while(m_norm>0.000000001 && iterations<m_MM_iterations)
	{
		//cout<<"norm: "<<norm(dgamma, M)<<"iterations: "<<iterations<<endl;
		//PRIMERA PARTE
		iterations++;
        
		// for (i=0;i<N;i++) for(j=0;j<N;j++)g[i][j]=f[i][j]; //esta linea en matlab hace falta, pero aqui no.
        
		for(i=0;i<m_P;i++)
			for(j=0;j<m_N;j++)
				g[i][j]=gamma[f[i][j]];
		
		//for (j=0;j<N;j++)
        //		aux[P-1][j]=g[P-1][j];
        
        memcpy(aux[m_P-1], g[m_P-1], sizeof(float)*m_N);
        
		for(i=m_P-2;i>=0;i--)
			for (j=0;j<m_N;j++)
				aux[i][j]=g[i][j]+aux[i+1][j];
        
		for(i=m_P-1;i>=0;i--)
			memcpy(g[i],aux[i],sizeof(float)*m_N);
        
		
		for(j=0;j<m_N;j++)
			g[m_P-1][j]=0;
        
		
		for(i=0;i<m_M-1;i++)
			for(j=0;j<m_N;j++)
				g[i][j]=(float)1/g[i][j]; //en este for anidado va a dar error porque la ultima fila de g es de zeros
		//como debe, y no estamos seguros d que se controle. Hay que ver si el N-1 esta bien puesto en el primer
		//for
		
		cumsum(g,m_P,m_N,g);
        
		//SEGUNDA PARTE
		for (j=0;j<m_N;j++)
		{
			for(i=0;i<m_M;i++)
			{
				aux_i= r[i][j] % m_P;
				aux_j= r[i][j] / m_P;
				r2[i][j]=g[aux_i][aux_j];
			}
		}
		
		for (i=0;i<m_M;i++)
		{
			sr2[i] = sum(r2[i],m_N);
			newgamma[i] = w[i] / sr2[i];
			dgamma[i] = newgamma[i] - gamma[i];
		}
        
		memcpy(gamma, newgamma, sizeof(float)*m_M);

        m_norm=norm(dgamma, m_M);        
	}
    


	float sumGamma=sum(gamma,m_M);
    for (int i=0;i<m_M;i++) m_inverse_weights[i]=gamma[i]/sumGamma;

    CalculateMostProbableInverseSolution();

    return true;
}

/*
 * It applies random shake_power times a perturbation over the given individual.
 */
void CPlackettLuceModel::ShakeByInsert(int * permutation, int shake_power)
{
	int i,j;
	for (int iter=0;iter<shake_power;iter++)
	{
		//permute randomly genes in position i and j to scape the stackeness in both neighborhood.
		i = rand() % m_problem_size;
		j = rand() % m_problem_size;
        InsertAt2(permutation,i,j,m_problem_size);
	}
}

/*
 * Samples the plackett-luce model learnt by shaking the most probable solution.
 */
void CPlackettLuceModel::SampleByShake(int * permutation, int * most_probable)
{
    memcpy(permutation, most_probable,sizeof(int)*m_problem_size);
    ShakeByInsert(permutation,m_problem_size/20);
}

/*
 * Samples the plackett-luce model learnt.
 */
inline void CPlackettLuceModel::Sample(int * permutation)
{
	int i,j,k , remaining,pos;
	//int located[ProblemSize];
	float randVal,acumul,acumul_aux;
	for (i=0;i<m_problem_size;i++) { permutation[i]=-1;m_located[i]=0;}
	for (i=0;i<m_problem_size;i++)
	{

		//get ready data
		remaining=m_problem_size-i;
		acumul_aux=0;
		pos=0;
		for (j=0;j<m_problem_size;j++)
		{
			if (m_located[j]==0)//permutation[j]==-1)
			{
                distribution[pos]=m_weights[j];
				acumul_aux+=m_weights[j];
   				objects[pos]=j;
				pos++;
			}
		}

        randVal=drand48()*acumul_aux;
        //randVal=(float)rand()/(((float)RAND_MAX/acumul_aux));
        acumul=distribution[0];
        for(k=1;(k<remaining && acumul<randVal);k++)
        {
            acumul += distribution[k];
        }
        k--;
        
        permutation[i]=objects[k];
        m_located[objects[k]]=1;
    }
}


/*
 * Samples the straight plackett-luce model fixing the items that exceed from a probability threshold.
 */
void CPlackettLuceModel::SampleWithThreshold(int * permutation)
{
	int i,j,k , remaining,pos;
    float newThreshold;
	float randVal,acumul,acumul_aux;
	for (i=0;i<m_problem_size;i++) { permutation[i]=-1;m_located[i]=0;}
	for (i=0;i<m_problem_size;i++)
	{
		//get ready data
		remaining=m_problem_size-i;
		acumul_aux=0;
		pos=0;
		for (j=0;j<m_problem_size;j++)
		{
			if (m_located[j]==0)//permutation[j]==-1)
			{                
                distribution[pos]=m_weights[j];
				acumul_aux+=m_weights[j];
   				objects[pos]=j;
				pos++;
			}
		}
        newThreshold=m_sampling_threshold*acumul_aux;

        for (k=0;k<remaining && distribution[k]<newThreshold;k++);
        
        if (k==remaining)
        {
            //randVal=(float)rand()/(((float)RAND_MAX/acumul_aux));
            randVal=drand48()*acumul_aux;
            acumul=distribution[0];
            for(k=1;(k<remaining && acumul<randVal);k++)
            {
                acumul += distribution[k];
            }
        
            k--;
        }

        permutation[i]=objects[k];
        m_located[objects[k]]=1;
    }
}

/*
 * Samples the inverse plackett-luce model fixing the items that exceed from a probability threshold.
 */
void CPlackettLuceModel::SampleInverseWithThreshold(int *permutation)
{
	int i,j,k , remaining,pos;
    float newThreshold;
	float randVal,acumul,acumul_aux;
	for (i=0;i<m_problem_size;i++) { permutation[i]=-1;m_located[i]=0;}
	for (i=0;i<m_problem_size;i++)
	{
        
		//get ready data
		remaining=m_problem_size-i;
		acumul_aux=0;
		pos=0;
		for (j=0;j<m_problem_size;j++)
		{
			if (m_located[j]==0)//permutation[j]==-1)
			{
                distribution[pos]=m_inverse_weights[j];
				acumul_aux+=m_inverse_weights[j];
   				objects[pos]=j;
				pos++;
			}
		}
        newThreshold=m_sampling_threshold*acumul_aux;
        
        for (k=0;k<remaining && distribution[k]<newThreshold;k++);
        
        if (k==remaining)
        {

            //randVal=(float)rand()/(((float)RAND_MAX/acumul_aux));
            randVal=drand48()*acumul_aux;
            acumul=distribution[0];
            for(k=1;(k<remaining && acumul<randVal);k++)
            {
                acumul += distribution[k];
            }
        
            k--;
        }
                permutation[i]=objects[k];
        m_located[objects[k]]=1;

    }
}

/*
 * Samples the inverse plackett-luce model learnt.
 */
void CPlackettLuceModel::SampleInverse(int * permutation)
{
	int i,j,k , remaining,pos;
    //	int m_located[ProblemSize];
	float randVal,acumul,acumul_aux;
	for (i=0;i<m_problem_size;i++) { permutation[i]=-1;m_located[i]=0;}
	for (i=0;i<m_problem_size;i++)
	{
        
		//get ready data
		remaining=m_problem_size-i;
		acumul_aux=0;
		pos=0;
		for (j=0;j<m_problem_size;j++)
		{
			if (m_located[j]==0)//permutation[j]==-1)
			{
                distribution[pos]=m_inverse_weights[j];
				acumul_aux+=m_inverse_weights[j];
   				objects[pos]=j;
				pos++;
			}
		}
        
//		randVal=(float)rand()/(((float)RAND_MAX/acumul_aux));
         randVal=drand48()*acumul_aux;
        
		acumul=distribution[0];
		for(k=1;(k<remaining && acumul<randVal);k++)
		{
			acumul += distribution[k];
  		}
        
		k--;
        
      
        permutation[i]=objects[k];
        m_located[objects[k]]=1;
    }
}

/*
 * This method updates the best solution of the optimization processes if the most probable solutions
 * are better than those solutions obtained by the population.
 */
void CPlackettLuceModel::PostProcesses()
{
    if (BEST->Value()<m_mostProbable->Value())
    {
        BEST->SetGenes(m_mostProbable->Genes());
        BEST->SetValue(m_mostProbable->Value());
        CONVERGENCE_ITERATION=GEN_NUM;
        CONVERGENCE_EVALUATIONS=EVALUATIONS;
        
    }
    
   /* if (m_twomodels_learning)
    {
        if (BEST->Value()<m_inverse_mostProbable->Value())
        {
            BEST->SetGenes(m_inverse_mostProbable->Genes());
            BEST->SetValue(m_inverse_mostProbable->Value());
            CONVERGENCE_ITERATION=GEN_NUM;
            CONVERGENCE_EVALUATIONS=EVALUATIONS;
        }
    }*/
}

/*
 * Returns the probability of the 'indiv' individual for Plackett-Luce model.
 */
double CPlackettLuceModel::Probability(int * indiv,bool straight)
{
    double probability=1;
    double aux,aux2;
    int i,j;
    if (straight)
    {
        for (i=0;i<m_problem_size;i++)
        {
            aux=m_weights[indiv[i]];
            aux2=0;
            for (j=i;j<m_problem_size;j++)
                aux2=aux2+m_weights[indiv[j]];
            probability=probability*(aux/aux2);
        }
    }
    else
    {
        int genes[m_problem_size];
        for (i=0;i<m_problem_size;i++)
            genes[i]=indiv[m_problem_size-1-i];
        
        for (i=0;i<m_problem_size;i++)
        {
            aux=m_inverse_weights[genes[i]];
            aux2=0;
            for (j=i;j<m_problem_size;j++)
                aux2=aux2+m_inverse_weights[genes[j]];
            probability=probability*(aux/aux2);
        }
    }
	return probability;
}


/*
 * Returns the most probable solution given the vector of weights
 */
void CPlackettLuceModel::CalculateMostProbableSolution()
{
    RandomKeysDescending(m_aux1,m_weights,m_problem_size);
    //Invert
    for(int i=0; i<m_problem_size; i++) m_mostProbable->Genes()[m_aux1[i]]=i;
    
    m_mostProbable->SetValue(MIN_INTEGER);
}

/*
 * Returns the most probable solution given the vector of weights
 */
void CPlackettLuceModel::CalculateMostProbableInverseSolution()
{
    RandomKeysDescending(m_aux1,m_inverse_weights,m_problem_size);
    Invert(m_aux1,m_problem_size,m_aux2);
    for (int i=0;i<m_problem_size;i++)
        m_inverse_mostProbable->Genes()[i]=m_aux2[IND_SIZE-i-1];
    m_inverse_mostProbable->SetValue(MIN_INTEGER);
    
}


/*
 * Returns the most probable solution given the vector of weights
 */
CIndividual * CPlackettLuceModel::GetMostProbableSolution()
{
    return m_mostProbable;
}

/*
 * Returns the most probable solution given the vector of weights
 */
CIndividual * CPlackettLuceModel::GetMostProbableInverseSolution()
{
    return m_inverse_mostProbable;
}





