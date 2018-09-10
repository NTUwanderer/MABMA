//
//  MallowsModel.cpp
//  RankingEDA
//
//  Created by Josu Ceberio Uribe on 1/23/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "MallowsModel.h"

#include "EDA.h"
#include "Problem.h"
#include "SortingTools.h"
#include "PermutationTools.h"
#include "Tools.h"
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <list>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;


/*
 * Class constructor.
 */
CMallowsModel::CMallowsModel()
{
    m_problem_size=IND_SIZE;
    m_consensus_ranking=new int[m_problem_size];
    
    m_sample_size=SEL_SIZE;
    m_samples= new int*[m_sample_size];
    
    m_frequency_matrix= new int*[m_problem_size];
    for (int i=0;i<m_problem_size;i++)
        m_frequency_matrix[i]=new int[m_problem_size];
    
    m_psi_parameters= new float[m_problem_size-1];
    
    //reservar memoria para las probabilidades de Vs.
    m_VProbs=new float*[m_problem_size-1];
    for (int i=0;i<m_problem_size-1;i++)
    {
        m_VProbs[i]=new float[m_problem_size];
    }
    
    
    //fix upper bound for theta.
    if (m_problem_size<=100)
        m_upper_theta_bound=4.7;
    else if (m_problem_size<=200)
        m_upper_theta_bound=5.5;
    else if (m_problem_size<=250)
        m_upper_theta_bound=4.4;
    else
        m_upper_theta_bound=6;
    
    m_lower_theta_bound=m_upper_theta_bound*0.25;
    
    //sampling stage auxiliary data structures
    aux= new int[m_problem_size];
    aux_v= new int[m_problem_size]; //last position always 0.
    
    //learning stage auxiliary data structures
    VjsMean= new float[m_problem_size-1];
    VjsNonMean= new int[m_problem_size-1];
    Vjs = new int[m_problem_size-1];
    invertedB = new int[m_problem_size];
    composition = new int[m_problem_size];
}

/*
 * Class destructor.
 */
CMallowsModel::~CMallowsModel()
{
    delete [] m_consensus_ranking;
    delete [] m_samples;
    
    for (int i=0;i<m_problem_size-1;i++)
    {
        delete [] m_VProbs[i];
    }
    delete [] m_VProbs;
    delete [] m_psi_parameters;
    
    //sampling stage auxiliary data structures
    delete [] aux;
    delete [] aux_v;
    
    //learning stage auxiliary data structures
    delete [] VjsMean;
    delete [] VjsNonMean;
    delete [] Vjs;
    delete [] invertedB;
    delete [] composition;
}

/******************************  SIMULATING MODEL  *******************************/
/*
 * Simulates a new individual sampling the probabilistic model.
 */
CIndividual * CMallowsModel::Simulate()
{
    //1.- Create New individual.
    //cout << "1. Create a new individual."<<endl;
    CIndividual * individual = new CIndividual(m_problem_size);
    int* genes=individual->Genes();
    
    //2.- Sample permutation from Vs.
    //cout << "2. Sample permutation."<<endl;
    Sample(genes);
    
    return individual;
}

/*
 * Given the consensus ranking, it samples a new individual.
 */
void CMallowsModel::Sample(int * permutation)
{
    int i, index, limit;
    float randVal, acumul;
    
    for (i=0;i < m_problem_size;i++) aux[i]=-1;
    
    //generate samples and calcualte likelihood
    //cout << "2.1. Sample a Vj."<<endl;
    
    for(i=0;i<m_problem_size-1;i++)
    {
        //muestreo de las n-1 posiciones del vector V que definen la permutacion
        randVal=(float)rand()/((float)RAND_MAX+1);
        //      randVal=drand48();
        acumul=m_VProbs[i][0];
        index=0;
        limit=m_problem_size-1-i;
        for(;(index<limit && acumul<randVal);index++)
            acumul += m_VProbs[i][index+1];
        aux_v[i]=index;
    }
    aux_v[m_problem_size-1]=0;
    
    GeneratePermuFromV(aux_v,aux,m_problem_size);
    
    //Compose(aux,m_consensus_ranking,permutation,m_problem_size);
    for(int i=0;i<m_problem_size;i++)
        permutation[i]=aux[m_consensus_ranking[i]];
}


bool CMallowsModel::Learn(CPopulation * population)
{
    //0.- Preprocess the individuals to fit in the population.
    for(int i=0;i<m_sample_size;i++)
        m_samples[i] = population->m_individuals[i]->Genes();
    
    //1.- Calculate the consensus ranking/permutation.
    CalculateConsensusRanking(m_samples, m_sample_size, m_consensus_ranking);
    
    //2.- Calculate theta parameters
    m_theta_parameter= CalculateThetaParameter(m_consensus_ranking, m_samples, m_sample_size, m_upper_theta_bound);
    
    //3.- Calculate psi constant.
    CalculatePsiConstants(m_theta_parameter, m_psi_parameters);
    
    //4.- Calculate Vj probabilities matrix.
    float upper, lower, checkSum;
    int j,r;
    for(j=0;j<m_problem_size-1;j++)
    {
        checkSum=0;
        for(r=0;r<m_problem_size-j;r++)
        { //v(j,i) proba dql v_j tome valor i
            upper=exp((-1) * r * m_theta_parameter);
            lower=m_psi_parameters[j];
            m_VProbs[j][r] =  upper/lower;
            checkSum+=m_VProbs[j][r];
        }
    }
    return true;
}

/*
 * Calculates de consensus permutation from the given population cases.
 */
void CMallowsModel::CalculateConsensusRanking(int** samples, int sample_size, int* consensus_ranking)
{
    int i, j;
    
    //Initialize matrix to zeros.
    for (i=0;i<m_problem_size;i++)
    {
        for (j=0;j<m_problem_size;j++)
            m_frequency_matrix[i][j]=0;
    }
    
    int geneValue, genePosition;
    
    //Fill frecuency matrix reviewing individuals in the population.
    //para cada permutación muestreada:
    for (i = 0; i < sample_size; i++)
    {
        for (j=0;j<m_problem_size;j++)
        {
            geneValue=j;
            genePosition=samples[i][j];
            m_frequency_matrix[genePosition][geneValue]++;
        }
    }
    
    //cout<<"Calculate consensus vector: "<<endl;
    //Calculate consensus vector.
    float * consensusVector=new float[m_problem_size];
    int job,position;
    float positionMean;
    for (job=0;job<m_problem_size;job++)
    {
        positionMean=0;
        for (position=0;position<m_problem_size;position++)
        {
            positionMean=positionMean+position*m_frequency_matrix[position][job];
        }
        //positionMean=positionMean/sample_size;
        consensusVector[job]=positionMean;
    }
    
    RandomKeys(consensus_ranking, consensusVector,m_problem_size);
    
    delete[] consensusVector;
}


/*
 * Calculates the spread theta parameter from the ConsensusRanking and the individuals in the population.
 */
float CMallowsModel::CalculateThetaParameter(int*consensus_ranking, int** samples, int samples_num, float upper_bound)
{
    float initialTheta = 0.001;
    
    int i,j;
    for (i=0;i<m_problem_size-1;i++) {VjsNonMean[i]=0;VjsMean[i]=0;}
    
    //Calculate Vjmean vector from the population
    //cout << "Calculate Vjsmean vector"<<endl;
    Invert(consensus_ranking,m_problem_size,invertedB);
    
    for (i=0;i<samples_num;i++)
    {
        int* individua=samples[i];
        
        Compose(individua,invertedB,composition,m_problem_size);
        
        vVector(Vjs, composition, m_problem_size);
        
        for (j=0;j<m_problem_size-1;j++)
        {
            VjsNonMean[j]=VjsNonMean[j]+Vjs[j];
        }
    }
    
    float valuee;
    for (i=0;i<m_problem_size-1;i++)
    {
        valuee=(float)VjsNonMean[i]/samples_num;
        VjsMean[i]=valuee;
    }
    
    //calculate array of thetas.
    
    float theta=NewtonRaphsonMethod(initialTheta, VjsMean, upper_bound);
    
    return theta;
}

/*
 * Calculates the total Psi normalization constant from the ThetaParameter and psi-s vector.
 */
void CMallowsModel::CalculatePsiConstants(float theta, float* psi)
{
    //en el param psi se dejan los valores de las ctes d normalización psi_j
    //returns psiTotal: la multiplicacion de los n-1 psi_j
    int n=m_problem_size;
    //calculate psi
    int i, j;
    for(i=0;i< n - 1;i++)
    {
        j=i+1;//el indice en el paper es desde =1 .. n-1
        psi[i] = (1 - exp((-1)*(n-j+1)*(theta)))/(1 - exp((-1)*theta));
    }
}

/*
 * Generates a permutation from a v vector.
 */
void CMallowsModel::GeneratePermuFromV(int*v,int*permu,int n)
{
    
    int i,j;
    
    permu[v[0]]=0;
    for (i=1;i<n;i++){
        j=0;
        while (v[i]>=0)
        {
            v[i] = v[i] - (permu[j]==-1);
            j++;
        }
        permu[j-1]=i;
    }
    
    
    /*
     permu[0]=n-1;
     int at,i;
     for (int element=n-2;element>=0;element--)
     {
     at=v[element];
     
     for (i=n-2;i>=at;i--)
     {
     permu[i+1]=permu[i];
     }
     permu[at]=element;
     
     //InsertAt(permu,element,at,n);
     }
     */

}

/*
 * Calculates the probability of a given individuals in the current model.
 */
float CMallowsModel::Probability(int * individual)
{
    Invert(individual, m_problem_size, invertedB);
    Compose(m_consensus_ranking,invertedB,composition,m_problem_size);
    vVector(aux_v,composition,m_problem_size);
    
    float dist=0;//Kendall(permu,ConsensusRanking,ProblemSize);
    float Psi=1;
    for (int i=0;i<m_problem_size-1;i++)
    {
        dist+=aux_v[i];
        Psi=Psi*m_psi_parameters[i];
    }
    
    float probability=exp(-m_theta_parameter*dist)/Psi;
    
    return probability;
}


/*
 * Estimates theta parameters by means of running newton iterative algorithm.
 */
float CMallowsModel::NewtonRaphsonMethod(float initialGuess, float* Vjs, float upper_bound)
{
    float xacc=0.0001;
    float theta= rtsafe(initialGuess, upper_bound, xacc, Vjs);
    return theta;
}

/*
 * Calculates Newton algorithms f and fdev functions
 */
void CMallowsModel::funcd(float theta, float *ff ,float *ffdev,  float* Vjs)
{
    int n = m_problem_size;
    *ff = f(theta, n,Vjs);
    *ffdev =fdev(theta,n);
}

/*
 * Newton - Rapshon execution algorithm.
 */
float CMallowsModel::rtsafe(float x1, float x2, float xacc, float* Vjs)
{
    int j;
    float df,dx,dxold,f,fh,fl;
    float temp,xh,xl,rts;
    
    funcd(x1,&fl,&df,Vjs);
    funcd(x2,&fh,&df,Vjs);
    //if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    //  cout<<"Root must be bracketed in rtsafe"<<endl;
    if (fl == 0.0)
        return x1;
    if (fh == 0.0)
        return x2;
    if (fl < 0.0)
    {
        xl=x1;
        xh=x2;
    }
    else
    {
        xh=x1;
        xl=x2;
    }
    rts=0.5*(x1+x2);
    rts=x1;//<-fijamos un valor inicial para el theta.
    dxold=fabs(x2-x1);
    dx=dxold;
    funcd(rts,&f,&df,Vjs);
    for (j=1;j<=1000;j++)
    {
        //Initialize the guess for root, the “stepsize before last,” and the last step.
        //Loop over allowed iterations.
        if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) //Bisect if Newton out of range,
            || (fabs(2.0*f) > fabs(dxold*df)))
        {   //or not decreasing fast enough.
            dxold=dx;
            dx=0.5*(xh-xl);
            rts=xl+dx;
            if (xl == rts)
                return rts;
        }
        else
        {
            dxold=dx;
            dx=f/df;
            temp=rts;
            rts -= dx;
            if (temp == rts)
                return rts;
            //Change in root is negligible. Newton step acceptable. Take it.
        }
        //cout<<"DX: "<<dx<<endl;
        if (fabs(dx) < xacc){
            return rts;
        }
        funcd(rts,&f,&df,Vjs); //The one new function evaluation per iteration.
        //Orient the search so that f(xl) < 0.
        //Convergence criterion.
        
        if (f < 0.0)    //Maintain the bracket on the root.
            xl=rts;
        else
            xh=rts;
        
        //cout<<"rts: "<<rts<<" xl: "<<xl<<" xh: "<<xh<<endl;
        
    }
    cout<<"Maximum number of iterations exceeded in rtsafe"<<endl;

    return m_upper_theta_bound;
}

/*
 * Theta parameter estimation function.
 */
inline float CMallowsModel::f(float theta, int n, float* VjMeans)
{
    float oper1,oper2;
    float oper=( 1 /(float)( exp(theta) -1 ) );
    float results=0;
    for (int j=1;j<n;j++)
    {
        oper1=( (n - j + 1) /(float)(exp((n - j + 1)*theta) - 1 ) );
        oper2=VjMeans[j-1];
        results=results + oper -oper1 -oper2;
    }
    return results;
}

/*
 * Theta parameter estimation function derivation.
 */
inline float CMallowsModel::fdev(float theta, int n)
{
    float oper1= (-1)*exp(theta) / pow(exp(theta)-1, 2);
    float oper=0;
    for (int j=1;j<n;j++)
    {
        oper= oper+ oper1 +  ( pow(n-j+1,2) * exp(theta*(n-j+1)) ) /(float)pow(exp(theta*(n-j+1))-1,2);
    }
    //cout<<"fdev: "<<oper<<endl;
    return oper;
}
