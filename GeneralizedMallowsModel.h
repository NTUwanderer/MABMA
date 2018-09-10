/*
 *  GeneralizedMallowsModel.h
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 4/26/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */


#ifndef GENERALIZEDMALLOWSMODEL_H__
#define GENERALIZEDMALLOWSMODEL_H__

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>


class CGeneralizedMallowsModel : public CRankingModel
{
	
public:
	
    /*
     * The constructor.
     */
	CGeneralizedMallowsModel();
	
    /*
     * The destructor.
     */
    virtual ~CGeneralizedMallowsModel();
	
    /*
     * Given a population of samples, it learns a Generalized Mallows model from the data.
     */
    bool Learn(CPopulation * population);
    
    /*
     * Creates a new individual sampling the generalized mallows probabilistic model.
     */
    CIndividual * Simulate();
    
	/*
     * Calculates ThetaParameters average value.
     */
	float GetThetasAverage();
	
	/*
	 * Calculates the probability of the individual given the probabilistic model.
	 */
	float Probability(int * individual);
    
    /*
     * Calculates the post process of the GM model. Restarts etc etc
     */
    void PostProcesses(CPopulation * population);
	
private:
    
    /*
     * Problem size.
     */
    int m_problem_size;
    
    /*
     * Sample size.
     */
    int m_sample_size;
	
	/*
     * Consensus ranking.
     */
	int* m_consensus_ranking;
    
    /*
     * Theta parameters.
     */
    float * m_theta_parameters;
    
    /*
     * Psi parameters.
     */
    float * m_psi_parameters;
    
    /*
     * Matrix of V probabilities.
     */
	float ** m_VProbs;
    
    /*
     * The set of samples used to learn the model
     */
    int ** m_samples;
    
    /*
     * Frequency matrix employed when calculating the consensus ranking.
     */
    int ** m_frequency_matrix;
    
    /*
     * Defines upper theta bound.
     */
    float m_upper_theta_bound;
    
    /*
     * Defines lower theta bound.
     */
    float m_lower_theta_bound;
    
    /*
     * Auxiliary data structures for sampling stage.
     */
    int * aux;
    int * aux_v;
    
    /*
     * Auxiliary data structures for learning stage.
     */
    float * VjsMean;
	int * VjsNonMean;
	int * Vjs;
	int * invertedB;
	int * composition;
    
    
	/*
     * Calculates de consensus permutation from the given population cases.
     */
	void CalculateConsensusRanking(int** samples, int sample_size, int* consensus_ranking);
	
    /*
     * Calculates the spread theta parameter from the ConsensusRanking and the individuals in the population.
     */
    void CalculateThetaParameters(int*consensus_ranking, int** samples, int samples_num, float * theta_parameters, float upper_bound);
    
    /*
     * Calculates the total Psi normalization constant from the ThetaParameter and psi-s vector.
     */
    void CalculatePsiConstants(float* thetas, float* psi);
    
	/*
     * Given the consensus ranking, it samples a new individual given the s
     */
    inline void Sample(int * permutation);
    
    /*
     * Generates a permutation from a given vector of Vs.
     */
    void GeneratePermuFromV(int*v,int*permu,int n);
    
    /*******************************  NEWTON RAPSHON ALGORITHM  ********************************/
	
	/*
	 * Estimates theta parameters by means of running newton iterative algorithm.
	 */
	float NewtonRaphsonMethod(float initialGuess, float* Vjs, int j, float upper_bound);
	
	/*
	 * Theta parameter estimation function.
	 */
	inline float f(float theta, int n, float* VjMeans, int j);
	
	/*
	 * Theta parameter estimation function derivation.
	 */
	inline float fdev(float theta, int n, int j);
	
	/*
	 * Calculates Newton algorithms f and fdev functions.
	 */
	void funcd(float theta, float*ff, float*ffdev, float* VjMeans, int j);
	
	/*
	 * Newton - Rapshon execution algorithm.
	 */
	float rtsafe(float x1, float x2, float xacc, float* VjMeans, int j);
	
};
#endif
