//
//  MallowsModel.h
//  RankingEDA
//
//  Created by Josu Ceberio Uribe on 1/23/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDA__MallowsModel__
#define __RankingEDA__MallowsModel__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>


class CMallowsModel : public CRankingModel
{
	
public:
	
    /*
     * The constructor.
     */
	CMallowsModel();
	
    /*
     * The destructor.
     */
    virtual ~CMallowsModel();
	
    /*
     * Given a population of samples, it learns a Generalized Mallows model from the data.
     */
    bool Learn(CPopulation * population);
    
    /*
     * Creates a new individual sampling the generalized mallows probabilistic model.
     */
    CIndividual * Simulate();
    
	/*
	 * Calculates the probability of the individual given the probabilistic model.
	 */
	float Probability(int * individual);
	
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
    float m_theta_parameter;
    
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
    float CalculateThetaParameter(int*consensus_ranking, int** samples, int samples_num, float upper_bound);
    
    /*
     * Calculates the total Psi normalization constant from the ThetaParameter and psi-s vector.
     */
    void CalculatePsiConstants(float theta, float* psis);
    
	/*
     * Given the consensus ranking, it samples a new individual given the s
     */
    void Sample(int * permutation);
	
	/*
	 * Generates a permutation from a v vector.
	 */
	void GeneratePermuFromV(int * v, int * permu, int n);
	
	
	/*******************************  NEWTON-RAPSHON ALGORITHM for THETAS ESTIMATION ********************************/
	
	/*
	 * Estimates theta parameters by means of running newton iterative algorithm.
	 */
	float NewtonRaphsonMethod(float initialGuess, float* Vjs, float upper_bound);
	
	/*
	 * Theta parameter estimation function.
	 */
	inline float f(float theta, int n, float* VjMeans);
	
	/*
	 * Theta parameter estimation function derivation.
	 */
	inline float fdev(float theta, int n);
	
	/*
	 * Calculates Newton algorithms f and fdev functions.
	 */
	void funcd(float theta, float*ff, float*ffdev, float* VjMeans);
	
	/*
	 * Newton - Rapshon execution algorithm.
	 */
	float rtsafe(float x1, float x2, float xacc, float* VjMeans);
	
};

#endif /* defined(__RankingEDA__MallowsModel__) */
