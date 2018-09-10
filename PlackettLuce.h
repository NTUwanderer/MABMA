/*
 *  PlackettLuce.h
 *  MixtureMallowsEDA
 *
 *  Created by Josu Ceberio Uribe on 7/16/12.
 *  Copyright 2012 University of the Basque Country. All rights reserved.
 *
 */


#ifndef PLACKETTLUCE_H__
#define PLACKETTLUCE_H__

#include <list>
#include <vector>
#include "RankingModel.h"
#include "Population.h"
#include "PlackettLuceInference.h"
#include "Individual.h"

class CPlackettLuceModel : public CRankingModel
{
public:
    
//    CPlackettLuceInference * m_pl_inference;
    
	/*
	 * The constructor.
	 */
	CPlackettLuceModel();
    
    CPlackettLuceModel(int, int);
	
    /*
	 * The destructor.
	 */
    virtual ~CPlackettLuceModel();
    
     /*
	 * Initialler
	 */   
    void init(int, int);
	
    /*
     * Given a population of samples, it learns a Plackett-Luce model from the data.
     */
    bool Learn(CPopulation * population);
    
	/*
	 * Simulates a new individual sampling the probabilistic model.
	 */
    void Sample(int * permutation);
    
	/*
     * Calculates ThetaParameters average value.
     */
	double GetWeightsAverage();
	
	/*
	 * Calculates the probability of the individual given the probabilistic model.
	 */
	double Probability(int * individual, bool straight);
    
    /*
     * This method updates the best solution of the optimization processes if the most probable solutions
     * are better than those solutions obtained by the population.
     */
    void PostProcesses();
    
private:
    
    /*
     * It applies random shake_power times a perturbation over the given individual.
     */
    void ShakeByInsert(int * permutation, int shake_power);
    
    /*
	 * It estimates the Plackett-Luce model from the given cases with the MM algorithm.
	 */
	bool LearnMM(CPopulation * population, int sel_total);

    bool LearnMMInverse(CPopulation * population, int sel_total);

    void LearnBayesianInference(int ** data,int sel_total);

    
    void SampleInverse(int * permutation);
    
    void SampleWithThreshold(int * permutation);
    
    void SampleInverseWithThreshold(int * permutation);

    /*
     * Samples the plackett-luce model learnt by shaking the most probable solution.
     */
    void SampleByShake(int * permutation, int * most_probable);

    
    CIndividual * SimulateInverse();
    
    /*
     * Calculates the most probable solution given the vector of weights
     */
    void CalculateMostProbableSolution();
    
    /*
     * Calculates the most probable solution given the inverse vector of weights
     */
    void CalculateMostProbableInverseSolution();
    
    /*
     * Returns the most probable solution given the vector of weights
     */
    CIndividual * GetMostProbableSolution();
	
    
    /*
     * Returns the most probable solution given the vector of weights
     */
    CIndividual * GetMostProbableInverseSolution();
	
	
	/*
     * Problem size.
     */
    int m_problem_size;

    /*
     * Number of samples.
     */
    int m_samples_size;
    
    /*
     * Number of maximum iterations to perform by the MM algorithm.
     */
    int m_MM_iterations;
    
    /*
     * Denotes if the sampling is carried out with threshold fixing.
     */
    bool m_twomodels_learning;
    
    /*
     * Denotes if the sampling is carried out with threshold fixing.
     */
    bool m_sampling_by_threshold;
    
    /*
     * Probability threshold to fix an item when sampling.
     */
    float m_sampling_threshold;
    
	/*
     * Parameter vector v
     */
	float * m_weights;
    
    /*
     * Inverse parameter vector v;
     */
    float * m_inverse_weights;
	
    /* 
     * Euclidean norm of the weights.
     */
    double m_norm;
    
    /*
     * Most probable individual according to the straight model.
     */
    CIndividual * m_mostProbable;
    
    /*
     * Most probable individual according to the inverse model.
     */
    CIndividual * m_inverse_mostProbable;
    
    /*
     * Auxiliary data.
     */
    int * m_aux1;
    int * m_aux2;
    
    // Number of items fixed in the template.
    int m_k;
    
    //Parameters for learning
    int m_M;
	int m_N; //number of rankings in the population.
	int m_P;
    int ** f;
    int ** r;
    float ** r2;
    float ** g;
    float * dgamma;
	float * gamma;
    float * w;
	int * pp;
    float **aux;
    float * newgamma;
    float * sr2;
    int * m_last;
    int m_estimation_type;
    
    //parameters for sampling
    float * distribution;
	int * objects;
    
    //forward-backward sampling arrays
    int * m_located;
    int * m_p1;
    int * m_p2;

};
#endif
