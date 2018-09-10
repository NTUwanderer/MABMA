//
//  PlackettLuceInference.h
//  PlackettLuceEDA
//
//  Created by Josu Ceberio Uribe on 10/16/12.
//  Copyright (c) 2012 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __PlackettLuceEDA__PlackettLuceInference__
#define __PlackettLuceEDA__PlackettLuceInference__

#include <iostream>
#include "Gamma.h"

class CPlackettLuceInference
{
public:
    
    /*
	 * The constructor.
	 */
    CPlackettLuceInference(int problem_size, int sel_size, int iterations);
	
    /*
	 * The destructor.
	 */
    virtual ~CPlackettLuceInference();

    double InferPL_Truncated(int ** data, float * weights, float * variances, bool standardInitialization);

    void update_factor (int * w, Gamma * q, Gamma * messages_from_factor, int k, double r);
    
    void NormalizeMarginals(Gamma*  marginals);
    
    float log_evidence_scale (int * w, Gamma * q, Gamma * messages_from_factor, int k);
    
private:
    
    //Variables of the method update factor.
    double * m_rates;
    double * m_shapes;
    double * m_means;
    double * m_a;
    double * m_c;
    double * m_d;
    
    //Inference method variables    
    double * m_log_qmass;
    double * m_aux_means;
    int m_K;
    int m_N;
    int m_topf;
    Gamma * m_marginals;
    //vector<vector<vector<Gamma***>>> m_messages_fact2var;
    Gamma *** m_messages_fact2var;
    
    //General parameters
    int m_iterations;

};
#endif /* defined(__PlackettLuceEDA__PlackettLuceInference__) */
