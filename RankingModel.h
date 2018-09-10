//
//  RankingModel.h
//  RankingEDA
//
//  Created by Josu Ceberio Uribe on 1/23/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDA__RankingModel__
#define __RankingEDA__RankingModel__

#include <iostream>
#include "Population.h"
#include "Individual.h"


class CRankingModel
{
public:
    
	/*
	 * The constructor.
	 */
	CRankingModel();
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the kendall model.
	 */
	virtual ~CRankingModel();
	
	/*
	 * Virtual learning function.
	 */
	virtual bool Learn(CPopulation * population);
    
    /*
	 * Virtual sampling function.
	 */
	virtual void Sample(int * genes);
    
    /*
     * Virtual post-processes function.
     */
    virtual void PostProcesses(class CPopulation * population);
    
private:
};

#endif /* defined(__RankingEDA__RankingModel__) */
