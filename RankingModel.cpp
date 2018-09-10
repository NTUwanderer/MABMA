//
//  RankingModel.cpp
//  RankingEDA
//
//  Created by Josu Ceberio Uribe on 1/23/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "RankingModel.h"
#include "Population.h"

/*
 * The constructor.
 */
CRankingModel::CRankingModel()
{
	
}

/*
 * The destructor. It frees the memory allocated at the construction of the kendall model.
 */
CRankingModel::~CRankingModel()
{
	//cout<<"Deleting distance model"<<endl;
}

/*
 * Virtual learning function.
 */
bool CRankingModel::Learn(CPopulation * population)
{
    return true;
}


/*
 * Virtual sampling function.
 */
void CRankingModel::Sample(int * genes)
{
    
}


/*
 * Virtual post-processes function.
 */
void CRankingModel::PostProcesses(CPopulation * population)
{
	
}