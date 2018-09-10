/*
 *  PBP.h
 *  MixtureMallowsEDA
 *
 *  Created by Josu Ceberio Uribe on 6/25/12.
 *  Copyright 2012 University of the Basque Country. All rights reserved.
 *
 */
#ifndef __PBP_H__
#define __PBP_H__

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdio.h>

using namespace std;
using std::stringstream;
using std::string;

class PBP
{
public:
		
	/*
	 * The constructor.
	 */
	PBP();
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the kendall model.
	 */
	virtual ~PBP();
	
	/*
	 * Virtual evaluation function.
	 */
	virtual int Evaluate(int*genes)=0;

    /*
	 * Virtual instance reading function.
	 */
	virtual int Read(string filename)=0;
    
private:
};
#endif
