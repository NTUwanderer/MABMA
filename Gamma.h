//
//  Gamma.h
//  PlackettLuceEDA
//
//  Created by Josu Ceberio Uribe on 10/16/12.
//  Copyright (c) 2012 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __PlackettLuceEDA__Gamma__
#define __PlackettLuceEDA__Gamma__

#include <iostream>

class Gamma
{
	
public:
	
    /// Rate parameter for the distribution
    double Rate;
    
    /// Shape parameter for the distribution
    double Shape;
    
    double Point;
    
    // Class constructor
    Gamma();
    
    //Class constructor.
    Gamma(double shape, double rate);
	
    // The destructor.
    virtual ~Gamma();
    
    void SetShapeAndRate(double Shp, double Rat);
    
    /*
     * Asks whether the instance is a point mass
     */
    bool IsPointMass();
    
    /*
     * Gets the expected value E(x) - calculated as shape/rate
     */
    float GetMean();
	   
    /*
     * Sets this Gamma instance to be a uniform distribution
     */
    void SetToRatio(Gamma * numerator, Gamma * denominator);
    
    /*
     * Sets the parameters to represent the product of two Gammas.
     */
    void SetToProduct(Gamma * a, Gamma * b);
    
    void Print2();
    
    void Copy(Gamma * gam);
    
    void InitializeByExpectationAndVariance(float expectation,float variance);
private:
	
};


#endif /* defined(__PlackettLuceEDA__Gamma__) */

