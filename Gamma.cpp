//
//  Gamma.cpp
//  PlackettLuceEDA
//
//  Created by Josu Ceberio Uribe on 10/16/12.
//  Copyright (c) 2012 Josu Ceberio Uribe. All rights reserved.
//

#include "Gamma.h"
#include <limits>
using namespace std;
/*
 *Class constructor.
 */
Gamma::Gamma()
{
	Shape = 1;
    Rate = 0;
}

/*
 *Class constructor.
 */
Gamma::Gamma(double shape, double rate)
{
	Shape = shape;
    Rate = rate;

}


/*
 * Class destructor.
 */
Gamma::~Gamma()
{
	
}


void Gamma::SetShapeAndRate(double Shp, double Rat)
{
    Shape=Shp;
    Rate=Rat;
}

/*
 * Asks whether the instance is a point mass
 */
bool Gamma::IsPointMass()
{
    return (Shape == numeric_limits<float>::infinity( ));
}

/*
 * Gets the expected value E(x) - calculated as shape/rate
 */
float Gamma::GetMean()
{
    if (IsPointMass())
    {
        Point=Rate;
        return Point;
    }
    else if (Rate == 0.0)
    {
        return numeric_limits<float>::infinity( );
   
    }
    else
    {
        if (Shape/Rate<0) {
            //cout<<"neg"<<endl;
            //Print2();
            //exit(1);
            //return 0;
        }
        return Shape / Rate;
    }
}

/*
 * Sets the parameters to represent the product of two Gammas.
 */
void Gamma::SetToProduct(Gamma * a, Gamma * b)
{
   // cout<<"sartu da"<<endl;
    if (a->IsPointMass())
    {
        //cout<<"hemen2"<<endl;
        if (b->IsPointMass() && a->Point!=b->Point) {
            //throw new AllZeroException();
        }
        Point = a->Point;
        return;
    }
    else if (b->IsPointMass())
    {
      //  cout<<"hemen3"<<endl;
        Point = b->Point;
        return;
    }
    else
    {
       // cout<<"hemen"<<endl;
        Shape = a->Shape + b->Shape - 1;
        Rate = a->Rate + b->Rate;
    }
}

/* 
 * Sets the parameters to represent the ratio of two Gammas.
 */
void Gamma::SetToRatio(Gamma * numerator, Gamma * denominator)
{
    if (numerator->IsPointMass())
    {
        if (denominator->IsPointMass()) {
            if (numerator->Point==denominator->Point)
            {
                //SetToUniform();
                Shape=1;
                Rate=0;
            }
            else
            {
        //        throw new DivideByZeroException();
            }
        } else {
            Point = numerator->Point;
        }
    }
    else if (denominator->IsPointMass())
    {
        // throw new DivideByZeroException();
    }
    else
    {
        Shape = numerator->Shape - denominator->Shape + 1;
        Rate = numerator->Rate - denominator->Rate;
    }

}

void Gamma::Print2()
{
    cout<<"Shape: "<<Shape;
    cout<<"     Rate: "<<Rate<<endl;
    
}

void Gamma::Copy(Gamma * gam)
{
    Shape=gam->Shape;
    Rate=gam->Rate;
    Point=gam->Point;
}

void Gamma::InitializeByExpectationAndVariance(float expectation,float variance)
{
    Rate=expectation/variance;
    Shape=Rate*Rate*variance;
}
