/***************************************************************************
                          myrand.cc  -  description
                             -------------------
    begin                : Sep 24 2001
    copyright            : (C) 2001 by Tian-Li Yu
    email                : tianliyu@cc.ee.ntu.edu.tw
 ***************************************************************************/

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "myrand.h"

#ifdef PI
#undef PI
#endif

#define PI 3.14159265
#define N 624

const time_t MyRand::timeseed = time(NULL);
bool MyRand::alreadySranded = false;

MyRand::MyRand ()
{

    unsigned long init_key[N];
	
	if(alreadySranded == false)
	{
		srand ((unsigned long) timeseed);
		alreadySranded = true;
	}

    for (int i = 0; i < N; i++)
        init_key[i] = (unsigned long) timeseed * rand ();

    init_by_array (init_key, N);
}


MyRand::~MyRand ()
{
}


bool MyRand::flip (double prob)
{
    return (uniform () < prob);
}


/** From [0,1) */
double
MyRand::uniform ()
{
    return genrand_res53 ();
}


/** From [a,b) */
double
MyRand::uniform (double a, double b)
{
    return uniform () * (b - a) + a;
}


/** Int From [a,b] */
int
MyRand::uniformInt (int a, int b)
{
    return (a + (int) (uniform () * (b - a + 1)));
}


/** Generate a random array of size num, from [a,b] */
/** num <= b-a+1 */
void
MyRand::uniformArray (int *array, int num, int a, int b)
{
	int l = b - a + 1;
    int *base = new int[l];
    int i;
    int r;
	int j;

    // for (i = 0; i < l; i++)
        // base[i] = a + i;

    for (i = 0; i < num; i++) {
		if(i%l == 0)
			for (j = 0; j < l; j++)
				base[j] = a + j;	
        r = uniformInt (0, b - a - (i%l));
        array[i] = base[r];
        base[r] = base[b - a - (i%l)];
		
    }

    delete[]base;
}
