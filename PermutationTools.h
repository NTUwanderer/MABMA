/*
 *  Tools.h
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 9/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */
#include <string.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include "Population.h"
#include <stdlib.h>
#include <stdio.h>

using std::stringstream;
using std::string;


/*
 * It determines if the given int sequecen if it is indeed a permutation or not.
 */
bool isPermutation(int * permutation, int size);

/*
 *	Calculates V_j-s vector.
 */
void vVector(int*v, int*permutation, int n);

/*
 * Inverts a permutation.
 */
void Invert(int*permu, int n, int* inverted);

/*
 * Implements the compose of 2 permutations of size n.
 */
void Compose(int*s1, int*s2, int*res, int n);

/*
 * Generates a permutation from a v vector.
 */
void GeneratePermuFromV(int*v,int*permu,int n);

/*
 * Generates a random permutation of size 'n' in the given array.
 */
void GenerateRandomPermutation(int size, int * permutation);

int * vVector2Permutations(int* permutationA, int*permutationB, int size);

/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int* permutationA, int*permutationB, int size);
