/*
 *  Tools.cpp
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 9/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */
#include "PermutationTools.h"
#include "Tools.h"
#include "SortingTools.h"
#include "Problem.h"
#include <string.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using std::stringstream;
using std::string;
using std::cerr;
using std::cout;
using std::endl;


/*
 * It determines if the given int sequecen if it is indeed a permutation or not.
 */
bool isPermutation(int * permutation, int size)
{
	int * flags=new int[size];
	for (int i=0;i<size;i++) flags[i]=1;
	
	for (int i=0;i<size;i++)
	{
		int value=permutation[i];
		flags[value]=0;
	}
	
	int result,sum=0;
	for(int i=0;i<size;i++)
		sum+=flags[i];
	if (sum==0) result=true;
	else result=false;
	delete [] flags;
	return result;
}

/*
 * Generates a permutation from a v vector.
 */
void GeneratePermuFromV(int*v,int*permu,int n)
{
	int at;
	permu[0]=n-1;
	for (int j=n-2;j>=0;j--)
	{
		at=v[j];
		InsertAt(permu,j,at,n);
	}
}

/*
 *	Calculates V_j-s vector.
 */
void vVector(int*v, int*permutation, int n)
{
	for(int i=0;i<n-1;i++)
		v[i]=0;
	
	for (int i = n-2; i >= 0; i--)
		for (int j = i+1; j < n; j++)
			if(permutation[i] > permutation[j])
				v[i]++;
}

/*
 * Inverts a permutation.
 */
void Invert(int*permu, int n, int* inverted)
{
	for(int i=0; i<n; i++) inverted[permu[i]]=i;
}

/*
 * Implements the compose of 2 permutations of size n.
 */
void Compose(int*s1, int*s2, int*res, int n)
{
	for(int i=0;i<n;i++)
		res[i]=s1[s2[i]];
}

/*
 * Generates a random permutation of size 'n' in the given array.
 */
void GenerateRandomPermutation(int n, int * permutation) 
{
	for (int i = 0; i < n; ++i) 
	{
		int j = rand() % (i + 1);
		permutation[i] = permutation[j];
		permutation[j] = i;
	}
}



/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int* permutationA, int*permutationB, int size)
{
	int distance = 0;
	int * v= new int[size-1];
	int * composition= new int[size-1];
	int * invertedB= new int[size];
	
	Invert(permutationB, size, invertedB);	
	Compose(permutationA,invertedB,composition,size); 
	vVector(v,composition,size);
	for (int i = 0; i < size-1; i++) distance += v[i];
	
	delete[] composition;
	delete[] invertedB;
	delete[] v;
	return distance;
}

int * vVector2Permutations(int* permutationA, int*permutationB, int size)
{
	int * v= new int[size-1];
	int * composition= new int[size-1];
	int * invertedB= new int[size];
	
	Invert(permutationB, size, invertedB);	
	Compose(permutationA,invertedB,composition,size); 
	vVector(v,composition,size);

	delete [] composition;
	delete [] invertedB;
	
	return v;
}
