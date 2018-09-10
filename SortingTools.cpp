/*
 *  SortingTools.cpp
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 10/5/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "SortingTools.h"
#include "Variables.h"

//using std::cerr;
//using std::cout;
//using std::endl;
int Partition(int * a, double * criteriaValues, int left, int right, int pivotIndex, SORTING_CRITERIA criteria)
{
	if (criteria==ASCENDING)
	{
		int aux;
		double pivot =criteriaValues[a[pivotIndex]];
		do
		{
			while (criteriaValues[a[left]] < pivot) left++;
			while (criteriaValues[a[right]] > pivot) right--;
			if (left < right && criteriaValues[a[left]] != criteriaValues[a[right]])
			{
				//swap
				aux=a[left];
				a[left]=a[right];
				a[right]=aux;
			}
			else
			{
				return right;
			}
		}
		while (left < right);
	}
	else 
	{
		int aux;
		double pivot =criteriaValues[a[pivotIndex]];
		do
		{
			while (criteriaValues[a[left]] > pivot) left++;
			while (criteriaValues[a[right]] < pivot) right--;
			if (left < right && criteriaValues[a[left]] != criteriaValues[a[right]])
			{
				//swap
				aux=a[left];
				a[left]=a[right];
				a[right]=aux;
			}
			else
			{
				return right;
			}
		}
		while (left < right);
	}

    return right;
}

void QuickSort(int * a, double * criteriaValues, int left, int right, SORTING_CRITERIA criteria)
{
	//WARNING!!! THIS CODE DOES NOT WORK CORRENTLY WHEN CRITERIA VALUES ARE DUPLICATED.
	
//	cout<<"----------------------------------------------"<<endl;
//	PrintArray(a,20,"Permu: ");
//	cout<<"Left index: "<<left<<" Right index: "<<right<<" Left pos: "<<a[left]<<" Right pos: "<<a[right]<<endl;
	//PrintArrayDouble(criteriaValues,20,"Criteria: ");
	if (left < right)
	{
		int pivot = (left + right) / 2; // middle
		int pivotNew = Partition(a, criteriaValues, left, right, pivot, criteria);
		QuickSort(a, criteriaValues, left, pivotNew - 1, criteria);
		QuickSort(a, criteriaValues, pivotNew + 1, right, criteria);
	}
//	cout<<"----------------------------------------------"<<endl;
}

void RandomKeys( int * a, float * criteriaValues, int ProblemSize)
{
	bool fixedValues[ProblemSize];
	float criteria, min;
	int i, j;
	for (i=0;i<ProblemSize;i++){
		fixedValues[i]=false;
		a[i]=0;
	}
	int minPos=0;
	for (i=0;i<ProblemSize;i++)
	{
		min=MAX_INTEGER;
		for (j=0;j<ProblemSize;j++)
		{
			criteria=criteriaValues[j];
			//cout<<"j: "<<j<<" min: "<<min<<" criteriaValue: "<<criteria<<" fixedValue: "<<fixedValues[j]<<endl;
			if (!fixedValues[j] && min>criteria )
			{
				min=criteria;
				minPos=j;
			//	cout<<"********* new min: "<<min<<" position: "<<j<<endl;
			}
		}

		fixedValues[minPos]=true;
		a[minPos]=i;
		//cout<<"++++++++ fixed pos: "<<minPos<< " job assigned: "<<i<<endl;
		//cout<<"------"<<endl;
	}
}

void RandomKeysDescending( int * a, float * criteriaValues, int ProblemSize)
{
	bool fixedValues[ProblemSize];
	double criteria, max;
	int i, j;
	for (i=0;i<ProblemSize;i++){
		fixedValues[i]=false;
		a[i]=0;
	}
	int maxPos=0;
	for (i=0;i<ProblemSize;i++)
	{
		max=MIN_INTEGER;
		for (j=0;j<ProblemSize;j++)
		{
			criteria=criteriaValues[j];
			//cout<<"j: "<<j<<" min: "<<min<<" criteriaValue: "<<criteria<<" fixedValue: "<<fixedValues[j]<<endl;
			if (!fixedValues[j] && max<criteria )
			{
				max=criteria;
				maxPos=j;
                //	cout<<"********* new min: "<<min<<" position: "<<j<<endl;
			}
		}
        
		fixedValues[maxPos]=true;
		a[maxPos]=i;
		//cout<<"++++++++ fixed pos: "<<minPos<< " job assigned: "<<i<<endl;
		//cout<<"------"<<endl;
	}
}
