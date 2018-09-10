/*
 *  SortingTools.h
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 10/5/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

enum SORTING_CRITERIA {ASCENDING, DESCENDING};

int Partition(int * a, double * criteriaValues, int left, int right, int pivotIndex, SORTING_CRITERIA criteria);

void QuickSort(int * a, double * criteriaValues, int left, int right, SORTING_CRITERIA criteria);

void RandomKeys( int * a, double * criteriaValues, int ProblemSize);

void RandomKeys( int * a, float * criteriaValues, int ProblemSize);

void RandomKeysDescending( int * a, float * criteriaValues, int ProblemSize);
