/*
 *  Tools.h
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include "Variables.h"
using std::istream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;



void PrintInFile(int* array, int length, string text, ofstream& file);

/*
 * Determines if a given string contains a certain substring.
 */
bool strContains(const string inputStr, const string searchStr);

/*
 * Parses the string to a char array.
 */
char* StringToChar(string value);

/*
 * Prints in standard output 'length' integer elements of a given array.
 */
void PrintArray(int* array, int length, string text);

/*
 * Prints in standard output 'length' float elements of a given array.
 */
void PrintArray(float* array, int length, string text);

/*
 * Prints in the standard output the 'size' elements of the given array..
 */
void PrintArray(int * array, int size);

/*
 * Prints in the standard output given matrix.
 */
void PrintMatrix(int** matrix, int length, int length2, string text);

/*
 * Prints in standard output 'length' double elements of a given array.
 */
void PrintArrayDouble(double* array, int length, string text);


/*
 * Prints in standard output 'lengthxlength' double elements of a given matrix.
 */
void PrintMatrixDouble(double** matrix, int length, int length2, string text);

/*
 * Prints the 'length' double elements of the given array in the ofstream 'file'.
 */
void PrintArrayDoubleFile(double* array, int length, string text,ofstream file);

/*
 * Calculates the average of first num elements of the ints array given.
 */
double Mean(int* array, int num);

/*
 * Calculates the variance of first num elements of the int array given.
 */
double Variance(int* array, int num);

/*
 * Calculates the average of first num elements of the doubles array given.
 */
double Mean(double* array, int num);

/*
 * Calculates the variance of first num elements of the doubles array given.
 */
double Variance(double* array, int num);

/*
 * Returns the transpose given matrix.
 */
double** TransposeMatrix(double** matrix, int size);

/*
 * Copies the first 'size' elements of the given array 'a' to the array 'result'
 */
void CloneArray(int * a, int * result, int size);

/*
 * Inserts a given value in a position of the array an the values after that position are displaced 1 position.
 */

void InsertAt(int*array, int element, int at, int n);


/*
 * Inserts a given value in a position of the array an the values after that position are displaced 1 position.
 */
void InsertAt2 (int * array, int i, int j, int n);

/*
 * Removes a position in a given array of doubles and swap the solutions on the left.
 */

void RemoveAt(double*array,int positionIndex, int n);

/*
 * Removes a position in a given array of ints and swap the solutions on the left.
 */
void RemoveAt(int*array,int positionIndex, int n);

/*
 * Returns the max value of the given array.
 */
int Max(int*array, int n);

/*
 * Returns the min value of the given array.
 */

int Min(int*array, int n);

/*
 * Returns the max value's position in the given array of ints.
 */
int MaxPosition(int*array, int n);

/*
 * Returns the max value's position in the given array of doubles.
 */
int MaxPositionDouble(double*array, int n);

/*
 * Sorts the doubles array in the descending sequence.
 */
void DescendingSequence(double * array, int* descendingSequence, int size);

/*
 * Sorts the doubles array in the ascending sequence.
 */
void AscendingSequence(int * array, int* indexes,int* resultingSequence, int size);
