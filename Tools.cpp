/*
 *  Tools.cpp
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "Tools.h"


/*
 * Prints in a given file 'length' integer elements of a given array.
 */
void PrintInFile(int* array, int length, string text, ofstream& file)
{
	file<<text;
	for (int i=0;i<length;i++){
		file<<array[i]<<" ";
	}
	file<<" "<<endl;
}

/*
 * Determines if a given string contains a certain substring.
 */
bool strContains(const string inputStr, const string searchStr)
{
	size_t contains;
	
	contains = inputStr.find(searchStr);
	
	if(contains != string::npos)
		return true;
	else
		return false;
}

/*
 * Parses the string to a char array.
 */
char* StringToChar(string value)
{
	// stringstream ss;
	char* charValue = const_cast<char*>(value.c_str());
	// ss << value;
	// ss >> charValue;
	return charValue;
}


/*
 * Prints in standard output 'length' integer elements of a given array.
 */
void PrintArray(int* array, int length, string text)
{
	cout<<text;
	for (int i=0;i<length;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}


/*
 * Prints in the standard output the 'size' elements of the given array..
 */
void PrintArray(int * array, int size)
{
	for (int i=0;i<size;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}


/*
 * Prints in standard output 'length' double elements of a given array.
 */
void PrintArray(float* array, int length, string text)
{
	cout<<text;
	//printf("%s",text);
	//printf("\n ");
	for (int i=0;i<length;i++){
		//	cout<<array[i]<<" ";
		printf(" %3.4f,",array[i]);
	}
	printf("\n ");
	//cout<<" "<<endl;
}


/*
 * Prints in the standard output given matrix.
 */
void PrintMatrix(int** matrix, int length, int length2, string text)
{
	cout<<text;
	for (int i=0;i<length;i++)
	{
		cout<<""<<endl;
		for (int j=0;j<length2;j++)
		{
			cout<<matrix[i][j]<<" ";
		}
	}
	cout<<" "<<endl;
}

/*
 * Prints in standard output 'length' double elements of a given array.
 */
void PrintArrayDouble(double* array, int length, string text)
{
	cout<<text;
	//printf("%s",text);
	//printf("\n ");
	for (int i=0;i<length;i++){
		//	cout<<array[i]<<" ";
		printf(" %3.4f ",array[i]);
	}
	printf("\n ");
	//cout<<" "<<endl;
}

/*
 * Prints the 'length' double elements of the given array in the ofstream 'file'.
 */
void PrintArrayDoubleFile(double* array, int length, string text,ofstream file)
{
	file<<text;
	for (int i=0;i<length-1;i++){
		file<<array[i]<<";";
	}
	file<<array[length-1]<<endl;
}


/*
 * Prints in standard output 'lengthxlength' double elements of a given matrix.
 */
void PrintMatrixDouble(double** matrix, int length, int length2, string text)
{
	cout<<text<<endl;
	//printf(text);
	for (int i=0;i<length;i++)
	{
		
		//	cout<<""<<endl;
		for (int j=0;j<length2;j++)
		{
			//cout<<matrix[i][j]<<" ";
			printf("%3.4f ",matrix[i][j]);
		}
		printf("\n");
	}
	//cout<<" "<<endl;
}

/*
 * Calculates the average of first num elements of the ints array given.
 */
double Mean(int* array, int num)
{
    double result=0;
    for(int i=0;i<num;++i)
        result+=array[i];
    return result/num;
};

/*
 * Calculates the variance of first num elements of the int array given.
 */
double Variance(int* array, int num)
{
	double avg=Mean(array,num);
	double var=0;
	for (int i = 0; i < num; ++i )
	{
		var += ((array [i] - avg) * (array [i] - avg)) ;
	}
	var /= (num - 1) ;
	return var;
};

/*
 * Calculates the average of first num elements of the doubles array given.
 */
double Mean(double* array, int num)
{
    double result=0;
    for(int i=0;i<num;++i)
        result+=array[i];
    return result/num;
};

/*
 * Calculates the variance of first num elements of the doubles array given.
 */
double Variance(double* array, int num)
{
	double avg=Mean(array,num);
	double var=0;
	for (int i = 0; i < num; ++i )
	{
		var += ((array [i] - avg) * (array [i] - avg)) ;
	}
	var /= (num - 1) ;
	return var;
};

/*
 * Returns the transpose given matrix.
 */
double** TransposeMatrix(double** matrix, int size)
{
	double** transposed= new double*[size];
	for (int i=0;i<size;i++)
		transposed[i]= new double[size];
	
	for (int i=0;i<size;i++)
	{
		for (int j=0;j<size;j++)
		{
			transposed[i][j]=matrix[j][i];
		}
	}
	return transposed;
}

/*
 * Copies the first 'size' elements of the given array 'a' to the array 'result'
 */
void CloneArray(int * a, int * result, int size)
{
	for (int i=0;i<size;i++) result[i]=a[i];
}


/*
 * Inserts a given value in a position of the array an the values after that position are displaced 1 position.
 */
void InsertAt(int*array,int element, int at, int n)
{
	
	for (int i=n-2;i>=at;i--)
	{
		array[i+1]=array[i];
	}
	array[at]=element;
}

/*
 * This method moves the value in position i to the position j.
 */
void InsertAt2 (int * array, int i, int j, int n)
{
	if (i!=j)
	{
		int res[n];
		int val=array[i];
		if (i<j)
		{
			memcpy(res,array,sizeof(int)*i);
			
			for (int k=i+1;k<=j;k++)
				res[k-1]=array[k];
			
			res[j]=val;
			
			for (int k=j+1;k<n;k++)
				res[k]=array[k];
		}
		else if (i>j)
		{
			memcpy(res,array,sizeof(int)*j);
			
			res[j]=val;
			
			for (int k=j;k<i;k++)
				res[k+1]=array[k];
			
			for (int k=i+1;k<n;k++)
				res[k]=array[k];
		}
		memcpy(array,res,sizeof(int)*n);
	}
}


/*
 * Returns the max value of the given array.
 */
int Max(int*array, int n)
{
	int max_value=0;
	for(int i=0;i<n;i++)
	{
		if (array[i]>max_value)
			max_value=array[i];
	}
	return max_value;
}

/*
 * Returns the min value of the given array.
 */
int Min(int*array, int n)
{
	int min_value=MAX_INTEGER;
	for(int i=0;i<n;i++)
	{
		if (array[i]<min_value)
			min_value=array[i];
	}
	return min_value;
}

/*
 * Returns the max value's position in the given array of ints.
 */
int MaxPosition(int*array, int n)
{
	int max_value=MIN_INTEGER;
	int max_value_position=0;
	for(int i=0;i<n;i++)
	{
		if (array[i]>max_value)
		{
			max_value=array[i];
			max_value_position=i;
		}
	}
	return max_value_position;
}

/*
 * Returns the max value's position in the given array of doubles.
 */
int MaxPositionDouble(double*array, int n)
{
	int max_value=MIN_INTEGER;
	double max_value_position=0;
	for(int i=0;i<n;i++)
	{
		if (array[i]>max_value)
		{
			max_value=array[i];
			max_value_position=i;
		}
	}
	return max_value_position;
}

/*
 * Returns the min value's position in the given array.
 */
int MinPosition(int*array, int n)
{
	int min_value=MAX_INTEGER;
	int min_value_position=0;
	for(int i=0;i<n;i++)
	{
		if (array[i]<min_value)
		{
			min_value=array[i];
			min_value_position=i;
		}
	}
	return min_value_position;
}

/*
 * Sorts the doubles array in the descending sequence.
 */
void DescendingSequence(double * array, int* descendingSequence, int size)
{
	for (int i=0;i<size;i++)
	{
		int pos=MaxPositionDouble(array, size);
		descendingSequence[i]=pos;
		array[pos]=0;
	}
}

/*
 * Sorts the doubles array in the ascending sequence.
 */
void AscendingSequence(int * array, int* indexes, int* resultingSequence, int size)
{
	int * aux= new int[size];
	//memcpy(aux, array, sizeof(int)*size);
	for (int i=0;i<size;i++)
	{
		int pos=MinPosition(array, size);
		resultingSequence[i]=indexes[pos];
		aux[i]=array[pos];
		array[pos]=MAX_INTEGER;
	
	}
	memcpy(array,aux, sizeof(int)*size);
	delete[] aux;
}

/*
 * Removes a position in a given array of doubles and swap the solutions on the left.
 */
void RemoveAt(double*array,int positionIndex, int n)
{
	for (int i=positionIndex;i<n-1;i++)
	{
		array[i]=array[i+1];
	}
//	array[n-1]=0;//removable. It might be wrong.
}


/*
 * Removes a position in a given array of ints and swap the solutions on the left.
 */
void RemoveAt(int*array,int positionIndex, int n)
{
	for (int i=positionIndex;i<n-1;i++)
	{
		array[i]=array[i+1];
	}
	array[n-1]=0;
}
