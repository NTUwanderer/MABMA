// Individual.h: interface for the CIndibiduoa class.
//
// Description:
//		The class which handles the individual's information:
//		its genes and its value.
//
// Author:	Ramon Etxeberria
// Date:	1999-09-25
//
// Note:
//		The problem's information must be loaded before any
//		individual is created.


#ifndef _INDIVIDUAL_
#define _INDIVIDUAL_

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>

using namespace std;
using std::istream;
using std::ostream;

class CIndividual
{
public:

	// The constructor. The constructed individual has
	// all zeroes as its genes.
	CIndividual(int length);

	// The destructor. It frees the memory allocated at
	// the construction of the individual.
	virtual ~CIndividual();

	// It returns an array with the genes of the individual.
	int * Genes();

	// It returns the value of the individual.
	int Value();
	
	void SetValue(int metrica);
	void SetGenes(int * genes);
	void SetGenes(int* genes, int fitness);

	// Prints in standard output 'length' integer elements of a given array.
	void PrintGenes();
	
	// Clone the individual.
	CIndividual * Clone();

	// Input and output of an individual: its value and genes.
	friend ostream & operator<<(ostream & os,CIndividual * & individual);
	friend istream & operator>>(istream & is,CIndividual * & individual);

       
private:

	// The variable storing the individual's value. 
	int m_value;

	 // The genes of the individual.
	int * m_genes;
	
	//The size of the individual
	int m_size;
	
	//Auxiliray array.
	int * aux_genes;

};

#endif
