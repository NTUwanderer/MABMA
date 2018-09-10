// Individual.cpp: implementation of the CIndividual class.
//

#include "Individual.h"
#include "Problem.h"
#include "EDA.h"
#include "Variables.h"
#include <limits.h>

// The object storing the values of all the individuals
// that have been created during the execution of the
// program. 

CIndividual::CIndividual(int length): m_value(MIN_INTEGER)
{
	m_size = length;
	m_genes = new int[m_size];
	aux_genes=new int[m_size];
}

CIndividual::~CIndividual()
{
	delete [] m_genes;
	delete [] aux_genes;
	m_genes=NULL;
}

int CIndividual::Value()
{
	// The INT_MIN value indicates that the value of
	// the individual has not been calculated yet. In
	// that case it is calculated and stored before
	// it is returned.
	if(m_value==MIN_INTEGER)
	{
        
        if (MODEL_TYPE ==0 || MODEL_TYPE==1 ) //Mallows or Generalized Mallows model
        {
            for(int i=0; i<m_size; i++) aux_genes[m_genes[i]]=i;
		
            //evaluate the solution
            m_value=Metric(aux_genes, m_size);
        }
        else
        {
            m_value=Metric(m_genes, m_size);
        }
	}
	return m_value;
}


int * CIndividual::Genes()
{
	return m_genes;
}

ostream & operator<<(ostream & os,CIndividual * & individual)
{
	// os << individual->Value() << ",";

        os << individual->m_genes[0];
	for(int i=1;i<individual->m_size;i++)
		os << "," << individual->m_genes[i];
        os << ".";

	return os;
}

istream & operator>>(istream & is,CIndividual * & individual)
{
  char k; //to avoid intermediate characters such as ,.

  is >> individual->m_value;
  is >> individual->m_genes[0];
  for(int i=1;i<individual->m_size;i++)
    is >> k >> individual->m_genes[i];
  is >> k;

  return is;
}

void CIndividual::SetValue(int metrica)
{
	m_value=metrica;
}

void CIndividual::SetGenes(int * genes)
{
	memcpy(m_genes, genes, sizeof(int)*m_size);
	m_value=MIN_INTEGER;
}

void CIndividual::SetGenes(int* genes, int fitness)
{
    memcpy(m_genes, genes, sizeof(int)*m_size);
    m_value = fitness;
}

/*
 * Prints in the standard output genes.
 */
void CIndividual::PrintGenes()
{
	for (int i=0;i<m_size;i++){
		cout<<m_genes[i]<<" ";
	}
	cout<<" "<<endl;
}

CIndividual * CIndividual::Clone()
{
	CIndividual * ind = new CIndividual(m_size);
	ind->SetGenes(m_genes);
	ind->SetValue(m_value);
	return ind;
}



