/***************************************************************************
 *   Copyright (C) 2004 by Tian-Li Yu                                      *
 *   tianliyu@cc.ee.ntu.edu.tw                                             *
 ***************************************************************************/

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "global.h"
#include "chromosome.h"

Chromosome::Chromosome ()
{
    length = 0;
    gene.clear();
    evaluated = false;
    ispt = false;
    pertable = NULL;
}


Chromosome::Chromosome (int n_length)
{
    gene.clear();
    init (n_length, NULL);
}

Chromosome::Chromosome(int n_length, IEvaluator *eva)
{
    gene.clear();
    init(n_length, eva);
}


Chromosome::~Chromosome ()
{
    if (pertable != NULL)
        delete[]pertable;

}


void Chromosome::init (int n_length)
{
    length = n_length;

    gene.resize(n_length);
    evaluated = false;
    ispt = false;
    evaluator = NULL;
    pertable = new int[length];
    for (int i = 0; i < length; i++){
        pertable[i] = 0;
    }
}


void Chromosome::init (int n_length, IEvaluator *eva)
{
    length = n_length;


    gene.resize(n_length);
    evaluated = false;
    ispt = false;
    evaluator = eva;
    pertable = new int[length];
    for (int i = 0; i < length; i++){
        pertable[i] = 0;
    }
}

int Chromosome::getVal (int index) const
{
    if (index < 0 || index > length)
        outputErrMsg ("Index overrange in Chromosome::operator[]");

    return gene[index];
}


void Chromosome::setVal (int index, int val)
{
    if (index < 0 || index > length)
        outputErrMsg ("Index overrange in Chromosome::operator[]");

    gene.at(index) = val;
    evaluated = false;
    ispt = false;
}


double Chromosome::getFitness ()
{
    if (evaluated)
        return fitness;
    else
        return (fitness = evaluate ());
}


bool Chromosome::isEvaluated () const
{
    return evaluated;
}


double Chromosome::evaluate ()
{
    if(!evaluator)
        outputErrMsg ("Must have evaluator in Chromosome");

    evaluated = true;
    return evaluator->calc(this);

}


Chromosome & Chromosome::operator= (const Chromosome & c)
{
    this->clone(c);

    return *this;
}

void Chromosome::clone(const Chromosome& c)
{
    if (length != c.length) {
        length = c.length;
        this->gene.clear();
        init (length);
    }

    evaluated = c.evaluated;
    fitness = c.fitness;
    evaluator = c.evaluator;

    this->gene = c.gene;
}


void Chromosome::printf () const
{
    dump();
}


int Chromosome::getLength () const
{
    return length;
}


double Chromosome::getMaxFitness () const
{
    // For OneMax
    return ((double)length-1e-6);
}

void Chromosome::dump() const
{
    for(int i = 0; i < length; i++)
    {
        std::cout << std::setw(5) << gene[i];
    }
    std::cout << std::endl;
}

void Chromosome::setGenes(const int* t)
{
    if (length == 0)
        outputErrMsg ("Chromosome length can't be 0");
    for(int i = 0; i < length; i++)
    {
        setVal(i, t[i]);
    }
}






int Chromosome::twoDis( Chromosome & s) 
{
    double ran = myRand.uniform();
    if (ran>R1){
        return edgeDis(s);
    }
    else {
        return nodeDis2(s);
    }

} // random choose distance


int Chromosome::edgeDis(Chromosome & s)
{
    if (!ispt)
        buildpt();
    if (!s.ptVal())
        s.buildpt();

    int res = 0;

    int tval = 0;
    for (int i = 0; i < length; i++){
        tval = s.getVal(i);

        if (s.getVal((i + 1) % length) != gene[(pertable[tval] + 1) % length] && s.getVal((i + 1) % length) != gene[(pertable[tval] + length - 1) % length])
            res++;
        if (s.getVal((i + length - 1) % length) != gene[(pertable[tval] + 1) % length] && s.getVal((i + length - 1) % length) != gene[(pertable[tval] + length - 1) % length])
            res++;
    }

    return res;
} // edgeDis

int Chromosome::nodeDis(Chromosome & s)
{
    if (!ispt)
        buildpt();
    if (!s.ptVal())
        s.buildpt();

    int res = 0;
    int tval = 0;
    int val = 0;
    for (int i = 0; i < length; i++){
        tval = s.getVal(i);
        val = pertable[tval] - i;
        if (val < 0)
            val = val * -1;
        res += val;
    }

    return res;
} // nodeDis

int Chromosome::nodeDis2(Chromosome & s)
{
    if (!ispt)
        buildpt();
    if (!s.ptVal())
        s.buildpt();

    int res = 0;
    for (int i = 0; i < length; i++){
        if (gene[i] != s.getVal(i))
            res++;
    }

    return res;
} // nodeDis2


int Chromosome::orderDis(Chromosome & s)
{
    if (!ispt)
        buildpt();
    if (!s.ptVal())
        s.buildpt();

    int res = 0;
    for (int i = 0; i < length; i++){
        for (int j = i + 1; j < length; j++){
            if (pertable[s.getVal(i)]>pertable[s.getVal(j)]){
                res++;
            }
        }
    }


    return res;
} // orderDis



int Chromosome::orderDis2(Chromosome & s)
{
    if (!ispt)
        buildpt();
    if (!s.ptVal())
        s.buildpt();

    int res = 0;
    for (int i = 0; i < length - 1; i++){

        if (pertable[s.getVal(i)] >= length - 1)
            res++;
        if (s.getVal(i + 1) != gene[pertable[s.getVal(i)] + 1])
            res++;

    }


    return res;
} // orderDis2

void Chromosome::buildpt()
{
    for (int i = 0; i < length; i++){
        pertable[getVal(i)] = i;
    }

    return ;
} // build permutation table

bool Chromosome::ptVal()
{

    return ispt;
} // Is permutation table built
