#ifndef _CHROMOSOME_H
#define _CHROMOSOME_H

#include <vector>
#include "global.h"
#include "evaluator.h"


class Chromosome
{
    public:
        Chromosome ();
        Chromosome (int n_ell);
        Chromosome (int n_ell, IEvaluator *evaluator);

        ~Chromosome ();

        Chromosome& operator= (const Chromosome & c);

        void init (int n_ell);
        void init (int n_ell, IEvaluator*);

        int getVal (int index) const;
        void setVal (int index, int val);

        double getFitness ();

        /** real evaluator */
        double evaluate ();

        bool isEvaluated () const;

        void printf () const;

        int getLength () const;

        double getMaxFitness () const;

        void dump() const;

        void clone(const Chromosome& );

        void setGenes(const int*); //set genes with an int array
        int twoDis( Chromosome & s) ;
        int edgeDis(Chromosome & s);
        int nodeDis(Chromosome & s);
        int nodeDis2(Chromosome & s);
        int orderDis(Chromosome & s);
        int orderDis2(Chromosome & s);
        void buildpt();
        bool ptVal();

    protected:
        int * pertable;
        bool ispt;
        std::vector<int> gene;
        int length;
        double fitness;
        bool evaluated;
        IEvaluator *evaluator;
};
#endif
