#ifndef _MODEL_ORDER_H
#define _MODEL_ORDER_H

#include "model.h"
#include "myrand.h"
#include "Array2D.h"


class Model_Order : virtual public IModel
{
public:
    Model_Order();
    Model_Order(int n_ell, Chromosome* pop, int size);
    void init(int n_ell, Chromosome* pop, int size);
    virtual ~Model_Order();
    virtual void build();
    virtual Chromosome resample();
    virtual Chromosome resampleWithTemplate(const Chromosome*);
    virtual void update();

    virtual void dump();
protected:
    
	void clear();
    
    
	Array2D<double>& histogram;
	int *firstPosition;
	
    double epsilon;
private:
    int N_SIZE;
    MyRand random;
};



#endif
