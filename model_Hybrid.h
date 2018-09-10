#ifndef _MODEL_HYBRID_H
#define _MODEL_HYBRID_H
#include <string>
#include <vector>
#include "model.h"
#include "model_EHBSA.h"
#include "model_NHBSA.h"
#include "myrand.h"

class Model_Hybrid : public Model_EHBSA, public Model_NHBSA
{
public:
	Model_Hybrid();
	Model_Hybrid(double r1, double r2, int e_ell, Chromosome* pop, int n_size);
    void init(int n_ell, Chromosome* _pop, int p_size);
    virtual void build();
    virtual Chromosome resample();
    virtual Chromosome resampleWithTemplate(const Chromosome*);

    virtual void update();

	virtual void dump();
protected:
    void clearHM();
    Chromosome bitwiseResample(const Chromosome*);
    Chromosome weightResample(const Chromosome*);

private:
	MyRand random;
	double r1, r2;
	int problemSize;
};

#endif
