#ifndef _MODEL_PLACKETTLUCE_H_
#define _MODEL_PLACKETTLUCE_H_

#include "model.h"
#include "PlackettLuce.h"
#include "Population.h"

class Model_PlackettLuce : public IModel
{
public:
    Model_PlackettLuce();
    Model_PlackettLuce(int n_ell, Chromosome* pop, int size);
    ~Model_PlackettLuce();
    void init(int n_ell, Chromosome* pop, int size);
    void build();
    Chromosome resample();
    Chromosome resampleWithTemplate(const Chromosome*);
    void update();

    void dump();
protected:
    void init();
    int problemSize;
    std::string MODEL_NAME;
private:
    MyRand random;
    CPlackettLuceModel* plModel;
    CPopulation* m_population;
};
#endif
