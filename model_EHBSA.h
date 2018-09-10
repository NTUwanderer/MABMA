#ifndef _MODEL_EHBSA_H
#define _MODEL_EHBSA_H
#include <string>
#include <vector>
#include "model.h"
#include "myrand.h"
#include "Array2D.h"


class Model_EHBSA : virtual public IModel
{
public:
    Model_EHBSA();
    Model_EHBSA(int n_ell, Chromosome* _pop, int p_size);
    virtual ~Model_EHBSA();
    void init(int n_ell, Chromosome* population, int n_size);
    virtual void build();
    virtual Chromosome resample();
    virtual Chromosome resampleWithTemplate(const Chromosome*);
    virtual void update();


    virtual void dump();
protected:
    void initEHM();
    void clearEHM();

    Array2D<double>& EdgeHistogramMatrix;
    int *firstPositionMatrix;
    double totalEdges;
    MyRand random;

    double epsilon;
private:
    int N_SIZE;
};


class Model_EHBSAwoTemplate : public Model_EHBSA
{
public:
    Model_EHBSAwoTemplate(): Model_EHBSA() { };
    Model_EHBSAwoTemplate(int n_ell, Chromosome* _pop, int p_size) : Model_EHBSA(n_ell, _pop, p_size) { MODEL_NAME = "EHBSAwoTemplate"; };		
    Chromosome resample();

    // const std::string MODEL_NAME = "EHBSAwoTemplate";
private:

};


class Model_EHBSAwCircleResample : public Model_EHBSA
{
public:
    Model_EHBSAwCircleResample(): Model_EHBSA() {};
    Model_EHBSAwCircleResample(int n_ell, Chromosome* _pop, int p_size):Model_EHBSA(n_ell, _pop, p_size) { MODEL_NAME = "EHBSAwCircleResample"; };		
    Chromosome resample();
protected:
    // const std::string MODEL_NAME = "EHBSAwCircleResample";
};

#endif
