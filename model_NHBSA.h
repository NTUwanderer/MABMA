#ifndef _MODEL_NHBSA_H
#define _MODEL_NHBSA_H
#include <string>
#include <vector>
#include "model.h"
#include "myrand.h"
#include "Array2D.h"

class Chromosome;
class IEvaluator;


class Model_NHBSA : virtual public IModel
{
public:
    Model_NHBSA();
    Model_NHBSA(int n_ell, Chromosome* pop, int n_size);
    void init(int n_ell, Chromosome* pop, int n_size);
    virtual ~Model_NHBSA();
    virtual void build();
    virtual Chromosome resample();
    virtual Chromosome resampleWithTemplate(const Chromosome*);
    virtual void update();


    virtual void dump();
protected:
    void initNHM();
    void clearNHM();

    //double **NodeHistogramMatrix;
    Array2D<double>& NodeHistogramMatrix;
    double totalHistogram;
    MyRand random;

    double epsilon;
private:
};

class Model_NHBSAwoTemplate : public Model_NHBSA
{
public:
    Model_NHBSAwoTemplate(): Model_NHBSA() {};
    Model_NHBSAwoTemplate(int n_ell, Chromosome* population, int size) : Model_NHBSA(n_ell, population, size) { MODEL_NAME = "NHBSAwoTemplate"; };	
    Chromosome resample();

    // const std::string MODEL_NAME = "NHBSAwoTemplate";
private:

};

#endif
