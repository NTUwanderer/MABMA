#ifndef _MODEL_H
#define _MODEL_H
#include <string>
#include <vector>
#include "myrand.h"

class Chromosome;
class IEvaluator;

class IModel
{
public:
    virtual ~IModel() {};
    virtual void init(int n_ell, Chromosome* , int n_size) = 0;
    virtual void build() = 0;
    virtual Chromosome resample() = 0;
    virtual Chromosome resampleWithTemplate(const Chromosome*) = 0;
    virtual void update() = 0;
    const char* getModelName(){return MODEL_NAME.c_str();}

    virtual void dump() = 0;
protected:
    Chromosome* population;
    int n_size;
    int problemSize;
    std::string MODEL_NAME;
private:

};

#endif
