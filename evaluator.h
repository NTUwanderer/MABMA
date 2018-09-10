#ifndef _EVALUATOR_H
#define _EVALUATOR_H
#include <string>
#include <vector>

#include "Array2D.h"

class Chromosome;

class IEvaluator
{
public:
	bool init() {NFE = 0; return true;}
	virtual double calc(const Chromosome*) = 0;
	virtual double getBest() = 0;
	virtual int getProblemSize() = 0;
	virtual ~IEvaluator() {}
	
	unsigned int getNFE(){return NFE;}
protected:
	void increaseNFE() {NFE++;}
	unsigned int NFE;
};

/*------------------TSP PART-----------------*/
class TSP_Evaluator : public IEvaluator
{
public:
    TSP_Evaluator();
    TSP_Evaluator(const char* fileName);
    ~TSP_Evaluator();
    bool readData(const char*);
    void readBest(const char*);

    double calc(const Chromosome*);
    double getBest();
    int getProblemSize();

protected:
    Array2D<double> vDistanceMatrix;
	
    int problemSize;
    double best;
};


class TSPFixed_Evaluator : public TSP_Evaluator
{
public:
    TSPFixed_Evaluator() : TSP_Evaluator() {}
    TSPFixed_Evaluator(char* fileName);
    ~TSPFixed_Evaluator() {}

    double calc(const Chromosome*);
};

/*---------------Flat Evaluator--------------------*/
//fitness return 0 anyway
class Flat_Evaluator : public IEvaluator
{
public:
	Flat_Evaluator():problemSize(51) { IEvaluator::init(); }
	Flat_Evaluator(int a): problemSize(a) { IEvaluator::init(); }
	~Flat_Evaluator() {}
 	double calc(const Chromosome* a){ this->increaseNFE(); return 1.0;}
	double getBest() { return 0.0;}
	int getProblemSize() { return this->problemSize;}
private:
	int problemSize;
};

/*---------------FSSP Evaluator--------------------*/
class FSSP_Evaluator : public IEvaluator
{
public:
    FSSP_Evaluator();
	FSSP_Evaluator(const char* fileName);
	~FSSP_Evaluator();
	bool readData(const char*);

	double calc(const Chromosome*);
	double getBest();
	int getProblemSize();

protected:
	void calcDistanceMatrix();
	
private:
	int problemSize;
	double best;
	
	int job_num;
	int machine_num;
	bool isValid;
	int** job_times;
	int* due_times;
	
	long initial_seed;
	int current_job;
	
	int* flowTable;
};

/*---------------LOP Evaluator--------------------*/
class LOP_Evaluator : public IEvaluator
{
public:
    LOP_Evaluator();
	LOP_Evaluator(const char* fileName);
	~LOP_Evaluator();
	bool readData(const char*);
	void readBest(const char*);
	
	double calc(const Chromosome*);
	double getBest();
	int getProblemSize();
	

protected:
	
private:
	int problemSize;
	double best;
	
	bool isReverse;
	
	Array2D<double> costMatrix;
};

/*---------------QAP Evaluator--------------------*/
class QAP_Evaluator : public IEvaluator
{
public:
    QAP_Evaluator();
	QAP_Evaluator(const char* fileName);
	~QAP_Evaluator();
	bool readData(const char*);
	void readBest(const char*);
	
	double calc(const Chromosome*);
	double getBest();
	int getProblemSize();
	
	void dump();
protected:
	
private:
	int problemSize;
	double best;

	bool isReverse;
	
	Array2D<double> distanceMatrix;
	Array2D<double> flowMatrix;
};


/*---------------CVRP Evaluator--------------------*/
class CVRP_Evaluator : public IEvaluator
{
public:
    CVRP_Evaluator();
	CVRP_Evaluator(const char* fileName);
	~CVRP_Evaluator();
    void init(const char*);
	bool readData(const char*);
	void readBest(const char*);
	
	double calc(const Chromosome*);
	double getBest();
	int getProblemSize();
	
	void dump();
protected:
	
private:
	int problemSize; //Nc + Nv
	double best;
	
	Array2D<double> cost;
    int capacity;
    int numTrucks;
    int numCustomer;
	int* demand;
};
#endif
