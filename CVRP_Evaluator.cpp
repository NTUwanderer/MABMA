#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iostream>
#include <stack>
#include <algorithm>
#include "evaluator.h"
#include "global.h"

using namespace std;

#define for_(i, a, b) for (int i=(a);i<(b);++i)
#define sqr(a) (a)*(a)

CVRP_Evaluator::CVRP_Evaluator()
{
    init(NULL);
}

CVRP_Evaluator::CVRP_Evaluator(const char* fileName)
{
    init(fileName);
}

CVRP_Evaluator::~CVRP_Evaluator()
{

}

//problem size = Nv - 1 + Nc - 1 = Nv + Nc -2
void CVRP_Evaluator::init(const char* fileName)
{
    
    problemSize = 0;
    numCustomer = 0;
    best = 0.0;
    demand = NULL;
    capacity = 0;
    numTrucks = 0;
    
    if(fileName == NULL) return;
    
    readData(fileName);
    readBest(fileName);
}

bool CVRP_Evaluator::readData(const char* filename)
{
  string line;
  ifstream file(filename);
  
  if(!file)
    return false;
    
  while (file) {
    getline(file, line);
    
    if (line.empty()) continue;
    
    string kw, x;
    istringstream ss(line);
    
    ss >> kw;
    
    if (kw == "DIMENSION") {
      int n;
      ss >> x >> n;
      
      numCustomer = n - 1;
      cost.init(n, n);
      demand = new int[n];
    }
    else if (kw == "CAPACITY") {
      int c;
      ss >> x >> c;
      capacity = c;
    }
    else if (kw == "TRUCKS") {
      int t;
      ss >> x >> t;
      numTrucks = t;
      problemSize = numTrucks + numCustomer - 1;
    }
    else if (kw == "NODE_COORD_SECTION") {
      int n = numCustomer + 1;
      
      vector<int> x(n), y(n);
      
      int a, b, c;
      for (;;) {
        getline(file, line);
        
        if (line[0] == 'D')
          break;
        
        istringstream ss(line);
        
        ss >> a >> b >> c;
        
        x[a-1] = b;
        y[a-1] = c;
      }
      
      for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            double dist2 = sqr(x[i]-x[j]) + sqr(y[i]-y[j]);
            double dist = sqrt(dist2);
            cost[i][j] = floor(dist + 0.5);
            cost[j][i] = floor(dist + 0.5);
        }
      }
      
      // DEMAND_SECTION
      
      for (;;) {
        getline(file, line);
        
        if (line[0] == 'D')
          break;
        
        istringstream ss(line);
        int v, d;
        ss >> v >> d;
        demand[v-1] = d;
      }
      
      return true;
    }
  }
  return false;
}

void CVRP_Evaluator::readBest(const char* filename)
{
    string line;
    string fileName(filename);
    
    int pos = fileName.find("cvrp/") + 5; //at "/"
  	fileName.insert(pos, "sol/opt-");
    fileName = fileName.substr(0, fileName.length() - 4);
    
    ifstream file(fileName.c_str());
    if (!file) {
      cout << "file not found: " << fileName << endl;
      return;
    }
    
    while (file) 
    {
        getline(file, line);
        if (line.empty()) continue;
        if(line[0] == 'c')
        {
            string s;
            int b;
            
            istringstream ss(line);
            ss >> s >> b;
            this->best = b;
            break;
        }
    }
    
}


/**
    0,1,...,Nc-1 : indices for customers
    Nc, Nc+1, ..., Nc+Nv-1 for vehicles
    
    As Nc + 1 vertices, Vo is for depot,
    Vi+1 = Ci
*/	
#define K 99 
double CVRP_Evaluator::calc(const Chromosome* c)
{
    stack<int> route;  
    int* weight = new int[numTrucks];
    double* length = new double[numTrucks];
    
    int vehicle = 0;
    
    while(!route.empty())
        route.pop();
    for(int i = 0; i <= problemSize; i++)
    {
        int index = problemSize;
        if(i < problemSize)
            index = c->getVal(i);
        if(index < numCustomer)
            route.push(index + 1);
        else 
        {
            
            if(route.empty())
            {
                length[vehicle] = 0.0;
                weight[vehicle] = 0;
            } // if
            else
            {
                int w = 0;
                double l = 0.0;
                int top = 0;
                while(!route.empty())
                {
                    int next = route.top();
                    route.pop();                 
                    l += cost[next][top];                 
                    w += demand[next]; 
                    top = next;            
                } // while
                l += cost[0][top];                
                
                length[vehicle] = l;
   
                weight[vehicle] = w;   
                vehicle++;                
            } // else
        } // else
    } // for
    
    double fitness = 0.0;
    for (int i = 0; i < vehicle; i++)
    {
        double exceed = (double)(weight[i] - capacity) / ((double)capacity);
        double penalty = K * length[i] * exceed;
        double fitness_i = length[i] + max(0.0, penalty);
        fitness += fitness_i;
    }
    return fitness;
}

double CVRP_Evaluator::getBest()
{
    return this->best;
}

int CVRP_Evaluator::getProblemSize()
{
    return this->problemSize;
}

void CVRP_Evaluator::dump()
{

}
