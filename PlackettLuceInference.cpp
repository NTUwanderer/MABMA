//
//  PlackettLuceInference.cpp
//  PlackettLuceEDA
//
//  Created by Josu Ceberio Uribe on 10/16/12.
//  Copyright (c) 2012 Josu Ceberio Uribe. All rights reserved.
//

#include "PlackettLuceInference.h"
#include "Tools.h"
#include <math.h>
#include <list>
using std::list;

/*
 * The constructor method
 */
CPlackettLuceInference::CPlackettLuceInference(int problem_size, int sel_size, int iterations)
{
    m_N = sel_size;
    m_K = problem_size;
    m_topf = m_K;

    //update factor method variables
    m_rates= new double[sel_size];
    m_shapes= new double[sel_size];
    m_means= new double[sel_size];
    
    m_a= new double[sel_size];
    m_c= new double[sel_size];
    m_d= new double[sel_size];
    
    m_log_qmass = new double[m_K];
    m_aux_means= new double[m_K];
    
    m_iterations=iterations;
    
    //Inference structure initializations
    m_marginals= new Gamma[m_K];
    for (int i=0;i<m_K;i++)
    {
        m_marginals[i].SetShapeAndRate(3.0, 2.0); //3.0 and 2.0
    }
    
    m_messages_fact2var= new Gamma**[m_N];
    for (int n=0;n<m_N;n++)
    {
        m_messages_fact2var[n]=new Gamma*[m_topf];
        //vector<vector<Gamma*>> aux2;
        for (int k=0;k<m_topf;k++)
        {
            m_messages_fact2var[n][k]= new Gamma[m_K-k];
        }
    }

}

/*
 * The destructor method
 */
CPlackettLuceInference::~CPlackettLuceInference()
{
    
    //update factor method variables
    delete [] m_rates;
    delete [] m_shapes;
    delete [] m_means;
    
    delete [] m_a;
    delete [] m_c;
    delete [] m_d;
    
    //inference method variables
    delete [] m_log_qmass;
    delete [] m_aux_means;
    
    delete [] m_marginals;
    
    /*
    for (int n=0;n<m_N;n++)
    {
        for (int k=0;k<m_topf;k++)
        {
            delete [] m_messages_fact2var[n][k];
        }
        delete [] m_messages_fact2var[n];
    }
    delete [] m_messages_fact2var;
     */
}

void PrintArrayOfGammas(Gamma * array, int size)
{
    for(int i=0;i<size;i++)
    {
        cout<<"("<<array[i].Shape<<", "<<array[i].Rate<<")    ";
    }
    cout<<endl;
}

double logMassUGam(Gamma * gamma)
{
    double a= gamma->Shape;
    double b= gamma->Rate;
    double alogb= a*log(b);
    double gamlna= lgamma(a);//estos logaritmos pueden ser incorrectos!!

//    cout<<"gam: "<<gamlna<<endl;
//    gamlna=log(a);
//    cout<<"gam2: "<<gamlna<<endl;
    double result = gamlna-alogb;
    return result;
}


void GammaDiv(Gamma* numerator, Gamma * denominator, Gamma * result)
{
    result->SetToRatio(numerator, denominator);
}

void GammaMul(Gamma* a, Gamma * b, Gamma * result)
{
    result->SetToProduct(a, b);

}
/// Finds the projection of cGam(x|a,b)+dGam(x|a+1,b) onto a gamma
/// Note that this is not strictly a KL-projection because it uses
/// (x,x^2) as features rather than (x,ln(x))
void MomentsProjection(double a, double b, double c, double d, Gamma * result)
{
    double ev = ((c * a) + (d * (a + 1.0))) / (b * (c + d));
    double ev2 = ((c * a * (a + 1.0)) + (d * (a + 1.0) * (a + 2.0))) / (b * b * (c + d));
    double aprime = (ev * ev) / (ev2 - (ev * ev));
    double bprime = aprime / ev;
    result->SetShapeAndRate(aprime, bprime);
}


void CPlackettLuceInference::update_factor (int * w, Gamma * q, Gamma * messages_from_factor, int k, double r)
{
    if (abs((long long)r)>0.00001)
    {
        /// Calculate the messages into the factor. In general, messages are indexed
        /// by variable indices in rank order
        //printf "ord length: %d, marg len: %d messages len: %d, K: %d\n" w.Length q.Length messages_from_factor.Length K
        //cout<<"Calculating messages into the factor..."<<endl;
        Gamma * aux= new Gamma();
        Gamma * aux1= new Gamma();
        Gamma * aux2= new Gamma();
        Gamma * aux3= new Gamma();
        Gamma * messages_to_factor = new Gamma[k];
        for (int l=0;l<k;l++)
        {
            GammaDiv(&q[w[m_K-k+l]],&messages_from_factor[l],aux);
            messages_to_factor[l].Copy(aux);
            
        }
        /// Get the g terms which are   (m:a->wl ^2) * m:wl->a
        /// The ^2 comes from 1-alpha = 2 when alpha-divergence = -1
       // cout<<"g terms.."<<endl;

        Gamma * g = new Gamma[k];
        for (int l=0;l<k;l++)
        {
            GammaMul(&messages_to_factor[l], &messages_from_factor[l],aux1);
            GammaMul(aux1, &messages_from_factor[l],aux2);
            g[l].Copy(aux2);
        }

        //  cout<<"Initializations for projections..."<<endl;
        double sumMeans=0;

        for (int l=0;l<k;l++)
        {
            m_rates[l] = g[l].Rate;
            m_shapes[l] = g[l].Shape;
            m_means[l] = m_shapes[l]/m_rates[l];
            sumMeans += m_means[l];
        }
        
        double mult=m_rates[0]/ (m_shapes[0] - 1.0);

       //  cout<<"k: "<<k<<" m_rate[0]: "<<m_rates[0]<<" m_shapes[0]: "<<m_shapes[0]<<" mult: "<<mult<<endl;
       //if (mult<0)
       // {
       // exit(1);
       // }
        /// The projections
        //cout<<"Calculating projections..."<<endl;
        Gamma * proj = new Gamma[k];
        for (int l=0;l<k;l++)
        {
            if (l==0)
            {
                m_a[l]=m_shapes[l]-1.0;
                m_c[l]=mult*(sumMeans-m_means[0]);
                m_d[l]=1.0;
            }
            else
            {
                m_a[l]=m_shapes[l];
                m_c[l]=1.0+ (mult*(sumMeans-m_means[0]-m_means[l]));
                m_d[l]=mult*m_means[l];
            }
            MomentsProjection(m_a[l],m_rates[l],m_c[l],m_d[l],aux);
            proj[l].Copy(aux);
        }
        
        /// Update the messages and the marginals
        for (int l=0;l<k;l++)
        {
            int wl=w[m_K-k+l];
            GammaMul(&q[wl], &q[wl],aux3);
            GammaDiv(aux3,&proj[l],aux);
            q[wl].Copy(aux);
        
            GammaDiv(&q[wl], &messages_to_factor[l],aux1);
            messages_from_factor[l].Copy(aux1);
        
        }
        
        delete aux;
        delete aux1;
        delete aux2;
        delete aux3;
        //deleting some structures
        delete [] proj;
        delete [] g;
        delete [] messages_to_factor;
       // exit(1);
    }
  //  else
  //  {
  //      // With 0 exponent, no update is possible or necessary
  //  }
}

float CPlackettLuceInference::log_evidence_scale (int * w, Gamma * q, Gamma * messages_from_factor, int k)
{

    int wk = w[m_K-1];
    // double lba = log(q[wk].Rate) - log(q[wk].Shape);
    
    double * loge_arr= new double[k];
    Gamma * g= new Gamma();
    for (int l=0;l<k;l++)
    {
        int wl=w[m_K-k+l];
        double add=2.0; if (l==0) add=1.0;
        double logmult= 0.0;
        if (l!=0)
            logmult= log(q[wk].Rate) - log(q[wk].Shape);
        double shp_num = q[wl].Shape - messages_from_factor[l].Shape + add;
        double rat_num = q[wl].Rate - messages_from_factor[l].Rate;
        g->Shape=shp_num;
        g->Rate=rat_num;
        double log_e = logmult + logMassUGam(g) - logMassUGam(&q[wl]);
        
        loge_arr[l]=log_e;
    }
    delete g;
    double sume = 0;
    for (int le=0;le<k;le++)
        sume+=exp(loge_arr[le]);

    
    double log_ev_scale = -log(sume);

    delete [] loge_arr;
    
    return log_ev_scale;    
}




void CPlackettLuceInference::NormalizeMarginals(Gamma* marginals)
{
 //   double * means = new double[marginals->size()];
    for (int i=0;i<m_K;i++)
        m_aux_means[i]=marginals[i].GetMean();
    
    double sumMeans = 0;
    for (int i=0;i<m_K;i++) sumMeans+=m_aux_means[i];
    
    
    
    for (int i=0;i<m_K;i++)
    {
        marginals[i].Rate=sumMeans*marginals[i].Rate;
    }
    
    //delete [] means;
}

double CPlackettLuceInference::InferPL_Truncated(int ** data, float * weights, float * variances, bool standardInitialization)
{
    //some initializations
//  PrintArray(weights,m_K,"in: ");
    if (standardInitialization==true)
    {
        for (int i=0;i<m_K;i++)
        {
            m_marginals[i].SetShapeAndRate(3.0,2.0);
        }
    }
    else
    {
      
        for (int i=0;i<m_K;i++)
        {
            m_marginals[i].InitializeByExpectationAndVariance(weights[i],weights[i]*0.001);
            
        }
      //  cout<<"out"<<endl;
      //  for (int i=0;i<m_K;i++)  cout<<m_marginals[i].GetMean()<<"   "<<endl;
     //   exit(1);
    }
    

    
    //Initialize from PL factors to variables
    //cout<<"Initializing factors..."<<endl;
    for (int n=0;n<m_N;n++)
    {
        for (int k=0;k<m_topf;k++)
        {
            for (int i=0;i<m_K-k;i++)
                m_messages_fact2var[n][k][i].SetShapeAndRate(1, 0);
        }
    }

    //The schedule
    //cout<<"The Schedule..."<<endl;
    for (int pass=0 ; pass<1 ; pass++)
    {
        //Update all messages
        for (int i=0;i<m_N;i++)
        {
            for (int k=0;k<m_topf;k++)
            {
                update_factor(data[i], m_marginals, m_messages_fact2var[i][k], m_K-k, 1.0);
            }
        }
    }

    double log_evidence=0;
    /*
    for (int i=0;i<m_K;i++) m_log_qmass[i] = logMassUGam(&m_marginals[i]);
   
    double sum_log_qmass = 0;
    for (int i=0;i<m_K;i++) sum_log_qmass+=m_log_qmass[i];
        
    //Now get the evidence scales
    //cout<<"Get evidence scales..."<<endl;
    double sum_log_scales=0.0;
    for (int i=0;i<m_N;i++)
    {
        for (int k=0;k<m_topf;k++)
        {
            sum_log_scales+= log_evidence_scale(data[i], m_marginals, m_messages_fact2var[i][k],m_K-k);
        }
    }
    

    log_evidence=sum_log_qmass + sum_log_scales;
     */
    //Normalize marginals
    NormalizeMarginals(m_marginals);
   
    for (int i=0;i<m_K;i++)
    {
        weights[i]=m_marginals[i].GetMean();
        variances[i]=weights[i]/m_marginals[i].Rate;
  //      cout<<"expectation: "<<weights[i]<<" ,  "<<weights[i]/m_marginals[i].Rate<<endl;
    }
 //   PrintArray(weights,m_K,"out: ");
   
    return log_evidence;

}
