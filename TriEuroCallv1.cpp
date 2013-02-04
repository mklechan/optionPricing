#include <math.h>
#include <algorithm>
#include "Hmwk4.hpp"
#include <vector>
#include <cmath>
#include <ql/quantlib.hpp>
using namespace std;
using namespace QuantLib;

Real TriEuroCallv1 (  Real S,  Real X,  Real r, Real Vol, Real t, int N) {
    Real delta_t = t/N;
    Real DF = exp(-r*(delta_t));
    Real Vol_sqr=pow(Vol,2.0);
    Real u = exp(Vol*sqrt(3.0*delta_t)); 
    Real d = 1.0/u;
    Real Pu = ( r*delta_t*(1 - d) + pow(r*delta_t,2.0) + pow(Vol,2.0)*delta_t ) / ( (u - d) * (u - 1) );
    Real Pd = ( r*delta_t*(1 - u) + pow(r*delta_t,2.0) + pow(Vol,2.0)*delta_t ) / ( (u - d) * (1 - d) );
    Real Pm = 1 - Pu - Pd;

    vector< vector<Real> > Tree;       
    vector<Real> SVec;
    SVec.push_back(S);
    for (int step=1;step<=N;++step){ 
	Tree.push_back(SVec);
	SVec.insert(SVec.begin(),SVec[0]*d); 
	SVec.push_back(SVec[SVec.size()-1]*u);
    }


    int m = SVec.size();
    vector<Real> values_next = vector<Real>(m);       
    for (int i=0; i<m; ++i) values_next[i] = max(0.0, SVec[i]-X);
    vector<Real> values;
    for (int step=N-1; step>=0; --step) {
	   m = Tree[step].size();
	   values = vector<Real> (m);
	   for (int i=0; i<m; ++i) {
	      values[i] = (Pu*values_next[i+2]+Pm*values_next[i+1] + Pd*values_next[i])*DF;
	   }
	   values_next=values;
    }
    return values[0];
}
