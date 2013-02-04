#include <cmath>             
#include <algorithm>             
#include <vector>
#include <ql/quantlib.hpp>
#include "Hmwk4.hpp"
using namespace std;
using namespace QuantLib;
Real BinAmPut (Real S, Real X, Real Rf, Real Vol, Real T, int N){

   Real delta = T/N;
   Real R = exp(Rf*(delta));            
   Real Rinv = 1.0/R;                    
   Real u = Ud(delta, Vol);
   Real uu = u*u;
   Real d = 1.0/u;
   Real p = Pd(Rf, delta, Vol);
   Real q = 1.0-p;
   vector<Real> prices(N+1);       
   prices[0] = S*pow(d, N);  
   for (int i=1; i<=N; ++i) prices[i] = uu*prices[i-1];

   vector<Real> put_values(N+1);       
   for (int i=0; i<=N; ++i) put_values[i] = max(0.0, (X-prices[i])); 

   for (int step=N-1; step>=0; --step) {
      for (int i=0; i<=step; ++i) {
	     put_values[i] = (p*put_values[i+1]+q*put_values[i])*Rinv;
	     prices[i] = d*prices[i+1];
	     put_values[i] = max(put_values[i],(X-prices[i]));
      }
   }
   return put_values[0];
}
