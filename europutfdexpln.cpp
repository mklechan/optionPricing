/*
 * europutfdexpln.cpp
 *
 *  Created on: May 16, 2011
 *      Author: mklechan
 */

#include <cmath>
#include <vector>
#include <fstream>
using namespace std;

double EuroPutFDExplicitLN( double S,  double K,  double r,  double sigma,  double time,  int no_S_steps,  int no_t_steps) {
   int N=no_t_steps;
   double delta_t = time/N;
   double sigma_sqr = pow(sigma,2);
   int M=no_S_steps; if ((no_S_steps%2)==1) { ++M; }
   double X = log(S);
//   double Klog = log(K);
   double vega = r - sigma_sqr/2;
//   double delta_X = sigma*sqrt(4*delta_t);
   double delta_X = sigma*sqrt(8*delta_t);
   vector<double> X_values(M+1);
   for (int m=0;m<=M;m++) {
	   if (m < M/2 )  { X_values[m] = X - (M/2-m)*delta_X;}
	   else if (m == M/2 ) { X_values[m] = X;}
	   else {X_values[m] = X + (m-M/2)*delta_X;}
   }

   double Pu = delta_t*(sigma_sqr/(2*pow(delta_X,2.0)) + vega/(2*delta_X));
   double Pm = 1-delta_t*sigma_sqr/pow(delta_X,2.0) - r*delta_t;
   double Pd = delta_t*(sigma_sqr/(2*pow(delta_X,2.0)) - vega/(2*delta_X));
   vector<double> f_next(M+1);
   for (int m=0;m<=M;++m) { f_next[m]=max(0.0,K - exp(X_values[m])); };
   double f[M+1];
   for (int t=N-1;t>=0;--t) {
      f[0]=K;
      for (int m=1;m<M;++m) {
         f[m]=Pu*f_next[m-1]+Pm*f_next[m]+Pd*f_next[m+1];
      };
      f[M] = 0;
      for (int m=0;m<=M;++m) { f_next[m] = f[m]; };
   };

    ofstream myfile; myfile.precision(15); myfile.open("explicit.csv");
    for (int i=0;i<=M;i++) {myfile << exp(X_values[i]) << ",";} myfile << endl;
    for (int i=0;i<=M;i++) {myfile << exp(-r*time)*f[i] << ",";} myfile << endl;
    myfile.close();
   return exp(-r*time)*f[M/2+1];
};
