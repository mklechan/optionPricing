/*
 * europutfdimp.cpp
 *
 *  Created on: May 16, 2011
 *      Author: mklechan
 */
#include <cmath>
#include <ql/quantlib.hpp>
#include <algorithm>

using namespace std;
using namespace QuantLib;

double EuroPutFDImplicit( double S, double K, double r, double sigma, double time, int no_S_steps, int no_t_steps) {
    double sigma_sqr = sigma*sigma;
    int M=no_S_steps + (no_S_steps%2); 
    
    double delta_S = 2.0*S/M;
    Array S_values(M+1);
    for (int m=0;m<=M;m++) { S_values[m] = m*delta_S; };
    int N=no_t_steps;
    double delta_t = time/N;

    Matrix A(M+1,M+1);
    A[0][0] = 1.0;
    for (int j=1;j<M;++j) {
        A[j][j-1] = 0.5*j*delta_t*(r-sigma_sqr*j);
        A[j][j]  = 1.0 + delta_t*(r+sigma_sqr*j*j);
        A[j][j+1] = 0.5*j*delta_t*(-r-sigma_sqr*j);
    };
    A[M][M]=1.0;
    Array B(M+1);
    for (int m=0;m<=M;++m){ B[m] = max(0.0,K-S_values[m]); };
    Array F(M+1);
    F = inverse(A)*B;
    for(int t=N-1;t>0;--t) {
        B = F;
        F = inverse(A)*B;
    };
    return F[M/2];
};

