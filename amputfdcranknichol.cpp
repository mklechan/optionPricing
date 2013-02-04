/*
 * europutfdcranknichol.cpp
 *
 *  Created on: May 16, 2011
 *      Author: mklechan
 */

#include <cmath>
#include <ql/quantlib.hpp>
#include <algorithm>
#include <fstream>
using namespace std;
using namespace QuantLib;

double AmPutFDCrankNichLN( double S, double K, double r, double sigma, double time, int no_S_steps, int no_t_steps) {
    int N=no_t_steps;
    double delta_t = time/N;
    double sigma_sqr = pow(sigma,2.0);
    int M=no_S_steps + (no_S_steps%2);
    double X = log(S);
    double delta_X = sigma*sqrt(4.0*delta_t);
    double vega = r - sigma_sqr/2;
    Array X_values(M+1);
    for (int m=0;m<=M;m++) {
    	   if (m < M/2 )  { X_values[m] = X + (M/2-m)*delta_X; }
    	   else if (m == M/2 ) { X_values[m] = X;}
    	   else {X_values[m] = X - (m-M/2)*delta_X;}
    }
    cout << X_values << endl;
    double Pu = -0.25*delta_t*(sigma_sqr/pow(delta_X,2.0) + vega/delta_X);
    double Pm = 1 + delta_t*sigma_sqr/(2.0*pow(delta_X,2.0)) + r*delta_t/2.0;
    double Pd = -0.25*delta_t*(sigma_sqr/pow(delta_X,2.0) - vega/delta_X);

    Matrix A(M+1,M+1);
    A[0][0] = 1.0;
    A[0][1] = -1.0;
    for (int j=1;j<M;++j) {
        A[j][j-1] = Pu;
        A[j][j]  = Pm;
        A[j][j+1] = Pd;
    };
    A[M][M]=-1.0;
    A[M][M-1]=1.0;
    //cout << A << endl;
    Array B(M+1);
    double Cup; double Cdown; double Cm;
    for (int m=1;m<M;++m){
    	 Cup = max(0.0,K-exp(X_values[m-1]));
    	 Cm = max(0.0,K-exp(X_values[m]));
    	 Cdown = max(0.0,K-exp(X_values[m+1]));
    	 B[m]= -Pu*Cup - (Pm-2)*Cm - Pd*Cdown;
    };
    B[0]=0; B[M]=-(exp(X_values[M])-exp(X_values[M-1]));

    Array F(M+1);
    F = inverse(A)*B;
    for(int t=N-1;t>0;--t) {
        B = F;
        B[0]=0; B[M]=-(exp(X_values[M])-exp(X_values[M-1]));
        F = inverse(A)*B;
        for (int m=1;m<M;++m) { 
            F[m] = max(F[m], K-exp(X_values[m]));
        };
    };
    cout << endl;
    ofstream myfile; myfile.precision(15); myfile.open("amputcranknich.csv");
    for (int i=0;i<=M;i++) {myfile << exp(X_values[i]) << ",";} myfile << endl;
    for (int i=0;i<=M;i++) {myfile << exp(-r*time)*F[i] << ",";} myfile << endl;
    myfile.close();
    return exp(-r*time)*F[M/2+1];
};
