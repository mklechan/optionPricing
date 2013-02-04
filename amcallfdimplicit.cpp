#include <cmath>
#include <ql/quantlib.hpp>
#include <algorithm>
#include <fstream>

using namespace std;
using namespace QuantLib;

double AmCallFDImplicit( double S, double K, double r, double sigma, double time, int no_S_steps, int no_t_steps) {
    double sigma_sqr = pow(sigma,2);
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
	A[j][j]   = 1.0 + delta_t*(r+sigma_sqr*j*j);  
	A[j][j+1] = 0.5*j*delta_t*(-r-sigma_sqr*j);   
    };
    A[M][M]=1.0;
    Array B(M+1);
    for (int m=0;m<=M;++m){ B[m] = max(0.0,S_values[m]-K); };
    Array F=inverse(A)*B;
    for(int t=N-1;t>0;--t) {
	B = F;
	F = inverse(A)*B; 
	for (int m=1;m<M;++m) {	
	    F[m] = max(F[m], S_values[m]-K);
	};
    };
    ofstream myfile; myfile.precision(15); myfile.open("amcallimplicit.csv");
    for (int i=0;i<=M;i++) {myfile << S_values[i] << ",";} myfile << endl;
    for (int i=0;i<=M;i++) {myfile << F[i] << ",";} myfile << endl;
    myfile.close();

    return F[M/2];
};
