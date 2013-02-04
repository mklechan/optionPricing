#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
using namespace std;

double AmCallFDExplicit( double S,  double K, double r, double sigma, double time, int no_S_steps, int no_t_steps) {
    double sigma_sqr = pow(sigma,2);
    int M = no_S_steps + (no_S_steps%2); 
    double delta_S = 2.0*S/M;
    vector<double> S_values(M+1,0.0);
    for (int m=1; m<=M; m++) { S_values[m] = m*delta_S; };
    int N = no_t_steps;
    double delta_t = time/N;
    
    vector<double> a(M,0.0);
    vector<double> b(M,0.0);
    vector<double> c(M,0.0);
    double r1=1.0/(1.0+r*delta_t);
    double r2=delta_t/(1.0+r*delta_t);
    for (int j=1;j<M;j++){
	a[j] = r2*0.5*j*(-r+sigma_sqr*j);
	b[j] = r1*(1.0-sigma_sqr*j*j*delta_t);
	c[j] = r2*0.5*j*(r+sigma_sqr*j);
    };
    vector<double> f_next(M+1,0.0);
    for (int m=0;m<=M;++m) { f_next[m]=max(0.0,S_values[m]-K); };
    vector<double> f(M+1,0.0);
    for (int t=N-1;t>=0;--t) {
	f[0]=0;
	for (int m=1;m<M;++m) {
	    f[m]=a[m]*f_next[m-1]+b[m]*f_next[m]+c[m]*f_next[m+1];
	    f[m] = max(f[m],S_values[m]-K);  
	};
	f[M] = S_values[M]-K;
	for (int m=0;m<=M;++m) { f_next[m] = f[m]; };
    };

    ofstream myfile; myfile.precision(15); myfile.open("amcallexplicit.csv");
    for (int i=0;i<=M;i++) {myfile << S_values[i] << ",";} myfile << endl;
    for (int i=0;i<=M;i++) {myfile << f[i] << ",";} myfile << endl;
    myfile.close();

    double C = f[M/2];
    return C;
};
