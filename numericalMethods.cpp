/*
 * numericalHelper.cpp
 *
 *  Created on: May 22, 2011
 *      Author: mklechan
 */


#include <iostream>
#include <ql/quantlib.hpp>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <Eigen>

using namespace std;
using namespace QuantLib;
using Eigen::MatrixXd;
using Eigen::VectorXd;


double NormCdf (double X) {
   double d1 =	0.0498673470;
   double d2 =	0.0211410061;
   double d3 =	0.0032776263;
   double d4 = 	0.0000380036;
   double d5 =	0.0000488906;
   double d6 =	0.0000053830;
   double result;
   if (X >= 0) {
	  result = 1 - 0.5 * pow((1 + d1*X + d2*pow(X,2.0) + d3*pow(X,3.0) \
			  + d4*pow(X,4.0) + d5*pow(X,5.0) + d6*pow(X, 6.0)) , -16.0);
   }
   else {
	  result = 1 - NormCdf(-X);
   }
   return result;

}

void csvDump(double* arrayOut, int size, const char* name) {
	ofstream myfile;
	myfile.precision(32);
	myfile.open (name);
	for (int i=0; i<size; i++)
	{
	myfile << arrayOut[i] << ",";
	}
	myfile.close();
}

double meanEigen(const Eigen::VectorXd& Y) {
    return Y.sum()/Y.size();
}

double mean(int arraysize, double* Y) {
	double accum = 0;
    for (int i=0; i<arraysize; i++) { accum += Y[i]; }
    return accum/Real(arraysize);
}

double variance(int arraysize, double* Y) {
	double meanY = mean(arraysize, Y); double accum = 0;
    for (int i=0; i<arraysize; i++) { accum += pow((Y[i]-meanY),2.0); }
    return accum/double(arraysize);
}

double stdDev(int arraysize, double* Y) {
    double varianceY = variance(arraysize, Y);
    return sqrt(varianceY);
}


double correlation(int arraysize, double* X, double* Y) {
    //compute covariance
    double X_mean = mean(arraysize, X); double Y_mean = mean(arraysize, Y); double accum = 0;
    for (int i=0; i<arraysize; i++) {
        accum += (X[i]-X_mean) * (Y[i]-Y_mean);
    }
    double covariance = accum / double(arraysize);
    double sigmaX = stdDev(arraysize, X); double sigmaY = stdDev(arraysize, Y);
    double correlation = covariance / (sigmaX * sigmaY);
    return correlation;
}

VectorXd PMNormGen(int num){

	VectorXd NormArray(num);
	double V1,V2,W;
    BigInteger seed=SeedGenerator::instance().get();
    MersenneTwisterUniformRng unifMt(seed);
    for (int i=0; i<num; i++) {
       do {
	       V1 = 2.0*unifMt.next().value-1;
	       V2 = 2.0*unifMt.next().value-1;
	       W = V1*V1+V2*V2;

       } while ( W > 1 || W==0);
   	   NormArray(i) = V1*sqrt(-2.0*log(W)/W);
  	   NormArray(i+1) = V2*sqrt(-2.0*log(W)/W);
  	   i++;
    }
    return NormArray;
}
