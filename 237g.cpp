//============================================================================
// Name        : 237g.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <ql/quantlib.hpp>
#include <fstream>
#include <boost/math/constants/constants.hpp>
#include <boost/timer.hpp>

using namespace std;
using namespace QuantLib;
const double pi = boost::math::constants::pi<double>();


double* qlibMTUnif(int num){

    	double * UiArray = new double [num];
        BigInteger seed=SeedGenerator::instance().get();
        MersenneTwisterUniformRng unifMt(seed);

    	for (int i=0; i<num; i++) {
    		UiArray[i]= unifMt.next().value;
    	}
        return UiArray;

}

long XnGen(long X) {

	long y;
	long m = pow(2.0,31)-1;
	y = long(pow(7.0,5)*X) % m;
	return y;
}

double* myLGMUniform(int num) {
	long m = pow(2.0,31)-1;
	long seed=SeedGenerator::instance().get();
	long XnArray [num];
	double * UiArray = new double [num];
	XnArray[0]=XnGen(seed);
	for (int i=1; i<num; i++) {
		XnArray[i]=XnGen(XnArray[i-1]);
	}
	for (int i=0; i< num - 1; i++) {
		UiArray[i]=double(XnArray[i])/double(m);
	}
   // for (int i=0; i<100; i++) {
    //	cout << UiArray[i] << endl;
  //  }
    return UiArray;
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

double* BMNormGen(int num){
	double* NormArray = new double [num];
	double* MTUni1 = qlibMTUnif(num/2);
	double* MTUni2 = qlibMTUnif(num/2);

	for (int i=0; i<num; i++) {
		NormArray[i] = sqrt(-2.0*log(MTUni1[i/2])) * cos(2.0*pi*MTUni2[i/2]);
		NormArray[i+1] = sqrt(-2.0*log(MTUni1[i/2])) * sin(2.0*pi*MTUni2[i/2]);
		i++;
	}
    return NormArray;
}

double* PMNormGen(int num, int &size){
	double* NormArray = new double [num];
	double* MTUni1 = qlibMTUnif(num/2);
	double* MTUni2 = qlibMTUnif(num/2);
    int X = 0;
	for (int i=0; i<num/2; i++) {
	    double V1 = 2.0*MTUni1[i]-1;
	    double V2 = 2.0*MTUni2[i]-1;
	    double W = V1*V1+V2*V2;
	    if (W < 1) {
	    	NormArray[X] = V1*sqrt(-2.0*log(W)/W);
	    	NormArray[X+1] = V2*sqrt(-2.0*log(W)/W);
	    	X+=2;
	    }
	}
	size = X;
    return NormArray;
}

double* PMNormGenForce(int num){

	double* NormArray = new double [num];
	double V1,V2,W;
    BigInteger seed=SeedGenerator::instance().get();
    MersenneTwisterUniformRng unifMt(seed);
    for (int i=0; i<num; i++) {
       do {
	       V1 = 2.0*unifMt.next().value-1;
	       V2 = 2.0*unifMt.next().value-1;
	       W = V1*V1+V2*V2;

       } while ( W > 1 || W==0);
   	   NormArray[i] = V1*sqrt(-2.0*log(W)/W);
  	   NormArray[i+1] = V2*sqrt(-2.0*log(W)/W);
    }
    return NormArray;
}

int main() {


	// Problem 1
	//Generate 10000 numbers with LGM vs MT from Quantlib
	cout.precision(32);
    int NumRnd = 10000;
    double* lgmUi;
    double* MTUi;
	lgmUi = myLGMUniform(NumRnd);
	csvDump(lgmUi, NumRnd, "myLGMUni1.csv");
	MTUi = qlibMTUnif(NumRnd);
	csvDump(MTUi, NumRnd, "MTUni1.csv");

    //Problem 2
	//Create r.v. with Xi=-1 p=0.3  Xi=0 p=0.5 Xi=1 p=0.2
    double Bern1[NumRnd];
    for(int i=0; i<NumRnd; i++){
    	if (lgmUi[i] < 0.2) {
    		Bern1[i]=-1;
    	}
    	else if (lgmUi[i] < 0.6) {
    		Bern1[i]=0;
    	}
    	else Bern1[i]=1;
    }
    csvDump(Bern1, NumRnd, "Bern1.csv");

    //Problem 3
    // Create 1000 binomial r.v's with n=44 p=0.64
    double Binom1[1000];
    for(int i=0; i<1000; i++) {
    	double* lgm44 = myLGMUniform(44);
    	int BinSuccess = 0;
    	for(int x=0; x<44; x++){
    		if (lgm44[x] < 0.64) {
    		BinSuccess++;
    		}
    	}
    	Binom1[i]=BinSuccess;
    }
    csvDump(Binom1, 1000, "Binom1n44p64.csv");
    //Problem 3b Probability that a Binomial(44, 0.64) is at least 40
    int Bin40 = 0;
    for(int i=0; i<1000; i++) {
       if (Binom1[i] > 40) {
    	   Bin40++;
       }
    }
    cout << "Number of < 40s " << Bin40 << endl;
    cout << "The probability of a Binomial(44, 0.64) greater than 40 is: " << double(Bin40/1000) << endl;


    //Problem 4
    // Generate 10000 exponential r.v's with lambda=1.5
    double lambda = 1.5;
    double expRnd[10000];
    for (int i=0; i<10000; i++) {
    	expRnd[i]=-lambda*log(1-lgmUi[i]);
    }
    csvDump(expRnd, 10000, "expRnd.csv");

    //Problem5
    //Generate 5000 Uniform(0,1) and then gen. 5000 Normal(0,1) with Box Mueller
    int NumRnd2 = 5000;
    boost::timer t0;
    double* NormBM = BMNormGen(NumRnd2);
    cout << "Cputime of Box Muller Normal RNG NumRnd2: " << t0.elapsed_max() << endl;
    csvDump(NormBM, NumRnd2, "BMNorm.csv");
    //Generate X number of Normal RN Polar Marsaglia
    int NVSize = 0;
    double* NormPM = PMNormGen(NumRnd2,NVSize);
    cout << "Elements in PM Normvec is: " << NVSize << endl;
    csvDump(NormPM, NVSize, "PMNorm.csv");

    //Problem5f
    // Must make NumRnd2 Normal RN
    boost::timer t1;
    double* NormPM5k = PMNormGenForce(NumRnd2);
    cout << "Cputime of Polar Marsaglia Normal RNG NumRnd2: " << t1.elapsed_max() << endl;
    csvDump(NormPM5k, NumRnd2, "PMNormForce.csv");

    //Problem 6 compare in Matlab (use histfit)
    //Problem 7 Use Quantlib
    double qlibNormVec[NumRnd2];
    BigInteger seed=12324;
    MersenneTwisterUniformRng   unifMt(seed);
    BoxMullerGaussianRng<MersenneTwisterUniformRng> bmGauss(unifMt);
    for (int i=0; i<NumRnd2; i++) {
    	qlibNormVec[i]=bmGauss.next().value;
    }
    csvDump(qlibNormVec, NumRnd2, "qlibNormVec.csv");

    //Problem 8 compare all 3 normal vectors in Matlab gui
    //Problem 9 generate 10000 Bernoulli, Normal, Expon, Uni in Excel
	return 0;

}
