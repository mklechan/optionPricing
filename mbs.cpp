/*
 * mbs.cpp
 *
 *  Created on: May 29, 2011
 *      Author: mklechan
 */

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <time.h>
#include <ql/quantlib.hpp>
#include <boost/math/distributions.hpp>
using namespace std;
using namespace QuantLib;

double* PMNormGen(int num){

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

// main function
double mbsFunc(double kappa, double Rbar) {

   cout << "Processing MBS with kappa = " << kappa << "and Rbar = " << Rbar << endl;
   //CIR process variables
   //double kappa = 0.4;
   double R0 = 0.01;
   double sigma = 0.5;
   //double Rbar = 0.06;
   double T = 30;
   double dt = 1.0/360.0;
   int num_sims = 1000;
   // Mortage Variables
   double WAC = 0.08;
   double MBS_P0;
   double MBS_i;
   double PVtminus1; double Dt; double Ct; double r10;
   double RIt; double BUt; double SGt; double SYt; double CPRt; double IPt; double TPPt;

   double DaysInMonth = 30;
   double PV0 = 100000;  //principle of loan
   double r = WAC/12.0;
   double N = T*12;
   bool R10_Analytic_MC = 1;

   double Rt;
   double Rtplus1;
   double Ravg;
   double negIntCnt = 0;

     // Problem 1a
     Array RtArray(int((N+120)*DaysInMonth));
     double accum = 0;
     for (int i=1;i<=num_sims;i++) {
        MBS_i= 0;
        Rt = R0;
        PVtminus1 = PV0;
        RtArray[0] = R0;
        double* NormVec = PMNormGen((N+120)*DaysInMonth);
        for (int k=1;k<=(N+120)*DaysInMonth;k++) {
           Rtplus1 = Rt + (kappa*(Rbar-Rt)*dt) + (sigma*sqrt(Rt)*sqrt(dt)*NormVec[k-1]);
           if (Rtplus1<0) negIntCnt++;
           Rt = fabs(Rtplus1);
           RtArray[k] = Rt;
        }
        for (int j=1;j<=N;j++) {
           double accumR = 0;
           for (int k=1;k<=j*DaysInMonth;k++) {
              accumR = accumR + RtArray[k];
           }
           Ravg = accumR/(j*DaysInMonth);
           Dt = exp(-Ravg*j/12);
           // Find r10
           if (R10_Analytic_MC==1) {
              accumR = 0;
              for (int k=int(1+(DaysInMonth*(j-1))); k<=(int(DaysInMonth*(j-1))+(10*12*DaysInMonth)); k++) {
                 accumR = accumR + RtArray[k];
              }
              r10 = accumR/(10*12*DaysInMonth);
           }
           else {
              r10 = (RtArray[int((DaysInMonth*(j-1)))]*exp(-kappa*10)) + (Rbar*(1-exp(-kappa*10)));
           }
           // compute CPRt
           RIt = 0.28 + 0.14*atan(-8.57+430*(WAC-r10));
           BUt = 0.3 + (0.7*PVtminus1/PV0);
           SGt = j/30 < 1 ? j/30 : 1;
           switch (j % 12) {
            case 1 : SYt = 0.94; break;
            case 2 : SYt = 0.76; break;
            case 3 : SYt = 0.74; break;
            case 4 : SYt = 0.95; break;
            case 5 : SYt = 0.98; break;
            case 6 : SYt = 0.92; break;
            case 7 : SYt = 0.98; break;
            case 8 : SYt = 1.10; break;
            case 9 : SYt = 1.18; break;
            case 10 : SYt = 1.22; break;
            case 11 : SYt = 1.23; break;
            case 0 : SYt = 0.98; break;
         }
         CPRt = RIt*BUt*SGt*SYt;
         IPt = PVtminus1*r;
         TPPt = (PVtminus1*r*((1/(1-pow(1+r, -N+j-1)))-1)) + ((PVtminus1) - (PVtminus1*r*((1/(1-pow(1+r, -N+j-1)))-1)))*(1-pow(1-CPRt, 1.0/12.0));
         Ct = IPt + TPPt;
         MBS_i = MBS_i + Dt*Ct;
         PVtminus1 = PVtminus1 - TPPt;
   }
         //writeFile2 << MBS_i << endl;
         accum = accum + MBS_i;
  }
   MBS_P0 = accum/num_sims;
   cout << "MBS Price = " << MBS_P0 << endl;
   cout << "Times negative interest rates seen: " << negIntCnt << endl;

   return MBS_P0;
}
