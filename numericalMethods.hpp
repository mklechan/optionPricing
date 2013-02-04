/*
 * numericalMethods.hpp
 *
 *  Created on: May 22, 2011
 *      Author: mklechan
 */

#ifndef NUMERICALMETHODS_HPP_
#define NUMERICALMETHODS_HPP_

#include <Eigen>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

VectorXd PMNormGen(int num);
void csvDump(double* arrayOut, int size, const char* name);
double mean(int arraysize, double* Y);
double variance(int arraysize, double* Y);
double stdDev(int arraysize, double* Y);
double correlation(int arraysize, double* X, double* Y);
double NormCdf(double X);
#endif /* NUMERICALMETHODS_HPP_ */
