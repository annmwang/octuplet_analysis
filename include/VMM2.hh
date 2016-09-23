//#ifndef VMM2_HH
#define VMM2_HH
#include <iostream>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include "TMath.h"
#include "TProfile.h"
#include  <vector>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TVector3.h"
#include "TGraph.h"
#include <fstream>	

using namespace std;

//////////////////////////////////////////////////// 
// VMM2 geometry and reconstruction module for running on VMM2 (2016) data. 
// Data obtained from Paolo Giromini, June 2016. 
// This code based on beam.C and associated classes. 
// NOTE: It is assumed that data in the root file is ordered chronologically by the branch Ntrig  
// (see Run2.root for branch structure) 
////////////////////////////////////////////////////

//global variables
std::vector <double> param;
std::vector <double> xpos;
std::vector <double> zpos;

//the array of parameters
double *cx, *mx, *cy, *my;
double *xx[4] = {cx, mx, cy, my}; 

//get z, assign z coordinate to board (mmfe8) number
double getz (const int mmfe8) { 
  if (mmfe8 == 0)
    return 0.0;
  if (mmfe8 == 1)
    return 11.2;
  if (mmfe8 == 2)
    return 32.4;
  if (mmfe8 == 3)
    return 43.6;
  if (mmfe8 == 4)
    return 113.6;
  if (mmfe8 == 5)
    return 124.8;
  if (mmfe8 == 6)
    return 146.0;
  if (mmfe8 == 7)
    return 157.2;
  return -1;
}


//get alpha, the board's strip angle measured clockwise from the vertical. Angle measured in radians.
double getalpha (const int i) { 
  if (i == 0)
    return 0;
  if (i == 1)
    return 0;
  if (i == 2)
    return -0.0261799;
  if (i == 3)
    return 0.0261799;
  if (i == 4)
    return -0.0261799;
  if (i == 5)
    return 0.0261799;
  if (i == 6)
    return 0;
  if (i == 7)
    return 0;
  return -1;
}


//get translation in y
double gettrans (const int i) { 
  if (i == 0)
    return 0;
  else if (i == 1)
    return 0;
  else if (i == 2)
    return 17.9;
  else if (i == 3) 
    return 17.9; 
  else if (i == 4)
    return 17.9;  
  else if (i == 5)
    return 17.9; 
  else if (i == 6)
    return 0;
  else if (i == 7)
    return 0;
  else
    return -1;
} 

//function to be minimized
double minfunc(const double *xx) {
  double sum = 0;
  // this currently assumes 8 points, need to fix
  for (int i=0; i < xpos.size() ; i++) {
    double magfx = xx[0] + xx[1]*zpos[i];
    //cout << "magfx= " << magfx << endl;
    double magfy = xx[2] + xx[3]*zpos[i]+gettrans(i);
    //cout << "magfy= " << magfy << endl;
    double a = getalpha(i);
    //cout << "a= " << a << endl;
    double dx = TMath::Tan(a)*magfy;
    //cout << "dx= " << dx << endl;
    double d = (magfx + dx - xpos[i])/sigma;
    //cout << "d= " << d << endl;
    double dsquared = TMath::Power(d, 2);
    //cout << "dsquared= " << dsquared << endl;
    sum = sum + dsquared;
    //cout << "sum= " << sum << endl;
  }
  double totsum = sum;
  //cout << "totsum= " << totsum << endl;
  return totsum;
}


// double minfunc(const double *xx) {
//   double sum = 0;
//   for (int i=0; i < xpos.size() ; i++) {
//     double magfx = xx[0] + xx[1]*zpos[i];
//     //cout << "magfx= " << magfx << endl;
//     double magfy = xx[2] + xx[3]*zpos[i]+gettrans(i);
//     //cout << "magfy= " << magfy << endl;
//     double a = getalpha(i);
//     //cout << "a= " << a << endl;
//     double a1 = TMath::Cos(a); //a is measured in radians
//     //cout << "a1= " << a1 << endl;
//     double a2 = TMath::Sin(-a);
//     //cout << "a2= " << a2 << endl;
//     double sigma = 0.4; //TODO
//     TVector3 fxy(magfx,magfy,0);
//     TVector3 xpo(xpos[i],0,0);
//     TVector3 uperp(a1, a2, 0);
//     double d1 = fxy.Dot(uperp);
//     //cout << "d1= " << d1 << endl;
//     double d2 = xpo.Dot(uperp);
//     //cout << "d2= " << d2 << endl;
//     double d = (d1-d2)/sigma;
//     //cout << "d= " << d << endl;
//     double dsquared = TMath::Power(d, 2);
//     //cout << "dsquared= " << dsquared << endl;
//     sum = sum + dsquared;
//     //cout << "sum= " << sum << endl;
//   }
//   double totsum = sum;
//   //cout << "totsum= " << totsum << endl;
//   return totsum;
// }

//getx
double getx(int mmfe8, double chword) {
	double x;
	if (mmfe8 ==1 | mmfe8 == 7)
	{
		cout << "BOARD: " << mmfe8 << " channel " << chword << endl;
		x = (chword - 1)*.4 + .1; 
		cout << "POS: " << x << endl;
	}
	if (mmfe8==2 | mmfe8 == 4)
	{
		cout << "BOARD: " << mmfe8 << " channel " << chword << endl;
		x = (chword - 1)*.4 + .1 - 5.24; 
	}
	if (mmfe8==3 | mmfe8 ==5)
	{
		cout << "BOARD: " << mmfe8 << " channel " << chword << endl;
		x = 204.6 - ((chword - 1)*.4 + .1 - 5.24); 
	}
	if (mmfe8==0 | mmfe8==6)
	{		
		cout << "BOARD: " << mmfe8 << " channel " << chword << endl;
		x = 204.6 - ((chword - 1)*.4 + .1); 
	}
	return x;
} 

int NumericalMinimization(const char * minName = "Minuit2",
			  const char *algoName = "" ,
			  int randomSeed = -1)
{
  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(minName, algoName);

  // set tolerance , etc...
  min->SetMaxFunctionCalls(100000000); // for Minuit/Minuit2
  min->SetMaxIterations(10000);  // for GSL
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);

  // create function wrapper for minmizer
  // a IMultiGenFunction type
  ROOT::Math::Functor f(&minfunc,4);
  double step[4] = {0.01,0.01,0.01,0.01}; //TODO
  //  double step[4] = {0.01,0.01,0.01,0.01}; //TODO
  // starting point

  double variable[4] = {100, 0, 100, 0}; //TODO
  if (randomSeed >= 0) {
    TRandom2 r(randomSeed);
    variable[0] = r.Uniform(-20,20);
    variable[1] = r.Uniform(-20,20);
    variable[2] = r.Uniform(-20,20);
    variable[3] = r.Uniform(-20,20);
  }

  min->SetFunction(f);

  // Set the free variables to be minimized
  min->SetVariable(0,"cx",variable[0], step[0]); //TODO When we ran tests on fake data these were limited to -204.6 to 409.2
  min->SetVariable(1,"mx",variable[1], step[1]);
  min->SetVariable(2, "cy", variable[2], step[2]); //TODO When we ran tests on fake data these were limited to -204.6 to 409.2
  min->SetVariable(3, "my", variable[3], step[3]);

  // do the minimization
  min->Minimize();

  const double *xs = min->X();
  std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] << "," << xs[3] << "): " << min->MinValue()  << std::endl;
  for (int j = 0; j<4; j++)
    param.push_back(xs[j]);
  ofstream myfile;
  myfile.open ("dump1_vec.dat");
  myfile << "param" << "\n";
  for(std::vector<double>::const_iterator l = param.begin(); l != param.end(); ++l)
    myfile << *l << "\n";
  myfile.close();
  std::cout << "Parameters: f(" << param[0] << "," << param[1] << "," << param[2] << "," << param[3] << ") " << std::endl;
  std::cout << "MINVALUE: " << min->MinValue() << " f(xs): " << f(xs) << endl;
  // expected minimum is 0
  // if ( min->MinValue()  < 1.E-2  && f(xs) < 1.E-2)
  //   std::cout << "Minimizer " << minName << " - " << algoName
  // 	      << "   converged to the right minimum" << std::endl;
  // else {
  //   std::cout << "Minimizer " << minName << " - " << algoName
  // 	      << "   failed to converge !!!" << std::endl;
  //   Error("NumericalMinimization","fail to converge");
  // }
  return -1;
  //minfunc->Draw();

}


//void
//VMM2::Loop ();  //loops over VMM2_data tree and reconstructs line of best fit through each trigger's hits

//main function
//int 
//main (const int argc, const char *const *const argv) ;
