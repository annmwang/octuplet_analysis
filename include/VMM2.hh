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
std::vector <double> xpos;
std::vector <double> zpos;
std::vector <double> param;

//the array of parameters
double *cx, *mx, *cy, *my;
double *xx[4] = {cx, mx, cy, my}; 

//get z, assign z coordinate to board (mmfe8) number
double getz (const int mmfe8) { 
  if (mmfe8 == 0)
    return 0;
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
  if (i == 1)
    return 0;
  if (i == 2)
    return -17.9;
  if (i == 3) 
    return -17.9; 
  if (i == 4)
    return -17.9;  
  if (i == 5)
    return -17.9; 
  if (i == 6)
    return 0;
  if (i == 7)
    return 0;
  return -1;
} 

//function to be minimized
double getx(int mmfe8, double chword) {
	double x;
	if (mmfe8 ==1 | mmfe8 == 7)
	{
		cout << "BOARD: " << mmfe8 << " channel " << chword << endl;
		x = (chword)*.4 +.1; 
		cout << "POS: " << x << endl;
	}
	if (mmfe8==2 | mmfe8 == 4)
	{
		cout << "BOARD: " << mmfe8 << " channel " << chword << endl;
		//x = (chword)*.4+.1 + 4.768456; 
		x = (chword)*.4+.1 + 5.24; 
	}
	if (mmfe8==3 | mmfe8 ==5)
	{
		cout << "BOARD: " << mmfe8 << " channel " << chword << endl;
		x = 204.6 - (chword)*.4+.1 - 5.24; 
		//x = 204.6 - (chword)*.4+.1 - 4.768456; 
	}
	if (mmfe8==0 | mmfe8==6)
	{		
		cout << "BOARD: " << mmfe8 << " channel " << chword << endl;
		x = 204.6 - (chword)*.4 +.1; 
	}
	return x;
} 

//void
//VMM2::Loop ();  //loops over VMM2_data tree and reconstructs line of best fit through each trigger's hits

//main function
//int 
//main (const int argc, const char *const *const argv) ;
