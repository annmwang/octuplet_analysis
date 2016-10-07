///
///  \file   RunAlignment.C
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Oct
///

#include <iostream>
#include <stdlib.h>

#include "include/ClusterTree.hh"
#include "include/GeoOctuplet.hh"
#include "include/SimpleTrackFitter.hh"

using namespace std;

double EvaluateMetric(const double* param);

ClusterTree* g_tree;
int g_N_event;
int g_RunNum;

int main(int argc, char* argv[]){

  char inputFileName[400];
  char outputFileName[400];
  
  if ( argc < 7 ){
    cout << "Error at Input: please specify input ClusterTree file, ";
    cout << "output .root file and run number.";
    cout << "Example:   ./MakeClusterTree.x -i input.root -o output.root -r RunNum" << endl;
    return 0;
  }

  bool b_input = false;
  bool b_out   = false;
  bool b_run   = false;
  for (int i=1;i<argc-1;i++){
    if (strncmp(argv[i],"-i",2)==0){
      sscanf(argv[i+1],"%s", inputFileName);
      b_input = true;
    }
    if (strncmp(argv[i],"-o",2)==0){
      sscanf(argv[i+1],"%s", outputFileName);
      b_out = true;
    }
    if (strncmp(argv[i],"-r",2)==0){
      g_RunNum = atoi(argv[i+1]);
      b_run = true;
    }
  }

  if(!b_input){
    cout << "Error at Input: please specify input file (-i flag)" << endl;
    return 0;
  }

  if(!b_out){
    cout << "Error at Input: please specify output file (-o flag)" << endl;
    return 0;
  }

  if(!b_run){
    cout << "Error at Input: please specify run number (-r flag)" << endl;
    return 0;
  }

  TFile* f = new TFile(inputFileName, "READ");
  if(!f){
    cout << "Error: unable to open input file " << inputFileName << endl;
    return false;
  }
  TTree* T = (TTree*) f->Get("ClusterTree");
  if(!T){
    cout << "Error: cannot find tree ClusterTree in " << inputFileName << endl;
    return false;
  }

  g_tree = (ClusterTree*) new ClusterTree(T);
  g_N_event = T->GetEntries();

  ROOT::Math::Minimizer* minimizer =
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
  minimizer->SetMaxFunctionCalls(10000000);
  minimizer->SetMaxIterations(100000);
  minimizer->SetTolerance(0.01);
  minimizer->SetPrintLevel(2);
  
  ROOT::Math::Functor* functor = new ROOT::Math::Functor(&EvaluateMetric, 21);
  minimizer->SetFunction(*functor);
  int ivar = 0;
  for(int i = 0; i < 7; i++){
    minimizer->SetVariable(ivar, Form("tranX_%d",i+1), 0., 0.001);
    ivar++;
  }
  for(int i = 0; i < 7; i++){
    minimizer->SetVariable(ivar, Form("tranY_%d",i+1), 0., 0.001);
    ivar++;
  }
  for(int i = 0; i < 7; i++){
    minimizer->SetVariable(ivar, Form("tranZ_%d",i+1), 0., 0.001);
    ivar++;
  }

  minimizer->Minimize();

  double out_tranX[7];
  double out_tranX_err[7];
  double out_tranY[7];
  double out_tranY_err[7];
  double out_tranZ[7];
  double out_tranZ_err[7];
  
  const double* param = minimizer->X();
  const double* error = minimizer->Errors();
  for(int i = 0; i < 7; i++){
    out_tranX[i] = param[i];
    out_tranX_err[i] = error[i];
  }
  for(int i = 0; i < 7; i++){
    out_tranY[i] = param[7+i];
    out_tranY_err[i] = error[7+i];
  }
  for(int i = 0; i < 7; i++){
    out_tranZ[i] = param[14+i];
    out_tranZ_err[i] = error[14+i];
  }
  
  TFile* fout = new TFile(outputFileName, "RECREATE");
  fout->cd();

  TTree* out_tree = new TTree("AlignmentTree",
			      "AlignmentTree");
  for(int i = 0; i < 7; i++){
    out_tree->Branch(Form("tranX_%d",i+1), &out_tranX[i]);
    out_tree->Branch(Form("tranX_err_%d",i+1), &out_tranX_err[i]);
  }
  for(int i = 0; i < 7; i++){
    out_tree->Branch(Form("tranY_%d",i+1), &out_tranY[i]);
    out_tree->Branch(Form("tranY_err_%d",i+1), &out_tranY_err[i]);
  }
  for(int i = 0; i < 7; i++){
    out_tree->Branch(Form("tranZ_%d",i+1), &out_tranZ[i]);
    out_tree->Branch(Form("tranZ_err_%d",i+1), &out_tranZ_err[i]);
  }

  out_tree->Fill();
  fout->cd();
  out_tree->Write();
  fout->Close();   
}

double EvaluateMetric(const double* param){
  if(!g_tree)
    return 0.;
  
  GeoOctuplet GEO;
  GEO.SetRunNumber(g_RunNum);
  // apply parameters to geometry
  for(int i = 0; i < 7; i++)
    GEO.TranslateX(param[i],i+1);
  for(int i = 0; i < 7; i++)
    GEO.TranslateY(param[i+7],i+1);
  for(int i = 0; i < 7; i++)
    GEO.TranslateZ(param[i+14],i+1);

  double res, sum_res2 = 0;

  SimpleTrackFitter FITTER;
  for(int evt = 0; evt < g_N_event; evt++){
    g_tree->GetEntry(evt);

    if(evt%10000 == 0)
      cout << "Processing event # " << evt << " | " << g_N_event << endl;

    MMClusterList clusters;
    for(int i = 0; i < g_tree->N_clus; i++){
      MMHit clus(g_tree->clus_MMFE8->at(i),
		 0, g_tree->clus_Channel->at(i));
      clusters.AddCluster(MMCluster(clus));
    }

    MMTrack track = FITTER.Fit(clusters, GEO);
    
    for(int i = 0; i < g_tree->N_clus; i++){
      res = GEO.GetResidualX(clusters[i], track);
      sum_res2 += res*res;
    }
  } // end event loop

  for(int i = 0; i < 7; i++)
    cout << "tranX " << i+1 << ": " << param[i] << endl;
  for(int i = 0; i < 7; i++)
    cout << "tranY " << i+1 << ": " << param[i+7] << endl;
  for(int i = 0; i < 7; i++)
    cout << "tranZ " << i+1 << ": " << param[i+14] << endl; 
  cout << "myFCN = " << sum_res2 << endl;
  return sum_res2;
}
