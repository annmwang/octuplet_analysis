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

int g_N_track;
vector<MMCluster*> g_cluster;
vector<MMClusterList*> g_clusters;
vector<MMClusterList*> g_clusters_all;

SimpleTrackFitter* g_FITTER;

bool g_tranX;
bool g_tranY;
bool g_tranZ;
bool g_rotX;
bool g_rotY;
bool g_rotZ;

int main(int argc, char* argv[]){

  char inputFileName[400];
  char outputFileName[400];
  
  if ( argc < 7 ){
    cout << "Error at Input: please specify input ClusterTree file, ";
    cout << "output .root file and run number.";
    cout << "Example:   ./RunAlignment.x -i input.root -o output.root -r RunNum" << endl;
    cout << "Other options:" << endl;
    cout << "   --tranX : consider translations in X" << endl;
    cout << "   --tranY : consider translations in Y" << endl;
    cout << "   --tranZ : consider translations in Z" << endl;
    cout << "   --rotX : consider translations around X axis" << endl;
    cout << "   --rotY : consider translations around Y axis" << endl;
    cout << "   --rotZ : consider translations around Z axis" << endl;
    return 0;
  }

  bool b_input = false;
  bool b_out   = false;
  bool b_run   = false;

  g_tranX = false;
  g_tranY = false;
  g_tranZ = false;
  g_rotX = false;
  g_rotY = false;
  g_rotZ = false;
  for (int i=1;i<argc;i++){
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

    if (strncmp(argv[i],"--tranX",7)==0){
      g_tranX = true;
    }
    if (strncmp(argv[i],"--tranY",7)==0){
      g_tranY = true;
    }
    if (strncmp(argv[i],"--tranZ",7)==0){
      g_tranZ = true;
    }
    if (strncmp(argv[i],"--rotX",6)==0){
      g_rotX = true;
    }
    if (strncmp(argv[i],"--rotY",6)==0){
      g_rotY = true;
    }
    if (strncmp(argv[i],"--rotZ",6)==0){
      g_rotZ = true;
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

  if(!g_tranX && !g_tranY && !g_tranZ &&
     !g_rotX && !g_rotY && !g_rotZ)
    g_tranX = true;

  TFile* f = new TFile(inputFileName, "READ");
  if(!f){
    cout << "Error: unable to open input file " << inputFileName << endl;
    return 0;
  }
  TTree* T = (TTree*) f->Get("ClusterTree");
  if(!T){
    cout << "Error: cannot find tree ClusterTree in " << inputFileName << endl;
    return 0;
  }

  GeoOctuplet GEO;
  GEO.SetRunNumber(g_RunNum);
  g_FITTER = new SimpleTrackFitter();

  g_tree = (ClusterTree*) new ClusterTree(T);
  g_N_event = T->GetEntries();

  g_N_track = 0;
  g_clusters.clear();
  g_cluster.clear();
  g_clusters_all.clear();
  for(int i = 0; i < g_N_event; i++){
    g_tree->GetEntry(i);
    int Nc = g_tree->N_clus;
    if(Nc < 6)
      continue;

    MMClusterList all_clusters;

    for(int i = 0; i < Nc; i++){
      MMHit clus(g_tree->clus_MMFE8->at(i),
                 0, g_tree->clus_Channel->at(i), g_RunNum);
      all_clusters.AddCluster(MMCluster(clus));
    }

    MMClusterList* clusters_tmp = new MMClusterList();
    g_clusters_all.push_back(clusters_tmp);
    for (auto clus: all_clusters)
      g_clusters_all.back()->AddCluster(*clus);
  }

  cout << g_clusters_all.size() << " tracks total" << endl;

  g_FITTER = new SimpleTrackFitter();

  ROOT::Math::Minimizer* minimizer =
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
  minimizer->SetMaxFunctionCalls(10000000);
  minimizer->SetMaxIterations(100000);
  minimizer->SetTolerance(0.01);
  minimizer->SetPrintLevel(2);
  
  ROOT::Math::Functor* functor = new ROOT::Math::Functor(&EvaluateMetric, 42);
  minimizer->SetFunction(*functor);
  int ivar = 0;
  for(int i = 0; i < 7; i++){
    minimizer->SetVariable(ivar, Form("tranX_%d",i+1), 0., 0.01);
    ivar++;
  }
  for(int i = 0; i < 7; i++){
    minimizer->SetVariable(ivar, Form("tranY_%d",i+1), 0., 0.1);
    ivar++;
  }
  for(int i = 0; i < 7; i++){
    minimizer->SetVariable(ivar, Form("tranZ_%d",i+1), 0., 0.1);
    ivar++;
  }
  for(int i = 0; i < 7; i++){
    minimizer->SetVariable(ivar, Form("rotX_%d",i+1), 0., 0.001);
    ivar++;
  }
  for(int i = 0; i < 7; i++){
    minimizer->SetVariable(ivar, Form("rotY_%d",i+1), 0., 0.001);
    ivar++;
  }
  for(int i = 0; i < 7; i++){
    minimizer->SetVariable(ivar, Form("rotZ_%d",i+1), 0., 0.001);
    ivar++;
  }

  if(!g_tranX)
    for(int i = 0; i < 7; i++)
      minimizer->FixVariable(i);
  if(!g_tranY)
    for(int i = 0; i < 7; i++)
      minimizer->FixVariable(7+i);
  if(!g_tranZ)
    for(int i = 0; i < 7; i++)
      minimizer->FixVariable(14+i);
  if(!g_rotX)
    for(int i = 0; i < 7; i++)
      minimizer->FixVariable(21+i);
  if(!g_rotY)
    for(int i = 0; i < 7; i++)
      minimizer->FixVariable(28+i);
  if(!g_rotZ)
    for(int i = 0; i < 7; i++)
      minimizer->FixVariable(35+i);

  // no transY on X-planes
  if(g_tranY)
    for(int i = 0; i < 7; i++)
      if (i == 0 || i == 5 || i == 6)
        minimizer->FixVariable(7+i);

  // tie Z-position of Y-planes together
  if(g_tranZ)
    for(int i = 0; i < 7; i++)
      if (i == 2 || i == 4)
        minimizer->FixVariable(14+i);

  minimizer->Minimize();

  double out_tranX[7];
  double out_tranX_err[7];
  double out_tranY[7];
  double out_tranY_err[7];
  double out_tranZ[7];
  double out_tranZ_err[7];

  double out_rotX[7];
  double out_rotX_err[7];
  double out_rotY[7];
  double out_rotY_err[7];
  double out_rotZ[7];
  double out_rotZ_err[7];
  
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
    if (i == 2 || i == 4){
      out_tranZ[i] = param[14+i-1];
      out_tranZ_err[i] = error[14+i-1];
    }
    else {
      out_tranZ[i] = param[14+i];
      out_tranZ_err[i] = error[14+i];
    }
  }
  for(int i = 0; i < 7; i++){
    out_rotX[i] = param[21+i];
    out_rotX_err[i] = error[21+i];
  }
  for(int i = 0; i < 7; i++){
    out_rotY[i] = param[28+i];
    out_rotY_err[i] = error[28+i];
  }
  for(int i = 0; i < 7; i++){
    out_rotZ[i] = param[35+i];
    out_rotZ_err[i] = error[35+i];
  }

  TFile* fout = new TFile(outputFileName, "RECREATE");
  fout->cd();

  TTree* out_tree = new TTree("AlignmentTree",
			      "AlignmentTree");
  
  cout << "Final Parameters:" << endl;
  for(int i = 0; i < 7; i++){
    out_tree->Branch(Form("tranX_%d",i+1), &out_tranX[i]);
    out_tree->Branch(Form("tranX_err_%d",i+1), &out_tranX_err[i]);
    if(g_tranX)
      cout << "tran X [mm] " << i+1 << ": " << out_tranX[i] << " +/- " << out_tranX_err[i] << endl;
  }
  for(int i = 0; i < 7; i++){
    out_tree->Branch(Form("tranY_%d",i+1), &out_tranY[i]);
    out_tree->Branch(Form("tranY_err_%d",i+1), &out_tranY_err[i]);
    if(g_tranY)
      cout << "tran Y [mm] " << i+1 << ": " << out_tranY[i] << " +/- " << out_tranY_err[i] << endl;
  }
  for(int i = 0; i < 7; i++){
    out_tree->Branch(Form("tranZ_%d",i+1), &out_tranZ[i]);
    out_tree->Branch(Form("tranZ_err_%d",i+1), &out_tranZ_err[i]);
    if(g_tranZ)
      cout << "tran Z [mm] " << i+1 << ": " << out_tranZ[i] << " +/- " << out_tranZ_err[i] << endl;
  }
  for(int i = 0; i < 7; i++){
    out_tree->Branch(Form("rotX_%d",i+1), &out_rotX[i]);
    out_tree->Branch(Form("rotX_err_%d",i+1), &out_rotX_err[i]);
    if(g_rotX)
      cout << "rot X [rad] " << i+1 << ": " << out_rotX[i] << " +/- " << out_rotX_err[i] << endl;
  }
  for(int i = 0; i < 7; i++){
    out_tree->Branch(Form("rotY_%d",i+1), &out_rotY[i]);
    out_tree->Branch(Form("rotY_err_%d",i+1), &out_rotY_err[i]);
    if(g_rotY)
      cout << "rot Y [rad] " << i+1 << ": " << out_rotY[i] << " +/- " << out_rotY_err[i] << endl;
  }
  for(int i = 0; i < 7; i++){
    out_tree->Branch(Form("rotZ_%d",i+1), &out_rotZ[i]);
    out_tree->Branch(Form("rotZ_err_%d",i+1), &out_rotZ_err[i]);
    if(g_rotZ)
      cout << "rot Z [rad] " << i+1 << ": " << out_rotZ[i] << " +/- " << out_rotZ_err[i] << endl;
  }

  out_tree->Fill();
  fout->cd();
  out_tree->Write();
  fout->Close();   
}

double EvaluateMetric(const double* param){
  GeoOctuplet GEO;
  GEO.SetRunNumber(g_RunNum);
  // apply parameters to geometry
  if(g_tranX)
    for(int i = 0; i < 7; i++)
      GEO.TranslateX(param[i],i+1);
  if(g_tranY)
    for(int i = 0; i < 7; i++)
      GEO.TranslateY(param[i+7],i+1);
  if(g_tranZ)
    for(int i = 0; i < 7; i++){
      if (i == 2 || i == 4)
        GEO.TranslateZ(param[i+14-1] ,i+1);
      else
        GEO.TranslateZ(param[i+14] ,i+1);
    }
  if(g_rotX)
    for(int i = 0; i < 7; i++)
      GEO.RotateX(param[i+21],i+1);
  if(g_rotY)
    for(int i = 0; i < 7; i++)
      GEO.RotateY(param[i+28],i+1);
  if(g_rotZ)
    for(int i = 0; i < 7; i++)
      GEO.RotateZ(param[i+35],i+1);

  double res, sum_res2 = 0;
  int evts = (int)(g_clusters_all.size());
  for(int evt = 0; evt < evts; evt++){
    if(evt%10000 == 0)
      cout << "Processing event # " << evt << " | " << evts << endl;
    MMTrack track = g_FITTER->Fit(*g_clusters_all[evt], GEO);
    res = GEO.GetQuadraticSumOfResidualsX(*g_clusters_all[evt], track, false);
    sum_res2 += res;
  }
  sum_res2 = sum_res2 / evts;

  if(g_tranX)
    for(int i = 0; i < 7; i++)
      cout << "tranX " << i+1 << ": " << param[i] << endl;
  if(g_tranY)
    for(int i = 0; i < 7; i++)
      cout << "tranY " << i+1 << ": " << param[i+7] << endl;
  if(g_tranZ)
    for(int i = 0; i < 7; i++)
      cout << "tranZ " << i+1 << ": " << param[i+14] << endl; 
  if(g_rotX)
    for(int i = 0; i < 7; i++)
      cout << "rotX " << i+1 << ": " << param[i+21] << endl;
  if(g_rotY)
    for(int i = 0; i < 7; i++)
      cout << "rotY " << i+1 << ": " << param[i+28] << endl;
  if(g_rotZ)
    for(int i = 0; i < 7; i++)
      cout << "rotZ " << i+1 << ": " << param[i+35] << endl;

  cout << "myFCN = " << sum_res2 << endl;
  return sum_res2;
}
