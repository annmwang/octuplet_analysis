///
///  \file   MakeClusterTree.C
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Oct
///

#include <iostream>
#include <stdlib.h>

#include "include/PDOToCharge.hh"
#include "include/TDOToTime.hh"
#include "include/MMDataAnalysis.hh"
#include "include/MMPacmanAlgo.hh"
#include "include/GeoOctuplet.hh"
#include "include/SimpleTrackFitter.hh"

using namespace std;

int main(int argc, char* argv[]){

  char inputFileName[400];
  char outputFileName[400];
  char PDOFileName[400];
  char TDOFileName[400];
  
  if ( argc < 5 ){
    cout << "Error at Input: please specify input/output .root files ";
    cout << " and (optional) PDO/TDO calibration files" << endl;
    cout << "Example:   ./MakeClusterTree.x -i input.root -o output.root" << endl;
    cout << "Example:   ./MakeClusterTree.x -i input.root -o output.root";
    cout << " -p PDOcalib.root -t TDOcalib.root" << endl;
    return 0;
  }

  bool b_input = false;
  bool b_out   = false;
  bool b_pdo   = false;
  bool b_tdo   = false;
  for (int i=1;i<argc-1;i++){
    if (strncmp(argv[i],"-i",2)==0){
      sscanf(argv[i+1],"%s", inputFileName);
      b_input = true;
    }
    if (strncmp(argv[i],"-o",2)==0){
      sscanf(argv[i+1],"%s", outputFileName);
      b_out = true;
    }
    if (strncmp(argv[i],"-p",2)==0){
      sscanf(argv[i+1],"%s", PDOFileName);
      b_pdo = true;
    }
    if (strncmp(argv[i],"-t",2)==0){
      sscanf(argv[i+1],"%s", TDOFileName);
      b_tdo = true;
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

  // board ID's
  map<int,int> ib;
  vector<int> iboards;
  
  PDOToCharge* PDOCalibrator;
  if(b_pdo)
    PDOCalibrator = new PDOToCharge(PDOFileName);
  else
    PDOCalibrator = new PDOToCharge();

  TDOToTime* TDOCalibrator;
  if(b_tdo)
    TDOCalibrator = new TDOToTime(TDOFileName);
  else
    TDOCalibrator = new TDOToTime();

  MMPacmanAlgo* PACMAN = new MMPacmanAlgo();

  SimpleTrackFitter* FITTER = new SimpleTrackFitter();

  GeoOctuplet* GEOMETRY = new GeoOctuplet();

  MMDataAnalysis* DATA;
  TFile* f = new TFile(inputFileName, "READ");
  if(!f){
    cout << "Error: unable to open input file " << inputFileName << endl;
    return false;
  }
  TTree* T = (TTree*) f->Get("COMB_data");
  if(!T){
    cout << "Error: cannot find tree COMB_data in " << inputFileName << endl;
    return false;
  }

  DATA = (MMDataAnalysis*) new MMDataAnalysis(T);

  int Nevent = DATA->GetNEntries();
  
  TFile* fout = new TFile(outputFileName, "RECREATE");
  fout->cd();
  int N_clus;
  vector<int> clus_MMFE8;
  vector<int> clus_Index;
  vector<double> clus_Charge;
  vector<double> clus_Time;
  vector<double> clus_Channel;

  TTree* cluster_tree = new TTree("ClusterTree",
				  "ClusterTree");
  cluster_tree->Branch("N_clus", &N_clus);
  cluster_tree->Branch("clus_MMFE8", &clus_MMFE8);
  cluster_tree->Branch("clus_Index", &clus_Index);
  cluster_tree->Branch("clus_Charge", &clus_Charge);
  cluster_tree->Branch("clus_Time", &clus_Time);
  cluster_tree->Branch("clus_Channel", &clus_Channel);

  for(int evt = 0; evt < Nevent; evt++){
       DATA->GetEntry(evt);
    if(evt%10000 == 0)
      cout << "Processing event # " << evt << " | " << Nevent << endl;

    if(GEOMETRY->RunNumber() < 0){
      GEOMETRY->SetRunNumber(DATA->RunNum);

      iboards = GEOMETRY->MMFE8list();
      for(int i = 0; i < iboards.size(); i++)
	ib[iboards[i]] = i;
    }

    if(!DATA->sc_EventHits.IsGoodEvent())
      continue;
    
    // Calibrate PDO -> Charge
    PDOCalibrator->Calibrate(DATA->mm_EventHits);
    // Calibrate TDO -> Time
    TDOCalibrator->Calibrate(DATA->mm_EventHits);
  
    // initialize PACMAN info for this event
    PACMAN->SetEventTrigBCID(DATA->mm_trig_BCID);
    PACMAN->SetEventPadTime(0); // add this

    int Nboards = DATA->mm_EventHits.GetNBoards();
    
    vector<MMClusterList> all_clusters;
    for(int i = 0; i < Nboards; i++){
      if(DATA->mm_EventHits[i].GetNHits() == 0)
	continue;
      
      MMClusterList clusters = PACMAN->Cluster(DATA->mm_EventHits[i]);
      if(clusters.GetNCluster() > 0)
	all_clusters.push_back(clusters);
    }

    N_clus = 0;
    clus_MMFE8.clear();
    clus_Index.clear();
    clus_Charge.clear();
    clus_Time.clear();
    clus_Channel.clear();

    MMClusterList fit_clusters;

    int Ncl = all_clusters.size();
    // write clusters to tree
    for(int i = 0; i < Ncl; i++){
      // add highest charge cluster from each board;
      int Nc = all_clusters[i].GetNCluster();
      for(int c = 0; c < Nc; c++){
	if(all_clusters[i][c].GetNHits() > 1){
	  fit_clusters.AddCluster(all_clusters[i][c]);
	  N_clus++;
	  clus_MMFE8.push_back(all_clusters[i][c].MMFE8());
	  clus_Index.push_back(ib[all_clusters[i][c].MMFE8()]);
	  clus_Charge.push_back(all_clusters[i][c].Charge());
	  clus_Time.push_back(all_clusters[i][c].Time());
	  clus_Channel.push_back(all_clusters[i][c].Channel());
	  break;
	}
      }
    }

    if(N_clus < 5)
      continue;
    
    MMTrack track_all = FITTER->Fit(fit_clusters, *GEOMETRY);
    double sumresX2 = 0.;
    for(int i = 0; i < N_clus; i++){
      double resX = GEOMETRY->GetResidualX(fit_clusters[i], track_all);
      sumresX2 += resX*resX;
    }

    if( sumresX2 > double(N_clus-2)*0.1)
      continue;

    if(N_clus >= 5)
      cluster_tree->Fill(); 

  } // end event loop

  fout->cd();
  cluster_tree->Write();
  fout->Close();
    
}
