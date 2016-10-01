///
///  \file   RunMMAnalysisTemplate.C
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///

#include "TH1D.h"
#include "TH2D.h"
#include <iostream>

#include "include/MMPlot.hh"
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
    cout << "Example:   ./RunMMAnalysisTemplate.x -i input.root -o output.root" << endl;
    cout << "Example:   ./RunMMAnalysisTemplate.x -i input.root -o output.root";
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

  // PDO calibration object
  PDOToCharge* PDOCalibrator;
  if(b_pdo)
    PDOCalibrator = new PDOToCharge(PDOFileName);
  else
    PDOCalibrator = new PDOToCharge();

  // TDO calibration object
  TDOToTime* TDOCalibrator;
  if(b_tdo)
    TDOCalibrator = new TDOToTime(TDOFileName);
  else
    TDOCalibrator = new TDOToTime();

  // clustering algorithm object
  MMPacmanAlgo* PACMAN = new MMPacmanAlgo();

  // Octuplet geometry object
  GeoOctuplet* GEOMETRY = new GeoOctuplet();

  // track fitting object
  SimpleTrackFitter* FITTER = new SimpleTrackFitter();

  // data object
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

  // open output file
  TFile* fout = new TFile(outputFileName, "RECREATE");
  // set style for plotting
  MMPlot();

  for(int evt = 0; evt < Nevent; evt++){
       DATA->GetEntry(evt);
    if(evt%10000 == 0)
      cout << "Processing event # " << evt << " | " << Nevent << endl;

    // tell the octuplet geometry object
    // the run number if not already done
    if(GEOMETRY->RunNumber() < 0)
      GEOMETRY->SetRunNumber(DATA->RunNum);

    // DATA->sc_EventHits (SCEventHits class) is the collection
    // of scintillator hits (SCHit class) for the event

    // check scint hit requirements for good event
    if(!DATA->sc_EventHits.IsGoodEvent())
      continue;

    // DATA->mm_EventHits (MMEventHits class) is the collection
    // of MM hits (MMHit class) for the event
    
    // Calibrate PDO -> Charge
    PDOCalibrator->Calibrate(DATA->mm_EventHits);
    // Calibrate TDO -> Time
    TDOCalibrator->Calibrate(DATA->mm_EventHits);
  
    // initialize PACMAN info for this event
    PACMAN->SetEventTrigBCID(DATA->mm_trig_BCID);

    // how many duplicate hits in the event
    // (number of hits with at least 1 dup)
    int Ndup_evt = DATA->mm_EventHits.GetNDuplicates();
    
    // Loop through boards with hits in event
    int Nboard = DATA->mm_EventHits.GetNBoards();
    for(int i = 0; i < Nboard; i++){
      int Ndup_board = DATA->mm_EventHits[i].GetNDuplicates();
      
      // Loop through hits on the board
      int Nhit = DATA->mm_EventHits[i].GetNHits();
      for(int j = 0; j < Nhit; j++){
	double PDO  = DATA->mm_EventHits[i][j].PDO();
	double TDO  = DATA->mm_EventHits[i][j].TDO();
	double chan = DATA->mm_EventHits[i][j].Channel();
	// ... plus other methods - see include/MMHit.hh

	// can also grab a reference to a hit for easier
	// typing:
	const MMLinkedHit& hit = DATA->mm_EventHits[i][j];
	
	// how many duplicates are there for this hit?
	int Ndup_hit = hit.GetNHits();
      }
    }

    // a vector of lists of clusters
    vector<MMClusterList> all_clusters;
    for(int i = 0; i < Nboard; i++){
      // ignore empty boards
      if(DATA->mm_EventHits[i].GetNHits() == 0)
	continue;

      // cluster a board's collection of hits into clusters
      MMClusterList clusters = PACMAN->Cluster(DATA->mm_EventHits[i]);

      // if there's at least one cluster
      // in list save it 
      if(clusters.GetNCluster() > 0)
	all_clusters.push_back(clusters);
    }

    // how many non-empty cluster lists?
    int Ncl = all_clusters.size();
    for(int i = 0; i < Ncl; i++){
      // how many clusters in a list?
      int Nc = all_clusters[i].GetNCluster();
      for(int j = 0; j < Nc; j++){
	// get a reference to cluster [i][j]
	const MMCluster& clus = all_clusters[i][j];
	double Q  = clus.Charge();
	double T  = clus.Time();
	int Nhole = clus.NHoles();
	int Ndup_clus = clus.GetNDuplicates();
	// ....
      }
    }

    // track fitting
    MMClusterList fit_clusters;
    for(int i = 0; i < Ncl; i++){
      // add highest charge cluster from each board
      // (sorted by charge in MMClusterList)
      if(all_clusters[i].GetNCluster() > 0)
	fit_clusters.AddCluster(all_clusters[i][0]);
    }

    // save some event displays for those with
    // at least 1 cluster on each board
    if(fit_clusters.GetNCluster() < 8)
      continue;

    // find best-fit track using list of clusters and
    // an octuplet geometry object
    MMTrack track = FITTER->Fit(fit_clusters, *GEOMETRY);

    // 2D plot with U/V error bars
    TCanvas* can = Plot_Track2D(Form("track2D_%d",DATA->mm_EventNum), track, *GEOMETRY, &fit_clusters); 
    fout->cd();
    can->Write();
    delete can;

    // 2D plot with U/V track shift
    TCanvas* canY = Plot_Track2DY(Form("track2DY_%d",DATA->mm_EventNum), track, *GEOMETRY, &fit_clusters); 
    fout->cd();
    canY->Write();
    delete canY;

    // 3D plot
    TCanvas* can3D = Plot_Track3D(Form("track3D_%d",DATA->mm_EventNum), track, *GEOMETRY, &fit_clusters); 
    fout->cd();
    can3D->Write();
    delete can3D;
  }
  
  fout->Close();
    
}
