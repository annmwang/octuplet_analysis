//
///  \file   HitEffAnalysis.C
///
///  \author Ann Miao Wang
///          (anwang@cern.ch)
///
///  \date   2016 Nov
///
///  \note   Studying hit efficiency of the octuplet

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include <iostream>
#include <stdlib.h>
#include "TROOT.h"

//#include "include/AtlasStyle.hh"
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
  char AlignFileName[400];
  
  if ( argc < 5 ){
    cout << "Error at Input: please specify input/output .root files ";
    cout << " and (optional) PDO/TDO calibration files" << endl;
    cout << "Example:   ./HitEffAnalysis.x -i input.root -o output.root" << endl;
    cout << "Example:   ./HitEffAnalysis.x -i input.root -o output.root";
    cout << " -p PDOcalib.root -t TDOcalib.root" << endl;
    cout << "Other options:" << endl;
    cout << "   -p PDOcalib.root : PDO calibration file" << endl;
    cout << "   -t TDOcalib.root : TDO calibration file" << endl;
    cout << "   -a alignment.root : alignment file" << endl;
    return 0;
  }

  bool b_input = false;
  bool b_out   = false;
  bool b_pdo   = false;
  bool b_tdo   = false;
  bool b_align   = false;
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
    if (strncmp(argv[i],"-a",2)==0){
      sscanf(argv[i+1],"%s", AlignFileName);
      b_align = true;
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

  GeoOctuplet* GEOMETRY = new GeoOctuplet();
  if(b_align)
    GEOMETRY->SetAlignment(AlignFileName);

  SimpleTrackFitter* FITTER = new SimpleTrackFitter();

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

  // histograms for track analysis
  vector<TH1D*> board_track_x;
  vector<TH1D*> board_track_y;
  vector<TH1D*> board_itrack_resX;
  vector<TH1D*> board_otrack_resX;
  vector<TH1D*> track_deg_top21;

  vector<TH1D*> track_sumresX2;
  vector<TH1D*> track_FitProb;

  // counting events with single clusters
  int single_cluster[] = {0,0,0,0,0,0,0,0};

  double x_res = 0.18; // for now

  // Hits hits hits
  int mm_evt = 0;
  TEfficiency* pEff;
  TH1D* passedhits;
  TH1D* totalhits;

  vector<TH1D*> x_hits;

  int five_cluster = 0;
  int fit_five_cluster = 0;
  int cut_five_cluster = 0;

  int gothit[8] = {0};
  int exphit[8] = {0};
  vector < vector<Double_t> > hiteff;
  vector < vector<Double_t> > hitefferrors;
  int HT = 5; // boards hit thresh

  for(int i = 0; i < 6; i++){ 
    track_deg_top21.push_back(new TH1D(Form("track_deg_top21_bot%d",i),
				       Form("track_deg_top21_bot%d",i),
				       100,-25.,25.));
  }
  for(int i = 0; i < 8; i++){
    

    board_track_x.push_back(new TH1D(Form("b_t_X_%d",i),
				     Form("b_t_X_%d",i),
				     1001,-50.5,250.5));
    board_track_y.push_back(new TH1D(Form("b_t_Y_%d",i),
				     Form("b_t_Y_%d",i),
				     1001,-50.5,250.5));
    board_itrack_resX.push_back(new TH1D(Form("b_it_resX_%d",i),
					 Form("b_it_resX_%d",i),
					 1001,-5.,5.));
    board_otrack_resX.push_back(new TH1D(Form("b_ot_resX_%d",i),
					 Form("b_ot_resX_%d",i),
					 1001,-5.,5.));
    track_sumresX2.push_back(new TH1D(Form("t_sumresX2_%d",i),
              Form("t_sumresX2_%d",i),
              1024,0.,5.));
    track_FitProb.push_back(new TH1D(Form("t_FitProb_%d",i),
              Form("t_FitProb_%d",i),
              1024,0.,1.));
    x_hits.push_back(new TH1D(Form("x_hits_%d",i),
              Form("x_hits_%d",i),
              200,-50.5,250.5));

  }

  passedhits = new TH1D("passed_hits","passed_hits",8, -0.5,7.5);
  totalhits = new TH1D("total_hits","total_hits",8, -0.5,7.5);
  TFile* fout = new TFile(outputFileName, "RECREATE");
  fout->mkdir("event_displays");
  MMPlot();

  for(int evt = 0; evt < Nevent; evt++){
    DATA->GetEntry(evt);
    mm_evt = DATA->mm_EventNum;

    if(evt%(Nevent/10) == 0) 
      //cout << "Processing event # " << evt << " | " << Nevent << endl;

    if(GEOMETRY->RunNumber() < 0){
      GEOMETRY->SetRunNumber(DATA->RunNum);

      iboards = GEOMETRY->MMFE8list();
      for(int i = 0; i < iboards.size(); i++)
      	ib[iboards[i]] = i;
    }
    //cout << "Event: " << mm_evt << endl;
    if(!DATA->sc_EventHits.IsGoodEvent())
      continue;
    //cout << "passec scint" << endl;

    // Calibrate PDO -> Charge
    PDOCalibrator->Calibrate(DATA->mm_EventHits);
    // Calibrate TDO -> Time
    TDOCalibrator->Calibrate(DATA->mm_EventHits);

    // initialize PACMAN info for this event
    PACMAN->SetEventTrigBCID(DATA->mm_trig_BCID);
    PACMAN->SetEventPadTime(0); // add this

    // histogram group A - no BCID cuts, no "strict" calibration cuts
    // book histograms for MM hits
    int Nboard = DATA->mm_EventHits.GetNBoards();
    for(int i = 0; i < Nboard; i++){
      int Nhit = DATA->mm_EventHits[i].GetNHits();
      for(int j = 0; j < Nhit; j++){
    	const MMLinkedHit& hit = DATA->mm_EventHits[i][j];
    	// hit has duplicate
      	int Ndup = hit.GetNHits();
      	if(Ndup > 1){
      	  const MMLinkedHit* phit = &hit;
      	  while(phit){
    	    phit = phit->GetNext();
      	  }
      	}
      }
    }

    // Calculate scint information
    int NScHits = DATA->sc_EventHits.GetNHits();
    
    pair < SCHit*,SCHit* > topPair;
    pair < SCHit*,SCHit* > botPair;
    topPair = DATA->sc_EventHits.GetTopPair()[0]; 
    botPair = DATA->sc_EventHits.GetBotPair()[0]; 

    int shifted_top = topPair.first->Channel() - 16 - 1; // shifted top scintillator!
    double deltaX = -20.*(shifted_top - botPair.first->Channel())+8.; 
    double deltaZ = 274.3;
    double sc_theta = atan(deltaX/deltaZ); 
    double sc_theta_deg = atan(deltaX/deltaZ)/TMath::Pi()*180.; 

    // histogram group B
    // clustering applies BCID cuts, strips away hits with values of zero for MMFE8,PDO,TDO,Q,calib
    int TotalMMHits = 0;

    vector<MMClusterList> all_clusters;
    for(int i = 0; i < Nboard; i++){
      if(DATA->mm_EventHits[i].GetNHits() == 0)
    	continue;
      
      MMClusterList clusters = PACMAN->Cluster(DATA->mm_EventHits[i]);
      TotalMMHits += PACMAN->GetGoodHits();
      if (clusters.GetNCluster() == 1) 
        single_cluster[ib[DATA->mm_EventHits[i].MMFE8()]]++;
      if(clusters.GetNCluster() > 0)
      	all_clusters.push_back(clusters);
    }
    if (TotalMMHits < 1)
      continue;
    int Ncl = all_clusters.size();
    // vector<MMClusterList> pruned_clusters;
    // // prune cluster lists using "slope roads"
    // for(int i = 0; i < Nboard; i++){
    //   MMClusterList temp_clusters;
    //   GeoPlane plane = GEOMETRY->Get(i);
    //   int Nc = all_clusters[i].GetNCluster();
    //   for(int c = 0; c < Nc; c++) {
    //     if(all_clusters[i][c].GetNHits() > 0){
    //       double x_pos = 102.3+plane.LocalXatYbegin(all_clusters[i][c].Channel());
    //       if ((x_pos < 115.) && botPair.first->Channel() == 4 && topPair.first->Channel() == 21)
    //         temp_clusters.AddCluster(all_clusters[i][c]);
    //     }
    //   }
    //   pruned_clusters.push_back(temp_clusters);
    // }





    MMClusterList aggr_clusters;
    for(int i = 0; i < Ncl; i++){
      // add together all clusters
      int Nc = all_clusters[i].GetNCluster();
      for(int c = 0; c < Nc; c++) {
        if(all_clusters[i][c].GetNHits() > 0){
          aggr_clusters.AddCluster(all_clusters[i][c]);
        }
      }
    }
    for (int i = 0; i < 8; i++){
      //if (botPair.first->Channel() == 4 && topPair.first->Channel() == 20) {
        GeoPlane plane = GEOMETRY->Get(i);
        for (int j = 0; j < aggr_clusters.GetNCluster(); j++) {
          if (ib[aggr_clusters.Get(j).MMFE8()] == i) {
            double x_pos = 102.3+plane.LocalXatYbegin(aggr_clusters.Get(j).Channel());
            //cout << "Board: " << i << " Position: " << x_pos << endl;
            x_hits[i]->Fill(102.3+plane.LocalXatYbegin(aggr_clusters.Get(j).Channel()));
          }
        }
      //}
    }

    for(int i = 0; i < Ncl; i++){
      int Nc = all_clusters[i].GetNCluster();
    }

/*    // track fitting
    MMClusterList fit_clusters;
    MMClusterList fit_clusters_a; // first five boards only
    for(int i = 0; i < Ncl; i++){
      // add highest charge cluster from each board;
      int Nc = all_clusters[i].GetNCluster();
      for(int c = 0; c < Nc; c++) {
      	if(all_clusters[i][c].GetNHits() > 0){
      	  fit_clusters.AddCluster(all_clusters[i][c]);
          if (ib[all_clusters[i][c].MMFE8()] < 5) {
            fit_clusters_a.AddCluster(all_clusters[i][c]);
          }
      	  break;
        }
      }
    }
    int Nclus_all_a = fit_clusters_a.GetNCluster();
    */

    // track fitting with all possible clusters
    
    // here we save all possible cluster combinations
    int cluster_n[8] = {0,0,0,0,0,0,0,0}; // iterator through clusters
    int listsize = 0; // keep track of how many boards have != 0 clusters

    vector<MMClusterList> fit_clusters_a; 
    for (int i = 0; i < all_clusters.size(); i++) {
      if (all_clusters[i].GetNCluster() > 0){
        cluster_n[i] = 0;
        listsize++; 
      }
      else
        cluster_n[i] = -1; // don't loop through if no clusters        
    }

    if (listsize < 1) // if nothing, stop
      continue;
    
    // cout << "SIZE: " << listsize << endl;
    
    // add all possible combinations of clusters
    while (cluster_n[0] != all_clusters[0].GetNCluster()){
      MMClusterList fit_clusters;
      for (int j = 0; j < listsize; j++){
        fit_clusters.AddCluster(all_clusters[j][cluster_n[j]]);
        // cout << "BOARD: " << j << " CLUSTER # " << cluster_n[j] << endl;
      }
      fit_clusters_a.push_back(fit_clusters);
      cluster_n[listsize-1]++;
      for (int i = listsize-1; (i > 0) && (cluster_n[i] == (all_clusters[i].GetNCluster())); i--) {
        cluster_n[i] = 0;
        cluster_n[i-1]++;
      }
    }

    // how many combinations?
    int Nfit = fit_clusters_a.size();
    cout << "COMBINATIONS: " << Nfit << endl;

    // let's not explode
    if (Nfit > 10000)
      continue;

    int Nclus_all_a = fit_clusters_a[0].GetNCluster();
    
    // For comparing with Paolo
    //cout << "nev 3 " << mm_evt << " " << Nclus_all_a << " " <<endl;
    //if (mm_evt == 20000)
    //  break;

    int Nplane = GEOMETRY->GetNPlanes();

    // here we only fit boards if we have HT or more clusters
    

    if (Nclus_all_a < HT) 
      continue;
    five_cluster++;

    // if (Nclus_all < HT) 
    //   continue;
    // saving clusters 
    /*vector < MMClusterList> hit_clusters;
    vector <int> denom_ints; // save the boards for which we have eff. denominator candidates
    //MMClusterList hit_clusters_temp;
    int other_hits = 0;

   for (int i = 0; i < Nplane; i++){ // calculating hit efficiency for board i
      MMClusterList hit_clusters_temp;
      for (int j = 0; j < Nplane; j++){ //loop through all other planes
      	if (i == j) continue;
      	bool had_cluster = false;
      	for (int c = 0; c < Nclus_all; c++){ // look for cluster on that board
      	  int b = ib[fit_clusters[c].MMFE8()];
      	  if (b == j) 
      	    had_cluster = true;
      	}
      	if (had_cluster == true) {
      	  other_hits++;
      	  had_cluster = false;
      	}
      }
      if (other_hits >= HT){
      	for (int c = 0; c < Nclus_all; c++){
  	      int b = ib[fit_clusters[c].MMFE8()]; //taking from a list of the highest charge clusters
      	  if (b != i) {
      	    hit_clusters_temp.AddCluster(fit_clusters[c]);
      	  }
      	}
      	hit_clusters.push_back(hit_clusters_temp); // save this fucking piece of shit
      	denom_ints.push_back(i);
      }
      other_hits = 0;
    }*/
      
    // picking the cluster combination that gives you the best track
  
    vector <double> list_sumres;
    vector <MMTrack> list_tracks;
    int imin_res = -1;
    double min_res = 10000000000.;
    double sumresX2;
    for (int k = 0; k < fit_clusters_a.size(); k++) {
      MMTrack track_all = FITTER->Fit(fit_clusters_a[k], *GEOMETRY);
      // Calculate the chi2 for the track, using the "other board's highest charge clusters"
      sumresX2 = 0.;
      for(int c = 0; c < HT; c++){
        const MMCluster& clus = fit_clusters_a[k][c];
        // fill residuals
        double resX = GEOMETRY->GetResidualX(clus, track_all);
        sumresX2 += resX*resX;
      }
      list_sumres.push_back(sumresX2);
      list_tracks.push_back(track_all);
      if (sumresX2 < min_res){
        min_res = sumresX2;
        imin_res = k;
      }
    }
    // Uncomment this if you want to see what cluster it picks
/*    cout << "MINRES: " << min_res << " " << imin_res << endl;
    for (int i = 0; i < fit_clusters_a[imin_res].GetNCluster(); i++){
      int boardi = ib[fit_clusters_a[imin_res].Get(i).MMFE8()];
      cout << "BOARD: " << boardi << endl;
      GeoPlane plane = GEOMETRY->Get(boardi);
      double tempx = 102.3+plane.LocalXatYbegin(fit_clusters_a[imin_res].Get(i).Channel());
      cout << "XPOS: " << tempx << endl;
    } */     
//    if (mm_evt > 600)
//      break;
    
    MMTrack track_all = list_tracks[imin_res]; 
    
    // track angles
    double track_theta = atan(track_all.SlopeX());
    double track_theta_deg = atan(track_all.SlopeX())/TMath::Pi()*180.;
    for (int i = 0; i < 6; i++) { 
      if (botPair.first->Channel() == i && topPair.first->Channel() == 21) {
  	    cout << botPair.first->Channel() << " SC angle: " << sc_theta_deg << " Track angle: " << track_theta_deg << endl;
  	    track_deg_top21[i]->Fill(track_theta_deg);
      }
    }


    cout << "SC angle: " << sc_theta_deg << " Track angle: " << track_theta_deg << endl;

    double fitProb = TMath::Prob(min_res/pow(x_res,2),HT-4);
    track_sumresX2[0]->Fill(min_res/pow(x_res,2));
    track_FitProb[0]->Fill(fitProb); //assuming HT-1 degrees of freedom

    // Fiducial Cut, projected onto board 7 
    std::vector<double> vx;
    TVector3 p;
    GeoPlane plane = GEOMETRY->Get(7);
    p = plane.Intersection(track_all);


    // if ((p.X() > 215.) || (p.X() < -15.)) continue;
    // if ((p.Y() > 233.) || (p.Y() < 2.)) continue;

    int fid_cut = 1;
    if (((p.X() > 200.) || (p.X() < 0.)) || (p.Y() > 217.9) || (p.Y() < 17.9))
      fid_cut = 0;
    //cout << "nev " << mm_evt << " fitprob=   " << fitProb << " acc= " << fid_cut << endl;

    if (fitProb < pow(10,-5.)) continue;
    fit_five_cluster++;

    if (fid_cut == 0.)
      continue;
    cut_five_cluster++;

    // past denominator selection

    for (int k = 5; k < Nplane; k++){
      exphit[k]++;
      // save all clusters for target hit eff board
      MMClusterList target_clusters;
      for(int i = 0; i < Ncl; i++){
        int Nc = all_clusters[i].GetNCluster();
        for(int j = 0; j < Nc; j++){
        if(all_clusters[i][j].GetNHits() > 0){
          if (ib[all_clusters[i][j].MMFE8()] == k)
            target_clusters.AddCluster(all_clusters[i][j]);
          }
        }
      }

      double min_res = 999999.;
      // calculate the residuals for ALL clusters in the hit eff board
      for(int c = 0; c < target_clusters.GetNCluster(); c++){
        const MMCluster& clus = target_clusters[c];
        // fill residuals
        double resX = GEOMETRY->GetResidualX(clus, track_all);
        if (fabs(resX) < fabs(min_res)) {
          min_res = resX;
        }
      }
      // save the cluster with the smallest residual
      if (min_res < 999999.)
        board_itrack_resX[k]->Fill(min_res);
    }


/*    for (int k = 0; k < denom_ints.size(); k++){

      // ID of the board we're interested in measuring the hit efficiency of
      int BID = denom_ints[k];

      // fit the track
      MMTrack track_all = FITTER->Fit(hit_clusters[k], *GEOMETRY);
      
      // save all clusters for target hit eff board
      MMClusterList target_clusters;
      for(int i = 0; i < Ncl; i++){
    	  int Nc = all_clusters[i].GetNCluster();
      	for(int j = 0; j < Nc; j++){
    	  if(all_clusters[i][j].GetNHits() > 0){
    	    if (ib[all_clusters[i][j].MMFE8()] == BID)
    	      target_clusters.AddCluster(all_clusters[i][j]);
      	  }
      	}
      }

      // Calculate the chi2 for the track, using the "other board's highest charge clusters"
      double sumresX2 = 0.;
      for(int c = 0; c < HT; c++){
      	const MMCluster& clus = hit_clusters[k][c];
      	// fill residuals
      	double resX = GEOMETRY->GetResidualX(clus, track_all);
      	sumresX2 += resX*resX;
      }

      double fitProb = TMath::Prob(sumresX2/pow(x_res,2),HT-4);
      track_sumresX2[BID]->Fill(sumresX2/pow(x_res,2));
      // cout << "Event: " << evt<< endl;
      // cout << "Filled hist! " << endl;
      track_FitProb[BID]->Fill(fitProb); //assuming HT-1 degrees of freedom

      // cut to remove bad fits
      if (fitProb > pow(10,-5.)) continue;
      exphit[BID]++;

      double min_res = 999999.;
      // calculate the residuals for ALL clusters in the hit eff board
      for(int c = 0; c < target_clusters.GetNCluster(); c++){
    	  const MMCluster& clus = target_clusters[c];
      	// fill residuals
      	double resX = GEOMETRY->GetResidualX(clus, track_all);
      	if (fabs(resX) < fabs(min_res)) {
      	  min_res = resX;
      	}
      }
      // save the cluster with the smallest residual
      if (min_res < 999999.)
        board_itrack_resX[BID]->Fill(min_res);
    }*/
  }
  
  /*vector < TF1 > sidebands;
  for (int i = 5; i < 8; i++) { 
    TF1 *f1 = new TF1("f1","pol0",-5.,-2.5);
    TF1 *f2 = new TF1("f2","pol0",2.5,5.);
    board_itrack_resX[i]->Fit("f1","R");
    board_itrack_resX[i]->Fit("f2","R");
    double f = 0.5*(f1->GetParameter(0) + f2->GetParameter(0));
    double cfit = board_itrack_resX[i]->GetNbinsX()*f;
    cout << "Exp events: " << exphit[i] << endl;
    cout << "Subtracted events: " << cfit << endl;
    gothit[i] = board_itrack_resX[i]->GetEntries()-cfit;
    cout << "Full events " << board_itrack_resX[i]->GetEntries() << endl;
    passedhits->SetBinContent(i,gothit[i]);
    cout << "GOT HIT: " << gothit[i] << endl;
    totalhits->SetBinContent(i,exphit[i]);
  }
  AtlasStyle* test = new AtlasStyle();
  test->SetAtlasStyle();
  //  gStyle->SetErrorX(0.0);
  pEff = new TEfficiency(*passedhits,*totalhits);
  pEff->SetStatisticOption(TEfficiency::kFNormal);
  pEff->SetMarkerStyle(4);

  
  for (int i = 0; i < 8; i++){
    double eff = double(gothit[i])/double(exphit[i]);
  }
  
  cout << "/////////////////////////////////" << endl;
  cout << "EVENT SUMMARY: " << endl;
  cout << "Event candidates: " << five_cluster <<" After fit: " << fit_five_cluster << " After cuts: " <<cut_five_cluster <<  endl;
  cout << "/////////////////////////////////" << endl;
  for (int i =0; i<8; i++){
      cout << "Board " << i << ": Exp Hits " << exphit[i] << endl;
      cout << "Board " << i << ": Got Hits " << gothit[i] << endl;
      cout << "Board " << i << ": Hit % " << double(gothit[i])/double(exphit[i]) << endl;
   }
  fout->cd();
  fout->mkdir("histograms");
  for(int i = 0; i < 8; i++){
    fout->cd("histograms");
    board_track_x[i]->Write();
    board_track_y[i]->Write();
    board_itrack_resX[i]->Write();
    board_otrack_resX[i]->Write();
  }

  for(int i = 0; i < 6; i++){
    track_deg_top21[i]->Write();
  }
  
  for(int i = 0; i < 8; i++){
    fout->cd("histograms");
    track_sumresX2[i]->Write();
    track_FitProb[i]->Write();
    x_hits[i]->Write();
  }
//  x_hits->Write();
  passedhits->Write();
  totalhits->Write();
  pEff->SetDirectory(gDirectory);
  pEff->Write();

  TCanvas* c1 = new TCanvas("c1","",600,400);
  TGraphAsymmErrors* gr = pEff->CreateGraph("");
  gr->SetMarkerColor(kBlue-7);

  gr->SetPointEXlow(0,0.0);
  gr->SetPointEXhigh(0,0.0);
  gr->SetPointEXlow(1,0.0);
  gr->SetPointEXhigh(1,0.0);
  gr->SetPointEXlow(2,0.0);
  gr->SetPointEXhigh(2,0.0);
  gr->Draw("AP");
  gr->GetXaxis()->SetTitle("Board Number");
  gr->GetYaxis()->SetTitle("Hit Efficiency");
  c1->Print("HitEff_5.pdf");*/

  fout->Close();
    
}
