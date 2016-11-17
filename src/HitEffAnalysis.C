//
///  \file   HitEffAnalysis_v2.C
///
///  \author Ann
///          (anwang@cern.ch)
///
///  \date   2016 Nov
///
///  \note   Studying hit efficiency

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TEfficiency.h"
#include <iostream>
#include <stdlib.h>

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

  vector<TH1D*> track_sumresX2;
  vector<TH1D*> track_FitProb;

  // counting events with single clusters
  int single_cluster[] = {0,0,0,0,0,0,0,0};

  double x_res = 0.18; // for now

  // gothit/exphit = hit efficiency for board assuming 7 other boards hit
  
  // Hits hits hits
  TEfficiency* pEff;
  TH1D* passedhits;
  TH1D* totalhits;

  int gothit[8] = {0};
  int exphit[8] = {0};
  vector < vector<Double_t> > hiteff;
  vector < vector<Double_t> > hitefferrors;
  int HT = 5; // boards hit thresh
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
  }

  passedhits = new TH1D("passed hits","passed hits",8, -0.5,8.5);
  totalhits = new TH1D("total hits","total hits",8, -0.5,8.5);

  TFile* fout = new TFile(outputFileName, "RECREATE");
  fout->mkdir("event_displays");
  MMPlot();

  for(int evt = 0; evt < Nevent; evt++){
    DATA->GetEntry(evt);
    if(evt%(Nevent/10) == 0) 
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

    // histogram group B
    // clustering applies BCID cuts, strips away hits with values of zero for MMFE8,PDO,TDO,Q,calib
    vector<MMClusterList> all_clusters;
    for(int i = 0; i < Nboard; i++){
      if(DATA->mm_EventHits[i].GetNHits() == 0)
    	continue;
      
      MMClusterList clusters = PACMAN->Cluster(DATA->mm_EventHits[i]);
      if (clusters.GetNCluster() == 1) 
        single_cluster[ib[DATA->mm_EventHits[i].MMFE8()]]++;
      if(clusters.GetNCluster() > 0)
      	all_clusters.push_back(clusters);
    }

    int Ncl = all_clusters.size();

    // track fitting
    MMClusterList fit_clusters;
    for(int i = 0; i < Ncl; i++){
      // add highest charge cluster from each board;
      int Nc = all_clusters[i].GetNCluster();
      for(int c = 0; c < Nc; c++) {
      	if(all_clusters[i][c].GetNHits() > 0){
      	  fit_clusters.AddCluster(all_clusters[i][c]);
      	  break;
        }
    	}
    }

    int Nclus_all = fit_clusters.GetNCluster();
    int Nplane = GEOMETRY->GetNPlanes();

    // here we only fit boards if we have 7 or more clusters
    
    if (Nclus_all < HT) 
      continue;

    // saving clusters 
    vector < MMClusterList> hit_clusters;
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
    }

    for (int k = 0; k < denom_ints.size(); k++){

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
    }
  }
  
  vector < TF1 > sidebands;
  for (int i = 0; i < 8; i++) { 
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

  pEff = new TEfficiency(*passedhits,*totalhits);
  
  for (int i = 0; i < 8; i++){
    double eff = double(gothit[i])/double(exphit[i]);
  }
  
  // TCanvas *canHitEff = new TCanvas("Hit Eff","Hit Efficiencies",700,500);
  // TGraph *gr[8];
  // //  TGraphErrors *gr[8];
  // TMultiGraph *mg = new TMultiGraph();
  // Int_t Nc = 10;
  // Double_t fC_cut[10] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  // //  Double_t fC_cut_errors[10] = {0.};
  // TLegend *leg = new TLegend(0.7,0.2,0.8,0.4);
  // for (int i = 0; i < 8; i++){
  //   gr[i] = new TGraph(Nc,fC_cut,&(hiteff[i][0]));
  //   //    gr[i] = new TGraphErrors(Nc,fC_cut,&(hiteff[i][0]),fC_cut_errors,&(hitefferrors[i][0]));
  //   gr[i]->SetMarkerColor(i+2);
  //   //    gr[i]->SetMarkerStyle(21);
  //   mg->Add(gr[i]);
  //   leg->AddEntry(gr[i],Form("Board %d",i),"p");
  // }
  // mg->GetXaxis()->SetTitle("Cluster Charge Min (fC)");
  // mg->GetYaxis()->SetTitle("Efficiency");
  // mg->SetMaximum(1.);
  // mg->SetMinimum(0.);
  // mg->Draw("APL");
  // leg->Draw();
  // canHitEff->Print(Form("HitEfficiency_%d.pdf",HT));
  // delete canHitEff;
  cout << "/////////////////////////////////" << endl;
  cout << "EVENT SUMMARY: " << endl;
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

  for(int i = 0; i < 8; i++){
    fout->cd("histograms");
    track_sumresX2[i]->Write();
    track_FitProb[i]->Write();
  }
  passedhits->Write();
  totalhits->Write();
  pEff->Write();

  TCanvas* c1 = new TCanvas("c1","",600,400);
  pEff->Draw("AP");
  c1->Print("HitEff.pdf");

  fout->Close();
    
}
