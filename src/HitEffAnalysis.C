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
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "TMath.h"
#include "TColor.h"
#include "TLegend.h"

#include "include/MMPlot.hh"
#include "include/PDOToCharge.hh"
#include "include/TDOToTime.hh"
#include "include/MMDataAnalysis.hh"
#include "include/MMPacmanAlgo.hh"
#include "include/GeoOctuplet.hh"
#include "include/SimpleTrackFitter.hh"
#include "include/HighQTrackFitter.hh"
#include "include/CombinatoricTrackFitter.hh"
#include "include/CombinatoricChi2TrackFitter.hh"
#include "include/ScintillatorClusterFilterer.hh"
#include <sys/stat.h>
#include <sys/param.h>
#include <unistd.h>
using namespace std;

int main(int argc, char* argv[]){
  int runnum;

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

  //MMPacmanAlgo* PACMAN = new MMPacmanAlgo(5,5.,0.5);
  MMPacmanAlgo* PACMAN = new MMPacmanAlgo(5,2.,0.5);

  GeoOctuplet* GEOMETRY = new GeoOctuplet();
  if(b_align)
    GEOMETRY->SetAlignment(AlignFileName);

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

  int max_EventNum =  T->GetMaximum("mm_EventNum");

  DATA = (MMDataAnalysis*) new MMDataAnalysis(T);

  int Nevent = DATA->GetNEntries();

  ScintillatorClusterFilterer* FILTER = new ScintillatorClusterFilterer();
  SimpleTrackFitter* FITTER = new SimpleTrackFitter();

  // Hits!

  int HT = 7; //hit threshold
  int NEventsHT = 0;
  TH1D* passedhits;
  TH1D* totalhits;
  TH1D* hitEff;
  TH1D* track_angles;
  vector<TH1D*> hDeltaX;
  TH1D* nBoardsHit;
  vector<TH1D*> nBoardsHit_anglebin;
  vector<TH1D*> targetBoardClustM;

  vector<TH1D*> passedhits_anglebin;
  vector<TH1D*> totalhits_anglebin;
  vector<TH1D*> hitEff_anglebin;
  vector<vector<TH1D*>> hDeltaX_anglebin;

  TH2D* nBoards_vs_track_angle;


  double minCut [8] = {5.,5.,5.,5.,5.,5.,5.,5.};
  double track_angle_bins[7] = {0.,5.,10.,15.,20.,25.,90.};
  //double track_angle_bins[7] = {0.,2.5,5.,7.5,10.,12.5,90.};

  int nangles = sizeof(track_angle_bins)/sizeof(*track_angle_bins);

  passedhits = new TH1D("passed_hits","passed_hits",8, -0.5,7.5);
  totalhits = new TH1D("total_hits","total_hits",8, -0.5,7.5);
  hitEff = new TH1D("hitEff","hitEff",8, -0.5,7.5);
  nBoards_vs_track_angle = new TH2D("nBoards_vs_track_angle","nBoards_vs_track_angle",
    5, 3.5,8.5,
    50,0.,5.);


  for (int i = 0; i < 8; i++){
      hDeltaX.push_back(new TH1D(Form("hDeltaX_%d",i),
           Form("hDeltaX_%d",i),
           100,-15,15.));
      targetBoardClustM.push_back(new TH1D(Form("targetBoardClustM_%d",i),
        Form("targetBoardClusM_%d",i),
        4,-0.5,4.5));
      }
  for (int i = 0; i < nangles-1; i++){
      passedhits_anglebin.push_back(new TH1D(Form("passedhits_anglebin_%d",i),
        Form("passedhits_anglebin_%d",i),
        8,-0.5,7.5)) ;
      totalhits_anglebin.push_back(new TH1D(Form("totalhits_anglebin_%d",i),
        Form("totalhits_anglebin_%d",i),
        8,-0.5,7.5)) ;
      hitEff_anglebin.push_back(new TH1D(Form("hitEff_anglebin_%d",i),
        Form("hitEff_anglebin_%d",i),
        8,-0.5,7.5)) ;
      vector<TH1D*> hDeltaX_anglebin_temp;
      for (int j = 0; j < 8; j++)
        hDeltaX_anglebin_temp.push_back(new TH1D(Form("hDeltaX_%d_%d",j,i),
        Form("hDeltaX_%d_%d",j,i),
        100,-15.,15.));
      hDeltaX_anglebin.push_back(hDeltaX_anglebin_temp);

      passedhits_anglebin[i]->Sumw2();
      totalhits_anglebin[i]->Sumw2();
      nBoardsHit_anglebin.push_back(new TH1D(Form("nBoardsHit_%d",i),
        Form("nBoardsHit_%d",i),
        9,-0.5,8.5));
  }
  //hDeltaX = new TH1D("hDeltaX","hDeltaX",250,-50.,50.);
  nBoardsHit = new TH1D("nBoardsHit","nBoardsHit",9,-0.5,8.5);
  //nBoardsHit = new TH1D("nBoardsHit","nBoardsHit",5,3.5,8.5);

  track_angles = new TH1D("track angles","track angles",50, -25.,25.);

  passedhits->Sumw2();
  totalhits->Sumw2();

  TFile* fout = new TFile(outputFileName, "RECREATE");
  fout->mkdir("event_displays");
  MMPlot();

  //Nevent /= 100;


  for(int evt = 0; evt < Nevent; evt++){
    DATA->GetEntry(evt);
    if (evt == 0)
      runnum = DATA->RunNum;
    if(evt%(Nevent/1000) == 0) 
      cout << "Processing event # " << evt << " | " << Nevent << endl;
    //    if (evt == 10000) break;
    //cout << "Nevent " << evt << endl;
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

    // initialize cluster making
    FILTER->SetRunNumber(DATA->RunNum);

    // book histograms for MM hits
    int Nboard = DATA->mm_EventHits.GetNBoards();
    
    vector<MMClusterList> all_clusters;
    for(int i = 0; i < Nboard; i++){
      if(DATA->mm_EventHits[i].GetNHits() == 0)
	    continue;
      
      MMClusterList clusters = PACMAN->Cluster(DATA->mm_EventHits[i]);
      if(clusters.GetNCluster() > 0)
        all_clusters.push_back(clusters);
    }
    
    int Ncl = all_clusters.size();
    for(int i = 0; i < Ncl; i++){
      int Nc = all_clusters[i].GetNCluster();
      for(int j = 0; j < Nc; j++){
	     const MMCluster& clus = all_clusters[i][j];
	     int b = ib[clus.MMFE8()];
      }
    }
    // track fitting
    MMClusterList fit_clusters;
    for(int i = 0; i < Ncl; i++){
      // add all clusters from each board
      int Nc = all_clusters[i].GetNCluster();
      for(int c = 0; c < Nc; c++)
	       // if(all_clusters[i][c].GetNHits() > 1)
	       fit_clusters.AddCluster(all_clusters[i][c]);
    }

    int Nclus_all = fit_clusters.GetNCluster();

    int nB = 0;

    bool BoardHit[8] = {false,false,false,false,false,false,false,false};
    int missing_board = -1;
    for (int c = 0; c < Nclus_all; c++){
        int b = ib[fit_clusters[c].MMFE8()]; //taking from a list of the highest charge clusters
        BoardHit[b] = true;
    }
    for (int i = 0; i < 8; i++){
      if (BoardHit[i] == true)
        nB++;
      else
        missing_board= i;
    }
    nBoardsHit->Fill(nB);
    if ((nB) < 4)
      continue;
    // checking nboard distribution as a function of angle
    pair < SCHit*,SCHit* > botPair;
    botPair = DATA->sc_EventHits.GetBotPair()[0];
    MMClusterList filtered_allclusters = FILTER->FilterClustersScint(fit_clusters, *GEOMETRY, botPair.first->Channel(), evt);

    MMTrack track = FITTER->Fit(filtered_allclusters, *GEOMETRY);

    // check if track is okay
    if(!track.IsFit() ||
      track.NX() < 2 ||
      track.NU()+track.NV() < 2 ||
      track.NX()+track.NU()+track.NV() < 5)
      continue;

    double precut_track_theta_deg = atan(track.SlopeX())/TMath::Pi()*180.;
    nBoards_vs_track_angle->Fill(nB,precut_track_theta_deg);
    int precut_anglebin = -1;
    for (int j = 0; j < nangles-1; j++)
      if (fabs(precut_track_theta_deg) >= track_angle_bins[j] 
          && fabs(precut_track_theta_deg) < track_angle_bins[j+1] )
            precut_anglebin = j;
    nBoardsHit_anglebin[precut_anglebin]->Fill(nB);

    //cout << "N clus: " << Nclus_all << endl;
    if(nB < HT)
      continue;
    //continue;
    //cout << "Got min clusters" << endl;
    NEventsHT++;

    // saving efficiency clusters for each hit eff calculation
    vector < MMClusterList> hit_clusters;
    vector <int> denom_ints;
    MMClusterList hit_clusters_temp;
    for (int i = 0; i < 8; i++){ // calculating hit efficiency for board i
      // if there's a missing board, we can only calculate the hit eff for that one board!
      if (missing_board > -1) 
        if (missing_board != i) // this should always just be 0 for HT = 8
          continue;
      hit_clusters_temp.Reset();
      for (int c = 0; c < Nclus_all; c++){
          int b = ib[fit_clusters[c].MMFE8()];
          if (b != i) { // if not target board, add cluster
            //cout << "BOARD: "<< b << endl;
            hit_clusters_temp.AddCluster(fit_clusters[c]);
          }
      }
      hit_clusters.push_back(hit_clusters_temp); // save this cluster for fitting
      denom_ints.push_back(i); // save which board we need to check
    }
    MMClusterList filtered_clusters;
    MMClusterList target_clusters;

    std::vector<double> vx;
    TVector3 p;
    bool outBound = false;
    //pair < SCHit*,SCHit* > botPair;
    int clustM[8] = {0,0,0,0,0,0,0,0};

    for (int k = 0; k < denom_ints.size(); k++){
      int BID = denom_ints[k]; // for board we're interested in for the hit eff
      // fit track
      filtered_clusters.Reset();
      //botPair = DATA->sc_EventHits.GetBotPair()[0];
      filtered_clusters = FILTER->FilterClustersScint(hit_clusters[k], *GEOMETRY, botPair.first->Channel(), evt);

      MMTrack track_all = FITTER->Fit(filtered_clusters, *GEOMETRY);

      // check if track is okay
      if(!track_all.IsFit() ||
        track_all.NX() < 2 ||
        track_all.NU()+track_all.NV() < 2 ||
        track_all.NX()+track_all.NU()+track_all.NV() < 4)
        continue;

      // track angles
      double track_theta = atan(track_all.SlopeX());
      double track_theta_deg = atan(track_all.SlopeX())/TMath::Pi()*180.;
      //track_angles->Fill(track_theta_deg);

      // save all clusters for target hit eff board
      target_clusters.Reset();
      for(int i = 0; i < Ncl; i++){
        int Nc = all_clusters[i].GetNCluster();
        for(int j = 0; j < Nc; j++){
          if(all_clusters[i][j].GetNHits() > 0){
            if (ib[all_clusters[i][j].MMFE8()] == BID)
              target_clusters.AddCluster(all_clusters[i][j]);
          }
        }
      }

      GeoPlane plane = GEOMETRY->Get(BID);
      p = plane.Intersection(track_all);

      // fiducial cut

      if ((p.X() > 200.) || (p.X() < 0.)) 
        continue;
      if ((p.Y() > 217.9) || (p.Y() < 17.9)) 
        continue;
      // if ((p.X() > 215.) || (p.X() < -15.)) 
      //   continue;
      // if ((p.Y() > 233.) || (p.Y() < 2.)) 
      //   continue;
      int anglebin = -1;
      for (int j = 0; j < nangles-1; j++)
        if (fabs(track_theta_deg) >= track_angle_bins[j] 
            && fabs(track_theta_deg) < track_angle_bins[j+1] )
            anglebin = j;
      //cout << "foudn angle bin: " << anglebin << endl;
      totalhits->Fill(BID);
      totalhits_anglebin[anglebin]->Fill(BID);

      bool targetBoardHit[8] = {false,false,false,false,false,false,false,false};

      for (int i = 0; i < target_clusters.GetNCluster(); i++){
        if (ib[target_clusters[i].MMFE8()] == BID){
          double deltaX = plane.GetResidualX(target_clusters[i].Channel(),track_all);
          hDeltaX[BID]->Fill(deltaX);
          hDeltaX_anglebin[anglebin][BID]->Fill(deltaX);
          if ((deltaX < minCut[k]) && targetBoardHit[BID] == false){
            //cout << "Filling hist" << endl;
            passedhits->Fill(BID);
            passedhits_anglebin[anglebin]->Fill(BID);
            targetBoardHit[BID] = true;
            clustM[BID]++;
            continue;
          }          
        }
      }
    }
      for (int k = 0; k  < 8; k++)
        targetBoardClustM[k]->Fill(clustM[k]);
      //cout << "made it past filling" << endl;

  }

  // make directory
  string dirname = "HitEfficiencyPlots_";
  cout << runnum << endl;
  dirname += to_string(runnum);
    
  mkdir(dirname.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  TCanvas* c1 = new TCanvas("c1","",600,400);
  gStyle->SetTextFont(42);
  int ngauss = 3; 
  for (int i = 0; i < 8; i++) { 
    //TF1 *f1 = new TF1("gaussian_1","gaus(0)",-2.4,-2.4);
    c1->cd();
    c1->Clear();
    c1->SetLogy();
    TF1 *f2 = nullptr;
    if (ngauss == 1){
      f2 = new TF1("gaussian and constant","pol0(0)+gaus(1)",-14.99,14.99);
      f2->SetParameter(0, 0.);
      f2->SetParameter(1, hDeltaX[i]->GetMaximum());
      f2->SetParameter(2, hDeltaX[i]->GetMean());
      f2->SetParameter(3, hDeltaX[i]->GetRMS()/2);     
    }
    else if (ngauss == 2){
      f2 = new TF1("double gaussian and constant","pol0(0)+gaus(1)+gaus(4)",-14.99,14.99);
      f2->SetParameter(0, 0.);
      f2->SetParLimits(0, 0.1, 5.);
      f2->SetParameter(1, hDeltaX[i]->GetMaximum());
      f2->SetParameter(2, hDeltaX[i]->GetMean());
      f2->SetParameter(3, hDeltaX[i]->GetRMS()/2);
      f2->SetParameter(4, hDeltaX[i]->GetMaximum()/2/1.5);
      f2->SetParameter(5, hDeltaX[i]->GetMean());
      f2->SetParameter(6, hDeltaX[i]->GetRMS()*1.5);       
    }
    else if (ngauss == 3){
      f2 = new TF1("triple gaussian and constant","pol0(0)+gaus(1)+gaus(4)+gaus(7)",-14.99,14.99);
      f2->SetParameter(0, 0.);
      //f2->SetParLimits(0, 0.1, 5.);
      f2->SetParLimits(0, 0., 100.);
      f2->SetParameter(1, hDeltaX[i]->GetMaximum());
      f2->SetParameter(2, hDeltaX[i]->GetMean());
      f2->SetParameter(3, hDeltaX[i]->GetRMS()/2);
      f2->SetParameter(4, hDeltaX[i]->GetMaximum()/2.5/1.5);
      f2->SetParameter(5, hDeltaX[i]->GetMean());
      f2->SetParameter(6, hDeltaX[i]->GetRMS()*2);      
      if (i!=7){
        f2->SetParameter(7, hDeltaX[i]->GetMaximum()/6/1.5);
        f2->SetParameter(8, hDeltaX[i]->GetMean());
        f2->SetParameter(9, hDeltaX[i]->GetRMS()*4);   
      }
      else{
        f2->SetParameter(7, hDeltaX[i]->GetMaximum()/100);
        f2->SetParameter(8, hDeltaX[i]->GetMean());
        f2->SetParameter(9, hDeltaX[i]->GetRMS()*6);
      }
    }
    hDeltaX[i]->Fit(f2,"RQN");
//     cout << "First 9 fit parameters: " << f2->GetParameter(0) << " " 
//     << f2->GetParameter(1) << " " << f2->GetParameter(2) << " " 
//          << f2->GetParameter(3) << " " << f2->GetParameter(4) << " " 
//          << f2->GetParameter(5) << " " << f2->GetParameter(6) << " "
//          << f2->GetParameter(7) << " " << f2->GetParameter(8) << " " 
//          << f2->GetParameter(9) << endl;
    TF1 *f1 = new TF1("pol0","pol0(0)",-14.9,14.9);
    f1->SetParameter(0,f2->GetParameter(0));
    double backgroundEvents = f1->Integral(-5.,5.);
    cout << "Nevents bkg: " << backgroundEvents << endl;
    double tempEvents = passedhits->GetBinContent(i+1);
    cout << "Nevents pre cut: " << tempEvents << endl;

    passedhits->SetBinContent(i+1,tempEvents-backgroundEvents);
    cout << "Nevents post cut: " << passedhits->GetBinContent(i+1) << endl;

    // plot beautification
    hDeltaX[i]->SetTitle(Form("x Residuals for Board %d",i));
    hDeltaX[i]->GetXaxis()->SetTitleOffset(1.3);
    hDeltaX[i]->GetYaxis()->SetTitleOffset(1.3);
    hDeltaX[i]->GetXaxis()->SetLabelSize(0.045);
    hDeltaX[i]->GetXaxis()->SetTitleSize(0.045);    
    hDeltaX[i]->GetYaxis()->SetLabelSize(0.045);
    hDeltaX[i]->GetYaxis()->SetTitleSize(0.045);
    hDeltaX[i]->GetXaxis()->SetTitle("x_{cluster} - x_{track, proj.} [mm]");
    hDeltaX[i]->GetYaxis()->SetTitle("Number of Clusters");

    //hDeltaX[i]->GetXaxis()->CenterTitle();
    hDeltaX[i]->Draw("pe");
    f2->SetLineStyle(7);
    f2->SetLineColor(kBlue);
    f1->SetLineStyle(1);
    f1->SetLineColor(kRed);
    f2->Draw("same");
    f1->Draw("same");
    TString temp = dirname;
    temp.Append(Form("/const_gauss_fit_%d.pdf",i));
    c1->Print(temp);

    // Do this bin by bin for angle bins
//     cout << "SIZE: " << sizeof(track_angle_bins)/sizeof(*track_angle_bins) << endl;
//     for (int j = 0; j < nangles-1; j++){
//       TF1* f3 = nullptr;
//       if (ngauss == 3){
//        f3 = new TF1("triple gaussian and constant","pol0(0)+gaus(1)+gaus(4)+gaus(7)",-14.99,14.99);
//         f3->SetParameter(0, 0.5);
//         f3->SetParLimits(0,0.,1000.);
//         f3->SetParameter(1, hDeltaX[i]->GetMaximum());
//         f3->SetParameter(2, hDeltaX[i]->GetMean());
//         f3->SetParameter(3, hDeltaX[i]->GetRMS()/2);
//         f3->SetParameter(4, hDeltaX[i]->GetMaximum()/2/1.5);
//         f3->SetParameter(5, hDeltaX[i]->GetMean());
//         f3->SetParameter(6, hDeltaX[i]->GetRMS()*1.5);      
//         f3->SetParameter(7, hDeltaX[i]->GetMaximum()/4/1.5);
//         f3->SetParameter(8, hDeltaX[i]->GetMean());
//         f3->SetParameter(9, hDeltaX[i]->GetRMS()*3);   
//       }
//       cout << "Fitting board " << i << " anglebin " << j << endl;
//       hDeltaX_anglebin[j][i]->Fit(f3, "RQ");
//       TF1 *f4 = new TF1("pol0","pol0(0)",-14.9,14.9);
//       f4->SetParameter(0,f3->GetParameter(0));
//       backgroundEvents = f4->Integral(-5.,5.);
//       cout << "background estimation: " << backgroundEvents << endl;
//       tempEvents = passedhits_anglebin[j]->GetBinContent(i+1);
//       passedhits_anglebin[j]->SetBinContent(i+1,tempEvents-backgroundEvents);
//     }
  }

  hitEff->Divide(passedhits,totalhits,1.,1.,"B");
  for (int i = 1; i < 9; i++){
    cout << "Efficiency for Board " << i-1 << ", " << hitEff->GetBinContent(i) << endl;
  }
  for (int i = 0; i < 5; i++){
    hitEff_anglebin[i]->Divide(passedhits_anglebin[i],totalhits_anglebin[i],1.,1.,"B");
  }

  cout << "/////////////////////////////////" << endl;
  cout << "EVENT SUMMARY: " << endl;
  cout << "Events: " << NEventsHT << endl;
  cout << "/////////////////////////////////" << endl;

  fout->cd();
  fout->mkdir("histograms");
  fout->cd("histograms");
    track_angles->Write();
//  x_hits->Write();
  passedhits->Write();
  totalhits->Write();
  for (int i = 0; i < 8; i++){
    hDeltaX[i]->Write();
    targetBoardClustM[i]->Write();
  }
  nBoardsHit->Write();
  for (int i = 0; i < nangles-1; i++)
    nBoardsHit_anglebin[i]->Write();
  nBoards_vs_track_angle->Write();
  hitEff->Write();
  TCanvas* c2 = new TCanvas("c2","",600,400);
  vector<int> colors = {kOrange+6,kPink+3,kGreen+2,kCyan-6,kAzure+7,kViolet-8};

  gStyle->SetErrorX(0.);
  c2->cd();
  c2->Clear();
  hitEff->SetMarkerColor(kBlue-7);
  hitEff->SetMarkerStyle(20);
  hitEff->SetMinimum(0.);
  hitEff->Draw("E1 P");
  hitEff->GetXaxis()->SetTitle("Board Number");
  //hitEff->GetXaxis()->CenterTitle();
  hitEff->GetYaxis()->SetTitle("Hit Efficiency");
  hitEff->SetMaximum(1.);
  hitEff->SetMinimum(0.);
  //hitEff->GetYaxis()->CenterTitle();
  hitEff->GetXaxis()->SetTitleOffset(1.3);
  hitEff->GetYaxis()->SetTitleOffset(1.3);
  hitEff->GetXaxis()->SetLabelSize(0.045);
  hitEff->GetXaxis()->SetTitleSize(0.045);    
  hitEff->GetYaxis()->SetLabelSize(0.045);
  hitEff->GetYaxis()->SetTitleSize(0.045);
  TString temp = dirname;
  temp.Append("/HitEff.pdf");
  c2->Print(temp);
  c2->Clear();


  //

  TLegend* leg = new TLegend(0.8,0.3,0.94,0.5);
  hitEff_anglebin[0]->SetMarkerStyle(20);
  hitEff_anglebin[0]->SetMinimum(0.3);
  hitEff_anglebin[0]->SetMaximum(1.);
  hitEff_anglebin[0]->GetXaxis()->SetTitle("Board Number");
  //hitEff->GetXaxis()->CenterTitle();
  hitEff_anglebin[0]->GetYaxis()->SetTitle("Hit Efficiency");
  //hitEff->GetYaxis()->CenterTitle();
  hitEff_anglebin[0]->GetXaxis()->SetTitleOffset(1.3);
  hitEff_anglebin[0]->GetYaxis()->SetTitleOffset(1.3);
  hitEff_anglebin[0]->GetXaxis()->SetLabelSize(0.045);
  hitEff_anglebin[0]->GetXaxis()->SetTitleSize(0.045);    
  hitEff_anglebin[0]->GetYaxis()->SetLabelSize(0.045);
  hitEff_anglebin[0]->GetYaxis()->SetTitleSize(0.045);
  hitEff_anglebin[0]->Draw("E1 P");
  for (int i = 0; i < nangles-3; i++){
    hitEff_anglebin[i]->SetMarkerColor(colors[i]);
    hitEff_anglebin[i]->SetMarkerStyle(20);
    hitEff_anglebin[i]->Draw("E1 PL same");
    leg->AddEntry(hitEff_anglebin[i],Form("%d#circ #leq #theta < %d#circ",int(track_angle_bins[i]),int(track_angle_bins[i+1])));
  }
  leg->Draw();
  temp = dirname;
  temp.Append("/HitEff_angles.pdf");
  c2->Print(temp);
  c2->Clear();

  fout->cd();
  fout->Close();
    
}
