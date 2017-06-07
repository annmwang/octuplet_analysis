//
///  \file   TPAnalysis.C
///
///  \author Ann Miao Wang
///          (anwang@cern.ch)
///
///  \date   2017 May
///
///  \note   Studying events from the trigger processor

#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
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

using namespace std;

int deltaBCID(int sciBCID, int sciph, int hitBCID, int offset){
  if ((sciBCID-hitBCID) < 0){
    return sciBCID+sciph/16.-hitBCID+offset;
  }
  else
    return sciBCID+sciph/16.-hitBCID+(offset-4096);
}

int main(int argc, char* argv[]){

  char inputFileName[400];
  char outputFileName[400];
  char PDOFileName[400];
  char TDOFileName[400];
  char AlignFileName[400];
  
  if ( argc < 5 ){
    cout << "Error at Input: please specify input/output .root files ";
    cout << " and (optional) PDO/TDO calibration files" << endl;
    cout << "Example:   ./TPAnalysis.x -i input.root -o output.root" << endl;
    cout << "Example:   ./TPAnalysis.x -i input.root -o output.root";
    cout << " -p PDOcalib.root -t TDOcalib.root" << endl;
    cout << "Other options:" << endl;
    cout << "   -p PDOcalib.root : PDO calibration file" << endl;
    cout << "   -t TDOcalib.root : TDO calibration file" << endl;
    cout << "   -a alignment.root : alignment file" << endl;
    return 0;
  }
  bool debug = false;
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

  // colors
  string pink = "\033[38;5;205m";
  string green = "\033[38;5;84m";
  string blue = "\033[38;5;27m";
  string end = "\033[0m";
  string warning = "\033[38;5;227;48;5;232m";

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

  MMPacmanAlgo* PACMAN = new MMPacmanAlgo(5,2.,0.5);
  cout << pink << "Using cluster seed 2 fC, 0.5 fC else" << end << endl;
  
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


  DATA = (MMDataAnalysis*) new MMDataAnalysis(T);

  int Nevent = DATA->GetNEntries();

  ScintillatorClusterFilterer* FILTER = new ScintillatorClusterFilterer();
  SimpleTrackFitter* FITTER = new SimpleTrackFitter();

  // Hits!

  int dB = 1944;
  if (debug)
    cout << "Assumed constant: " << dB << endl;

  int NEvents = 0;
  int NEventsGood = 0;
  int NEventsTrig = 0;
  int NEventsTrigGoodTP = 0;
  int NEventsTrigGoodTrack = 0;
  int NEventsTrigScintMatch = 0;
  
  TH1D* TPcands;
  TH1D* TPmatch;
  TH1D* TPscimatch;
  TH1D* effTP;
  TH1D* fracTPscimatch;

  TH1D* nomatch_timediff;
  TH1D* all_timediff;
  TH1D* packetBCID;
  
  TH1D* all_MMmatch;
  TH1D* nomatch_MMmatch;
  TH1D* match_MMmatch;


  TH1D* dEarlyBCID;
  TH1D* dEarlyBCIDph;
  TH1D* dAvgBCID;
  TH1D* dAvgBCIDph;
  TH1D* dMedBCID;
  TH1D* dMedBCIDph;


  TH1D* dEarlyBCIDfull;
  TH1D* dAvgBCIDfull;
  TH1D* earliestBCID;
  TH1D* avgBCID;

  // tp_sci + 1945 = tp_BCID[0]

  TPcands = new TH1D("Candidate events","Candidate events",200,1494.5*pow(10,6),1495.2*pow(10,6));
  TPmatch = new TH1D("Trigger candidate events with tp","Trigger candidate events with tp",200,1494.5*pow(10,6),1495.2*pow(10,6));
  effTP = new TH1D("Fraction of trigger cands with tp","Fraction of trigger cands with tp",200,1494.5*pow(10,6),1495.2*pow(10,6));
  TPscimatch = new TH1D("Events with tp and tpsci BCID match","Events with tp and tpsci BCID match",200,1494.5*pow(10,6),1495.2*pow(10,6));
  fracTPscimatch = new TH1D("Fraction of events with tp and tpsci BCID match","Fraction of events with tp and tpsci BCID match",200,1494.5*pow(10,6),1495.2*pow(10,6));

  nomatch_timediff = new TH1D("Timestamp difference for those without the scintillator match",
                              "Timestamp difference for those without the scintillator match",
                              100,-0.01,0.06);
  all_timediff = new TH1D("Timestamp difference for all candidate events with tp",
                              "Timestamp difference for all candidate events with tp",
                              100,-0.01,0.06);

  nomatch_MMmatch = new TH1D("Number of TP hits matched with MM for those without the scintillator match",
                              "Number of TP hits matched with MM for those without the scintillator match",
                              9,-0.5,8.5);
  match_MMmatch = new TH1D("Number of TP hits matched with MM for those with the scintillator match",
                           "Number of TP hits matched with MM for those with the scintillator match",
                           9,-0.5,8.5);
  all_MMmatch = new TH1D("Number of TP hits matched with MM for all candidate events with tp",
                              "Number of TP hits matched with MM for all candidate events with tp",
                              9,-0.5,8.5);

  packetBCID = new TH1D("Delta tp BCID","Delta tp BCID",10001, -5000,5000);

  dEarlyBCID = new TH1D("Delta Earliest BCID","Delta Earliest BCID",13, -6.5,6.5);
  dEarlyBCIDph = new TH1D("Delta Earliest BCID phase","Delta Earliest BCID phase",208, -6.5,6.5);

  dAvgBCID = new TH1D("Delta Avg BCID","Delta Avg BCID",13, -6.5,6.5);
  dAvgBCIDph = new TH1D("Delta Avg BCID phase","Delta Avg BCID phase",208, -6.5,6.5);

  dMedBCID = new TH1D("Delta Med BCID","Delta Med BCID",13, -6.5,6.5);
  dMedBCIDph = new TH1D("Delta Med BCID phase","Delta Med BCID phase",208, -6.5,6.5);

  earliestBCID = new TH1D("earliest BCID","earliest BCID",5000, -1.,10000.);
  avgBCID = new TH1D("avg BCID","avg BCID",5000, -1.,10000.);

  dEarlyBCIDfull = new TH1D("delta earliest BCID full","delta earliest BCID full",5000, -10000.,10000.);
  dAvgBCIDfull = new TH1D("delta avg BCID full","delta avg BCID full",5000, -10000.,10000.);

  TFile* fout = new TFile(outputFileName, "RECREATE");


  for(int evt = 0; evt < Nevent; evt++){
    DATA->GetEntry(evt);
    // if (evt > 100)
    //   break
    if (debug)
      cout << "Event " << evt << endl;
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
    NEvents++;

    if (DATA->mm_EventHits.EventNum() > 300001)
      break;
    
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

    int NX1 = 0;
    int NX2 = 0;
    int NUV = 0;

    bool Gothits = false;
    vector<MMClusterList> all_clusters;
    for(int i = 0; i < Nboard; i++){
      int ib = DATA->mm_EventHits[i].MMFE8Index();
      if (ib < 2)
        NX1++;
      else if (ib > 5)
        NX2++;
      else
        NUV++;
      if(DATA->mm_EventHits[i].GetNHits() == 0)
	      continue;
      else
	      Gothits = true;
      
      MMClusterList clusters = PACMAN->Cluster(DATA->mm_EventHits[i]);
      if(clusters.GetNCluster() > 0)
        all_clusters.push_back(clusters);
    }

    int Ncl = all_clusters.size();
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

    
    pair < SCHit*,SCHit* > botPair;
    botPair = DATA->sc_EventHits.GetBotPair()[0];
    MMClusterList filtered_allclusters = FILTER->FilterClustersScint(fit_clusters, *GEOMETRY, botPair.first->Channel(), evt);
    
    MMTrack track = FITTER->Fit(filtered_allclusters, *GEOMETRY);

    for (int i = 0; i < filtered_allclusters.size(); i++){
      const MMCluster& clus = filtered_allclusters[i];
    }

    if (!(NX1>0 && NX2> 0 && NUV > 0))
      continue;
//     if (!track.IsTrigCand()){
//       //    if (!track.IsFit()||
//       //!track.IsTrigCand()){
//       if (debug)
//         cout << pink << "Event " << evt << " didn't make the cut" << end << endl;
//       continue;
//     }
    NEventsGood++;
    TPcands->Fill(DATA->mm_EventHits.Time());

    // tp packets
    TPLinkedTrack tr = DATA->tp_EventTracks.Get();
    int Ntptrk = DATA->tp_EventTracks.GetNTrack();

    vector <TPTrack> candtracks;
    vector <int> nMatchHits;

    
    // Didn't find any tracks!
    if (Ntptrk == 0) {
      continue;
      cout << pink << "Event " << evt << " didn't have any trigger matches! " << end << endl;
    }

    if (debug)
      cout << blue << "Looking for best trigger match" << end << endl;


    
    int nTP;
    int nTrk = Ntptrk;
    while (nTrk > 0){
      vector <MMHit> tpHits;
      TPTrack tp_track = TPTrack();
      int nMatch = 0;
      nTP = tr.GetNHits();
      if (debug)
        cout << "hits: " << nTP << endl;
      for (int i = 0; i < nTP; i++){
        tp_track += tr.Get(i);
        MMHit planeHit(tr.Get(i).MMFE8(), tr.Get(i).VMM(), tr.Get(i).Channel(), DATA->RunNum);
        for (int j = 0; j < fit_clusters.size(); j++){
          if (fit_clusters[j].ContainsTP(planeHit))
            nMatch+= 1;
        }
      }
      tp_track.SetBCID(tr.BCID());
      candtracks.push_back(tp_track);
      nMatchHits.push_back(nMatch);
      
      nTrk--;
      if (nTrk > 0)
        tr = *tr.GetNext();
    }

    int imax;
    int maxMatch = -1;
    int maxN = -1;
    if (debug)
      cout << green << "Ntracks " << candtracks.size() << end << endl;
    for (int i = 0; i < candtracks.size(); i++){
      if (debug)
        cout << blue << "iteration: " << i << " nmatch " << nMatchHits[i] << end << endl;
      if (nMatchHits[i] > maxMatch){
        //if (nMatchHits[i] > maxMatch && candtracks[i].GetNHits()>maxN){
        maxMatch = nMatchHits[i];
        maxN = candtracks[i].GetNHits();
        imax = i;
        if (debug)
          cout << blue << "Looping through: " << i << " nmatched: " << nMatchHits[i] << " out of "<< maxN <<  end<< endl;
      }
    }

    if (debug)
      cout << blue << "Found best match with " << Ntptrk << " tracks, max hits "<< maxN << " match " << maxMatch << end << endl;    

    TPTrack bestTrack(candtracks[imax]);

    if (maxMatch < 1)
      continue;

    TPmatch->Fill(DATA->mm_EventHits.Time());
    NEventsTrig++;

    if (maxN < 4)
      continue;
    NEventsTrigGoodTP++;
    
    all_timediff->Fill(candtracks[imax].Time()-DATA->mm_EventHits.Time());
    all_MMmatch->Fill(maxMatch);

    if (!track.IsFit()||
        !track.IsTrigCand()){
      if (debug)
        cout << pink << "Event " << evt << " didn't make the cut" << end << endl;
      continue;
    }
    NEventsTrigGoodTrack++;

    int minBCID = 9999;
    double sumBCID = 0.;
    double medBCID;
    vector <int> tp_BCIDs;
    bool valid = false;
    if (Ntptrk > 0){
      if (maxN >= 4){
        for (int i = 0; i < nTP; i++){
          const TPHit thit = bestTrack.Get(i);
          if (thit.BCID() > -1){
            tp_BCIDs.push_back(thit.BCID());
            valid = true;
          }
          if (debug)
            cout << blue << thit.BCID() << end << endl;
        }
        
        if (valid){
          int maxBCID = *max_element(tp_BCIDs.begin(),tp_BCIDs.end());
          int lowBCID = *min_element(tp_BCIDs.begin(),tp_BCIDs.end());
          bool loopflag = (maxBCID - lowBCID) > 4000;
          for (int i = 0; i < nTP; i++){
            if ((tp_BCIDs[i] < 4000) && loopflag){
              tp_BCIDs[i] = tp_BCIDs[i] + 4096;
            }
            sumBCID += tp_BCIDs[i];
          }
            
          sumBCID = sumBCID/nTP;
          minBCID = *min_element(tp_BCIDs.begin(),tp_BCIDs.end());

          sort(tp_BCIDs.begin(), tp_BCIDs.end());
          double mid =tp_BCIDs.size() / 2;
          if ((int)mid%2 == 0)
            medBCID = (tp_BCIDs[mid]+tp_BCIDs[mid+1])/2;
          else
            medBCID = tp_BCIDs[mid+0.5];
        }
      }

      if (debug)
        cout << pink << "Avg BCID: " << sumBCID << end<< endl;

      int tpBCID = DATA->sc_EventHits.TPBCID();
      int tpPh = DATA->sc_EventHits.TPph();
      bool foundmatch = false;
      // if tp sci comes before tp bcid, offset ~ 1944
      // otherwise, offset is -(4096-1944)
      for (int i = 0; i < candtracks.size(); i++){
        packetBCID->Fill(deltaBCID(tpBCID, tpPh, candtracks[i].BCID(), dB));
        if (fabs(deltaBCID(tpBCID, tpPh, candtracks[i].BCID(), dB)) < 10){
          TPscimatch->Fill(DATA->mm_EventHits.Time());
          NEventsTrigScintMatch++;
          foundmatch = true;
          //cout << pink << "found match" << end << endl;
          break;
        }
      }
      if (deltaBCID(tpBCID, tpPh, candtracks[imax].BCID(), dB)< -2150){
        if (debug)
          cout << "Error flag cases" << endl;
        continue;
      }
      if (!foundmatch){
        nomatch_timediff->Fill(candtracks[imax].Time()-DATA->mm_EventHits.Time());
        nomatch_MMmatch->Fill(maxMatch);
      }
      else{
        match_MMmatch->Fill(maxMatch);
      }
//       all_timediff->Fill(candtracks[imax].Time()-DATA->mm_EventHits.Time());
//       all_MMmatch->Fill(maxMatch);
      
      if (valid){
        if ((tpBCID-minBCID) < 0){
          dEarlyBCID->Fill((tpBCID-minBCID+dB));
          dEarlyBCIDph->Fill(tpBCID+tpPh/16.-minBCID+dB);
          dAvgBCID->Fill(tpBCID-sumBCID+dB);
          dAvgBCIDph->Fill(tpBCID+tpPh/16.-sumBCID+dB);
          dMedBCID->Fill(tpBCID-medBCID+dB);
          dMedBCIDph->Fill(tpBCID+tpPh/16.-medBCID+dB);
        }
        else {
          dEarlyBCID->Fill((tpBCID-minBCID+(dB-4096)));
          dEarlyBCIDph->Fill(tpBCID+tpPh/16.-minBCID+(dB-4096));
          dAvgBCID->Fill((tpBCID-sumBCID+(dB-4096)));
          dAvgBCIDph->Fill(tpBCID+tpPh/16.-sumBCID+(dB-4096));
          dMedBCID->Fill((tpBCID-medBCID+(dB-4096)));
          dMedBCIDph->Fill(tpBCID+tpPh/16.-medBCID+(dB-4096));
        }
        if ((tpBCID-minBCID+dB) > 2000 && (tpBCID-minBCID+dB) < 2010){
          cout << tpBCID << endl;
          cout << minBCID << endl;
          break;
        }
        dEarlyBCIDfull->Fill((tpBCID-minBCID+dB));
        dAvgBCIDfull->Fill((tpBCID-sumBCID+dB));

        earliestBCID->Fill(minBCID);
        avgBCID->Fill(sumBCID);
      }
      //      break;
    }

  }

  gStyle->SetTextFont(42);
  cout << "/////////////////////////////////" << endl;
  cout << pink << "EVENT SUMMARY: " << end << endl;
  cout << "NEvents: " << NEvents << endl;
  cout << "NEvents with triggerable planes: " << NEventsGood << endl;
  cout << "NEvents with trigger packets: " << NEventsTrig << endl;
  cout << "NEvents with trigger packets that have 4+ hits: " << NEventsTrigGoodTP << endl;
  cout << "NEvents with trigger packets and with a triggerable track: " << NEventsTrigGoodTrack << endl;
  cout << "NEvents with trigger + scintillator match: " << NEventsTrigScintMatch << endl;
  cout << "/////////////////////////////////" << endl;

  fout->cd();
  fout->mkdir("histograms");
  fout->cd("histograms");
  fracTPscimatch->Divide(TPscimatch, TPcands, 1., 1., "B");
  effTP->Divide(TPmatch, TPcands, 1., 1., "B");
  TPcands->Write();
  TPmatch->Write();
  TPscimatch->Write();
  fracTPscimatch->Write();
  effTP->Write();
  nomatch_timediff->Write();
  all_timediff->Write();
  nomatch_MMmatch->Write();
  match_MMmatch->Write();
  all_MMmatch->Write();
  packetBCID->Write();

  earliestBCID->Write();
  avgBCID->Write();
  dEarlyBCID->Write();
  dEarlyBCIDph->Write();
  dAvgBCID->Write();
  dAvgBCIDph->Write();
  dMedBCID->Write();
  dMedBCIDph->Write();
  dEarlyBCIDfull->Write();
  dAvgBCIDfull->Write();
  
  TCanvas* c2 = new TCanvas("c2","",600,400);
  vector<int> colors = {kOrange+6,kPink+3,kGreen+2,kCyan-6,kAzure+7,kViolet-8};


  // gStyle->SetErrorX(0.);
  // c2->cd();
  // c2->Clear();
  // hitEff->SetMarkerColor(kBlue-7);
  // hitEff->SetMarkerStyle(20);
  // hitEff->SetMinimum(0.);
  // hitEff->Draw("E1 P");
  // hitEff->GetXaxis()->SetTitle("Board Number");
  // //hitEff->GetXaxis()->CenterTitle();
  // hitEff->GetYaxis()->SetTitle("Hit Efficiency");
  // //hitEff->GetYaxis()->CenterTitle();
  // hitEff->GetXaxis()->SetTitleOffset(1.3);
  // hitEff->GetYaxis()->SetTitleOffset(1.3);
  // hitEff->GetXaxis()->SetLabelSize(0.045);
  // hitEff->GetXaxis()->SetTitleSize(0.045);    
  // hitEff->GetYaxis()->SetLabelSize(0.045);
  // hitEff->GetYaxis()->SetTitleSize(0.045);
  // c2->Print("HitEfficiencyPlots/HitEff.pdf");
  // c2->Clear();


  // //

  // TLegend* leg = new TLegend(0.8,0.3,0.94,0.5);
  // hitEff_anglebin[0]->SetMarkerStyle(20);
  // hitEff_anglebin[0]->SetMinimum(0.3);
  // hitEff_anglebin[0]->SetMaximum(1.);
  // hitEff_anglebin[0]->GetXaxis()->SetTitle("Board Number");
  // //hitEff->GetXaxis()->CenterTitle();
  // hitEff_anglebin[0]->GetYaxis()->SetTitle("Hit Efficiency");
  // //hitEff->GetYaxis()->CenterTitle();
  // hitEff_anglebin[0]->GetXaxis()->SetTitleOffset(1.3);
  // hitEff_anglebin[0]->GetYaxis()->SetTitleOffset(1.3);
  // hitEff_anglebin[0]->GetXaxis()->SetLabelSize(0.045);
  // hitEff_anglebin[0]->GetXaxis()->SetTitleSize(0.045);    
  // hitEff_anglebin[0]->GetYaxis()->SetLabelSize(0.045);
  // hitEff_anglebin[0]->GetYaxis()->SetTitleSize(0.045);
  // hitEff_anglebin[0]->Draw("E1 P");
  // for (int i = 0; i < nangles-3; i++){
  //   hitEff_anglebin[i]->SetMarkerColor(colors[i]);
  //   hitEff_anglebin[i]->SetMarkerStyle(20);
  //   hitEff_anglebin[i]->Draw("E1 PL same");
  //   leg->AddEntry(hitEff_anglebin[i],Form("%d#circ #leq #theta < %d#circ",int(track_angle_bins[i]),int(track_angle_bins[i+1])));
  // }
  // leg->Draw();
  // c2->Print("HitEfficiencyPlots/HitEff_angles.pdf");
  // c2->Clear();

  fout->cd();
  fout->Close();
    
}
