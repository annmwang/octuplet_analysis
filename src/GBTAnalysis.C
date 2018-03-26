///  \file   GBTAnalysis.C
///
///  \author Ann Miao Wang
///          (anwang@cern.ch)
///
///  \date   2017 May
///
///  \note   Studying events from the trigger processor

#include "TH1D.h"
#include "TH2D.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include <iostream>
#include <stdlib.h>
#include "TMath.h"
#include "TColor.h"
#include "TRandom3.h"
#include "TLegend.h"
#include <chrono>

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

double channel_from_x_mid(double xpos, int board, GeoOctuplet* geo){
  int SignChannel = geo->Get(board).SignChannel();
  double xorigin  = geo->Get(board).Origin().X();
  double alpha    = geo->Get(board).StripAlpha();
  double channel  = ((xpos - xorigin - tan(alpha)*100) / (SignChannel*0.4)) + 256.5;
  return channel;
}

int deltaBCID(int sciBCID, int sciph, int hitBCID, int offset){
  if ((sciBCID-hitBCID) < 0){
    return sciBCID+sciph/16.-hitBCID+offset;
  }
  else
    return sciBCID+sciph/16.-hitBCID+(offset-4096);
}

void setstyle(){
  gROOT->SetBatch();
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPaintTextFormat(".2f");
  gStyle->SetTextFont(42);
  //  gStyle->SetOptFit(kTRUE);
}

bool triggerable(int xup, int xdown, int uv )
{
  if (xup > 0 && xdown > 0 && uv > 1)
    return true;
  return false;
}

double simulate(int art_win, float sig)
{
  // change these parameters
  int nToys = pow(10,5);
  TH1F *htdis_4;
  TH1F *htdis_3;
  TH1F *htdis;

  TH1F *h_rpairs_full;
  TH1F *h_rpairs_pass;

  TH1F *hmult;
  TH1F *hmult1;
  TRandom3 *rand3 = new TRandom3;
  cout << "Using ART window of "<< art_win << " BCs" << endl;

  
  double eff[8] = {0.72, 0.92, 0.91, 0.52,\
                   0.96, 0.77, 0.94, 0.86};


  // book histo
  hmult1= new TH1F("mult1", " ", 11 ,-0.5, 10.5); 
  hmult= new TH1F("mult", " ", 11 ,-0.5, 10.5); 
  htdis_3= new TH1F("tdis_3", " ", 40 ,-20., 20.);
  htdis_4= new TH1F("tdis_4", " ", 20 ,0., 20.);
  htdis= new TH1F("tdis", " ", 11 ,-.5, 10.5); 
  h_rpairs_full = new TH1F("Random_pairs","Random_pairs",40, -20, 20);
  h_rpairs_pass = new TH1F("Random_pairs_pass","Random_pairs_pass", 40, -20, 20);

  int bcart[8];

  double art[8];
  int art_min,art_max,dbc;
  
  double xpos;
  int hit[8];
  int n_mm = 0;
  int n_xup, n_xdown,n_uv;
  
  int nEvtTrig = 0; // events triggerable
  int nHitsTrig = 0;
  int nEvtPass = 0; // events which fall within BC window
  int nHitsPass = 0;

  for (int j = 0; j < nToys; j++) { // rndm loop
    // printf(" %d \n", j); 
    n_mm = 0;
    n_xup = 0;
    n_xdown = 0;
    n_uv = 0;
    art_min = 20.;
    art_max = -10.; 
    for (int l = 0; l < 8; l++) 
      hit[l]=0;
  
    // hit mm plane assuming certain efficiency
    for (int l = 0; l < 8; l++) {
      xpos = rand3->Rndm();
      if (xpos <= eff[l]) {
        hit[l] = 1;
        n_mm++;
        if(l<2) n_xdown++;
        if(l>5) n_xup++;
        if(l==2 ||l==3 ||l==5 || l==4) n_uv++;
      }
    }
    //    printf("mm hits %d %d %d %d \n", n_mm,n_xup,n_xdown,n_uv); 
    //if (triggerable(n_xup,n_xdown,n_uv)){
    if (!triggerable(n_xup,n_xdown,n_uv))
      continue;
    hmult->Fill(float(n_mm)); // multiplicity

    nEvtTrig++;
    nHitsTrig += n_mm;
    //printf("mm hits, TRG GOOD %d \n", n_mm);
    
    // if hit and triggerable, produce art signal with certain resolution
    for(int l=0; l<8; l++) {
      if ( hit[l] == 1 ) {
        art[l] = rand3->Gaus(162.5,sig);
        bcart[l] = art[l] / 25.;
        if ( bcart[l] < art_min ) 
          art_min = bcart[l];
      }
    }
    
    // save random pairs delta BCID
    for (int k=0; k<8; k++) {
      for (int kk=k+1; kk<8; kk++) {
        if (hit[k]==1 && hit[kk]==1) {
          h_rpairs_full->Fill(bcart[k]-bcart[kk]);
        }
      }
    }
    //  printf("here\n");
    //for (int h=0; h<8; h++) printf("%d ",hit[h]);
    // printf(" \n");
    
    // study effect of art window
    for ( int kk = 0; kk < 8; kk++) {
      if ( hit[kk] == 1 ) {
        //  printf("  plane %d \n",kk);
        dbc = bcart[kk] - art_min + 1;
        if ( dbc > art_win && kk < 2) {
          n_xdown = n_xdown - 1;
          n_mm = n_mm - 1;
          hit[kk] = 0;
        }
        
        if ( dbc > art_win && kk > 5 ) {
          n_xup = n_xup - 1;;
          n_mm = n_mm - 1;;
          hit[kk] = 0;
        } 
        
        if( dbc > art_win && (kk==2 || kk==3 || kk==4 || kk==5) ) {
          n_uv = n_uv - 1;
          n_mm = n_mm - 1;
          hit[kk] = 0;
        }
        
      }
    }
    // save random pairs with art window cutoff
    for ( int k = 0; k < 8; k++) {
      for(int kk = k+1; kk < 8; kk++) {
        if(hit[k]==1 && hit[kk]==1) {
          h_rpairs_pass->Fill(bcart[k] - bcart[kk]);
        }
      }
    }
    
      // printf("CUT  hits %d \n", n_mm); 
    if (!triggerable(n_xup, n_xdown, n_uv))
      continue;
    hmult1->Fill(float(n_mm));

    nEvtPass++;
    nHitsPass += n_mm;
    for (int kk = 0; kk < 8; kk++) {
      if(hit[kk] == 1 && (bcart[kk] - art_min + 1) > art_max) 
        art_max = bcart[kk] - art_min + 1;    
    }
    htdis->Fill(art_max);
    
    //  printf("END LOOP %d \n", j); 
  } // end rndm loop

     
  float trigEff = (float)nEvtPass/(float)nEvtTrig;
  cout << "Number of toys: " << nToys << endl;
  cout << "Number of triggerable events: " << nEvtTrig << endl;
  cout << "Number of triggerable events which pass the BC cut: " << nEvtPass << endl;
  cout << "Efficiency: " << trigEff << endl;
  cout << "Number of ART hits that fall in a triggerable event: " << nHitsTrig << endl;
  cout << "Number of ART hits that fall in a triggerable event that pass the BC cut: " << nHitsPass << endl;


  return trigEff;
}

void show_overflow(TH1D* hist, bool show_underflow, bool show_overflow){
  int nbins = hist->GetNbinsX();
  int underflow = hist->GetBinContent(0);
  int underflowerror = hist->GetBinError(0);
  int overflow = hist->GetBinContent(nbins+1);
  int overflowerror = hist->GetBinError(nbins+1);
  
  int firstbin = hist->GetBinContent(1);
  int firstbinerror = hist->GetBinError(1);
  
  int lastbin = hist->GetBinContent(nbins);
  int lastbinerror = hist->GetBinError(nbins);
  int newerror, newcontent;
  if (show_underflow && underflow != 0){
    newcontent = underflow + firstbin;
    if (firstbin == 0){
      newerror = underflowerror;
    }
    else {
      newerror = sqrt(underflowerror * underflowerror + firstbinerror * firstbinerror);
    }
    hist->SetBinContent(1,newcontent);
    hist->SetBinError(1,newerror);
    hist->SetBinContent(0,0);
    hist->SetBinError(0,0);
  }
  if (show_overflow && overflow != 0){
    newcontent = overflow + lastbin;
    if (lastbin == 0){
      newerror = overflowerror;
    }
    else {
      newerror = sqrt(overflowerror * overflowerror + lastbinerror * lastbinerror);
    }
    hist->SetBinContent(nbins,newcontent);
    hist->SetBinError(nbins,newerror);
    hist->SetBinContent(nbins+1,0);
    hist->SetBinError(nbins+1,0);
  }
}


int main(int argc, char* argv[]){

  char inputFileName[400];
  char outputFileName[400];
  char PDOFileName[400];
  char TDOFileName[400];
  char AlignFileName[400];
  TRandom3 *ran = new TRandom3;
  
  if ( argc < 5 ){
    cout << "Error at Input: please specify input/output .root files ";
    cout << " and (optional) PDO/TDO calibration files" << endl;
    cout << "Example:   ./GBTAnalysis.x -i input.root -o output.root -r runnumber" << endl;
    cout << "Example:   ./GBTAnalysis.x -i input.root -o output.root -r runnumber";
    cout << " -p PDOcalib.root -t TDOcalib.root" << endl;
    cout << "Other options:" << endl;
    cout << "   -p PDOcalib.root : PDO calibration file" << endl;
    cout << "   -t TDOcalib.root : TDO calibration file" << endl;
    cout << "   -a alignment.root : alignment file" << endl;
    return 0;
  }
  bool debug = false;
  bool sim = false;
  int ndisp = 0;

  bool b_input = false;
  bool b_out   = false;
  bool b_pdo   = false;
  bool b_tdo   = false;
  bool b_align   = false;
  
  int iRun = -1;

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
    if (strncmp(argv[i],"-r",2)==0){
      iRun = atoi(argv[i+1]);
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

  if (iRun == -1){
    cout << "Error at Input: please specify run number! (-r flag)" << endl;
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

  int NEvents = 0;
  int NEventsGood = 0;
  int NEventsTrig = 0;
  int NEventsTrigGoodTP = 0;

  int NTrig2U = 0;
  int NTrig2V = 0;
  int NTrig1U1V = 0;
  int NTrig1U2V = 0;
  int NTrig2U1V = 0;
  int NTrig2U2V = 0;

  
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

  TH1D* TP_theta;
  TH1D* MM_theta;

  TH1D* TP_phi;
  TH1D* MM_phi;

  // without extra roads
  TH1D* TP_xres_full;
  TH1D* TP_angres_full;

  TH1D* TP_xres;
  TH1D* TP_angres;

  TH1D* TP_yres;
  TH1D* TP_yres_1U1V;
  TH1D* TP_yres_1U1V_smallangle;
  TH1D* TP_yres_1U1V_adj;
  TH1D* TP_yres_1U1V_apart;
  TH1D* TP_yres_2U1V;
  TH1D* TP_yres_1U2V;
  TH1D* TP_yres_2U2V;
  TH2D* TP_yres_vs_thetay_2U2V;
  TH2D* TP_yres_vs_thetax_2U2V;
  
  TH1D* TP_yres_smallroad;
  TH1D* TP_yres_smallroad_1U1V;
  TH1D* TP_yres_smallroad_1U1V_smallangle;
  TH1D* TP_yres_smallroad_1U1V_adj;
  TH1D* TP_yres_smallroad_1U1V_apart;
  TH1D* TP_yres_smallroad_2U1V;
  TH1D* TP_yres_smallroad_1U2V;
  TH1D* TP_yres_smallroad_2U2V;

  TH1D* MM_y_1U1V;
  TH2D* MM_y_1U1V_ib;
  TH2D* MM_y_1U1V_time;
  TH2D* MM_y_1U1V_thetay;
  TH2D* MM_x_y;
  TH2D* MM_y_1U1V_x;
  TH1D* MM_x_1U1V;
  TH1D* MM_thetay;
  TH2D* MM_y_2U2V_thetay;

  TH1D* MM_y_1U1V_adj;
  TH1D* MM_y_1U1V_apart;

  TH1D* MM_y_2U1V;
  TH1D* MM_y_1U2V;
  TH1D* MM_y_2U2V;
  TH1D* TP_y_1U1V;
  TH1D* TP_y_1U1V_adj;
  TH1D* TP_y_1U1V_apart;
  TH1D* TP_y_2U1V;
  TH1D* TP_y_1U2V;
  TH1D* TP_y_2U2V;

  // no UV pairs res
  TH1D* TP_y_1U;
  TH1D* TP_y_2U;
  TH1D* TP_y_1V;
  TH1D* TP_y_2V;

  // no UV pairs res
  TH1D* MM_y_1U;
  TH1D* MM_y_2U;
  TH1D* MM_y_1V;
  TH1D* MM_y_2V;

  TH1D* MM_y;
  TH1D* TP_y;

  TH2D* TP_xres_angres_full;
  TH2D* TP_xres_angres;
  TH2D* TP_xres_angres_within5;
  TH2D* TP_xres_angres_outside5;

  TH1D* dEarlyBCID;
  TH1D* dEarlyBCIDph;
  TH1D* dAvgBCID;
  TH1D* dAvgBCIDph;
  TH1D* dMedBCID;
  TH1D* dMedBCIDph;
  TH1D* dLatestBCID;
  TH1D* dLatestBCIDph;

  // binned in # hits
  vector <TH1D*> AvgBCID_n; 

  TH1D* hARTwin;
  TH1D* hARTrpairs;
  TH1D* hARTphase;
  TH1D* hARTphase_avg;

  TH1D* passedTrig;
  TH1D* totalTrig;
  TH1D* trigEff;


  TH1D* dEarlyBCIDfull;
  TH1D* dAvgBCIDfull;
  TH1D* earliestBCID;
  TH1D* avgBCID;

  double t0, tf;

  if (iRun == 3522){
    t0 = 1494529235;
    tf = 1495127430;
  }
  else if (iRun == 3530){
    t0 = 1500495684;
    tf = 1500912122;
  }
  else if (iRun == 3528){
    t0 = 1499006282;
    tf = 1499438280;
  }
  else if (iRun == 3527){
    t0 = 1498515411;
    tf = 1499004849;
  }

  TPcands = new TH1D("tpcands","Candidate events with triggerable track",100,t0,tf);
  TPmatch = new TH1D("tpcands_wtp","Trigger candidate events with tp",100,t0,tf);
  TPcands->Sumw2();
  TPmatch->Sumw2();

  effTP = new TH1D("tpcandsfrac_wtp","Fraction of trigger cands with tp",100,t0,tf);
  TPscimatch = new TH1D("tpscimatch","Events with tp and tpsci BCID match",100,t0,tf);
  fracTPscimatch = new TH1D("tpscimatchfrac","Fraction of events with tp and tpsci BCID match", 100, t0,tf);

  nomatch_timediff = new TH1D("timediff_nomatch",
                              "Timestamp difference for those without the scintillator match",
                              100,-0.08,0.14);
  all_timediff = new TH1D("timediff_allcands_wtp",
                              "Timestamp difference for all candidate events with tp",
                              100,-0.08,0.14);

  nomatch_MMmatch = new TH1D("ntp_nomatch",
                              "Number of TP hits matched with MM for those without the scintillator match",
                              9,-0.5,8.5);
  match_MMmatch = new TH1D("ntp_scimatch",
                           "Number of TP hits matched with MM for those with the scintillator match",
                           9,-0.5,8.5);
  all_MMmatch = new TH1D("ntp_allmatch",
                              "Number of TP hits matched with MM for all candidate events with tp",
                              9,-0.5,8.5);

  TP_theta = new TH1D("TP_theta","TP theta",50, -30, 20);
  MM_theta = new TH1D("MM_theta","MM theta",50, -30, 20);
  TP_phi = new TH1D("TP_phi","TP phi",100, -100, 100);
  MM_phi = new TH1D("MM_phi","MM phi",100, -100, 100);

  // without cuts
  TP_xres_full = new TH1D("TP_xres_full","TP x residual full", 63, -10.5,10.5);
  TP_xres_full->StatOverflows(kTRUE);
  TP_angres_full = new TH1D("TP_angres_full","TP angular residual full",100, -100,100);
  TP_angres_full->StatOverflows(kTRUE);

  //  TP_xres = new TH1D("TP_xres","TP x residual", 63, -10.5,10.5);
  TP_xres = new TH1D("TP_xres","TP x residual", 200, -3,3);
  TP_xres->StatOverflows(kTRUE);
  //  TP_angres = new TH1D("TP_angres","TP angular residual",100, -100, 100);
  TP_angres = new TH1D("TP_angres","TP angular residual",162, -40.5, 40.5);
  TP_angres->StatOverflows(kTRUE);

  TP_y = new TH1D("TP_y","TP y distribution",150, -50.5,280.5);
  MM_y = new TH1D("MM_y","MM y distribution",150, -50.5,280.5);
  TP_y->StatOverflows(kTRUE);
  MM_y->StatOverflows(kTRUE);
//   TP_y = new TH1D("TP_y","TP y distribution",801, -300.5,500.5);
//   MM_y = new TH1D("MM_y","MM y distribution",801, -300.5,500.5);

  TP_yres = new TH1D("TP_yres","TP y residual",200, -350.5,350.5);
  //TP_yres = new TH1D("TP_yres","TP y residual",200, -20.5,20.5);
  TP_yres->StatOverflows(kTRUE);
  TP_yres_1U1V = new TH1D("TP_yres_1u1v","TP y residual 1u1v",200, -350.5,350.5);
  TP_yres_1U1V_smallangle = new TH1D("TP_yres_1u1v small angle","TP y residual 1u1v small angle",200, -350.5,350.5);
  TP_yres_1U1V_adj = new TH1D("TP_yres_1u1v adj","TP y residual 1u1v adj",200, -350.5,350.5);
  TP_yres_1U1V_apart = new TH1D("TP_yres_1u1v apart","TP y residual 1u1v apart",200, -350.5,350.5);
  TP_yres_1U2V = new TH1D("TP_yres_1u2v","TP y residual 1u2v",200, -350.5,350.5);
  TP_yres_2U1V = new TH1D("TP_yres_2u1v","TP y residual 2u1v",200, -350.5,350.5);
  TP_yres_2U2V = new TH1D("TP_yres_2u2v","TP y residual 2u2v",200, -350.5,350.5);
  TP_yres_vs_thetay_2U2V = new TH2D("TP_yres_v_thetay_2u2v","",100, -350.5, 350.5, 100, -TMath::Pi(), TMath::Pi());
  TP_yres_vs_thetax_2U2V = new TH2D("TP_yres_v_thetax_2u2v","",100, -350.5, 350.5, 100, -TMath::Pi(), TMath::Pi());

  TP_yres_smallroad = new TH1D("TP_yres_smallroad","TP y residual smallroad",200, -350.5,350.5);
  TP_yres_smallroad->StatOverflows(kTRUE);
  TP_yres_smallroad_1U1V = new TH1D("TP_yres_smallroad_1u1v","TP y residual smallroad 1u1v",200, -350.5,350.5);
  TP_yres_smallroad_1U1V_smallangle = new TH1D("TP_yres_smallroad_1u1v small angle","TP y residual smallroad 1u1v small angle",200, -350.5,350.5);
  TP_yres_smallroad_1U1V_adj = new TH1D("TP_yres_smallroad_1u1v_adj","TP y residual smallroad 1u1v adj",200, -350.5,350.5);
  TP_yres_smallroad_1U1V_apart = new TH1D("TP_yres_smallroad_1u1v_apart","TP y residual smallroad 1u1v apart",200, -350.5,350.5);
  TP_yres_smallroad_1U2V = new TH1D("TP_yres_smallroad_1u2v","TP y residual smallroad 1u2v",200, -350.5,350.5);
  TP_yres_smallroad_2U1V = new TH1D("TP_yres_smallroad_2u1v","TP y residual smallroad 2u1v",200, -350.5,350.5);
  TP_yres_smallroad_2U2V = new TH1D("TP_yres_smallroad_2u2v","TP y residual smallroad 2u2v",200, -350.5,350.5);

  TP_y_1U = new TH1D("TP_ydist_1u","TP_ydist_1u",801, -300.5,500.5);
  TP_y_2U = new TH1D("TP_ydist_2u","TP y distribution 2u",801, -300.5,500.5);
  TP_y_1V = new TH1D("TP_ydist_1v","TP y distribution 1v",801, -300.5,500.5);
  TP_y_2V = new TH1D("TP_ydist_2v","TP y distribution 2v",801, -300.5,500.5);

  MM_x_y = new TH2D("MM_x_y","MM x y",100, -300.5, 500.5, 100, -300.5, 500.5);
  MM_y_1U1V = new TH1D("MM_y_1u1v","MM y  1u1v",100, -300.5, 500.5);
  MM_y_1U1V_ib = new TH2D("MM_y_1u1v vs board","MM y  1u1v vs board",100, -300.5, 500.5, 8, -0.5, 7.5);
  MM_y_1U1V_time = new TH2D("time_1u1v","time 1u1v",100, -300.5, 500.5,200,1494.5*pow(10,6),1495.2*pow(10,6));

  MM_y_1U1V_thetay = new TH2D("MM_y_1u1v vs y track angle","MM y  1u1v vs y track angle",100, -300.5, 500.5, 40, -TMath::Pi(), TMath::Pi());
  MM_y_1U1V_x = new TH2D("MM_y_1u1v_vs_x","MM y  1u1v vs x",100, -300.5, 500.5, 100, -300.5, 500.5);
  MM_y_2U2V_thetay = new TH2D("MM_y_2u2v_vs_y_track_angle","MM y  2u2v vs y track angle",100, -300.5, 500.5, 40, -TMath::Pi(), TMath::Pi());

  MM_x_1U1V = new TH1D("MM x  1u1v","MM x  1u1v",100, -300.5, 500.5);
  MM_thetay= new TH1D("MM y track angle","MM y track angle", 40, -TMath::Pi(), TMath::Pi());

  MM_y_1U1V_adj = new TH1D("MM y  1u1v adj","MM y  1u1v adj",801, -300.5, 500.5);
  MM_y_1U1V_apart = new TH1D("MM y  1u1v apart","MM y  1u1v apart",801, -300.5, 500.5);
  MM_y_1U2V = new TH1D("MM y  1u2v","MM y  1u2v",801, -300.5, 500.5);
  MM_y_2U1V = new TH1D("MM y  2u1v","MM y  2u1v",801, -300.5, 500.5);
  MM_y_2U2V = new TH1D("MM y  2u2v","MM y  2u2v",801, -300.5, 500.5);

  TP_y_1U1V = new TH1D("TP y  1u1v","TP y  1u1v",801, -300.5, 500.5);
  TP_y_1U1V_adj = new TH1D("TP y  1u1v adj","TP y  1u1v adj",801, -300.5, 500.5);
  TP_y_1U1V_apart = new TH1D("TP y  1u1v apart","TP y  1u1v apart",801, -300.5, 500.5);
  TP_y_1U2V = new TH1D("TP y  1u2v","TP y  1u2v",801, -300.5, 500.5);
  TP_y_2U1V = new TH1D("TP y  2u1v","TP y  2u1v",801, -300.5, 500.5);
  TP_y_2U2V = new TH1D("TP y  2u2v","TP y  2u2v",801, -300.5, 500.5);

  MM_y_1U = new TH1D("MM_ydist_1u","MM y distribution 1u",801, -300.5,500.5);
  MM_y_2U = new TH1D("MM_ydist_2u","MM y distribution 2u",801, -300.5,500.5);
  MM_y_1V = new TH1D("MM_ydist_1v","MM y distribution 1v",801, -300.5,500.5);
  MM_y_2V = new TH1D("MM_ydist_2v","MM y distribution 2v",801, -300.5,500.5);

  TP_xres_angres = new TH2D("TP_xres_v_angres","TP x residual v. angular residual",63,-10.5,10.5,1000,-100,100);
  TP_xres_angres->StatOverflows(kTRUE);
  TP_xres_angres_full = new TH2D("TP_xres_v_angres_full","TP x residual v. angular residual full",483,-80.5,80.5,4000,-800,800);
  TP_xres_angres_full->StatOverflows(kTRUE);
  TP_xres_angres_within5 = new TH2D("TP_xres_v_angres_within_5_strips","TP x residual v. angular residual within 5 strips",63,-10.5,10.5,1000,-100,100);
  TP_xres_angres_outside5 = new TH2D("TP_xres_v_angres_outside_5_strips","TP x residual v. angular residual outside 5 strips",63,-10.5,10.5,1000,-100,100);

  packetBCID = new TH1D("Delta_tp_BCID","Delta tp BCID",10001, -5000,5000);

  dEarlyBCID = new TH1D("Delta_Earliest_BCID","Delta Earliest BCID",19, -9.5,9.5);
  dEarlyBCIDph = new TH1D("Delta_Earliest_BCID_phase","Delta Earliest BCID phase",304, -9.5,9.5);

  dAvgBCID = new TH1D("Delta_Avg_BCID","Delta Avg BCID",19, -9.5,9.5);
  dAvgBCIDph = new TH1D("Delta_Avg_BCID_phase","Delta Avg BCID phase",304, -9.5,9.5);

  dMedBCID = new TH1D("Delta_Med_BCID","Delta Med BCID",19, -9.5,9.5);
  dMedBCIDph = new TH1D("Delta_Med_BCID_phase","Delta Med BCID phase",304, -9.5,9.5);

  dLatestBCID = new TH1D("Delta_Latest_BCID","Delta Latest BCID",19, -9.5,9.5);
  dLatestBCIDph = new TH1D("Delta_Latest_BCID_phase","Delta Latest BCID phase",304, -9.5,9.5);

  hARTwin = new TH1D("ART_window","ART window",13, -0.5,12.5);
  hARTwin->Sumw2();
  hARTphase = new TH1D("ART_phase","ART phase",17, -8.5,8.5);
  hARTphase_avg = new TH1D("ART_phase_avg","ART phase avg",17, -8.5,8.5);
  hARTrpairs = new TH1D("ART_rpairs","ART rpairs",25, -12.5,12.5);
  passedTrig = new TH1D("passed_trig","passed_trig",8, 1.5,9.5);
  totalTrig = new TH1D("total_trig","total_trig",8, 1.5,9.5);
  trigEff = new TH1D("trigEff","trigEff",8, 1.5,9.5);

  earliestBCID = new TH1D("earliest_BCID","earliest BCID",5000, -1.,10000.);
  avgBCID = new TH1D("avg_BCID","avg BCID",5000, -1.,10000.);

  dEarlyBCIDfull = new TH1D("delta_earliest_BCID_full","delta earliest BCID full",5000, -10000.,10000.);
  dAvgBCIDfull = new TH1D("delta_avg_BCID_full","delta avg BCID full",5000, -10000.,10000.);

  for (int i = 4; i < 9; i++){
    TH1D* dummy = new TH1D(Form("Delta_Avg_BCID_for_%d_hits",i),Form("Delta Avg BCID for %d hits",i),19,-9.5,9.5);
    AvgBCID_n.push_back(dummy);
  }

  TFile* fout = new TFile(outputFileName, "RECREATE");
  fout->cd();
  fout->mkdir("event_displays");
  MMPlot();

  for(int evt = 0; evt < Nevent; evt++){
    DATA->GetEntry(evt);

//    if (evt > 1000)
//      break;

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

    
    // Calibrate PDO -> Charge
    PDOCalibrator->Calibrate(DATA->mm_EventHits);
    // Calibrate TDO -> Time
    TDOCalibrator->Calibrate(DATA->mm_EventHits);
  
    // initialize PACMAN info for this event
    //PACMAN->SetEventTrigBCID(DATA->mm_trig_BCID);
    //PACMAN->SetEventPadTime(0); // add this

    
    // initialize cluster making
    FILTER->SetRunNumber(DATA->RunNum);

    // set tbegin, tend

    // Set delta scintillator 
    int dB = DATA->tp_EventTracks.FetchSciOffset(DATA->RunNum);

    // book histograms for MM hits
    int Nboard = DATA->mm_EventHits.GetNBoards();

    bool Gothits = false;
    
    // cluster handling 
    vector<MMClusterList> all_clusters;
    for(int i = 0; i < Nboard; i++){
      if (DATA->mm_EventHits[i].GetNHits() == 0)
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

    
    pair < SCHit*,SCHit* > botPair;
    botPair = DATA->sc_EventHits.GetBotPair()[0];
    MMClusterList filtered_allclusters = FILTER->FilterClustersScint(fit_clusters, *GEOMETRY, botPair.first->Channel(), evt);
    
    MMTrack track = FITTER->Fit(filtered_allclusters, *GEOMETRY);
    
    // Fun event displays
    if (DATA->mm_EventNum == 2579 || DATA->mm_EventNum == 2625 || DATA->mm_EventNum == 3402){
      TCanvas* can;
      can = Plot_Track2D(Form("track2D_%05d_all", DATA->mm_EventNum), track, *GEOMETRY, &fit_clusters);
      can->Write();
      delete can;
    }

    //------------------------------------------------//
    // CUT: track quality and trigger candidate cut
    if (!track.IsFit()||
        !track.IsTrigCand()){
      if (debug)
        cout << pink << "Event " << evt << " didn't make the cut" << end << endl;
      continue;
    }

    NEventsGood++;
    TPcands->Fill(DATA->mm_EventHits.Time());
    //------------------------------------------------//

    MM_thetay->Fill(TMath::ATan(track.SlopeY()));
    DATA->tp_EventTracks.SetSciBCID(DATA->sc_EventHits.TPBCID(),DATA->sc_EventHits.TPph());
    DATA->tp_EventTracks.SetSciOffset(dB);


    //------------------------------------------------//
    // CUT: require at least one trigger in the time window
    int Ntptrk = DATA->tp_EventTracks.GetNTrack();
    
    vector <TPTrack> candtracks;
    vector <int> nMatchHits;
    
    // Didn't find any tracks!
    if (Ntptrk == 0) {
      if (debug)
        cout << pink << "Event " << evt << " didn't have any trigger matches! " << end << endl;
      continue;
    }
    //------------------------------------------------//


    // looking for best trigger match
    if (DATA->tp_EventTracks.Highlander(fit_clusters, true, 16)==nullptr)
      continue;

    TPmatch->Fill(DATA->mm_EventHits.Time());

    TPTrack bestTrack = *DATA->tp_EventTracks.Highlander(fit_clusters, true, 16);
    MMClusterList clusters_tp;
    for (auto art: bestTrack)
      clusters_tp.AddCluster(MMCluster(MMHit(art->MMFE8(), art->VMM(), art->VMMChannel(), DATA->RunNum)));
    MMTrack track_tp = FITTER->Fit(clusters_tp, *GEOMETRY, DATA->mm_EventNum);
    if (DATA->mm_EventNum == 2579 || DATA->mm_EventNum == 2625 || DATA->mm_EventNum == 3402){
      TCanvas* can;
      can = Plot_Track2D(Form("tptrack2D_%05d_all", DATA->mm_EventNum), track, *GEOMETRY, &clusters_tp);
      can->Write();
      delete can;
    }

    NEventsTrig++;

    int maxMatch = bestTrack.NMatch();

    //------------------------------------------------//
    // CUT: if vmm3, board 118 exists, skip for Runs 3525-3539
    bool badvmm = false;
    for (auto art: bestTrack){
      if (art->MMFE8Index()==0 && art->VMM()==3)
        badvmm = true;
    }
    if (badvmm && DATA->RunNum >= 3525 && DATA->RunNum <= 3539)
      continue;
    //------------------------------------------------//

    //------------------------------------------------//
    // CUT: require at least one hit in common with the full list of MM hits
    if (maxMatch < 1)
      continue;
    //------------------------------------------------//
    
    //------------------------------------------------//
    // CUT: get rid of error flag cases
    if (bestTrack.BCID() > 4095) {
      if (debug)
        cout << "Error flag cases" << endl;
      continue;
    }
    //------------------------------------------------//

    if (!bestTrack.IsTrigCand())
      continue;

    NEventsTrigGoodTP++;


    // RESOLUTIONS! //

    // TP x position resolution
    double dX,tpX,tpZ;
    tpX = DATA->tp_EventTracks.AverageX(bestTrack, *GEOMETRY);
    tpZ = DATA->tp_EventTracks.AverageZ(bestTrack, *GEOMETRY);
    double mmX = track.PointZ(tpZ).X();
    dX = tpX - mmX;
    if (debug)
      cout << blue << "dX: " << dX << end <<  endl;

    // TP angular resolution
    double mm_theta = TMath::ATan(track.SlopeX());
    double tp_theta = TMath::ATan(bestTrack.MxLocal());

    TP_theta->Fill(tp_theta/TMath::Pi()*180);
    MM_theta->Fill(mm_theta/TMath::Pi()*180);

    double dtheta = tp_theta-mm_theta;
    TP_xres_angres_full->Fill(dX,dtheta*1000);    
    TP_xres_full->Fill(dX);
    TP_angres_full->Fill(dtheta*1000);
    
    if (debug)
      cout << blue << "found angular res" << end << endl;


    //------------------------------------------------//
    // CUT: if x plane art hits not within 24*0.4 = 9.6 mm of the MM track
    // CUT: if uv plane art hits not within 37*0.4 = 14.8 mm of the MM track

    bool outofvmm = false;

    for (auto art: bestTrack){
      int ip = art->MMFE8Index();
      double artx_ch = art->Channel();
      double mmx_ch = channel_from_x_mid(track.PointZ(GEOMETRY->Get(ip).Origin().Z()).X(), ip, GEOMETRY);
      if ( (fabs(artx_ch-mmx_ch) > 24 && art->isX()))
        outofvmm = true;
//       if ( (fabs(artx_ch-mmx_ch) > 24 && art->isX()) || (fabs(artx_ch-mmx_ch) > 37  && (art->isU() || art->isV())))
//         outofvmm = true;
    }
    if (outofvmm){
      continue;
    }

    if (debug)
      cout << blue << "made out of vmm cuts" << end << endl;

    //------------------------------------------------//
    // TP x position resolution
    TP_xres->Fill(dX);

    TP_angres->Fill(dtheta*1000);
    TP_xres_angres->Fill(dX,dtheta*1000);
    //------------------------------------------------//

    // type 2 roads: UV @ small R
    bool outofsmallroad = false; // flag!

    for (auto art: bestTrack){
      int ip = art->MMFE8Index();
      double artx_ch = art->Channel();
      double mmx_ch = channel_from_x_mid(track.PointZ(GEOMETRY->Get(ip).Origin().Z()).X(), ip, GEOMETRY);
      if ( (fabs(artx_ch-mmx_ch) > 24 && art->isX()) || (fabs(artx_ch-mmx_ch) > 40  && (art->isU() || art->isV())))
        outofsmallroad = true;
    }

    //------------------------------------------------//
    // TP y resolution
    //------------------------------------------------//

    double avgX_U, avgZ_U;
    double avgX_V, avgZ_V;
    double MXG, MU, MV;
    double dY, dY_smallroad, dY_dum;
    double mmY, mmY_X;

    double A,B;
    A = (1/TMath::Sin(1.5/180.*TMath::Pi()));
    B = (1/TMath::Tan(1.5/180.*TMath::Pi()));

    // averages of U strips
    avgX_U = DATA->tp_EventTracks.AverageU(bestTrack, *GEOMETRY);
    avgZ_U = DATA->tp_EventTracks.AverageZU(bestTrack, *GEOMETRY);
    avgX_V = DATA->tp_EventTracks.AverageV(bestTrack, *GEOMETRY);
    avgZ_V = DATA->tp_EventTracks.AverageZV(bestTrack, *GEOMETRY);

    // slopes
    MXG = DATA->tp_EventTracks.MXG(bestTrack, *GEOMETRY, 0., track);
    MU = DATA->tp_EventTracks.MU(bestTrack, *GEOMETRY, 0., track);
    MV = DATA->tp_EventTracks.MV(bestTrack, *GEOMETRY, 0., track);

    // IP
    double IPx = track.PointY(0.).X();
    double IPy = track.PointY(0.).Y();
    double IPz = track.PointY(0.).Z();


    // counting UV config fractions
    if (bestTrack.NU() == 2 && bestTrack.NV() == 2)
      NTrig2U2V++;
    else if (bestTrack.NU() == 1 && bestTrack.NV() == 2)
      NTrig1U2V++;
    else if (bestTrack.NU() == 2 && bestTrack.NV() == 1)
      NTrig2U1V++;
    else if (bestTrack.NU() == 1 && bestTrack.NV() == 1)
      NTrig1U1V++;
    else if (bestTrack.NU() == 2)
      NTrig2U++;
    else if (bestTrack.NV() == 2)
      NTrig2V++;


    if (bestTrack.NU() == 0){
      double avgZ_VX = (avgZ_V+tpZ)/2;
      dY_dum = track_tp.PointZ( avgZ_VX ).Y();
      dY = -B*(tpX-avgX_V + (avgZ_V-tpZ)*track.SlopeX()) + 217.9;
      mmY = track.PointZ( avgZ_VX ).Y();
      mmY_X = track.PointZ( avgZ_VX ).X();
      if (bestTrack.NV() == 1){
        TP_y_1V->Fill(dY);
      }
      else{
        TP_y_2V->Fill(dY);
        MM_y_2V->Fill(mmY);
      }
    }
    else if (bestTrack.NV() == 0){
      double avgZ_UX = (avgZ_U+tpZ)/2;
      // old def
      dY_dum = track_tp.PointZ( avgZ_UX ).Y();
      dY = -B*(avgX_U-tpX + (tpZ-avgZ_U)*track.SlopeX()) + 217.9;
      mmY = track.PointZ( avgZ_UX ).Y();
      mmY_X = track.PointZ( avgZ_UX ).X();
      if (bestTrack.NU() == 1)
        TP_y_1U->Fill(dY);
      else{
        TP_y_2U->Fill(dY);
        MM_y_2U->Fill(mmY);
      }
    }
    else{
      double avgZ_UV = (avgZ_U+avgZ_V)/2;
      dY_dum = track_tp.PointZ( avgZ_UV ).Y();
      dY = -B*(avgX_U-avgX_V + (avgZ_V-avgZ_U)*bestTrack.MxLocal())/2 + 217.9;

      mmY = track.PointZ( avgZ_UV ).Y();
      mmY_X = track.PointZ( avgZ_UV ).X();
      if (bestTrack.NV() == 1 && bestTrack.NU() == 1){
        TP_yres_1U1V->Fill(dY-mmY);
        if (!outofsmallroad)
          TP_yres_smallroad_1U1V->Fill(dY-mmY);
        MM_y_1U1V->Fill(mmY);
        MM_x_1U1V->Fill(mmY_X);
        TP_y_1U1V->Fill(dY);


        int ib1 = -1;
        int ib2 = -1;
        vector <int> ib_x;
        for (auto art: bestTrack){
          if (art->isX())
            ib_x.push_back(art->MMFE8Index());
          if (art->isU() || art->isV()) {
            if (ib1 < 0)
              ib1 = art->MMFE8Index();
            else
              ib2 = art->MMFE8Index();
          }
        }
        if ( fabs(ib1-ib2) < 2 && ( (ib1 < 4 && ib2 < 4) || (ib1 >3 && ib2 > 3)) ){
          if (mmY >  20. && mmY < 40. && TMath::ATan(track.SlopeY()) > 0.6){
            if (ndisp < 10){
              ndisp++;
              TCanvas* can;
              can = Plot_Track2D(Form("track2D_%05d_all", DATA->mm_EventNum), track, *GEOMETRY, &fit_clusters);
              can->Write();
              delete can;
            }
          }
          for (int k = 0; k < ib_x.size(); k++){
            MM_y_1U1V_ib->Fill(mmY, ib_x[k]);
          }

          MM_y_1U1V_time->Fill(mmY, DATA->mm_EventHits.Time());
          MM_y_1U1V_thetay->Fill(mmY, TMath::ATan(track.SlopeY()));
          MM_y_1U1V_x->Fill(mmY, mmY_X);

          TP_yres_1U1V_adj->Fill(dY-mmY);
          if (!outofsmallroad)
            TP_yres_smallroad_1U1V_adj->Fill(dY-mmY);
          MM_y_1U1V_adj->Fill(mmY);
          TP_y_1U1V_adj->Fill(dY);
          if (fabs(TMath::ATan(track.SlopeX())) < 0.1){
            TP_yres_1U1V_smallangle->Fill(dY-mmY);
          if (!outofsmallroad)
            TP_yres_smallroad_1U1V_smallangle->Fill(dY-mmY);
          }
        }
        else{
          TP_yres_1U1V_apart->Fill(dY-mmY);
          if (!outofsmallroad)
            TP_yres_smallroad_1U1V_apart->Fill(dY-mmY);
          MM_y_1U1V_apart->Fill(mmY);
          TP_y_1U1V_apart->Fill(dY);
        }
      }

      if (bestTrack.NV() == 2 && bestTrack.NU() == 1){
        TP_yres_1U2V->Fill(dY-mmY);
        if (!outofsmallroad)
          TP_yres_smallroad_1U2V->Fill(dY-mmY);
        MM_y_1U2V->Fill(mmY);
        TP_y_1U2V->Fill(dY);
      }
      if (bestTrack.NV() == 1 && bestTrack.NU() == 2){
        TP_yres_2U1V->Fill(dY-mmY);
        if (!outofsmallroad)
          TP_yres_smallroad_2U1V->Fill(dY-mmY);
        MM_y_2U1V->Fill(mmY);
        TP_y_2U1V->Fill(dY);
      }
      if (bestTrack.NV() == 2 && bestTrack.NU() == 2){
        MM_y_2U2V_thetay->Fill(mmY, TMath::ATan(track.SlopeY()));
        TP_yres_2U2V->Fill(dY-mmY);
        if (!outofsmallroad)
          TP_yres_smallroad_2U2V->Fill(dY-mmY);
        MM_y_2U2V->Fill(mmY);
        TP_y_2U2V->Fill(dY);
        if (!outofsmallroad){
          TP_yres_vs_thetay_2U2V->Fill(dY-mmY,TMath::ATan(track.SlopeY()));
          TP_yres_vs_thetax_2U2V->Fill(dY-mmY,TMath::ATan(track.SlopeX()));
        }
      }
      MM_y->Fill(mmY);
      TP_y->Fill(dY);
    }
    MM_x_y->Fill(mmY, mmY_X);
    TP_y->Fill(dY);
    TP_yres->Fill(dY - mmY);
    if (!outofsmallroad)
      TP_yres_smallroad->Fill(dY - mmY);
    double tp_phi = -1;


    if (outofsmallroad){
      continue;
    }

    // ART PACKET ANALYSIS!

    double averageBCID = -100.;
    double earlyBCID = -100.;
    double lateBCID = -100.;

    // run by run correction to center
    double avgcorr, earlycorr, latecorr;
    
    if (DATA->RunNum==3522){
      avgcorr = 1.5+0.047;
      earlycorr = -0.375+0.0185;
      latecorr = 3.5;
    }
    else if (DATA->RunNum == 3530){
      avgcorr = 0.5-0.9+1.09;
      earlycorr = -0.375-0.9+0.02;
      latecorr = 3.5;
    }
    else if (DATA->RunNum == 3528){
      avgcorr = 0.5-0.7+1.05;
      earlycorr = -0.375-0.4+0.124;
      latecorr = 3.5;
    }
    else if (DATA->RunNum == 3527){
      avgcorr = 0.5-0.75+1.12;
      earlycorr = -0.375-0.5+0.041;
      latecorr = 3.5;
    }



    if (bestTrack.BCIDAverage() >= 0) {
      averageBCID = DATA->tp_EventTracks.deltaBCID(bestTrack.BCIDAverage())+avgcorr;
      earlyBCID = DATA->tp_EventTracks.deltaBCID(bestTrack.BCIDEarliest())+earlycorr;
      lateBCID = DATA->tp_EventTracks.deltaBCID(bestTrack.BCIDLatest())+latecorr;
      if (debug)
        cout << pink << "Avg BCID: " << averageBCID << end<< endl;
      dEarlyBCID->Fill(round(earlyBCID));
      dEarlyBCIDph->Fill(earlyBCID);
      dAvgBCID->Fill(round(averageBCID));
      if (bestTrack.GetNHits()>3)
        AvgBCID_n[bestTrack.GetNHits()-4]->Fill(round(averageBCID));
      dAvgBCIDph->Fill(averageBCID);
      dLatestBCID->Fill(round(lateBCID));
      dLatestBCIDph->Fill(lateBCID);
    }
    hARTwin->Fill(bestTrack.BCIDWindow());

    //TRandom3 *ran = new TRandom3;
    for (int i = 0; i < bestTrack.size(); i++){
      for (int j = i+1; j < bestTrack.size(); j++){
        if ((bestTrack[i].BCID() == -1) || (bestTrack[j].BCID() == -1) )
          continue;
        double randnum = ran->Uniform(1.);
        if (randnum >= 0.5)
          hARTrpairs->Fill(bestTrack[j].BCID()-bestTrack[i].BCID());
        else
          hARTrpairs->Fill(bestTrack[i].BCID()-bestTrack[j].BCID());
      }
    }

    for (int i = 2; i < 10; i++){
      if (bestTrack.BCIDWindow() <= i){
        passedTrig->Fill(i);
      }
    }
      

    if (bestTrack.BCIDEarliest() > 0 && bestTrack.BCIDAverage() > 0) {
      earliestBCID->Fill(bestTrack.BCIDEarliest());
      avgBCID->Fill(bestTrack.BCIDAverage());
      hARTphase->Fill(bestTrack.BCIDEarliestA()-bestTrack.BCIDEarliestB());
      hARTphase_avg->Fill(bestTrack.BCIDAverageA()-bestTrack.BCIDAverageB());
    }
    //------------------------------------------------//
    // CUT: if art hits not within 5*0.4 = 2.0 mm of the MM track
    bool outofcluster = false;
    for (auto art: bestTrack){
      int ip = art->MMFE8Index();
      double artx = GEOMETRY->Get(ip).LocalXatYbegin(art->Channel())+GEOMETRY->Get(ip).Origin().X();
      double mmx = track.PointZ(GEOMETRY->Get(ip).Origin().Z()).X();
      if (debug)
        cout << blue << "art: " << artx << " mmx " << mmx << " ip " << ip << " zpos " << GEOMETRY->Get(ip).Origin().Z() <<end << endl;
      if (fabs(artx-mmx) > 2.0 && art->isX())
        outofcluster = true;
    }
    if (outofcluster){
      TP_xres_angres_outside5->Fill(dX,dtheta*1000);
    }
    else
      TP_xres_angres_within5->Fill(dX,dtheta*1000);
    //------------------------------------------------//

  }


  // make directory
  string dirname = "GBTAnalysisPlots_";
  cout << iRun << endl;
  dirname += to_string(iRun);
    
  mkdir(dirname.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  chdir(dirname.c_str());                                                                                                                                                                                                             

  gStyle->SetErrorX(0.);
  gStyle->SetOptStat(0.);
  gStyle->SetTextFont(42);
  setstyle();
  cout << "/////////////////////////////////" << endl;
  cout << pink << "EVENT SUMMARY: " << end << endl;
  cout << "NEvents: " << NEvents << endl;
  cout << "NEvents with triggerable track: " << NEventsGood << endl;
  cout << "NEvents with trigger packets that have at least one hit and no errors with scint match: " << NEventsTrig << endl;
  cout << "NEvents with trigger packets that have at least one hit and no errors and have 4 hits with scint match: " << NEventsTrigGoodTP << endl;
  cout << "NTrigs with 2U2V: " << NTrig2U2V << endl;
  cout << "NTrigs with 1U2V: " << NTrig1U2V << endl;
  cout << "NTrigs with 2U1V: " << NTrig2U1V << endl;
  cout << "NTrigs with 1U1V: " << NTrig1U1V << endl;
  cout << "NTrigs with 2U: " << NTrig2U << endl;
  cout << "NTrigs with 2V: " << NTrig2V << endl;
  cout << "/////////////////////////////////" << endl;

  TCanvas* c1 = new TCanvas("c1","",800,800);
  c1->SetRightMargin(0.1);
  vector<int> colors = {kOrange+6,kPink+3,kGreen+2,kCyan-6,kAzure+7,kViolet-8};


  c1->cd();
  c1->Clear();

  fout->cd();
  fout->mkdir("histograms");
  fout->cd("histograms");

  // write histograms

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
  dLatestBCID->Write();
  dLatestBCIDph->Write();
  hARTwin->Write();
  hARTphase->Write();
  hARTphase_avg->Write();
  hARTrpairs->Write();
  TP_xres->Write();
  TP_angres->Write();
  TP_xres_full->Write();
  TP_angres_full->Write();
  TP_yres->Write();
  TP_yres_1U1V->Write();
  TP_yres_1U1V_smallangle->Write();
  TP_yres_1U1V_adj->Write();
  TP_yres_1U1V_apart->Write();
  TP_yres_1U2V->Write();
  TP_yres_2U1V->Write();
  TP_yres_2U2V->Write();
  TP_yres_smallroad->Write();
  TP_yres_smallroad_1U1V->Write();
  TP_yres_smallroad_1U1V_smallangle->Write();
  TP_yres_smallroad_1U1V_adj->Write();
  TP_yres_smallroad_1U1V_apart->Write();
  TP_yres_smallroad_1U2V->Write();
  TP_yres_smallroad_2U1V->Write();
  TP_yres_smallroad_2U2V->Write();
  MM_x_y->Write();
  MM_y_1U1V->Write();
  MM_x_1U1V->Write();
  MM_y_1U1V_ib->Write();
  MM_y_1U1V_thetay->Write();
  MM_y_1U1V_time->Write();
  MM_y_1U1V_x->Write();
  MM_y_2U2V_thetay->Write();
  MM_y_1U1V_adj->Write();
  MM_y_1U1V_apart->Write();
  MM_y_1U2V->Write();
  MM_y_2U1V->Write();
  MM_y_2U2V->Write();
  TP_y_1U1V->Write();
  TP_y_1U1V_adj->Write();
  TP_y_1U1V_apart->Write();
  TP_y_1U2V->Write();
  TP_y_2U1V->Write();
  TP_y_2U2V->Write();
  TP_yres_vs_thetax_2U2V->Write();
  TP_yres_vs_thetay_2U2V->Write();
  TP_theta->Write();
  MM_theta->Write();
  MM_thetay->Write();
  TP_phi->Write();
  MM_phi->Write();
  MM_y->Write();
  TP_y->Write();
  TP_y_1U->Write();
  TP_y_2U->Write();
  TP_y_1V->Write();
  TP_y_2V->Write();
  MM_y_1U->Write();
  MM_y_2U->Write();
  MM_y_1V->Write();
  MM_y_2V->Write();
  TP_xres_angres_full->Write();
  TP_xres_angres->Write();
  TP_xres_angres_within5->Write();
  TP_xres_angres_outside5->Write();
  dEarlyBCIDfull->Write();
  dAvgBCIDfull->Write();



  gStyle->SetLegendBorderSize(0.);
  MM_theta->SetLineColor(kRed-9);
  MM_theta->SetLineWidth(3);
  MM_theta->SetTitle("");
  TP_theta->SetLineColor(kCyan-6);
  TP_theta->SetLineWidth(3);
  TP_theta->SetTitle("");
  MM_theta->GetXaxis()->SetTitle("#theta_{track}");
  MM_theta->GetYaxis()->SetTitle("Events");
  //MM_theta->GetYaxis()->CenterTitle();
  MM_theta->GetXaxis()->SetTitleOffset(1.0);
  MM_theta->GetYaxis()->SetTitleOffset(1.8);
  MM_theta->GetXaxis()->SetLabelSize(0.045);
  MM_theta->GetXaxis()->SetTitleSize(0.045);    
  MM_theta->GetYaxis()->SetLabelSize(0.045);
  MM_theta->GetYaxis()->SetTitleSize(0.045);
  MM_theta->Draw("");
  TP_theta->Draw("same");
  TLegend* le = new TLegend(0.65,0.65,0.75,0.8);
  le->AddEntry(MM_theta, "MM");
  le->AddEntry(TP_theta,"TP");
  le->Draw();

  //c1->Print(Form("run_%d.pdf(",GEOMETRY->RunNumber()),"pdf");
  c1->Print("ang.pdf");
  c1->Clear();


  dEarlyBCID->SetLineColor(kBlue-7);
  dEarlyBCID->SetLineWidth(3);
  dEarlyBCID->SetTitle("");
  dEarlyBCID->GetXaxis()->SetTitle("Earliest BCID");
  dEarlyBCID->GetYaxis()->SetTitle("Events");
  //dEarlyBCID->GetYaxis()->CenterTitle();
  dEarlyBCID->GetXaxis()->SetTitleOffset(1.0);
  dEarlyBCID->GetYaxis()->SetTitleOffset(2.);
  dEarlyBCID->GetXaxis()->SetLabelSize(0.045);
  dEarlyBCID->GetXaxis()->SetTitleSize(0.045);    
  dEarlyBCID->GetYaxis()->SetLabelSize(0.045);
  dEarlyBCID->GetYaxis()->SetTitleSize(0.045);
  dEarlyBCID->Draw("");

  TLatex latex4;
  latex4.SetTextColor(kRed);
  latex4.SetTextSize(0.03);
  latex4.SetTextAlign(13);  //align at top   
  latex4.SetNDC();
  latex4.DrawLatex(0.25,0.8,Form("RMS = %3.1f BC",dEarlyBCID->GetRMS())); 

  //c1->Print(Form("run_%d.pdf",GEOMETRY->RunNumber()),"pdf");
  c1->Print("earliest_BCID.pdf");
  c1->Clear();


  // fit this piece of shit
  TF1* f2 = new TF1("gaussian","gaus(0)",-6,6);
  f2->SetParameter(0, dEarlyBCIDph->GetMaximum());
  f2->SetParameter(1, dEarlyBCIDph->GetMean());
  f2->SetParameter(2, dEarlyBCIDph->GetRMS()/2);

  dEarlyBCIDph->Fit(f2,"RQN");
  cout << "Early BCID fit parameters: " << f2->GetParameter(0) << " " << 
    f2->GetParameter(1) << " " << f2->GetParameter(2) << endl;


  dEarlyBCIDph->SetLineColor(kBlue-7);
  dEarlyBCIDph->SetLineWidth(3);
  dEarlyBCIDph->SetTitle("");
  dEarlyBCIDph->GetXaxis()->SetTitle("Earliest BCID (not rounded)");
  dEarlyBCIDph->GetYaxis()->SetTitle("Events");

  dEarlyBCIDph->GetXaxis()->SetTitleOffset(1.0);
  dEarlyBCIDph->GetYaxis()->SetTitleOffset(1.8);
  dEarlyBCIDph->GetXaxis()->SetLabelSize(0.045);
  dEarlyBCIDph->GetXaxis()->SetTitleSize(0.045);    
  dEarlyBCIDph->GetYaxis()->SetLabelSize(0.045);
  dEarlyBCIDph->GetYaxis()->SetTitleSize(0.045);
  dEarlyBCIDph->Draw("");
  f2->Draw("same");

  c1->Print("earliest_BCID_fine.pdf");
  c1->Clear();

  dAvgBCID->SetLineColor(colors[1]);
  dAvgBCID->SetLineWidth(3);
  dAvgBCID->SetTitle("");
  dAvgBCID->GetXaxis()->SetTitle("Average BCID");
  dAvgBCID->GetYaxis()->SetTitle("Events");
  //dAvgBCID->GetYaxis()->CenterTitle();
  dAvgBCID->GetXaxis()->SetTitleOffset(1.0);
  dAvgBCID->GetYaxis()->SetTitleOffset(2.);
  dAvgBCID->GetXaxis()->SetLabelSize(0.045);
  dAvgBCID->GetXaxis()->SetTitleSize(0.045);    
  dAvgBCID->GetYaxis()->SetLabelSize(0.045);
  dAvgBCID->GetYaxis()->SetTitleSize(0.045);
  dAvgBCID->Draw("");
  latex4.DrawLatex(0.25,0.8,Form("RMS = %3.1f BC",dAvgBCID->GetRMS())); 
  c1->Print("avg_BCID.pdf");
  //c1->Print(Form("run_%d.pdf",GEOMETRY->RunNumber()),"pdf");
  c1->Clear();

  dAvgBCIDph->Fit(f2,"RQN");
  cout << "Average BCID fit parameters: " << f2->GetParameter(0) << " " << 
    f2->GetParameter(1) << " " << f2->GetParameter(2) << endl;

  dAvgBCIDph->SetLineColor(colors[1]);
  dAvgBCIDph->SetLineWidth(3);
  dAvgBCIDph->SetTitle("");
  dAvgBCIDph->GetXaxis()->SetTitle("Average BCID (not rounded)");
  dAvgBCIDph->GetYaxis()->SetTitle("Events");
  //dAvgBCIDph->GetYaxis()->CenterTitle();
  dAvgBCIDph->GetXaxis()->SetTitleOffset(1.0);
  dAvgBCIDph->GetYaxis()->SetTitleOffset(1.8);
  dAvgBCIDph->GetXaxis()->SetLabelSize(0.045);
  dAvgBCIDph->GetXaxis()->SetTitleSize(0.045);    
  dAvgBCIDph->GetYaxis()->SetLabelSize(0.045);
  dAvgBCIDph->GetYaxis()->SetTitleSize(0.045);
  dAvgBCIDph->Draw("");
  f2->Draw("same");
  //c1->Print(Form("run_%d.pdf",GEOMETRY->RunNumber()),"pdf");
  c1->Print("avg_BCID_fine.pdf");
  c1->Clear();

  c1->SetLogy(1);


  // TP ANGRES
  show_overflow(TP_angres,true,true);
  TP_angres->SetLineColor(kBlue-7);
  TP_angres->SetLineWidth(3);
  TP_angres->SetMarkerStyle(8);
  TP_angres->SetMarkerSize(1);  
  TP_angres->SetTitle("");
  //  TP_angres->GetXaxis()->SetNdivisions(505);
  TP_angres->GetXaxis()->SetTitle("#theta_{TP} #minus #theta_{MM} (mrad)");
  TP_angres->GetXaxis()->SetTitle("#theta_{local} #minus #theta_{MM} (mrad)");
  TP_angres->GetYaxis()->SetTitle("Events");
  TP_angres->GetXaxis()->SetTitleOffset(1.0);
  TP_angres->GetYaxis()->SetTitleOffset(1.8);
  TP_angres->GetXaxis()->SetLabelSize(0.045);
  TP_angres->GetXaxis()->SetTitleSize(0.045);    
  TP_angres->GetYaxis()->SetLabelSize(0.045);
  TP_angres->GetYaxis()->SetTitleSize(0.045);
  TP_angres->Draw("EP");


  TF1* f3 = new TF1("gaussian","gaus(0)",-10,10.5);
  f3->SetParameter(0, TP_angres->GetMaximum());
  f3->SetParameter(1, TP_angres->GetMean());
  f3->SetParameter(2, TP_angres->GetRMS()/2);
  f3->SetParameter(6, TP_angres->GetMaximum()/100);
  f3->SetParameter(7, TP_angres->GetMean());
  f3->SetParameter(8, TP_angres->GetRMS()*20);
//   f3->SetParameter(9, 4);
//   f3->SetParLimits(9,1,6);
  TP_angres->Fit(f3,"RQN");
  f3->Draw("same");

  TLatex run;
  run.SetTextSize(0.04);
  run.SetTextAlign(13);  //align at top                                                                                                  
  run.SetNDC();
  //  run.DrawLatex(0.65,0.95,Form("Run %d",GEOMETRY->RunNumber()));
  run.DrawLatex(0.65,0.95,"Harvard CRTS");
  TLatex road;
  road.SetTextSize(0.04);
  road.SetTextAlign(13);  //align at top                                                                                                  
  road.SetNDC();
  //  road.DrawLatex(0.2,0.95,"1/4-VMM road");
  latex4.SetTextColor(kBlack);
  latex4.SetTextSize(0.03);
  latex4.SetTextAlign(13);  //align at top   
  latex4.SetNDC();
  latex4.DrawLatex(0.25,0.8,Form("RMS = %3.1f mrad",TP_angres->GetRMS())); 
  latex4.SetTextColor(kRed);
  //  latex4.DrawLatex(0.25,0.75, Form("#sigma_{core} = %3.2f mrad",f3->GetParameter(2))); 
  //  latex4.DrawLatex(0.25,0.7, Form("#sigma_{core} = %3.2f mrad",f3->GetParameter(5))); 
  c1->Print("TP_angres.pdf");
  c1->SetLogy(0);
  c1->Print("TP_angres_lin.pdf");
  c1->Clear();

  c1->SetLogy(1);

  TP_angres_full->SetLineColor(kBlue-7);
  TP_angres_full->SetLineWidth(3);
  TP_angres_full->SetTitle("");
  TP_angres_full->GetXaxis()->SetNdivisions(505);
  TP_angres_full->GetXaxis()->SetTitle("#theta_{local} #minus #theta_{MM} (mrad)");
  TP_angres_full->GetYaxis()->SetTitle("Events");
  TP_angres_full->GetXaxis()->SetTitleOffset(1.0);
  TP_angres_full->GetYaxis()->SetTitleOffset(1.8);
  TP_angres_full->GetXaxis()->SetLabelSize(0.045);
  TP_angres_full->GetXaxis()->SetTitleSize(0.045);    
  TP_angres_full->GetYaxis()->SetLabelSize(0.045);
  TP_angres_full->GetYaxis()->SetTitleSize(0.045);
  TP_angres_full->Draw("");
  run.DrawLatex(0.65,0.95,Form("Run %d",GEOMETRY->RunNumber()));
  road.DrawLatex(0.2,0.95,"1-VMM road");
  latex4.DrawLatex(0.25,0.8,Form("RMS = %3.1f mrad",TP_angres_full->GetRMS())); 
  //  latex4.DrawLatex(0.3,0.8,Form("RMS = %3.1f #pm %3.1f mrad",TP_angres_full->GetRMS(),TP_angres_full->GetRMSError())); 
  //c1->Print(Form("run_%d.pdf",GEOMETRY->RunNumber()),"pdf");
  c1->Print("TP_angres_full.pdf");
  c1->Clear();


  TP_xres->SetLineColor(kBlue-7);
  TP_xres->SetLineWidth(3);
  TP_xres->SetTitle("");
  TP_xres->GetXaxis()->SetTitle("x_{TP} #minus x_{MM} (mm)");
  TP_xres->GetYaxis()->SetTitle("Events");
  TP_xres->GetXaxis()->SetTitleOffset(1.0);
  TP_xres->GetYaxis()->SetTitleOffset(1.8);
  TP_xres->GetXaxis()->SetLabelSize(0.045);
  TP_xres->GetXaxis()->SetTitleSize(0.045);    
  TP_xres->GetYaxis()->SetLabelSize(0.045);
  TP_xres->GetYaxis()->SetTitleSize(0.045);
  TP_xres->Draw("");
  run.DrawLatex(0.65,0.95,Form("Run %d",GEOMETRY->RunNumber()));
  road.DrawLatex(0.2,0.95,"1/4-VMM road");
  TLatex latex3;
  latex3.SetTextColor(kRed);
  latex3.SetTextSize(0.03);
  latex3.SetTextAlign(13);  //align at top   
  latex3.SetNDC();
  latex3.DrawLatex(0.25,.8,Form("RMS = %3.2f mm",TP_xres->GetRMS())); 
  //  latex3.DrawLatex(0.3,.8,Form("RMS = %3.1f #pm %3.1f mm",TP_xres->GetRMS(),TP_xres->GetRMSError())); 
  //c1->Print(Form("run_%d.pdf",GEOMETRY->RunNumber()),"pdf");
  c1->Print("TP_xres.pdf");
  c1->Clear();



  TP_xres_full->SetLineColor(kBlue-7);
  TP_xres_full->SetLineWidth(3);
  TP_xres_full->SetTitle("");
  TP_xres_full->GetXaxis()->SetTitle("x_{TP} #minus x_{MM} (mm)");
  TP_xres_full->GetYaxis()->SetTitle("Events");
  TP_xres_full->GetXaxis()->SetTitleOffset(1.0);
  TP_xres_full->GetYaxis()->SetTitleOffset(1.8);
  TP_xres_full->GetXaxis()->SetLabelSize(0.045);
  TP_xres_full->GetXaxis()->SetTitleSize(0.045);    
  TP_xres_full->GetYaxis()->SetLabelSize(0.045);
  TP_xres_full->GetYaxis()->SetTitleSize(0.045);
  TP_xres_full->Draw("");
  run.DrawLatex(0.65,0.95,Form("Run %d",GEOMETRY->RunNumber()));
  road.DrawLatex(0.2,0.95,"1-VMM road");
  latex3.DrawLatex(0.25,.8,Form("RMS = %3.2f mm",TP_xres_full->GetRMS())); 
  //  latex3.DrawLatex(0.3,.8,Form("RMS = %3.1f #pm %3.1f mm",TP_xres_full->GetRMS(),TP_xres_full->GetRMSError())); 
  //c1->Print(Form("run_%d.pdf",GEOMETRY->RunNumber()),"pdf");
  c1->Print("TP_xres_full.pdf");
  c1->Clear();

  c1->SetLogy(0);
  gStyle->SetPadRightMargin(0.13);
  c1->SetRightMargin(0.2);

  TP_xres_angres_full->SetLineColor(kBlue-7);
  TP_xres_angres_full->SetLineWidth(3);
  TP_xres_angres_full->SetTitle("");
  TP_xres_angres_full->GetXaxis()->SetTitle("x_{TP} #minus x_{MM} (mm)");
  TP_xres_angres_full->GetYaxis()->SetTitle("#theta_{local} #minus #theta_{MM} (mrad)");
  TP_xres_angres_full->GetXaxis()->SetTitleOffset(1.0);
  TP_xres_angres_full->GetYaxis()->SetTitleOffset(1.8);
  TP_xres_angres_full->GetZaxis()->SetTitleOffset(1.8);
  TP_xres_angres_full->GetXaxis()->SetLabelSize(0.045);
  TP_xres_angres_full->GetXaxis()->SetTitleSize(0.045);    
  TP_xres_angres_full->GetYaxis()->SetLabelSize(0.045);
  TP_xres_angres_full->GetYaxis()->SetTitleSize(0.045);
  TP_xres_angres_full->GetZaxis()->SetLabelSize(0.045);
  TP_xres_angres_full->GetZaxis()->SetTitleSize(0.045);
  TP_xres_angres_full->GetZaxis()->SetTitle("Events");
   
  TP_xres_angres_full->Draw("colz");
  run.DrawLatex(0.7,0.95,Form("Run %d",GEOMETRY->RunNumber()));
  road.DrawLatex(0.2,0.95,"1-VMM road");
  TLatex latex5;
  latex5.SetTextColor(kRed);
  latex5.SetTextSize(0.03);
  latex5.SetTextAlign(13);  //align at top                                                                                                                            
  latex5.SetNDC();
  latex5.DrawLatex(0.25,0.85,Form("#Deltax-RMS = %3.2f mm, #Delta#theta-RMS = %3.1f mrad",TP_xres_angres_full->GetRMS(1),TP_xres_angres_full->GetRMS(2)));
  //c1->Print(Form("run_%d.pdf",GEOMETRY->RunNumber()),"pdf");
  c1->Print("TP_xres_angres_full.pdf");
  c1->Clear();

  TP_xres_angres->SetLineColor(kBlue-7);
  TP_xres_angres->SetLineWidth(3);
  TP_xres_angres->SetTitle("");
  TP_xres_angres->GetXaxis()->SetTitle("x_{TP} #minus x_{MM} (mm)");
  TP_xres_angres->GetYaxis()->SetTitle("#theta_{local} #minus #theta_{MM} (mrad)");
  TP_xres_angres->GetZaxis()->SetTitle("Events");
  TP_xres_angres->GetXaxis()->SetTitleOffset(1.0);
  TP_xres_angres->GetYaxis()->SetTitleOffset(1.8);
  TP_xres_angres->GetZaxis()->SetTitleOffset(1.8);
  TP_xres_angres->GetXaxis()->SetLabelSize(0.045);
  TP_xres_angres->GetXaxis()->SetTitleSize(0.045);    
  TP_xres_angres->GetYaxis()->SetLabelSize(0.045);
  TP_xres_angres->GetYaxis()->SetTitleSize(0.045);
  TP_xres_angres->GetZaxis()->SetLabelSize(0.045);
  TP_xres_angres->GetZaxis()->SetTitleSize(0.045);
  TP_xres_angres->Draw("colz");
  run.DrawLatex(0.7,0.95,Form("Run %d",GEOMETRY->RunNumber()));
  road.DrawLatex(0.2,0.95,"1/4-VMM road");
  TLatex latex6;
  latex6.SetTextColor(kRed);
  latex6.SetTextSize(0.03);
  latex6.SetTextAlign(13);  //align at top                                                                                                                            
  latex6.SetNDC();
  latex6.DrawLatex(0.25,0.85,Form("#Deltax-RMS = %3.2f mm, #Delta#theta-RMS = %3.1f mrad",TP_xres_angres->GetRMS(1),TP_xres_angres->GetRMS(2)));
  //c1->Print(Form("run_%d.pdf)",GEOMETRY->RunNumber()),"pdf");
  c1->Print("TP_xres_angres.pdf");
  c1->Clear();

  TP_xres_angres_within5->SetLineColor(kBlue-7);
  TP_xres_angres_within5->SetLineWidth(3);
  TP_xres_angres_within5->SetTitle("");
  TP_xres_angres_within5->GetXaxis()->SetTitle("x_{TP} #minus x_{MM} (mm)");
  TP_xres_angres_within5->GetYaxis()->SetTitle("#theta_{TP} #minus #theta_{MM} (mrad)");
  TP_xres_angres_within5->GetZaxis()->SetTitle("Events");
  TP_xres_angres_within5->GetXaxis()->SetTitleOffset(1.0);
  TP_xres_angres_within5->GetYaxis()->SetTitleOffset(1.8);
  TP_xres_angres_within5->GetXaxis()->SetLabelSize(0.045);
  TP_xres_angres_within5->GetXaxis()->SetTitleSize(0.045);    
  TP_xres_angres_within5->GetYaxis()->SetLabelSize(0.045);
  TP_xres_angres_within5->GetYaxis()->SetTitleSize(0.045);
  TP_xres_angres_within5->GetZaxis()->SetLabelSize(0.045);
  TP_xres_angres_within5->GetZaxis()->SetTitleSize(0.045);
  TP_xres_angres_within5->Draw("colz");
  run.DrawLatex(0.7,0.95,Form("Run %d",GEOMETRY->RunNumber()));
  TLatex latex7;
  latex7.SetTextColor(kRed);
  latex7.SetTextSize(0.03);
  latex7.SetTextAlign(13);  //align at top                                                                                                                            
  latex7.SetNDC();
  latex7.DrawLatex(0.25,0.85,Form("#Deltax-RMS = %3.2f mm, #Delta#theta-RMS = %3.1f mrad",TP_xres_angres_within5->GetRMS(1),TP_xres_angres_within5->GetRMS(2)));
  //  c1->Print("TP_xres_angres_within5.pdf");
  c1->Clear();

  TP_xres_angres_outside5->SetLineColor(kBlue-7);
  TP_xres_angres_outside5->SetLineWidth(3);
  TP_xres_angres_outside5->SetTitle("");
  TP_xres_angres_outside5->GetXaxis()->SetTitle("x_{TP} #minus x_{MM} (mm)");
  TP_xres_angres_outside5->GetYaxis()->SetTitle("#theta_{TP} #minus #theta_{MM} (mrad)");
  TP_xres_angres_outside5->GetZaxis()->SetTitle("Events");
  TP_xres_angres_outside5->GetXaxis()->SetTitleOffset(1.0);
  TP_xres_angres_outside5->GetYaxis()->SetTitleOffset(1.8);
  TP_xres_angres_outside5->GetXaxis()->SetLabelSize(0.045);
  TP_xres_angres_outside5->GetXaxis()->SetTitleSize(0.045);    
  TP_xres_angres_outside5->GetYaxis()->SetLabelSize(0.045);
  TP_xres_angres_outside5->GetYaxis()->SetTitleSize(0.045);
  TP_xres_angres_outside5->GetZaxis()->SetTitleSize(0.045);
  TP_xres_angres_outside5->Draw("colz");
  run.DrawLatex(0.7,0.95,Form("Run %d",GEOMETRY->RunNumber()));
  TLatex latex8;
  latex8.SetTextColor(kRed);
  latex8.SetTextSize(0.03);
  latex8.SetTextAlign(13);  //align at top                                                                                                                            
  latex8.SetNDC();
  latex8.DrawLatex(0.25,0.85,Form("#Deltax-RMS = %3.2f mm, #Delta#theta-RMS = %3.1f mrad",TP_xres_angres_outside5->GetRMS(1),TP_xres_angres_outside5->GetRMS(2)));
  c1->Print("TP_xres_angres_outside5.pdf");
  c1->Clear();

  vector <double> rms;
  vector <double> rmserr;
  vector <double> nhit;
  vector <double> nhiterr;
  for (int i = 0; i < AvgBCID_n.size();i++){
    AvgBCID_n[i]->Write();
    rms.push_back(AvgBCID_n[i]->GetRMS());
    rmserr.push_back(AvgBCID_n[i]->GetRMSError());
    nhit.push_back((double)(i+4));
    nhiterr.push_back(0.0);
  }
  TVectorD rms_d(rms.size(),&rms[0]);
  TVectorD nhit_d(nhit.size(),&nhit[0]);
  TVectorD rmserr_d(rmserr.size(),&rmserr[0]);
  TVectorD nhiterr_d(nhiterr.size(),&nhiterr[0]);
  TGraphErrors* AvgRMS = new TGraphErrors(nhit_d,rms_d,nhiterr_d,rmserr_d);

  AvgRMS->Write();

  c1->cd();

  AvgRMS->SetMarkerColor(kBlue-7);
  AvgRMS->SetMarkerStyle(8);
  AvgRMS->SetMarkerSize(2);
  AvgRMS->SetTitle("");
  AvgRMS->GetXaxis()->SetTitle("Number of ART hits");
  AvgRMS->GetYaxis()->SetTitle("RMS of Average BCIDs");
  //AvgRMS->GetYaxis()->CenterTitle();
  AvgRMS->GetXaxis()->SetTitleOffset(1.0);
  AvgRMS->GetYaxis()->SetTitleOffset(1.5);
  AvgRMS->GetXaxis()->SetLabelSize(0.045);
  AvgRMS->GetXaxis()->SetTitleSize(0.045);    
  AvgRMS->GetYaxis()->SetLabelSize(0.045);
  AvgRMS->GetYaxis()->SetTitleSize(0.045);
  AvgRMS->SetMinimum(0.0);
  AvgRMS->GetXaxis()->SetLimits(3.,9.);
  AvgRMS->Draw("AP");
  TLatex latex;
  latex.SetTextSize(0.06);
  latex.SetTextAlign(13);  //align at top                                                                                                  
  latex.DrawLatex(6.75,1.25,Form("Run %d",GEOMETRY->RunNumber()));
  //TF1* f1 = new TF1("f1","sqrt([0]/x+[1])",3.5,8.5);
  TF1* f1 = new TF1("f1","[0]/sqrt(x)+[1]",3.5,8.5);
  f1->SetParameter(0,1.);
  f1->SetParameter(1,0.7);
  f1->SetParLimits(1,0.65,0.8);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.0001);
  AvgRMS->Fit(f1,"RN");
  f1->SetLineStyle(7);
  f1->SetLineColor(kAzure+1);
  //  f1->Draw("same");
    

  c1->Print("avgrms.pdf");
  


  setstyle();
  TCanvas* c2 = new TCanvas("c2","",800,800);


  c2->cd();

  // sort of messed up estimate of tp eff
  passedTrig->Sumw2();
  totalTrig->Sumw2();
  trigEff->Divide(passedTrig,totalTrig,1.,1.,"B");

  c2->Clear();

  if (sim == true){
  float sig = 38.6;
  vector <double> trigEffs;
  vector <double> trigEffsErr;
  vector <double> dataTrigEffs;
  vector <double> dataTrigEffsErr;
  vector <double> bc_cuts;
  vector <double> bc_cuts_err;
  for (int i = 2; i < 10; i++){
    trigEffs.push_back(simulate(i,sig));
    trigEffsErr.push_back(0.);
    dataTrigEffs.push_back(trigEff->GetBinContent(i-1));
    cout << "trigEff " << i << " " << trigEff->GetBinContent(i-1) << endl;
    dataTrigEffsErr.push_back(trigEff->GetBinError(i-1));
    bc_cuts.push_back(i);
    bc_cuts_err.push_back(0.);
  }

  TVectorD trigEffs_d(trigEffs.size(),&trigEffs[0]);
  TVectorD trigEffsErr_d(trigEffsErr.size(),&trigEffsErr[0]);
  TVectorD bc_cuts_d(bc_cuts.size(),&bc_cuts[0]);
  TVectorD dataTrigEffs_d(dataTrigEffs.size(),&dataTrigEffs[0]);
  TVectorD dataTrigEffsErr_d(dataTrigEffsErr.size(),&dataTrigEffsErr[0]);
  TVectorD bc_cuts_err_d(bc_cuts_err.size(),&bc_cuts_err[0]);
  TGraphErrors* TriggerEffs = new TGraphErrors(bc_cuts_d,trigEffs_d,bc_cuts_err_d,trigEffsErr_d);
  TGraphErrors* TriggerDataEffs = new TGraphErrors(bc_cuts_d,dataTrigEffs_d,bc_cuts_err_d,dataTrigEffsErr_d);

  c1->cd();
  c1->Clear();

  TriggerEffs->SetMarkerColor(kBlue-7);
  TriggerEffs->SetMarkerStyle(8);
  TriggerEffs->SetMarkerSize(2);
  TriggerEffs->SetFillStyle(0);
  TriggerEffs->SetFillColor(0);
  TriggerDataEffs->SetLineColor(kMagenta-7);
  TriggerDataEffs->SetMarkerColor(kMagenta-7);
  TriggerDataEffs->SetMarkerStyle(8);
  TriggerDataEffs->SetMarkerSize(2);
  TriggerDataEffs->SetLineWidth(2);
  TriggerDataEffs->SetFillStyle(0);
  TriggerDataEffs->SetFillColor(0);
  TriggerEffs->SetTitle("");
  TriggerEffs->GetXaxis()->SetTitle("BC window");
  TriggerEffs->GetYaxis()->SetTitle("Trigger Efficiency");
  //TriggerEffs->GetYaxis()->CenterTitle();                                                                                                     
  TriggerEffs->GetXaxis()->SetTitleOffset(1.0);
  TriggerEffs->GetYaxis()->SetTitleOffset(1.5);
  TriggerEffs->GetXaxis()->SetLabelSize(0.045);
  TriggerEffs->GetXaxis()->SetTitleSize(0.045);
  TriggerEffs->GetYaxis()->SetLabelSize(0.045);
  TriggerEffs->GetYaxis()->SetTitleSize(0.045);
  TriggerEffs->SetMinimum(0.0);
  //TriggerEffs->GetXaxis()->SetLimits(3.,9.);
  TriggerEffs->Draw("AP");
  //TriggerDataEffs->Draw("PL same");
  TLatex latex2;
  latex2.SetTextSize(0.06);
  latex2.SetTextAlign(13);  //align at top
  //latex2.DrawLatex(6.75,.75,Form("#sigma=%3.2f",sig));
  TLegend* leg2 = new TLegend(0.65,0.3,0.94,0.5);
  leg2->AddEntry(TriggerEffs, Form("sim: #sigma=%3.1f",sig));
  leg2->AddEntry(TriggerDataEffs,"data");
  leg2->Draw();
  c1->Print("trigeffs.pdf");
  }
  c2->SetLogy();

  TLegend* leg = new TLegend(0.8,0.3,0.94,0.5);
  all_MMmatch->SetLineColor(kBlue-7);
  all_MMmatch->SetLineWidth(2);
  all_MMmatch->SetTitle("");
  all_MMmatch->GetXaxis()->SetTitle("Number of MM-TP matches");
  all_MMmatch->GetYaxis()->SetTitle("Events");
  all_MMmatch->SetMinimum(0.1);
  //all_MMmatch->GetYaxis()->CenterTitle();
  all_MMmatch->GetXaxis()->SetTitleOffset(1.0);
  all_MMmatch->GetYaxis()->SetTitleOffset(1.0);
  all_MMmatch->GetXaxis()->SetLabelSize(0.045);
  all_MMmatch->GetXaxis()->SetTitleSize(0.045);    
  all_MMmatch->GetYaxis()->SetLabelSize(0.045);
  all_MMmatch->GetYaxis()->SetTitleSize(0.045);
  all_MMmatch->Draw("");
  leg->AddEntry(all_MMmatch, "all trig cands");
  nomatch_MMmatch->SetLineColor(kCyan-6);
  nomatch_MMmatch->SetLineWidth(2);
  nomatch_MMmatch->Draw("same");
  leg->AddEntry(nomatch_MMmatch, "no scint match");
  match_MMmatch->SetLineColor(kMagenta-8);
  match_MMmatch->SetLineWidth(2);
  match_MMmatch->Draw("same");
  leg->AddEntry(match_MMmatch, "w/ scint match");
  leg->Draw();
  c2->Print("nmatch.pdf");
  c2->Clear();




  fout->cd();
  fout->Close();
    
}
