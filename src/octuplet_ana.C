/*
 =======================================================
analyze combined mm+scint data files
 =======================================================
 Author      :  P.GIROMINI - hu
 Date        : 5/10/2014
 Version     : 1.00
 Revision    :    
 =======================================================
// EDITED BY AMW Aug 2016
*/

#include <cstdlib>
#include "VMM2.hh"
#include "Octuplet.hh"
#include "scint.hh"
#include "PDOToCharge.hh"
#include "TDOToTime.hh"
#include "pacman.hh"
#include <new>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTimeStamp.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TFile.h"
#include "TAxis.h"
#include  <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include "TPad.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLine.h"
#include "TMath.h"

using namespace std;

int main()
{
  
#define NUMBOARD 8 
#define NUMVMM 8 
#define TDO_UPPER 440
#define MAX_TRIG_DELAY 40 //max diff in BCID between trigger + event
#define MIN_TRIG_DELAY 15 //min diff in BCID between trigger + event
#define PR 0 //print flag
  //initialize calibration
  PDOToCharge PDO2Charge("AllBoards_PDOcalib.root");
  
  //  TDOToTime TDO2Time("~/atlas/VMM2_Calibration/DATA/TDO_Aug3/AllBoards_TDOcalib.root");
  Octuplet oct;
  scint s;
  pacman p;
  
  // declare and construct histograms
  TH1F ** hTDC = new TH1F*[25];
  TH1F ** hATDC = new TH1F*[12];
  TH1F ** hDTDC = new TH1F*[12];
  TH2F ** hDVSATDC = new TH2F*[12];
  TH2F ** hSTDC = new TH2F*[14];
  TH2F * hTOPXY;
  TH2F * hBOTXY;
  TH2F * hTVSTOPX;
  TH2F * hTVSBOTX;
  TH2F * hMULT;
  TH2F * hSOLANG;
  TH1F ** hPad = new TH1F*[2];
  TH2F * hPadMult;
  TH1F * hPadT;
  TH1F * hPadat;
  TH1F * hPadD;
  TH2F *hpadvsx;
  TH2F *hpadvsy;
  TH2F *hPadvsTr;
  TH2F *hEventvsPDO;
  ///////////////////// micromega
  TH1F ** hPDO = new TH1F*[8];
  TH1F ** channel_hits = new TH1F*[8];
  TH1F ** hDeltaBCID = new TH1F*[8];
  TH2F ** hCHvsPDOraw = new TH2F*[8];
  TH2F ** hCHvsPDO = new TH2F*[8];
  TH1F ** hClusterCharge = new TH1F*[8];
  TH1F ** hClusterChargeMult2 = new TH1F*[8];
  TH1F ** hClusterMult = new TH1F*[8];
  TH1F ** hClusterHoles = new TH1F*[8];
  /////////////////////////////////
  
  ////////////////////
  hPadvsTr= new TH2F("padvstr","padvstr",150,150.,300., 200,150.,300.);
  hpadvsx= new TH2F("padvsx","padvsx",120,-120.,120, 100,0.,100.);
  hpadvsy= new TH2F("padvsy","padvsy",120,-120.,120, 100,0.,100.);
  hPadD=  new TH1F("padiff ","padiff ",200,-100.,100.);
  hPadat=  new TH1F("padat ","padat ",200,150.,350.);
  hPadT =  new TH1F("padT ","padT ",100,0.,100.);
  hPad[0]= new TH1F("padw ","padw ",200,150.,350.);
  hPad[1]= new TH1F("pade ","pade ",200,150.,350.);
  hSOLANG=new TH2F("SOLANG","SOLANG ",64,0.,6.4,30,0.7,1.);
  hMULT  = new TH2F("MULT","MULT ",10,0.,10.,10,0.,10.);
  hTDC[0]= new TH1F("BW1 ","BW1 ",150,150.,300.);
  hTDC[1]= new TH1F("BW2 ","BW2 ",150,150.,300.);
  hTDC[2]= new TH1F("BW3 ","BW3 ",150,150.,300.);
  hTDC[3]= new TH1F("BW4 ","BW4 ",150,150.,300.);
  hTDC[4]= new TH1F("BW5 ","BW5 ",150,150.,300.);
  hTDC[5]= new TH1F("BW6 ","BW6 ",150,150.,300.);
  hTDC[6]= new TH1F("BE1 "," BE1",150,150.,300.);
  hTDC[7]= new TH1F("BE2 ","BE2 ",150,150.,300.);
  hTDC[8]= new TH1F("BE3 ","BE3 ",150,150.,300.);
  hTDC[9]= new TH1F("BE4 ","BE4 ",150,150.,300.);
  hTDC[10]= new TH1F("BE5 ","BE5 ",150,150.,300.);
  hTDC[11]= new TH1F("BE6 ","BE6 ",150,150.,300.);
  hTDC[12]= new TH1F("TW1 ","TW1 ",150,150.,300.);
  hTDC[13]= new TH1F("TW2 ","TW2 ",150,150.,300.);
  hTDC[14]= new TH1F("TW3 ","TW3 ",150,150.,300.);
  hTDC[15]= new TH1F("TW4 ","TW4 ",150,150.,300.);
  hTDC[16]= new TH1F("TW5 ","TW5 ",150,150.,300);
  hTDC[17]= new TH1F("TW6 ","TW6 ",150,150.,300.);
  hTDC[18]= new TH1F("TE1 ","TE1 ",150,150.,300.);
  hTDC[19]= new TH1F("TE2 ","TE2 ",150,150.,300.);
  hTDC[20]= new TH1F("TE3 ","TE3 ",150,150.,300.);
  hTDC[21]= new TH1F("TE4 ","TE4 ",150,150.,300.);
  hTDC[22]= new TH1F("TE5 ","TE5 ",150,150.,300.);
  hTDC[23]= new TH1F("TE6 ","TE6 ",150,150.,300.);
  hTDC[24]= new TH1F("TB ","TB ",150,150.,300.);
  hATDC[0]= new TH1F("AB1 ","AB1 ",150,150.,300.);
  hATDC[1]= new TH1F("AB2 ","AB2 ",150,150.,300.);
  hATDC[2]= new TH1F("AB3 ","AB3 ",150,150.,300.);
  hATDC[3]= new TH1F("AB4 ","AB4 ",150,150.,300.);
  hATDC[4]= new TH1F("AB5 ","AB5 ",150,150.,300.);
  hATDC[5]= new TH1F("AB6 ","AB6 ",150,150.,300.);
  hATDC[6]= new TH1F("AT1 "," AT1",150,150.,300.);
  hATDC[7]= new TH1F("AT2 ","AT2 ",150,150.,300.);
  hATDC[8]= new TH1F("AT3 ","AT3 ",150,150.,300.);
  hATDC[9]= new TH1F("AT4 ","AT4 ",150,150.,300.);
  hATDC[10]= new TH1F("AT5 ","AT5 ",150,150.,300.);
  hATDC[11]= new TH1F("AT6 ","AT6 ",150,150.,300.);
  hDTDC[0]= new TH1F("DB1 ","DB1 ",100,-50.,50.);
  hDTDC[1]= new TH1F("DB2 ","DB2 ",100,-50.,50.);
  hDTDC[2]= new TH1F("DB3 ","DB3 ",100,-50.,50.);
  hDTDC[3]= new TH1F("DB4 ","DB4 ",100,-50.,50.);
  hDTDC[4]= new TH1F("DB5 ","DB5 ",100,-50.,50.);
  hDTDC[5]= new TH1F("DB6 ","DB6 ",100,-50.,50.);
  hDTDC[6]= new TH1F("DT1 "," DT1",100,-50.,50.);
  hDTDC[7]= new TH1F("DT2 ","DT2 ",100,-50.,50.);
  hDTDC[8]= new TH1F("DT3 ","DT3 ",100,-50.,50.);
  hDTDC[9]= new TH1F("DT4 ","DT4 ",100,-50.,50.);
  hDTDC[10]= new TH1F("DT5 ","DT5 ",100,-50.,50.);
  hDTDC[11]= new TH1F("DT6 ","DT6 ",100,-50.,50.);
  hSTDC[0]= new TH2F("SB1 ","SB1 ",150,150.,300.,150,150.,300.);
  hSTDC[1]= new TH2F("SB2 ","SB2 ",150,150.,300.,150,150.,300.);
  hSTDC[2]= new TH2F("SB3 ","SB3 ",150,150.,300.,150,150.,300.);
  hSTDC[3]= new TH2F("SB4 ","SB4 ",150,150.,300.,150.,150,300.);
  hSTDC[4]= new TH2F("SB5 ","SB5 ",150,150.,300.,150,150.,300.);
  hSTDC[5]= new TH2F("SB6 ","SB6 ",150,150.,300.,150,150.,300.);
  hSTDC[6]= new TH2F("ST1 ","ST1 ",150,150.,300.,150,150.,300.);
  hSTDC[7]= new TH2F("ST2 ","ST2 ",150,150.,300.,150,150.,300.);
  hSTDC[8]= new TH2F("ST3 ","ST3 ",150,150.,300.,150,150.,300.);
  hSTDC[9]= new TH2F("ST4 ","ST4 ",150,150.,300.,150,150.,300.);
  hSTDC[10]= new TH2F("ST5 ","ST5 ",150,150.,300.,150,150.,300.);
  hSTDC[11]= new TH2F("ST6 ","ST6 ",150,150.,300.,150,150.,300.);
  hSTDC[12]= new TH2F("TVB ","TVB ",150,150.,300.,150,150.,300.);
  hSTDC[13]= new TH2F("T_B ","T_B ",100,-50.,50.,150,150.,300.);
  hDVSATDC[0]= new TH2F("DAB1 ","DAB1 ",100,-50.,50.,150,150.,300.);
  hDVSATDC[1]= new TH2F("DAB2 ","DAB2 ",100,-50.,50.,150,150.,300.);
  hDVSATDC[2]= new TH2F("DAB3 ","DAB3 ",100,-50.,50.,150,150.,300.);
  hDVSATDC[3]= new TH2F("DAB4 ","DAB4 ",100,-50.,50.,150,150.,300.);
  hDVSATDC[4]= new TH2F("DAB5 ","DAB5 ",100,-50.,50.,150,150.,300.);
  hDVSATDC[5]= new TH2F("DAB6 ","DAB6 ",100,-50.,50.,150,150.,300.);
  hDVSATDC[6]= new TH2F("DAT1 ","DAT1 ",100,-50.,50.,150,150.,300.);
  hDVSATDC[7]= new TH2F("DAT2 ","DAT2 ",100,-50.,50.,150,150.,300.);
  hDVSATDC[8]= new TH2F("DAT3 ","DAT3 ",100,-50.,50.,150,150.,300.);
  hDVSATDC[9]= new TH2F("DAT4 ","DAT4 ",100,-50.,50.,150,150.,300.);
  hDVSATDC[10]= new TH2F("DAT5 ","DAT5 ",100,-50.,50.,150,150.,300.);
  hDVSATDC[11]= new TH2F("DAT6 ","DAT6 ",100,-50.,50.,150,150.,300.);
  hTOPXY= new TH2F("TOPXY ","TOPXY ",26,-130,130.,6,-60.,60.);
  hBOTXY= new TH2F("BOTXY ","BOTXY ",26,-130,130.,6,-60.,60.);
  hTVSBOTX= new TH2F("TVSBOTX ","TVSBOTX ",26,-130,130.,150,150.,300.);
  hTVSTOPX= new TH2F("TVSTOPX ","TVSTOPX ",26,-130,130.,100,150.,300.);
  hPadMult=new TH2F("padmult","padmult",10,0.,10.,10,0.,10.);
  hEventvsPDO = new TH2F("PDO hits binned in event number","PDO hits binned in event number",1000,0.,160000.,1000,-0.5,1600.5);
   //////////////// micromega histos
  for (int k=0; k < NUMBOARD; k++){
    hPDO[k] = new TH1F(Form("Board %d PDO",k), Form("Board %d PDO",k) ,500,-0.5,1000.5);
    channel_hits[k] = new TH1F(Form("Board %d Channel Hits",k), Form("Board %d Channel Hits",k) ,513,-0.5,512.5);
    hDeltaBCID[k] = new TH1F(Form("Board %d Delta BCID",k), Form("Board %d Delta BCID",k) ,100,-0.5,101.5);
    hCHvsPDOraw[k] = new TH2F(Form("Board %d Channel vs PDO counts",k), Form("Board %d Channel vs PDO counts",k) ,513,-0.5,512.5, 100,0.,1000.);
    hCHvsPDO[k] = new TH2F(Form("Board %d Channel vs PDO",k), Form("Board %d Channel vs PDO",k) ,513,-0.5,512.5, 100,-0.5,100.5);
    hClusterCharge[k] = new TH1F(Form("Board %d Charge of clusters",k), Form("Board %d Charge of clusters",k) ,50,0.,400.);
    hClusterChargeMult2[k] = new TH1F(Form("Board %d Charge of clusters with >1 multiplicity",k), Form("Board %d Charge of clusters with >1 multiplicity",k) ,50,0.,400.);
    hClusterMult[k] = new TH1F(Form("Board %d Strip Multiplicity of clusters",k), Form("Board %d Strip Multiplicity of clusters",k) ,21,-0.5,20.5);
    hClusterHoles[k] = new TH1F(Form("Board %d Hole Multiplicity of clusters",k), Form("Board %d Hole Multiplicity of clusters",k) ,21,-0.5,20.5);
  }
  /////////////////////////////////////////
  TCanvas *des = new TCanvas("des","des",635,490);
  TCanvas *prun = new TCanvas("prun","prun",635,490); 
  TCanvas *fixsl = new TCanvas("fixsl","fixsl",635,490);
  /////////////////////////////////
  int  nrun, nev, nhit, channel, cont,an_nev,sel_ev,ret_ev,ch_ev;
  int last_nev; //holds last event trig number
  int mm_nev, mm_hit,mm_trig,mm_vmm, mm_bcid,mm_board,mm_hits,mm_ch,mm_pdo,mm_tdo; 
  int sc_sec,sc_nsec,mm_sec,mm_nsec;
  int nch,mm_flag,time_flag;
  //   int numstrip;
  double ped, val,pdo_curr,tdo_curr,t0_off; 
  char mm[2],strg1[5],strg2[6],strg3[9],strg4[1];
  double pdo_cal[NUMBOARD][NUMVMM][64][2];
  double tdo_cal[NUMBOARD][NUMVMM][64][2];
  double tcorr[130];
  double eped,egain;
  double pi=3.1416;
  FILE* fileout;
  int display_on;
  // corrections to  equalize west and east TOF 
  double time_cut=20.;
  double atdc_corr[12];     //correction for W+E
  const int oversampling = 1;
  const double samplingtime = 0.250 + oversampling*0.500;

  // corrections for scint stuff
  
  double top_corr[12],bot_corr[12]; 
  bot_corr[0]=-0.6;
  bot_corr[1]=-1.5;
  bot_corr[2]=0.8;
  bot_corr[3]=-1.5;
  bot_corr[4]=1.4;
  bot_corr[5]=0.8;
  bot_corr[6]=-.8;
  bot_corr[7]=0.1;
  bot_corr[8]=-0.;
  bot_corr[9]=-.8;
  bot_corr[10]=1.4;
  bot_corr[11]=.5;
  top_corr[0]=-1.8;
  top_corr[1]=2.4;
  top_corr[2]=1.2;
  top_corr[3]=0.6;
  top_corr[4]=.5;
  top_corr[5]=-1.6; //used to be -3.6
  top_corr[6]=1.1;
  top_corr[7]=2.;
  top_corr[8]=-.1;
  top_corr[9]=-2.;
  top_corr[10]=.3;
  top_corr[11]=-1.4;
  double top_time=0.;
  double bot_time=0.;
  //  for(int i=0;i<12;i++) atdc_corr[i]=0;  
  //////////////////////////////////
  // corrections to equalize  e+w for t&b
  atdc_corr[0]=-2.8;
  atdc_corr[1]=-2.7;
  atdc_corr[2]=-2.6;
  atdc_corr[3]=-2.6;
  atdc_corr[4]=-2.6;
  atdc_corr[5]=-2.9;
  atdc_corr[6]=-1.4;
  atdc_corr[7]=-1.4;
  atdc_corr[8]=-1.4;
  atdc_corr[9]=-1.4;
  atdc_corr[10]=-1.4;
  atdc_corr[11]=-1.4;
  
  ////////////////////////////////////////////
  /// Set calibration constants
  // AW - here we need to have multiple gain/tdo calibration files
  
  for (int i=0 ; i < 64; i++) {
    for (int j=0; j < NUMVMM; j++) {
      for (int k=0; k < NUMBOARD; k++) {
	pdo_cal[k][j][i][0]=200.; //pedestal
	pdo_cal[k][j][i][1]=20.; //value
      }
    }
  }

  for (int i=0 ; i < 64; i++) {
    for (int j=0; j < NUMVMM; j++) {
      for (int k=0; k < NUMBOARD; k++) {
	tdo_cal[k][j][i][0]=12.;
	tdo_cal[k][j][i][1]=1.3;
      }
    }
  }

  //////////////////////////////////////////////////
  display_on=0;
  printf(" type 1 to see display, 2 to save it, 0 skip it it\n");
  scanf("%d",&display_on);
  
  ///////////////////////////////////////////  
  // open data file
  sel_ev=0;
  an_nev=0;
  ret_ev=0;
  ch_ev=0;
  
  ////////////////////////////////// test timestamp

  fileout = fopen ("data/comb_r3508_pruned.dat","r"); // combined data file
  
  if(fileout == NULL)
    {
      printf("ERROR: Failure opening DATA file\n");
      return 0;
    }
  
  vector<double> tdcinf(3);
  vector<double> mminf(11);
  
  // scan all events
  while (!feof(fileout)) {
    //  while(an_nev<50000) {

    //////////////////////////////////////////////////////////
    double tdc_padw[10],tdc_pade[10];
    int n_padw=0;
    int n_pade=0;
    double pad_at;
    vector< vector<double> > top_array;
    vector< vector<double> > bot_array;
    
    vector< vector<double> > mm_array[NUMBOARD]; 
    double tdcav[12],tdcdiff[12],topxy[2][6],botxy[2][6];
    double length,dtime,botxpos,botypos,topxpos,topypos;
    double mu_phi,mu_stheta,tdo_upcut; 

    mm_flag=0;
    last_nev = nev;
    
    fscanf (fileout, "%s %s %d %d %d %d %d %d %d %d %d %d", strg1,strg2,
	    &nrun, &nev,&mm_nev, &sc_sec,&mm_sec,&sc_nsec,
	    &mm_nsec, &mm_trig, &nhit, &mm_hit);

    if (last_nev == nev) {
      cout << "READING SAME EVENT "<<nev<<" TWICE! May have reached end of file" << endl;
      break;
    }

    an_nev++;

    if (an_nev % 10000 == 0) cout << an_nev << " Events Analysed" << endl;
    
    time_flag=1;
    
    // if the scint timestamp and micromega timestamp aren't similar, flag
    if( abs(sc_sec-mm_sec) >=2) time_flag=0;

    // set tdo upper cut
    //    tdo_upcut=999999;
    t0_off=19.6;
    
    ///////////////////////////////////
      printf ("EVENT HEADER: %d %d %d %d %d \n", nrun, nev,mm_nev, nhit,mm_hit);

    if (nev == 0 && mm_nev!=0) {
      cout << "MISMATCH OF TRIGGER NUMBERS, SCINT: " << nev << " MM: " << mm_nev << endl;
      break;
    }

    ///// read scintillators info: counter #, TDC
    for (int i=0 ; i < nhit; i++) {
      fscanf(fileout, "%d %d", &channel, &cont);
      tdcinf[0] = channel;
      tdcinf[1] = cont;
      if (i<nhit-1) {
	if (channel<12) {
	  //	  cout << "found bottom hit: " << channel << endl;
	  tdcinf[1]=tdcinf[1]+bot_corr[channel];
	  bot_array.push_back(tdcinf);
	  hTDC[channel]->Fill(cont);
	}
	// 1/2 middle counters
	if (channel==12) {
	  if (cont >=200 && cont <= 250){
	    //	  if (cont >=150 && cont <= 200){
	    tdc_padw[n_padw]=cont;
	    n_padw++;
	  }
	  hPad[0]->Fill(cont);
	}
	// the other middle counter
	if (channel==13) {
	  if (cont >=200 && cont <= 250) {
	    //	  if (cont >=150 && cont <= 200) {
	    tdc_pade[n_pade]=cont +2; //2 for a centering
	    n_pade++;
	  }
	  hPad[1]->Fill(cont);
	}  
	if (channel>15 && channel<28) {
	  tdcinf[1]=tdcinf[1]+top_corr[channel-16];
	  hTDC[channel-4]->Fill(cont);    
	  tdcinf[0]=tdcinf[0]-4;
	  top_array.push_back(tdcinf);
	}  
      }
    } // all scintillators read out

    ////////////////////////////////////////////////
    // read micromega if any

    int numprch[2];
    int bID = -1;
    int totHits = 0;
    int mmhits[NUMBOARD];
    numprch[0]=0;
    numprch[1]=0;
    for (int i=0 ; i < mm_hit; i++) {
      
      fscanf(fileout, "%s %d %d %d %d %d %d %d", mm, &mm_vmm,&mm_ch, &mm_pdo,&mm_tdo,&mm_bcid,&mm_board,&mm_hits);
      if (PR==1)
	printf(" board: %d vmm: %d ch:%d pdo:%d tdo:%d bcid:%d \n",mm_board,mm_vmm, mm_ch,mm_pdo,mm_tdo,mm_bcid);

      // // NOISY CHANNELS
      // if (  (mm_board == 111 && (mm_ch + mm_vmm*64)==342)
      // 	  | (mm_board == 116 && (mm_ch + mm_vmm*64)==25)
      // 	  | (mm_board == 101 && (mm_ch + mm_vmm*64)==450)
      // 	  | (mm_board == 109 && (mm_ch + mm_vmm*64)==185)
      // 	  | (mm_board == 112 && (mm_ch + mm_vmm*64)==322)
      // 	  | (mm_board == 112 && (mm_ch + mm_vmm*64)==324)
      // 	  | (mm_board == 107 && (mm_ch + mm_vmm*64)==512) ) continue;

      // // BCID cuts
      if (fabs(mm_trig-mm_bcid) > MAX_TRIG_DELAY) continue; 
      if (fabs(mm_trig-mm_bcid) < MIN_TRIG_DELAY) continue; 

      if (PR==1)
	printf("cutting on BCID %d, %d\n",MIN_TRIG_DELAY, MAX_TRIG_DELAY);

      bID = oct.getBoardIndex(mm_board,nrun);
      pdo_curr=float(mm_pdo);
      if (bID == -1)
	cout << "bID PROBLEM!"<< endl;
      hDeltaBCID[bID]->Fill(fabs(mm_trig-mm_bcid));
      hPDO[bID]->Fill(pdo_curr);
      hCHvsPDOraw[bID]->Fill((mm_ch + mm_vmm*64),mm_pdo);
      tdo_curr=float(mm_tdo);

      //double pc= pdo_cal[bID][mm_vmm][mm_ch][0]; 
      //double gc= pdo_cal[bID][mm_vmm][mm_ch][1];
      double pt= tdo_cal[bID][mm_vmm][mm_ch][0];
      double gt= tdo_cal[bID][mm_vmm][mm_ch][1];

      pdo_curr=PDO2Charge.GetCharge(mm_pdo, mm_board, mm_vmm, mm_ch);

      hEventvsPDO->Fill(nev,pdo_curr+200*bID);
      
      tdo_curr=(tdo_curr-pt)/gt;

      // old calibration format
      // pdo_curr=(pdo_curr-pc)/gc;
      // tdo_curr=(tdo_curr-pt)/gt;


      mminf[0]=mm_ch + mm_vmm*64;
      mminf[1]=pdo_curr;
      mminf[2]=tdo_curr;
      mminf[3]=mm_bcid;
      mm_array[bID].push_back(mminf);

    }
    fscanf(fileout, "%s",strg3); //end event readout
    
    // reorder micromega hits by increasing channel
    for (int k=0; k<NUMBOARD; k++){
      sort(mm_array[k].begin(), mm_array[k].end(), oct.compare);
      mmhits[k] = mm_array[k].size();
      totHits = mmhits[k] + totHits;
    }
   
    //    printf("micromega  hits  %d\n", totHits);      
   

    ///////////////////////////////////////////////////
    // start scintillator based event selection
    // study chamber later
    //========================= 
    int nbhit=bot_array.size();
    int nthit=top_array.size();

    sort(bot_array.begin(), bot_array.end(), oct.compare);
    sort(top_array.begin(), top_array.end(), oct.compare);
    if (PR==1)
      printf("reordered\n");
    //    cout << "nbhits: " << nbhit << " nthits: " << nthit << endl;
 
    //////////// look for W+E pairs
    int bothits = s.pairs(bot_array,160.,220.);
    //    int bothits = s.pairs(bot_array,120.,170.);
    //    cout << "Num bottom pairs: " << bothits << endl;
    int tophits= s.pairs(top_array,180.,240.);
    //    int tophits= s.pairs(top_array,130.,170.);
    //    cout << "Num top pairs: " << tophits << endl;

   //////////////////////////////////////////////////////////
   // keep events with only one top and  bottom  counter
   //////////////////////////////////////////////////////////

   // if(bothits>1 || tophits>1) printf(" nev%d bot %d top%d\n",nev,tophits,bothits);
   hMULT->Fill(bothits,tophits);
   if(bothits ==1 && tophits ==1 ) { /// to be removed
     ////////// paddle multiplicity
     hPadMult->Fill(n_padw,n_pade);
     if(n_padw==1 && n_pade==1) {
       sel_ev++;
       hPad[1]->Fill(tdc_pade[0]);
       hPad[0]->Fill(tdc_padw[0]);
       hPadD->Fill(tdc_pade[0]-tdc_padw[0]);
       
       pad_at=(tdc_pade[0]+tdc_padw[0])/2.;
       hPadat->Fill(pad_at);
       int outofbound=1;
       for(int i=0; i<nthit;i++)   {
	 for(int j=i+1; j<nthit;j++) {
	   if (top_array[i][2]==1.  && top_array[j][2]==1.) {
	     int ind=top_array[i][0];
	     if(ind<18) {
	       hSTDC[ind-6]->Fill(top_array[i][1],top_array[j][1]);
	       tdcav[ind-6]=(top_array[i][1]+top_array[j][1])/2.+atdc_corr[ind-6];
	       tdcdiff[ind-6]=(top_array[i][1]-top_array[j][1]);
	       // topxy[0][ind-12]=50-(ind-12)*20;
	       // topxy[1][ind-12]=-1*(-tdcdiff[ind-6]*3.144-25.);
	       // OLD COORDINATES
	       topxy[1][ind-12]=50-(ind-12)*20;
	       topxy[0][ind-12]=-tdcdiff[ind-6]*3.144;
	       if(abs(topxy[0][ind-12])>100.) outofbound=0;
	       topxpos=topxy[0][ind-12];
	       topypos= topxy[1][ind-12];   
	       //  tdcav[ind-6]= tdcav[ind-6]+AT_Corr[ind-6][0]*tdcdiff[ind-6]+
	       //  AT_Corr[ind-6][1]*tdcdiff[ind-6]*tdcdiff[ind-6];
	     }
	   }
	 }
       }
       
   
       for(int i=0; i<nbhit;i++)   {
	 for(int j=i+1; j<nbhit;j++) {
	   if(bot_array[i][2]==1. && bot_array[j][2]==1. ) {
	     int ind=bot_array[i][0];
	     if(ind<6) {	 
	       tdcdiff[ind]=(bot_array[i][1]-bot_array[j][1]);
	       // botxy[0][ind]=50-ind*20;
	       // botxy[1][ind]=-1*(-tdcdiff[ind]*3.144-25);
	       // OLD COORDINATES
	       botxy[1][ind]=50-ind*20;
	       botxy[0][ind]=-tdcdiff[ind]*3.144;
	       if(abs(botxy[0][ind])>100.) outofbound=0;
	       botxpos=botxy[0][ind];
	       botypos= botxy[1][ind];
	       length=sqrt(274.3*274.3+(botxpos-topxpos)*(botxpos-topxpos)+
			   (topypos-botypos)*(topypos-botypos))/100.;
	       dtime=2.*3.3*length;
	       mu_phi=TMath::ATan2(topypos-botypos,topxpos-botxpos);
	       if(mu_phi<0.) mu_phi=6.2832+mu_phi;
	       mu_stheta=2.47/length;
	       
	       //sqrt((botxpos-topxpos)*(botxpos-topxpos)+
	       //	   (topypos-botypos)*(topypos-botypos))/247.;
	       //////////   TDCAV CORRECTED FOR MUON TOF
	       tdcav[ind]=(bot_array[i][1]+bot_array[j][1])/2.+atdc_corr[ind]+dtime;
	       // tdcav[ind]= tdcav[ind]+AT_Corr[ind][0]*tdcdiff[ind]+
	       //   AT_Corr[ind][1]*tdcdiff[ind]*tdcdiff[ind];
	       hSTDC[ind]->Fill(bot_array[i][1]+dtime,bot_array[j][1]+dtime);
	     }
	   } 
	 }
       }

       for(int i=0; i<nbhit;i++) {
	 if(bot_array[i][2]==1.) {
	   int nhist=bot_array[i][0];
	   hTDC[nhist]->Fill(bot_array[i][1]);
	   if(nhist<6) {
	     
	     hATDC[nhist]->Fill(tdcav[nhist]);
	     hDTDC[nhist]->Fill(tdcdiff[nhist]);
	     hDVSATDC[nhist]->Fill(tdcdiff[nhist],tdcav[nhist]);
	     // if(outofbound==1) hBOTXY->Fill(botxy[0][nhist],botxy[1][nhist]);
	     bot_time=tdcav[nhist];
	   }
	 }
       }
     
       for(int i=0; i<nthit;i++) {
	 // printf("  top %f %f %f \n", top_array[i][0], top_array[i][1],top_array[i][2]);
	 if(top_array[i][2]==1.) {
	   int nhist=top_array[i][0];
	   if(nhist<18) {
	     hATDC[nhist-6]->Fill(tdcav[nhist-6]);  
	     hDTDC[nhist-6]->Fill(tdcdiff[nhist-6]);  
	     // if(outofbound==1) hTOPXY->Fill(topxy[0][nhist-12],topxy[1][nhist-12]);
	     hDVSATDC[nhist-6]->Fill(tdcdiff[nhist-6],tdcav[nhist-6]);
	     top_time=tdcav[nhist-6];
	   }
	   hTDC[nhist]->Fill(top_array[i][1]);
	   
	 }
       }
       double time_diff=(top_time-bot_time);
       double av_time=  (top_time+bot_time)/2.;
       // if(time_diff<0) av_time=av_time+0.719*time_diff;
       //if(time_diff>=0) av_time=av_time-0.579*time_diff;

       if(outofbound==1)  {  // add y pos cut
	 if(abs(top_time-bot_time)<=time_cut) { // used to be 8
	   
	   ret_ev++;
	   mm_flag=1;
	   hTDC[24]->Fill( av_time);
	   hPadT->Fill( pad_at-av_time);
	   hPadvsTr->Fill( av_time,pad_at);
	   hpadvsx->Fill(botxpos,pad_at-av_time);
	   hpadvsy->Fill(botypos,pad_at-av_time);
	   
	   hBOTXY->Fill(botxpos,botypos);
	   hTOPXY->Fill(topxpos,topypos); 
	   hSTDC[12]->Fill(top_time,bot_time);
	   hSTDC[13]->Fill(time_diff, av_time);
	   hSOLANG->Fill(mu_phi,mu_stheta);
	   hTVSTOPX->Fill(topxpos,top_time);
	   
	   hTVSBOTX->Fill(botxpos,bot_time);
	   
	   
	 } //end ot top-bottom dt
       } // end of outofbound
     } // end of pad requirement
   }      //end  if top/bott hit=1
   
 ////////////////////////////////////////////
 //////////////// micromega study
 ///////////////////////////////////////////

   int left=0;

   if(totHits>0 && mm_flag==1 && time_flag==1) {
     //   if(totHits>0 && mm_flag==1 && time_flag==1) {
     ch_ev++; 
     double totchal=0.;
     double totchar=0.;  
     double tcharge=0.;
     int nstrips=0;
     double art0=100.;
     double art1=100.;
     double art=100.;

     // count channel hits for each board
     for(int k=0; k<NUMBOARD; k++) {
       for(int i=0; i<mmhits[k]; i++) {
	 channel_hits[k]->Fill(mm_array[k][i][0]);
	 hCHvsPDO[k]->Fill(mm_array[k][i][0],mm_array[k][i][1]);
       }
     }

     //////////////////////////////////////////
     //////////////////////////////////////////
     ////////////////////// search for clusters
     
     //PACMAN CLUSTERING//

     // chomp forwards
     int nclus = 0;
     int clust_start, clust_beg, clust_mult;
     int clust_end = -1;
     double clust_cha,clust_time, clust_range;
     int irep=0;
     // cout << "having fun" << endl;
     //   for(int i=0; i<mmhits; i++) {
     for(int k=0; k<NUMBOARD; k++){
       //p.forward(mm_array,mmhits,5,TDO_UPPER,k,pad_at,4.,2.5); 
       //       p.forward(mm_array,mmhits,5,TDO_UPPER,k,pad_at,10.,2.5); //used to be 6.
       p.forward(mm_array,mmhits,5,TDO_UPPER,k,pad_at,15.,10.); //used to be 6.
     } // end cluster search
     
     //hnclus->Fill(nclus);
     //   hnclus->Fill(nclus+10*(1-left));
     
     // chomp backwards
     for (int k=0; k<NUMBOARD; k++) {
       p.backward(mm_array,mmhits,5,TDO_UPPER,k,pad_at,10.);
       //       p.backward(mm_array,mmhits,5,TDO_UPPER,k,pad_at,2.5);
       //       p.backward(mm_array,mmhits,5,TDO_UPPER,k,pad_at,2.5);
     }

     // calculate average strip for each cluster, saved in mm_array[k][i][9]
     double xav; 
     double run_clus=0.;
     int ind1=0;
     vector<double> chacont(2);
     vector< vector<double> > cluscha_ord[NUMBOARD];
   
     for(int k=0; k<NUMBOARD; k++) {
       run_clus = 0.;
       for(int i=0; i<mmhits[k]; i++) {
	 mm_array[k][i][9] = 0;
	 ind1=i;
	 if( mm_array[k][i][4] > run_clus) {
	   xav = 0.; 
	   run_clus = mm_array[k][i][4]; 
	   chacont[0]= mm_array[k][i][7]; // cluster total charge
	   chacont[1]= mm_array[k][i][4]; // cluster number
	   cluscha_ord[k].push_back(chacont);
	   for(int j=ind1; j<mmhits[k]; j++) {
	     if (mm_array[k][j][4] == run_clus)  xav = xav + mm_array[k][j][0] * mm_array[k][j][1]; //channels avg weighted by pdo
	   }
	   xav=xav/mm_array[k][i][7];
	   mm_array[k][i][9] = xav; 
	 }
       }
     }

     // reorder clusters according to most charge to least charge
     for (int k=0; k<NUMBOARD; k++){
       run_clus = 0.;
       sort(cluscha_ord[k].begin(), cluscha_ord[k].end(), oct.compare1);
       int num_clus = cluscha_ord[k].size(); // number of clusters for that board
       if (PR==1)
	 cout << "NUMBER OF CLUSTER: " << num_clus << "FOR BOARD " << k << endl;
       int irun;
       for(int j=0; j<num_clus; j++) { 
	 run_clus = cluscha_ord[k][j][1]; //cluster number
	 irun = 0;
	 for(int i=0; i<mmhits[k]; i++) {
	   if(irun == 0 && mm_array[k][i][4]==run_clus) {
	     mm_array[k][i][10]=j+1; //save what order the cluster is in the cluster ranking
	     irun++; 
	   }
	 }
       }
     }

     // count holes
     int is_hole, istep;
     int nBoardsHit = 0;
     int nXBoardsHit = 0;
     for (int k=0; k<NUMBOARD; k++){
       for(int i=0; i<mmhits[k]; i++) { 
	 if ( mm_array[k][i][10]==1. ) { //highest charge loop
	   nBoardsHit++;
	   if (k == (0||1||6||7)) nXBoardsHit++;
	   istep = mm_array[k][i][6]; //cluster end
	   is_hole = mm_array[k][istep][0]-mm_array[k][i][0]+1.-mm_array[k][i][5]; // number of holes
	   hClusterHoles[k]->Fill(is_hole);
	   hClusterCharge[k]->Fill(mm_array[k][i][7]);
	   if (mm_array[k][i][5]>1) {
	     hClusterChargeMult2[k]->Fill(mm_array[k][i][7]);
	   }
	   hClusterMult[k]->Fill(mm_array[k][i][5]);
	 }
       }
     }

     // event display of barycenter positions
     int go_on;
     int bloop = 0;
     int bxloop = 0;
     double board_x[NUMBOARD][2]; //to draw lines corresponding to chambers
     double board_z[NUMBOARD][2];
     for (int k=0; k<NUMBOARD; k++){
       board_x[k][0] = 0.;
       board_z[k][0] = getz(k);
       board_x[k][1] = 204.6; // max x-pos of board
       board_z[k][1] = getz(k);
     }

     double barycenter_x[nBoardsHit];
     double barycenter_y[nBoardsHit] ;
     double barycenter_z[nBoardsHit];
     double barycenterX_x[nXBoardsHit];
     double barycenterX_z[nXBoardsHit];
     int boards_hit[NUMBOARD] = {-1,-1,-1,-1,-1,-1,-1,-1};
     //    if(true) {
     if(display_on==1) {
       des->cd(0);
       for (int k=0; k < NUMBOARD; k++) {
	 if (PR==1)
	   cout << "BOARD: " << k << endl;
	 for (int i=0; i < mmhits[k]; i++) {
	   if (mm_array[k][i][10] == 1.) { //highest charge loop
	     barycenter_x[bloop] = getx(k, mm_array[k][i][9]);
	     barycenter_z[bloop] = getz(k);
	     cout << "X: " << barycenter_x[bloop] << "Z: " << barycenter_z[bloop] << endl;
	     xpos.push_back(getx(k,mm_array[k][i][9]));
	     zpos.push_back(getz(k));
	     boards_hit[k] = bloop; //save bloop
	     bloop++;
	   }
	 }
       }
       cout << "Beginning track reconstruction...\n" << endl;
       NumericalMinimization("Minuit2","",-1);
       if (PR == 1) {
	 for (int i=0; i < nBoardsHit; i++){
	   cout << "X: " << barycenter_x[i] << endl;
	   cout << "Z: " << barycenter_z[i] << endl;
	 }
       }
       TGraph *gr  = new TGraph(nBoardsHit, barycenter_x, barycenter_z);
       TLine *fitxz = new TLine(param[0],0.,param[0]+param[1]*getz(7),getz(7));

       // TGraph *grx  = new TGraph(nXBoardsHit, barycenterX_x, barycenterX_z);
       // 	 gr->Fit(fttr);
       // 	 double c0 = fttr->GetParameter(0);
       // 	 double c1 = fttr->GetParameter(1);
	 
       // 	 double delta_x = 0.;
       // 	 double ytemp = 0.;
       // 	 for (int j=2; j < 6; j++) { //only loop through u-v boards
       // 	   cout << "UV board: " << j << endl;
       // 	   if (boards_hit[j]>= 0){ //a u-v board was hit
       // 	     cout << "X:VALUE: " <<barycenter_x[boards_hit[j]] << " GUESSED: " << (barycenter_z[boards_hit[j]]-c0)/c1 << endl;
       // 	     delta_x = fabs(barycenter_x[boards_hit[j]]-(barycenter_z[boards_hit[j]]-c0)/c1);
       // 	     cout << "DELTA X: " << delta_x << endl;
       // 	     ytemp = gettrans(j) + 200.-delta_x/TMath::Tan(fabs(getalpha(j)));
       // 	     barycenter_y[boards_hit[j]] = ytemp;
       // 	     cout << "Y-VALUE: " << ytemp << endl;
       // 	   }
       // 	 }
       //}
       // if (uv board)
       // y = gettranslation in y + 20-(deltax/tan(1.5degrees))
       // if we have more than two uv points, then fit uv and extrapolate for y-positions of x boards
       
       gr->SetTitle("Highest charge barycenter cluster locations");
       gr->GetXaxis()->SetLimits(-0.1,204.6);
       gr->GetXaxis()->SetTitle("x position (mm)");
       gr->GetYaxis()->SetRangeUser(-0.1,158.);
       gr->GetYaxis()->SetTitle("z position (mm)");
       gr->SetMarkerColor(46);
       gr->SetMarkerStyle(20);
       gr->SetMarkerSize(1);
       gr->SetLineColor(kBlack);
       gr -> Draw("PA");
       TGraph *grb[NUMBOARD];
       for (int k=0; k<NUMBOARD; k++){
       	 grb[k] = new TGraph(2, board_x[k], board_z[k]);
       	 grb[k]->SetLineColor(kBlue);
       	 grb[k]->Draw("L same");
       }
       fitxz->Draw("");
   
       // clear arrays + vectors
       for (int i=0; i < nBoardsHit; i++){
	 barycenter_x[i] = 0;
	 barycenter_z[i] = 0;
       }
       param.clear();
       xpos.clear();
       zpos.clear();

       des->Update();
       printf("write 2 to print + continue, 0 to end\n");
       scanf("%d", &go_on);
       if (go_on==2) {
	 des->Print(Form("barycenter_plts/barycenter_%i.png",nev));
       }
       if (go_on == 0)
	 return 0;
     }
   
   }// end mm study
  } // end while
 fclose(fileout);
 printf("event analyzed %d ev sel %d ev cut %d ev with mm %d\n",an_nev,sel_ev,ret_ev,ch_ev);

  // write histo file
 TFile *fout=new TFile("Run3508_pruned_CLUSTER15_BCID15_40.root","recreate");
 fout->cd();
 for (int  i=0; i<25;i++) {
   hTDC[i]->Write();
 }
 for (int  i=0; i<12;i++) {
   hATDC[i]->Write();
   hDTDC[i]->Write();
   hDVSATDC[i]->Write();
 }
 //delete hSTDC[0];
 for (int  i = 0; i<14;i++) {
   hSTDC[i]->Write();
 }
 hTOPXY->Write();
 hSOLANG->Write();
 hBOTXY->Write();
 hTVSTOPX->Write();
 hTVSBOTX->Write();
 hMULT->Write();
 hPad[0]->Write();
 hPad[1]->Write();
 hPadMult->Write();
 hPadat->Write();
 hPadT->Write();
 hPadD->Write(); 
 hpadvsx->Write(); 
 hpadvsy->Write(); 
 hPadvsTr->Write();
 hEventvsPDO->Write();
 for (int k=0; k<NUMBOARD; k++) {
   hPDO[k]->Write();
   channel_hits[k]->Write();
   hDeltaBCID[k]->Write();
   hCHvsPDOraw[k]->Write();
   hCHvsPDO[k]->Write();
   hClusterCharge[k]->Write();
   hClusterChargeMult2[k]->Write();
   hClusterMult[k]->Write();
   hClusterHoles[k]->Write();
    }
 TCanvas *c1 = new TCanvas("c1","Root Canvas 1");
 c1->cd();
 for (int k=0; k<NUMBOARD; k++) {
   c1->SetLogy(1);
   hPDO[k]->Draw("");
   c1->Print(Form("Board_%d_PDO.pdf",k));
   c1->Clear();
   c1->SetLogy(0);
   channel_hits[k]->Draw();
   c1->Print(Form("Board_%d_ChannelHits.pdf",k));
   c1->Clear();
   c1->SetLogy(1);
   hDeltaBCID[k]->Draw("");
   c1->Print(Form("Board_%d_BCID.pdf",k));
   c1->Clear();
   c1->SetLogy(0);
   c1->SetLogz(1);
   hCHvsPDOraw[k]->Draw("colz");
   c1->Print(Form("Board_%d_PDOraw.pdf",k));
   c1->Clear();
   c1->SetLogy(0);
   c1->SetLogz(1);
   hCHvsPDO[k]->Draw("colz");
   c1->Print(Form("Board_%d_PDO.pdf",k));
   c1->Clear();
   c1->SetLogy(0);
   c1->SetLogz(0);
   hClusterCharge[k]->Draw();
   c1->Print(Form("Board_%d_ClusterCharge.pdf",k));
   c1->Clear();
   hClusterChargeMult2[k]->Draw();
   c1->Print(Form("Board_%d_ClusterChargeMult2.pdf",k));
   c1->Clear();
   hClusterMult[k]->Draw();
   c1->Print(Form("Board_%d_ClusterMult.pdf",k));
   c1->Clear();
   hClusterHoles[k]->Draw();
   c1->Print(Form("Board_%d_ClusterHoles.pdf",k));
   c1->Clear();
 }
 fout->Close();

 // delete pointers
 // for (int  i=0; i<25;i++) {
 //   delete hTDC[i];
 // }
 // for (int  i=0; i<12;i++) {
 //   delete hATDC[i],hDTDC[i],hDVSATDC[i]; 
 // }
 // for (int  i = 0; i<14;i++) {
 //   delete hSTDC[i];
 // }
 // delete hPad[0];
 // delete hPad[1];
 // for (int k = 0; k<NUMBOARD; k++) {
 //   delete hPDO[k];
 //   delete channel_hits[k];
 // }
 return 0;
}

    
