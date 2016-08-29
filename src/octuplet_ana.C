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
  
  //initialize calibration
  PDOToCharge PDO2Charge("/Users/sezata/atlas/VMM2_Calibration/DATA/PDO_Run3505/AllBoards_PDOcalib.root");
  
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
  ///////////////////// micromega
  TH1F ** channel_hits = new TH1F*[8];
  TH2F ** hCHvsPDO = new TH2F*[8];
  TH1F ** hClusterCharge = new TH1F*[8];
  TH1F ** hClusterMult = new TH1F*[8];
  TH1F ** hClusterHoles = new TH1F*[8];
  TH2F *hchthre;
  TH2F *hTCvNSL;
  TH2F *hTCvNSR;
  TH2F *hchvstl;
  TH2F *hchvstr;
  TH2F *hchvstl_raw;
  TH2F *hchvstr_raw;
  TH2F *hnclusvstcha;
  TH2F *houtlier;
  TH2F *hxdis;
  TH2F *hxdis_1;
  TH2F *hxdis_2;
  TH2F *hxdis_2br;
  TH2F *hxdis_ge3;
  TH2F *hxdis_ge3br;
  TH2F *hclustvscha_0;
  TH2F *hclustvscha_1;
  TH2F *hclustcha_01;
  TH2F *hclustcha_11;
  TH2F *hnstrpvscha_0;
  TH2F *hnstrpvscha_1;
  TH2F * hchalvsr;
  TH1F *hnclus;
  TH1F *hmiss_chal;
  TH1F *hmiss_char;
  TH2F *hclu_cha_nstrip;
  TH2F *hchavsxclu; 
  TH2F *hmax_cha;
  TH2F *htvsnst0;
  TH2F *htvsnst1;
  TH2F *havtvsx;
  TH2F *havtvsx23;
  TH2F *hart0;
  TH2F *hart1;
  TH1F *hart;
  TH1F *hartdiff0;
  TH1F *hartdiff1;
  TH2F *hclus_dis;
  TH2F *hearlytvsx;
  TH2F *hxdep;
  TH2F *hydep;
  TH1F *htslope;
  TH2F *htgvsy;
  TH2F *htgvsy_2strip;
  TH2F *hmtpc;
  TH2F *hmtpc_2strip;
  TH2F *hmtpc_1;
  TH2F *hmtpc_2;
  TH2F *hmtpc_3;
  TH2F *hmtpc_ck;
  TH2F *hsng_tvsch;
  TH2F *hsng_befcorr;
  TH2F *hdtime;
  TH1F *hzslope;
  TH1F *hpruslope;
  TH1F *hprutrue;
  TH2F *hxclvstpc2;
  TH2F *hxclvstpcge3;
  TH2F *hcludiffvsmu;
  TH2F *hztpc;
  TH2F *hztpc_2strip;
  TH2F *hftpc;
  TH2F *hfftpc;
  TH2F *hfftpc2;
  TH2F *hftpc_res;
  TH2F *hftpc2;
  TH2F *hftpc_res2;
  TH2F *hxpres;
  TH2F *hxpres2;
  TH2F *hresid;
  TH2F *hresidtrue;
  TH1F *hchid; 
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
   //////////////// micromega histos
  for (int k=0; k < NUMBOARD; k++){
    channel_hits[k] = new TH1F(Form("Board %d Channel Hits",k), Form("Board %d Channel Hits",k) ,513,-0.5,512.5);
    hCHvsPDO[k] = new TH2F(Form("Board %d Channel vs PDO",k), Form("Board %d Channel vs PDO",k) ,513,-0.5,512.5, 100,-0.5,50.5);
    hClusterCharge[k] = new TH1F(Form("Board %d Charge of clusters",k), Form("Board %d Charge of clusters",k) ,50,0.,400.);
    hClusterMult[k] = new TH1F(Form("Board %d StripMultiplicity of clusters",k), Form("Board %d Strip Multiplicity of clusters",k) ,21,-0.5,20.5);
    hClusterHoles[k] = new TH1F(Form("Board %d Hole Multiplicity of clusters",k), Form("Board %d Hole Multiplicity of clusters",k) ,21,-0.5,20.5);
  }
  hchthre= new TH2F("chthre ","chthre ",131,-1.,130.,180,-20.,160.);
  hTCvNSL= new TH2F("TCvNSL ","tcha vs nstrips ",30,0.,30.,360,-20.,600.); 
  hTCvNSR= new TH2F("TCvNSR ","tcha vs nstrips ",30,0.,30.,360,-20.,600.);
  hchvstl= new TH2F("chvstl ","cha vs time ",70,0.,140.,200,100.,900.); 
  hchvstr= new TH2F("chvstr ","cha vs time ",70,0.,140.,200,100.,900.);
  hchvstl_raw= new TH2F("chvstl_raw ","cha vs time ",70,0.,140.,200,100.,900.); 
  hchvstr_raw= new TH2F("chvstr_raw ","cha vs time ",70,0.,140.,200,100.,900.);
  hnclusvstcha= new TH2F("nclusvstcha ","nclust vs cha ",20,0.,20.,70,0.,140.);
  houtlier= new TH2F("outlier ","outlier ",100,0.,100.,300,-100.,500.);
  hclustvscha_0= new TH2F("tvsch0 ","tvscha0 ",250,0.,250.,400,-100.,700.);
  hclustvscha_1= new TH2F("tvsch1 ","tvscha1 ",250,0.,250.,400,-100.,700.);
  hclustcha_01= new TH2F("tch01 ","tvscha0 one ",250,0.,250.,400,-100.,700.);
  hclustcha_11= new TH2F("tch11 ","tvscha1 one ",250,0.,250.,400,-100.,700.);
  hxdis= new TH2F("xdis ","xdis ",130,0.,130.,200,0.,400.);
  hxdis_1= new TH2F("xdis_1 ","xdis ",130,0.,130.,200,0.,400.);
  hxdis_2= new TH2F("xdis_2 ","xdis ",130,0.,130.,200,0.,400.);
  hxdis_2br= new TH2F("xdis_2br","xdis ",130,0.,130.,200,0.,400.);
  hxdis_ge3= new TH2F("xdis_ge3 ","xdis ",130,0.,130.,200,0.,400.);
  hxdis_ge3br= new TH2F("xdis_ge3br","xdis ",130,0.,130.,200,0.,400.);
  hnclus= new TH1F("nclus ","nclus ",20,0.,20.);
  hmiss_chal= new TH1F("miss_chal "," miss chatge board 0 ",150,-50.,100.);
  hmiss_char= new TH1F("miss_char "," miss chatge board 1 ",150,-50.,100.);
  hnstrpvscha_0= new TH2F("nstrpvscha0 ","nstrip vs cha0 ",180,0.,180.,20,0.,20.);
  hnstrpvscha_1= new TH2F("nstrpvscha1 ","nstrip vs cha1 ",180,0.,180.,20,0.,20.);
  hchalvsr= new TH2F("chalvsr ","cha1 vs cha0 ",120,0.,240.,120,0.,240.);
  hclu_cha_nstrip= new TH2F("clu_cha_nstrip ","clust nstri vs cha ",20,0.,20.,120,0.,240.);
  hchavsxclu = new TH2F("chavsxclu ","clust x  vs cha ",140,0.,140.,150,0.,300.);
  hmax_cha= new TH2F("maxcha ","maxcha ",300,0.,300.,300,0.,300.);
  htvsnst0= new TH2F("tvsnst0 ","time vs mult ",20,0.,20.,300,-100.,500.);
  htvsnst1= new TH2F("tvsnst1 ","time vs mult ",20,0.,20.,300,-100.,500.);
  havtvsx = new TH2F("avtvsx ","time vs xdis ",130,0.,130.,300,100.,700.);
  havtvsx23 = new TH2F("avtvsx23 ","time vs xdis ",130,0.,130.,300,100.,700.);
  hart =new TH1F("art "," art ",120, 200.,800.);
  hart0 =new TH2F("art0 "," art0 ",200,0.,1000.,200,-0.,1000.);
  hart1 =new TH2F("art1 "," art0 ",200,0.,1000.,200,-0.,1000.);
  hartdiff0= new TH1F("artdiff0 "," artdiff ",220,-100.,900.);
  hartdiff1= new TH1F("artdiff1 "," artdiff ",220,-100.,900.);
  hclus_dis     =new TH2F("clus_dis"," clust dis ",260,-130.,130.,100,0.,10.);
  hearlytvsx     =new TH2F("earlytvsx"," clust early time dis ",130,-0.,130.,200,0.,1000.);
  hxdep  =new TH2F("xdep", "xdep",50,-100.,100.,130,0.,130.);
  hydep  =new TH2F("ydep", "xdep",6,-60.,60.,130,0.,130.);
  htslope  =new TH1F("tslope", "tslope",130,0.,130.);
  hzslope  =new TH1F("zslope", "x-z fit",128,0.,51.2);
  hzslope->SetMaximum(10.);
  hzslope->SetMinimum(-4.);
  hzslope->SetMarkerStyle(8);
  hprutrue  =new TH1F("prutrue", "fixed angle x-z fit",128,0.,51.2);
  hprutrue->SetMaximum(10.);
  hprutrue->SetMinimum(-4.);
  hprutrue->SetMarkerStyle(8);
  hpruslope  =new TH1F("pruslope", "x-z fit pruned",128,0.,51.2);
  hpruslope->SetMaximum(10.);
  hpruslope->SetMinimum(-4.);
  hpruslope->SetMarkerStyle(8);
  hmtpc=new TH2F("mtpc", "micro-TPC",6,-60.,60.,90,-90.,90.);
  hmtpc_ck=new TH2F("mtpc_ck", "z vs x",130, 0.,52.,200,-20.,20.);
  hmtpc_2strip=new TH2F("mtpc_2strip", "micro-TPC",6,-60.,60.,90,-90.,90.);
  htgvsy=new TH2F("tgvsy", "micro-TPC",6,-60.,60.,180,-15.,15.);
  htgvsy_2strip=new TH2F("tgvsy_2strip", "micro-TPC",6,-60.,60.,180,-15.,15.);
  hmtpc_1=new TH2F("mtpc_1", "micro-TPC",6,-60.,60.,130,0.,130);
  hmtpc_2=new TH2F("mtpc_2", "micro-TPC",6,-60.,60.,20,-0.5,19.5);
  hmtpc_3=new TH2F("mtpc_3", "micro-TPC",6,-60.,60.,20,-0.5,19.5);//20,-0.5,19.5);
  hsng_tvsch=new TH2F("sng_tvsch","sng_tvsch",100,0.,200.,400,0.,800.);
  hsng_befcorr=new TH2F("sng_befcorr","sng_befcorr",100,0.,200.,400,0.,800.);
  hdtime=new TH2F("dtime"," dtime", 80, -40., 40., 150,-150.,150.); 
  hxclvstpc2=new TH2F("xclvstpc2"," x vs xtps", 260, 0., 52., 400,0.,80.);
  hxclvstpcge3=new TH2F("xclvstpcge3"," x vs xtps", 260, 0., 52., 400, 0.,80.);
  hcludiffvsmu=new TH2F("cludiffvsmu"," x diff vs mu",11,-0.5,10.5, 300,-30.,30.);
  hztpc=new TH2F("ztpc"," x diff vs y", 6,-60., 60., 300,-30.,30.);
  hztpc_2strip=new TH2F("ztpc_2strip"," x diff vs y", 6,-60., 60., 300,-30.,30.);
  hftpc=new TH2F("ftpc"," x diff vs y", 6,-60., 60., 300,-30.,30.);
  hfftpc=new TH2F("fftpc"," x diff vs y", 6,-60., 60., 300,-30.,30.);
  hfftpc2=new TH2F("fftpc2"," x diff vs y", 6,-60., 60., 300,-30.,30.);
  hftpc_res=new TH2F("ftpc_res"," x diff vs y", 6,-60., 60., 300,-30.,30.);
  hftpc2=new TH2F("ftpc2"," x diff vs y", 6,-60., 60., 300,-30.,30.);
  hftpc_res2=new TH2F("ftpc_res2"," x diff vs y", 6,-60., 60., 300,-30.,30.);
  hxpres=new TH2F("xpres"," x res", 6,-60., 60., 400, -2.,2.);
  hxpres2=new TH2F("xpres2"," x res", 6,-60., 60., 400,-2.,2.);
  hresid=new TH2F("resid"," z residuals", 6,-60., 60., 100,-25.,25.);
  hresidtrue=new TH2F("residtrue"," z residuals", 6,-60., 60., 100,-25.,25.);
  hchid= new TH1F("chid","chi2",100,0.,100.);
  /////////////////////////////////////////
  //  TCanvas *des = new TCanvas("des","des");
  // TCanvas *prun = new TCanvas("prun","prun"); 
  // TCanvas *fixsl = new TCanvas("fixsl","fixsl");
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
  double time_cut=8.;
  double atdc_corr[12];     //correction for W+E
  double AT_Corr[12][2];    // flatteners of W+E vs W-E
  const int oversampling = 1;
  const double samplingtime = 0.250 + oversampling*0.500;

  // corrections to flatten e+w vs e-w  -questionables, skipped
  /////////////////////////////////////////////////////
  AT_Corr[0][0]=-0.033 ;
  AT_Corr[0][1]=-0.01 ;
  AT_Corr[1][0]=.05 ;
  AT_Corr[1][1]=-.0088 ;
  AT_Corr[2][0]=-.097 ;
  AT_Corr[2][1]=-.0094 ;
  AT_Corr[3][0]=-.053 ;
  AT_Corr[3][1]= -.01;
  AT_Corr[4][0]= -.048;
  AT_Corr[4][1]= -.0081;
  AT_Corr[5][0]= -.009;
  AT_Corr[5][1]=-.0089 ;
  AT_Corr[6][0]=-.0105 ;
  AT_Corr[6][1]=-.01 ;
  AT_Corr[7][0]=-.026 ;
  AT_Corr[7][1]= -.01;
  AT_Corr[8][0]= .089;
  AT_Corr[8][1]= -.0094;
  AT_Corr[9][0]= -.078;
  AT_Corr[9][1]=-.01 ;
  AT_Corr[10][0]= .0097;
  AT_Corr[10][1]= -.01;
  AT_Corr[11][0]=-0.017 ;
  AT_Corr[11][1]=-.012 ;
  ////////////////////////////////////////////////////////////
  
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
  top_corr[5]=-3.6;
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

  fileout = fopen ("data/comb_r3505.dat","r"); // combined data file
  
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
	  cout << "found bottom hit: " << channel << endl;
	  tdcinf[1]=tdcinf[1]+bot_corr[channel];
	  bot_array.push_back(tdcinf);
	  //	hTDC[channel]->Fill(cont);
	}
	// 1/2 middle counters
	if (channel==12) {
	  if (cont >=150 && cont <= 200){
	    tdc_padw[n_padw]=cont;
	    n_padw++;
	  }
	  // hPad[0]->Fill(cont);
	}
	// the other middle counter
	if (channel==13) {
	  if (cont >=150 && cont <= 200) {
	    tdc_pade[n_pade]=cont +2; //2 for a centering
	    n_pade++;
	  }
	  //	hPad[1]->Fill(cont);
	}  
	if (channel>15 && channel<28) {
	  tdcinf[1]=tdcinf[1]+top_corr[channel-16];
	  //	hTDC[channel-4]->Fill(cont);    
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
      printf(" board: %d ch:%d pdo:%d tdo:%d\n",mm_board,mm_ch,mm_pdo,mm_tdo);

      // NOISY CHANNEL
      //      if (mm_vmm==7 && mm_ch==43 && mm_board == 110) continue;
      bID = oct.getBoardIndex(mm_board,nrun);
      pdo_curr=float(mm_pdo);
      tdo_curr=float(mm_tdo);
        
      //double pc= pdo_cal[bID][mm_vmm][mm_ch][0]; 
      //double gc= pdo_cal[bID][mm_vmm][mm_ch][1];
      double pt= tdo_cal[bID][mm_vmm][mm_ch][0];
      double gt= tdo_cal[bID][mm_vmm][mm_ch][1];

      pdo_curr=PDO2Charge.GetCharge(mm_pdo, mm_board, mm_vmm, mm_ch);
      cout << "PDO! " << pdo_curr << endl;
      tdo_curr=(tdo_curr-pt)/gt;

      // old calibration format
      // pdo_curr=(pdo_curr-pc)/gc;
      // tdo_curr=(tdo_curr-pt)/gt;

      mminf[0]=mm_ch + mm_vmm*64;
      mminf[1]=pdo_curr;
      mminf[2]=tdo_curr;
      mm_array[bID].push_back(mminf);
      // if(mm_board==105) numprch[0]++;
      // if(mm_board==17) numprch[1]++;
    }
    fscanf(fileout, "%s",strg3); //end event readout
    cout << strg3 << endl;
    
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
    printf("reordered\n");
    cout << "nbhits: " << nbhit << " nthits: " << nthit << endl;
 
    //////////// look for W+E pairs
    int bothits = s.pairs(bot_array,120.,170.);
    cout << "Num bottom pairs: " << bothits << endl;
    int tophits= s.pairs(top_array,130.,170.);
    cout << "Num top pairs: " << tophits << endl;

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
	       length=sqrt(247.*247.+(botxpos-topxpos)*(botxpos-topxpos)+
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
	 if(abs(top_time-bot_time)<=time_cut) { // normallly <=4
	   
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
       p.forward(mm_array,mmhits,5,TDO_UPPER,k,pad_at,6.,2.5); 
     } // end cluster search
     
     //hnclus->Fill(nclus);
     //   hnclus->Fill(nclus+10*(1-left));
     
     // chomp backwards
     for (int k=0; k<NUMBOARD; k++) {
       p.backward(mm_array,mmhits,5,TDO_UPPER,k,pad_at,2.5);
       //       p.backward(mm_array,mmhits,5,TDO_UPPER,k,pad_at,2.5);
     }

     // calculate average strip for each cluster, saved in mm_array[k][i][8]
     double xav; 
     double run_clus=0.;
     int ind1=0;
     vector<double> chacont(2);
     vector< vector<double> > cluscha_ord[NUMBOARD];
   
     for(int k=0; k<NUMBOARD; k++) {
       for(int i=0; i<mmhits[k]; i++) {
	 mm_array[k][i][8] = 0;
	 ind1=i;
	 if( mm_array[k][i][3] > run_clus) {
	   xav = 0.; 
	   run_clus = mm_array[k][i][3]; 
	   chacont[0]= mm_array[k][i][6]; // cluster total charge
	   chacont[1]= mm_array[k][i][3]; // cluster number
	   cluscha_ord[k].push_back(chacont);
	   for(int j=ind1; j<mmhits[k]; j++) {
	     if (mm_array[k][j][3] == run_clus)  xav = xav + mm_array[k][j][0] * mm_array[k][j][1]; //channels avg weighted by pdo
	   }
	   xav=xav/mm_array[k][i][6];
	   mm_array[k][i][8] = xav; 
	 }
       }
     }
     
     // reorder clusters according to most charge to least charge
     for (int k=0; k<NUMBOARD; k++){
       sort(cluscha_ord[k].begin(), cluscha_ord[k].end(), oct.compare1);
       int num_clus = cluscha_ord[k].size(); // number of clusters for that board
       int irun;
       for(int j=0; j<num_clus; j++) { 
	 run_clus = cluscha_ord[k][j][1]; //cluster number
	 irun = 0;
	 for(int i=0; i<mmhits[k]; i++) {
	   if(irun == 0 && mm_array[k][i][3]==run_clus) {
	     mm_array[k][i][9]=j+1; //save what order the cluster is in the cluster ranking
	   irun++; 
	   }
	 }
       }
     }

     // count holes
     int is_hole, istep;
     for (int k=0; k<NUMBOARD; k++){
       for(int i=0; i<mmhits[k]; i++) { 
	 if ( mm_array[k][i][9]==1. ) { //highest charge loop
	   cout << "passed charge loop" << endl;
	   istep = mm_array[k][i][5]; //cluster end
	   is_hole = mm_array[k][istep][0]-mm_array[k][i][0]+1.-mm_array[k][i][4]; // number of holes
	   cout << "passed setting hole" << endl;
	   hClusterHoles[k]->Fill(is_hole);
	   hClusterCharge[k]->Fill(mm_array[k][i][6]);
	   hClusterMult[k]->Fill(mm_array[k][i][4]);
	   cout << "passed filling" << endl;
	 }
       }
     }

   int d_strip;
   double dt_strip;

   int ncntr_bot; // counts the bottom counter, indexed from 0 to 5
   // find track's angle from ybot
   ncntr_bot = botypos/20. + 2.5;
   // track_ang=(float(ncntr_bot)*5.+90.)/180.*pi;
//    m_av=TMath::Tan(track_ang);
//    track_ang=(float(ncntr_bot)*5.-2.5+90.)/180.*pi;
//    m_min=TMath::Tan(track_ang);
//    track_ang=(float(ncntr_bot)*5.+2.5+90.)/180.*pi;
//    m_max=TMath::Tan(track_ang);
//   ///////////////////////////
//    htslope->Reset();
//    hzslope->Reset();
//    for(int i=0; i<mmhits[0]; i++) { //loop 1
//      if ( mm_array[0][i][9] ==1. && mm_array[0][i][8]>=8. && mm_array[0][i][8]<=119. ) { // this line for a test loop 2 add fiducial
       
//        hmtpc_2->Fill(botypos,mm_array[0][i][4]); 
//        istep= mm_array[0][i][5];
//        is_hole=mm_array[0][istep][0]-mm_array[0][i][0]+1.-mm_array[0][i][4];
//        /* 
// 	 if(is_hole>=1) {     
// 	 printf("HolyEvent (%i,%i,%i) with %i holes in primary cluster!!!!\n",nev, mm_nev,is_hole);  
// 	 for(int l=0; l<mmhits;l++) {
// 	 printf("%5.1f %5.1f  %5.1f  %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f  \n",
// 	 mm_array[0][l][0], mm_array[0][l][1], mm_array[0][l][2],
// 	 mm_array[0][l][3], mm_array[0][l][4], mm_array[0][l][5], mm_array[0][l][6],
// 	 mm_array[0][l][7], mm_array[0][l][8], mm_array[0][l][9] );
// 	 }
	 
// 	 }
//       */ 
//       hmtpc_3->Fill(botypos,is_hole);
//       x_clu=0.2+mm_array[0][i][8]*0.4;
//       clu_mult=mm_array[0][i][4];
      
//       for(int j=i; j<istep+1; j++) {//loop3
// 	if(mm_array[0][j][3]>0) {// loop 4
// 	  ///////////// time profile

// 	  d_strip=mm_array[0][j][0]-mm_array[0][i][8] ;
// 	  dt_strip=mm_array[0][j][2]-mm_array[0][i][7] ;
// 	  if(botypos>20. && mm_array[0][i][4] >=3) hdtime->Fill(d_strip,dt_strip);
// 	  hsng_befcorr->Fill(mm_array[0][j][1],mm_array[0][j][2]);
// 	  ///////////////////////////////////////////
// 	  //correct for slewing
// 	  if(nrun == 3079 || nrun== 3092 || nrun==3104) {
// 	    corr=5.-0.084*mm_array[0][j][1];
// 	    if(mm_array[0][j][1]  < 10.) corr=corr-105.6*TMath::Exp(-0.444*mm_array[0][j][1]);
// 	  }
// 	  else if(nrun==3106) corr=8.-0.1306*mm_array[0][j][1];
// 	  else{
// 	    corr=3.564-0.0594*mm_array[0][j][1];
// 	    if(mm_array[0][j][1]  < 10.) corr=corr-74.9*TMath::Exp(-0.266*mm_array[0][j][1]);
// 	  }
	  
// 	  mm_array[0][j][2]= mm_array[0][j][2]+corr;
// 	  havtvsx->Fill(mm_array[0][j][0],  mm_array[0][j][2]-0.5*pad_at+120.);
// 	  if(botypos>20.) havtvsx23->Fill(mm_array[0][j][0],  mm_array[0][j][2]-0.5*pad_at+120.);
// 	  hsng_tvsch->Fill(mm_array[0][j][1],mm_array[0][j][2]-0.5*pad_at+120. );
// 	  x_mm=mm_array[0][j][0]*0.4+0.2;
	  
// 	  z_mm=t0_off-(mm_array[0][j][2]-0.5*pad_at+120.)/20.; //
	  
// 	  hmtpc_ck->Fill(x_mm,z_mm);
// 	  if(abs(d_strip)<=4.) {
// 	    zx[0]= x_mm;
// 	    zx[1]= z_mm;
// 	    clus_coord.push_back(zx);
// 	    truec_coord.push_back(zx);
// 	    htslope->Fill(mm_array[0][j][0],mm_array[0][j][2]-0.5*pad_at+120.  );
// 	    htslope->SetBinError(int(mm_array[0][j][0])+1, 5.);
// 	    hzslope->Fill(x_mm,z_mm);
// 	    hzslope->SetBinError(int(x_mm/0.4)+1,0.25);
	    
//    //double bincon=htslope->GetBinContent(int(mm_array[0][j][0])+1);
//    //double binerr=htslope->GetBinError(int(mm_array[0][j][0])+1);
// 	    if( i==istep) {
// 	      hmtpc_1->Fill(botypos,mm_array[0][i][0]);
// 	    }
// 	    npoint++;
// 	  }
// 	  // printf(" strip %5.1f charge %5.1f time %5.1f %5.1f %5.1f \n", mm_array[0][j][0],mm_array[0][j][1],mm_array[0][j][2],
//  //bincon,binerr);
// 	}//loop 4
// 	if(mm_array[0][j][0] <64 && mm_array[0][j][3]>0. && mm_array[0][j][2] >ftclu0) ftclu0= mm_array[0][j][2];
// 	if(mm_array[0][j][0] >=64 && mm_array[0][j][3]>0. && mm_array[0][j][2] >ftclu1) ftclu1= mm_array[0][j][2];
// 	if( mm_array[0][j][3]>0. && mm_array[0][j][2] >ftclu) ftclu= mm_array[0][j][2];
//       } //loop3
//   hearlytvsx->Fill(mm_array[0][i][8],ftclu);
//     }//loop2
//   } //loop1
//   hart0->Fill(ftclu0,art0);
//   hart1->Fill(ftclu1,art1);
//   if(ftclu0>110.) hartdiff0->Fill(art0-ftclu0);
//   if(ftclu1>110.) hartdiff1->Fill(art1-ftclu1);
//   ///////////////////////////////////////////////
  
//   ////////////////////////////fit tslope
//   double nang,ang,tgthe;
//   TF1 *fttr= new TF1("fttr", "pol1(0)");
//   TF1 *fttr1= new TF1("fttr1", "pol1(0)");
//   TF1 *fttr2= new TF1("fttr2", "pol1(0)");
//   fttr->SetParameter(0,400.);
//   fttr->SetParameter(1,5.);
//   fttr1->SetParameter(0,0.);
//   fttr1->SetParameter(1,1.);
//   fttr2->SetParameter(0,0.);
 
//   if(is_hole>=0) {
//     if(npoint>=2) {
//       ///////////// angle
//       htslope->Fit(fttr,"Q","",0.,130.);
//       /////////////////////////////////
      
//       double c1 = fttr->GetParameter(0);
//       double c2 = fttr->GetParameter(1);
//       //	printf("fit p1: %5.1f p2:%5.1f \n",c1,c2); 
//       double xx=float(npoint-1)*0.4;
//       double yy=-c2/20.*float(npoint-1);
//         tgthe=yy/xx;
// 	nang=TMath::ATan2(yy,xx)*180./pi;
//         if(nang>=0.) ang=90.-nang;
// 	else ang=-90.-nang;
// 	//printf("fit p1: %5.3f p2:%5.3f tgthe %5.3f ang %5.3f   : \n",c1,c2,tgthe,ang);
// 	////////////////////////////////////////////
// 	/////////////////////////// zcluster
// 	hzslope->Fit(fttr1,"Q","",0.,51.2);
// 	c1 = fttr1->GetParameter(0);
//         c2 = fttr1->GetParameter(1);
// 	chi2= fttr1->GetChisquare();
// 	xclu_tpc=(2.5-c1)/c2;
// 	//////////////////////////////////////
// 	/// residual check
// 	if(npoint>=3) {	
// 	  for(int l=0; l<npoint;l++) {
// 	    zfit=c1+c2*clus_coord[l][0];
// 	    deltaz=(clus_coord[l][1]-zfit)/0.25;
// 	    adeltaz=abs(deltaz);
// 	    clus_coord[l][2]=adeltaz;
// 	    hresid->Fill(botypos,deltaz);
// 	}
// 	}
// 	//////////////////////////////////////////////
// 	// event display
// 	if(npoint>=5 && botypos>=20. && display_on==1) {
// 	  des->cd(0);
// 	  hzslope->Draw();
// 	TLine* line = new TLine(0.0, 2.5, 51.2, 2.5);
// 	line->Draw("same");
// 	TLine* line1 = new TLine(x_clu, -4., x_clu, 10.);
//   	line1->Draw("same");   
// 	des->Update();
//         for(int k=0; k<npoint; k++) 
// 	  printf(" here1 %d %5.3f %5.3f %5.3f \n ",an_nev,clus_coord[k][0],clus_coord[k][1],clus_coord[k][2]);
//   printf(" write 1 to continue 0 to end\n");
//   scanf("%d",&goon);
//   if(goon==0) return 0;
// 	}
	
// 	//////////////////////////////////////////////////////// 
// 	/////////////////// pruning zslope
 
// 	sort(clus_coord.begin(),clus_coord.end(), compare2);
// 	int remove=0.;
// 	int ckp;
	
// 	hpruslope->Reset();
	
//     for(int kk=0; kk<npoint-3;kk++) {
//       remove=0;
//       if(remove==0 && clus_coord[kk][2]>4. )	{
// 	clus_coord[kk][0]= -1.;
//    clus_coord[kk][2]=100.;
//    remove++;
//       }   
      
//       if(remove>=1) {
// 	hpruslope->Reset();
//     for(int l=0; l<npoint; l++) {
//       if(clus_coord[l][0]>0.) {
// 	hpruslope->Fill(clus_coord[l][0],clus_coord[l][1]);
// 	hpruslope->SetBinError(int(clus_coord[l][0]/0.4)+1,0.25);
//       }   
//     }
//     fttr1->SetParameter(0,0.);
//      fttr1->SetParameter(1,1.);
//      hpruslope->Fit(fttr1,"Q","",0.,51.2);
//      c1 = fttr1->GetParameter(0);
//      c2 = fttr1->GetParameter(1);
//      chi2= fttr1->GetChisquare();
//      xclu_tpc=(2.5-c1)/c2;
     
//      for(int l=kk+1; l<npoint;l++) {
//        zfit=c1+c2*clus_coord[l][0];
//        deltaz=(clus_coord[l][1]-zfit)/0.25;
//        adeltaz=abs(deltaz);
//        clus_coord[l][2]=adeltaz;
//      }
//      sort(clus_coord.begin(),clus_coord.end(), compare2);
     
//       } //remove loop
      
//     } //kk loop
    
//     ////////// 
    
//     // event display 
//     if(npoint>=5 && botypos>=20. && display_on==1) {
//       prun->cd(0);
//       hpruslope->Draw();
//       TLine* line = new TLine(0., 2.5, 51.2, 2.5);
//       line->Draw("same");
//       TLine* line1 = new TLine(x_clu, -4., x_clu, 10.);
//       line1->Draw("same");     
//       prun->Update();  
// 	printf("reordered \n");
// 	for(int k=0; k<npoint; k++)  
// 	  printf("here 2 %d %5.3f %5.3f %5.3f \n ",an_nev,clus_coord[k][0],clus_coord[k][1],clus_coord[k][2]);
	
// 	printf(" write 1 to continue 0 to end\n");
// 	scanf("%d",&goon);
// 	if(goon==0) return 0;
//     }
    
//     /////////////////////////////////////////
//     //printf("fit p1: %5.3f p2:%5.3f x_tpc: %5.3f x_clu: %5.3f   : \n",c1,c2,xclu_tpc,x_clu); 
//     ///////////////////////////////////
// 	//////// fit with constrained angle
//     fttr2->FixParameter(1,m_av);
//     hzslope->Fit(fttr2,"Q","",0.,51.2);
//     c1 = fttr2->GetParameter(0);
//     c2 = fttr2->GetParameter(1);
//     xconst=(2.5-c1)/c2;
//  // event display
//     if(npoint>=5 && botypos>=20. && display_on==1) {
//       fixsl->cd(0);
//       hzslope->Draw();
//       TLine* line = new TLine(0.0, 2.5, 51.2, 2.5);
// 	line->Draw("same");
// 	TLine* line1 = new TLine(x_clu, -4., x_clu, 10.);
//   	line1->Draw("same");   
// 	fixsl->Update();
// 	scanf("%d",&goon);
// 	if (goon==2) {
// 	  des->Print(Form("des_%i.png",nev));
// 	  prun->Print(Form("prun_%i.png",nev));
// 	  fixsl->Print(Form("fixsl_%i.png",nev));
//     }
// 	if(goon==0) return 0;
//     }
//     ///////////////////////////////
//     /////////////RESIDUAL CHECK AND  PRUNING
//     if(npoint>=3) {	
// 	for(int l=0; l<npoint;l++) {
// 	  zfit=c1+c2*truec_coord[l][0];
//           deltaz=(truec_coord[l][1]-zfit)/0.25;
//           adeltaz=abs(deltaz);
//           truec_coord[l][2]=adeltaz;
// 	  hresidtrue->Fill(botypos,deltaz);
// 	}
//     }
    
//     sort(truec_coord.begin(),truec_coord.end(), compare2);
//     remove=0.;
//     ckp=npoint-remove;
//     if(ckp>3) {
//       for(int k=0; k<npoint; k++){
// 	if(truec_coord[k][2]>4. && npoint-remove>3) {
// 	  truec_coord[k][0]= -1.;       
// 	remove++;
// 	}
//       }
//     }
//     hprutrue->Reset();
//     if(remove>=1) {
      
//     for(int l=0; l<npoint; l++) {
//       if(truec_coord[l][0]>0.) {
// 	hprutrue->Fill(truec_coord[l][0],truec_coord[l][1]);
// 	hprutrue->SetBinError(int(truec_coord[l][0]/0.4)+1,0.25);
//       }   
//  }
//     hprutrue->Fit(fttr2,"Q","",0.,51.2);
//     c1 = fttr2->GetParameter(0);
//     c2 = fttr2->GetParameter(1);
// 	xconst=(2.5-c1)/c2;
//     } 
//     ////////////////////////////// ev display
    
//     ///////////////////////////////////////////////////
//     ////////////////////////////////////////////////////////
//     ///// fix angle to errors 
//     fttr2->FixParameter(1,m_min);
//          hzslope->Fit(fttr2,"Q","",0.,51.2);
// 	 c1 = fttr2->GetParameter(0);
// 	 c2 = fttr2->GetParameter(1);
// 	 xconst_min=(2.5-c1)/c2;
// 	 fttr2->FixParameter(1,m_max);
//          hzslope->Fit(fttr2,"Q","",0.,51.2);
//        c1 = fttr2->GetParameter(0);
//        c2 = fttr2->GetParameter(1);
//        xconst_max=(2.5-c1)/c2;
//        xpos_res=xconst_max-xconst_min;
// 	//////////////////////////// fill histos with results
//        if(botypos>20.)  hcludiffvsmu->Fill(clu_mult,x_clu-xclu_tpc);
//        if(npoint>=3) {
// 	  hchid->Fill(chi2);
// 	  if(botypos>20.)   hxclvstpcge3->Fill(x_clu,xclu_tpc);
// 	  hztpc->Fill(botypos,x_clu-xclu_tpc); 	  
// 	  hftpc->Fill(botypos,x_clu-xconst);
// 	  hfftpc->Fill(botypos,xclu_tpc-xconst);
// 	  hftpc_res->Fill(botypos,(x_clu-xconst)/xpos_res);
// 	  hxpres->Fill(botypos,xpos_res);
// 	  hmtpc->Fill(botypos,ang);
// 	  htgvsy->Fill(botypos,tgthe);
// 	}
// 	else {
// 	  if(botypos>20. )  hxclvstpc2->Fill(x_clu,xclu_tpc);
// 	  hftpc2->Fill(botypos,x_clu-xconst);
//           hfftpc2->Fill(botypos,xclu_tpc-xconst);
// 	  hftpc_res2->Fill(botypos,(x_clu-xconst)/xpos_res);
//      hxpres2->Fill(botypos,xpos_res);
//     hztpc_2strip->Fill(botypos,x_clu-xclu_tpc); 	  
//   hmtpc_2strip->Fill(botypos,ang);
// htgvsy_2strip->Fill(botypos,tgthe);
// 	}
   
//     }
//       }
//   //////////////////////////////////////
//   for(int i=0; i<mmhits[0]; i++) {
//     if ( mm_array[0][i][9] ==1.) {
//       //////////////////////////
//       istep= mm_array[0][i][5];
//       hydep->Fill(botypos,mm_array[0][i][8]);
//       hxdep->Fill(botxpos,mm_array[0][i][8]);
//   for(int j=0; j< mmhits[0]; j++) {
//      double xdist=mm_array[0][j][0]-mm_array[0][i][8];
//      hclus_dis->Fill(xdist,mm_array[0][j][1]/mm_array[0][i][6]);
//    }
//   ////////////////////////////////////////////////
//  if(mm_array[0][i][4] >=1) {
 
//       if(mm_array[0][i][8]<64.)  hclustvscha_0->Fill(mm_array[0][i][6], mm_array[0][i][7]-0.5*pad_at+120.);
//       if(mm_array[0][i][8]>=64.)  hclustvscha_1->Fill(mm_array[0][i][6], mm_array[0][i][7]-0.5*pad_at+120.);
//       if(mm_array[0][i][8]<64.)  hnstrpvscha_0->Fill(mm_array[0][i][6], mm_array[0][i][4]);
//       if(mm_array[0][i][8]>=64.)  hnstrpvscha_1->Fill(mm_array[0][i][6], mm_array[0][i][4]);
//       hchavsxclu->Fill(mm_array[0][i][8],mm_array[0][i][6]);
//       hxdis->Fill( mm_array[0][i][8],mm_array[0][i][6]);
//      if(mm_array[0][i][8]<64.)  htvsnst0->Fill( mm_array[0][i][4], mm_array[0][i][7]-0.5*pad_at+120.);
//       if(mm_array[0][i][8]>=64.)  htvsnst1->Fill( mm_array[0][i][4], mm_array[0][i][7]-0.5*pad_at+120.);
  
// }
//  ///////////////////////////////
//       if(mm_array[0][i][4] ==1) hxdis_1->Fill(mm_array[0][i][8],mm_array[0][i][6]);
//       if(mm_array[0][i][4] ==2) {
//       if(mm_array[0][istep][0]-mm_array[0][i][0]+1. ==mm_array[0][i][4])  hxdis_2->Fill(mm_array[0][i][8],mm_array[0][i][6]);
//       else hxdis_2br->Fill(mm_array[0][i][8],mm_array[0][i][6]);
//       }
//     if (mm_array[0][i][4] >2) {
//       if(mm_array[0][istep][0]-mm_array[0][i][0]+1. ==mm_array[0][i][4])  hxdis_ge3->Fill(mm_array[0][i][8],mm_array[0][i][6]);
//       else hxdis_ge3br->Fill(mm_array[0][i][8],mm_array[0][i][6]);
//     }
 
	
//     }
//   }

// ALL THE WAY TO HERE 
    ///////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
   //aw
 } // end mm study
  } // end while
 fclose(fileout);
 printf("event analyzed %d ev sel %d ev cut %d ev with mm %d\n",an_nev,sel_ev,ret_ev,ch_ev);

  // write histo file
 TFile *fout=new TFile("test.root","recreate");
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
 for (int k=0; k<NUMBOARD; k++) {
   channel_hits[k]->Write();
   hCHvsPDO[k]->Write();
   hClusterCharge[k]->Write();
   hClusterMult[k]->Write();
   hClusterHoles[k]->Write();
 }

// hchthre->Write(); 
// hTCvNSL->Write();
// hTCvNSR->Write();
// hchvstl->Write();
// hchvstr->Write();
// hchvstl_raw->Write();
// hchvstr_raw->Write();
// hnclusvstcha->Write();
// houtlier->Write();
//  hxdis->Write();
//  hxdis_1->Write();
//  hxdis_2->Write();
//  hxdis_2br->Write();
//  hxdis_ge3->Write();
//  hxdis_ge3br->Write();
// hclustvscha_0->Write();
// hclustvscha_1->Write();
// hclustcha_01->Write();
// hclustcha_11->Write();
// hnstrpvscha_0->Write();
// hnstrpvscha_1->Write();
// hchalvsr->Write();
 hnclus->Write();
// hmiss_chal->Write();
// hmiss_char->Write();
// hchavsxclu->Write();
// hmax_cha->Write();
// htvsnst0->Write();
// htvsnst1->Write();
// havtvsx->Write();
// havtvsx23->Write();
// hart0->Write();
// hart->Write();
// hart1->Write();
// hartdiff0->Write();
// hartdiff1->Write();
// hclus_dis->Write();
// hearlytvsx->Write();
// hxdep->Write();
// hydep->Write();
// hmtpc->Write();
// hmtpc_2strip->Write();
// hmtpc_1->Write();
// hmtpc_2->Write();
// hmtpc_3->Write();
// hmtpc_ck->Write();
// htgvsy->Write();
// htgvsy_2strip->Write();
// hsng_tvsch->Write();
// hsng_befcorr->Write();
// hdtime->Write();
// hzslope->Write();
// hztpc->Write();
// hztpc_2strip->Write();
// hxclvstpc2->Write();
// hxclvstpcge3->Write();
// hcludiffvsmu->Write();
// hftpc->Write();
// hftpc_res->Write();
// hftpc2->Write();
// hftpc_res2->Write();
// hfftpc->Write();
// hfftpc2->Write();
// hxpres->Write();
// hxpres2->Write();
// hresid->Write();
// hresidtrue->Write();
// hpruslope->Write();
// hprutrue->Write();
// hchid->Write();
 fout->Close();
 cout << "successfully closed" << endl;

 // delete pointers
 for (int  i=0; i<25;i++) {
   delete hTDC[i];
 }
 for (int  i=0; i<12;i++) {
   delete hATDC[i],hDTDC[i],hDVSATDC[i]; 
 }
 for (int  i = 0; i<14;i++) {
   delete hSTDC[i];
 }
 delete hPad[0];
 delete hPad[1];
 for (int k = 0; k<NUMBOARD; k++) {
   delete channel_hits[k];
 }
 return 0;
}

    
