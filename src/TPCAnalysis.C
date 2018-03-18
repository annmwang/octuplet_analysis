///
///  \file   TPCAnalysis.C
///
///  \author Tunaface
///          (tuna@cern.ch)
///
///  \date   2017 Dec
///

#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TFitResult.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <numeric>
#include <chrono>

#include "include/MMPlot.hh"
#include "include/PDOToCharge.hh"
#include "include/TDOToTime.hh"
#include "include/MMDataAnalysis.hh"
#include "include/MMClusterAlgo.hh"
#include "include/MMPacmanAlgo.hh"
#include "include/GeoOctuplet.hh"
#include "include/GeoPlane.hh"
#include "include/SimpleTrackFitter.hh"
#include "include/ScintillatorClusterFilterer.hh"

using namespace std;

// helper functions
void progress(double time_diff, int nprocessed, int ntotal);
double channel_from_x(double xpos, int board, GeoOctuplet* geo, int begin);
double theta(double slope);
std::tuple<double, double> fit(std::vector<double> xs, std::vector<double> zs);
std::tuple<double, double, double, double, double, double, double, double, double> fit_root(std::vector<double> xs, std::vector<double> zs, int invert);
double TPCUncertainty(double ztpc, double m, double b, double sigma2_b, double sigma2_m, double sigma_bm, int invert);
double AngleCorrection(int i, int j, double slope, GeoOctuplet* GEOMETRY);
int OkayForTPC(MMCluster* clus);
//std::map< std::tuple<int,int>, double > ChannelOffsets(std::string filename, int board);
std::vector< double > ChannelOffsets(std::string filename, int board, int this_vmm);
double OffsetByChannel(int vmm, int ch, std::map< std::tuple<int,int>, double > offsets);
double OffsetByBoard(int board);

int main(int argc, char* argv[]){

  char inputFileName[400];
  char outputFileName[400];
  char PDOFileName[400];
  char TDOFileName[400];
  char AlignFileName[400];
  
  if ( argc < 5 ){
    cout << "Error at Input: please specify input/output .root files ";
    cout << " and (optional) PDO/TDO calibration files" << endl;
    cout << "Example:   ./RunMMAnalysis.x -i input.root -o output.root" << endl;
    cout << "Example:   ./RunMMAnalysis.x -i input.root -o output.root";
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
  bool b_align = false;
  bool b_paolo = false;
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
    if (strncmp(argv[i+1],"-g",2)==0){
      b_paolo = true;
    }
  }

  if(!b_input) std::cout << "Error at Input: please specify  input file (-i flag)" << std::endl;
  if(!b_out)   std::cout << "Error at Input: please specify output file (-o flag)" << std::endl;
  if(!b_input || !b_out) return 0;

  // retrieve the channel offsets
  std::string fname = "/Users/alexandertuna/Downloads/octuplet_analysis/paolo_20180215/drift_offst.txt";
  std::vector< std::vector< std::vector<double> > > channel_offsets = {};
  //std::vector< std::map< std::tuple<int,int>, double > > channel_offsets = {};
  for (int bo = 0; bo < 8; bo++){
    channel_offsets.push_back( std::vector< std::vector<double> >() );
    for (int vm = 0; vm < 8; vm++)
      channel_offsets[bo].push_back( ChannelOffsets(fname, bo, vm) );
  }

  const int nboards = 8;
  int nboardshit = 0;
  int i = 0, ibo = 0, ich = 0, test = 0;
  //int jbo = 0;
  int debug = 0;

  // class defs
  PDOToCharge* PDOCalibrator;
  TDOToTime*   TDOCalibrator;
  MMPacmanAlgo*                PACMAN       = new MMPacmanAlgo(5,2.,0.5);
  GeoOctuplet*                 GEOMETRY     = new GeoOctuplet();
  SimpleTrackFitter*           FITTER       = new SimpleTrackFitter();
  ScintillatorClusterFilterer* FILTERER     = new ScintillatorClusterFilterer();
  MMDataAnalysis*    DATA;
  if(b_pdo) PDOCalibrator = new PDOToCharge(PDOFileName);
  else      PDOCalibrator = new PDOToCharge();
  if(b_tdo) TDOCalibrator = new TDOToTime(TDOFileName);
  else      TDOCalibrator = new TDOToTime();

  if(b_align)
    GEOMETRY->SetAlignment(AlignFileName);
  if(b_paolo)
    GEOMETRY->SetAlignmentPaolo(1);

  // zboard, post-alignment
  std::vector<double> zboard = {};
  for (int i = 0; i < GEOMETRY->GetNPlanes(); i++)
    zboard.push_back(GEOMETRY->Get(i).Origin().Z() + 2.7*(i%2==0 ? 1 : -1));

  TFile* f = new TFile(inputFileName, "READ");
  if(!f){
    std::cout << "Error: unable to open input file " << inputFileName << std::endl;
    return 0;
  }
  TTree* T = (TTree*) f->Get("COMB_data");
  if(!T){
    std::cout << "Error: cannot find tree COMB_data in " << inputFileName << std::endl;
    return 0;
  }

  DATA = (MMDataAnalysis*) new MMDataAnalysis(T);
  DATA->SetTP(0);
  int Nevent = DATA->GetNEntries();

  std::map< string, TH1D* > h1;
  std::map< string, TH2D* > h2;

  h2["track_res2_vs_y"]  = new TH2D("track_res2_vs_y",  " ;y;Sqrt(Sum of res. sq.); Tracks", 220, 0, 220,    100, 0.0,  5.0);
  h2["track_res2_vs_my"] = new TH2D("track_res2_vs_my", ";my;Sqrt(Sum of res. sq.); Tracks", 100, -4.0, 4.0, 100, 0.0,  5.0);
  h2["track_res2_vs_cy"] = new TH2D("track_res2_vs_cy", ";cy;Sqrt(Sum of res. sq.); Tracks", 100, -300, 500, 100, 0.0,  5.0);

  h2["scint_bot_vs_top"] = new TH2D("scint_bot_vs_top", ";SC bottom channel;SC top channel;Tracks", 8, -1.5, 6.5, 8, 14.5, 22.5);
  h2["scint_ch_vs_counts"] = new TH2D("scint_ch_vs_counts", ";SC CH;SC Counts;", 30, -1.5, 28.5, 400, 0, 400);

  h2["track_angle_x_scint_19"] = new TH2D("track_angle_x_scint_19", ";x angle;scint. bottom channel;Tracks, scint. top channel 19", 100, -30, 30, 8, -1.5, 6.5);
  h2["track_angle_x_scint_20"] = new TH2D("track_angle_x_scint_20", ";x angle;scint. bottom channel;Tracks, scint. top channel 20", 100, -30, 30, 8, -1.5, 6.5);
  h2["track_angle_x_scint_21"] = new TH2D("track_angle_x_scint_21", ";x angle;scint. bottom channel;Tracks, scint. top channel 21", 100, -30, 30, 8, -1.5, 6.5);
  h2["track_scint_0_vs_time"]  = new TH2D("track_scint_0_vs_time", ";Time [days];track angle, bottom scint. 0;Tracks", 1000, 0, 20, 100, -30, 30);
  h2["track_scint_1_vs_time"]  = new TH2D("track_scint_1_vs_time", ";Time [days];track angle, bottom scint. 1;Tracks", 1000, 0, 20, 100, -30, 30);
  h2["track_scint_2_vs_time"]  = new TH2D("track_scint_2_vs_time", ";Time [days];track angle, bottom scint. 2;Tracks", 1000, 0, 20, 100, -30, 30);
  h2["track_scint_3_vs_time"]  = new TH2D("track_scint_3_vs_time", ";Time [days];track angle, bottom scint. 3;Tracks", 1000, 0, 20, 100, -30, 30);
  h2["track_scint_4_vs_time"]  = new TH2D("track_scint_4_vs_time", ";Time [days];track angle, bottom scint. 4;Tracks", 1000, 0, 20, 100, -30, 30);
  h2["track_scint_5_vs_time"]  = new TH2D("track_scint_5_vs_time", ";Time [days];track angle, bottom scint. 5;Tracks", 1000, 0, 20, 100, -30, 30);

  h2["track_hits_v_board"]      = new TH2D("track_hits_v_board", ";Board;N(Clusters);Tracks", 8, -0.5, 7.5, 6, -0.5, 5.5);
  h2["track_hits_vs_time"]      = new TH2D("track_hits_vs_time",     ";Time [days];N(clusters in track);", 1000, 0, 30, 10, -0.5, 9.5);
  h2["track_hits_vs_time_fid"]  = new TH2D("track_hits_vs_time_fid", ";Time [days];N(clusters in track);", 1000, 0, 30, 10, -0.5, 9.5);
  h2["track_clusmult_vs_theta"] = new TH2D("track_clusmult_vs_theta", ";#theta#lower[0.5]{track} [degrees];hits in cluster", 100, -35, 35, 13, -0.5, 12.5);

  h2["track_param_x"] = new TH2D("track_param_x", ";x slope; x constant;Tracks",   100, -0.6, 0.6, 100,  -80, 280);
  h2["track_param_y"] = new TH2D("track_param_y", ";y slope; y constant;Tracks",   100,   -4,   4, 100, -300, 500);
  h2["track_angle_6"] = new TH2D("track_angle_6", ";#theta#lower[0.5]{track} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);
  h2["track_angle_7"] = new TH2D("track_angle_7", ";#theta#lower[0.5]{track} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);
  h2["track_angle_8"] = new TH2D("track_angle_8", ";#theta#lower[0.5]{track} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);

  h2["track_N1_board_vs_residual"]      = new TH2D("track_N1_board_vs_residual",      ";board;x_{cluster} - x_{track, proj.};Tracks", 8, -0.5, 7.5, 200, -5.0, 5.0);
  h2["track_N1_board_vs_residual_norm"] = new TH2D("track_N1_board_vs_residual_norm", ";board;x_{cluster} - x_{track, proj.};Tracks", 8, -0.5, 7.5, 200, -5.0, 5.0);
  h2["track_N1_board_vs_utpc"]          = new TH2D("track_N1_board_vs_utpc",          ";board;x_{utpc} - x_{track, proj.};   Tracks", 8, -0.5, 7.5, 200, -5.0, 5.0);
  h2["track_N1_board_vs_prob"]          = new TH2D("track_N1_board_vs_prob",          ";board;fit probability;   Tracks",             8, -0.5, 7.5, 200, -0.1, 1.1);

  h1["track_diff01_bary"] = new TH1D("track_diff01_bary", ";x_{bary,0} - x_{bary,1}; Tracks", 200, -5, 5);
  h1["track_diff67_bary"] = new TH1D("track_diff67_bary", ";x_{bary,6} - x_{bary,7}; Tracks", 200, -5, 5);
  h1["track_diff01_btpc"] = new TH1D("track_diff01_btpc", ";x_{bary,0} - x_{bary,1}; Tracks", 200, -5, 5);
  h1["track_diff67_btpc"] = new TH1D("track_diff67_btpc", ";x_{bary,6} - x_{bary,7}; Tracks", 200, -5, 5);
  h1["track_diff06_btpc"] = new TH1D("track_diff06_btpc", ";x_{bary,0} - x_{bary,6}; Tracks", 200, -5, 5);
  h1["track_diff17_btpc"] = new TH1D("track_diff17_btpc", ";x_{bary,1} - x_{bary,7}; Tracks", 200, -5, 5);
  h1["track_diff01_utpc"] = new TH1D("track_diff01_utpc", ";x_{utpc,0} - x_{utpc,1}; Tracks", 200, -5, 5);
  h1["track_diff67_utpc"] = new TH1D("track_diff67_utpc", ";x_{utpc,6} - x_{utpc,7}; Tracks", 200, -5, 5);
  h1["track_diff06_utpc"] = new TH1D("track_diff06_utpc", ";x_{utpc,0} - x_{utpc,6}; Tracks", 200, -5, 5);
  h1["track_diff17_utpc"] = new TH1D("track_diff17_utpc", ";x_{utpc,1} - x_{utpc,7}; Tracks", 200, -5, 5);
  h1["track_diff01_comb"] = new TH1D("track_diff01_comb", ";x_{comb,0} - x_{comb,1}; Tracks", 200, -5, 5);
  h1["track_diff67_comb"] = new TH1D("track_diff67_comb", ";x_{comb,6} - x_{comb,7}; Tracks", 200, -5, 5);

  h1["track_diff01_btpc_norm"] = new TH1D("track_diff01_btpc_norm", ";x_{bary,0} - x_{bary,1}; Tracks", 200, -5, 5);
  h1["track_diff67_btpc_norm"] = new TH1D("track_diff67_btpc_norm", ";x_{bary,6} - x_{bary,7}; Tracks", 200, -5, 5);
  h1["track_diff06_btpc_norm"] = new TH1D("track_diff06_btpc_norm", ";x_{bary,0} - x_{bary,6}; Tracks", 200, -5, 5);
  h1["track_diff17_btpc_norm"] = new TH1D("track_diff17_btpc_norm", ";x_{bary,1} - x_{bary,7}; Tracks", 200, -5, 5);
  h1["track_diff01_utpc_norm"] = new TH1D("track_diff01_utpc_norm", ";x_{utpc,0} - x_{utpc,1}; Tracks", 200, -5, 5);
  h1["track_diff67_utpc_norm"] = new TH1D("track_diff67_utpc_norm", ";x_{utpc,6} - x_{utpc,7}; Tracks", 200, -5, 5);
  h1["track_diff06_utpc_norm"] = new TH1D("track_diff06_utpc_norm", ";x_{utpc,0} - x_{utpc,6}; Tracks", 200, -5, 5);
  h1["track_diff17_utpc_norm"] = new TH1D("track_diff17_utpc_norm", ";x_{utpc,1} - x_{utpc,7}; Tracks", 200, -5, 5);
  h1["track_diff01_comb_norm"] = new TH1D("track_diff01_comb_norm", ";x_{comb,0} - x_{comb,1}; Tracks", 200, -5, 5);
  h1["track_diff67_comb_norm"] = new TH1D("track_diff67_comb_norm", ";x_{comb,6} - x_{comb,7}; Tracks", 200, -5, 5);

  h2["track_diff01_bary_theta"] = new TH2D("track_diff01_bary_theta", ";theta;x_{bary,0} - x_{bary,1}; Tracks", 100, -30, 30, 2000, -50, 50);
  h2["track_diff67_bary_theta"] = new TH2D("track_diff67_bary_theta", ";theta;x_{bary,6} - x_{bary,7}; Tracks", 100, -30, 30, 2000, -50, 50);
  h2["track_diff01_btpc_theta"] = new TH2D("track_diff01_btpc_theta", ";theta;x_{bary,0} - x_{bary,1}; Tracks", 100, -30, 30, 2000, -50, 50);
  h2["track_diff67_btpc_theta"] = new TH2D("track_diff67_btpc_theta", ";theta;x_{bary,6} - x_{bary,7}; Tracks", 100, -30, 30, 2000, -50, 50);
  h2["track_diff01_utpc_theta"] = new TH2D("track_diff01_utpc_theta", ";theta;x_{utpc,0} - x_{utpc,1}; Tracks", 100, -30, 30, 2000, -50, 50);
  h2["track_diff67_utpc_theta"] = new TH2D("track_diff67_utpc_theta", ";theta;x_{utpc,6} - x_{utpc,7}; Tracks", 100, -30, 30, 2000, -50, 50);
  h2["track_diff01_comb_theta"] = new TH2D("track_diff01_comb_theta", ";theta;x_{comb,0} - x_{comb,1}; Tracks", 100, -30, 30, 2000, -50, 50);
  h2["track_diff67_comb_theta"] = new TH2D("track_diff67_comb_theta", ";theta;x_{comb,6} - x_{comb,7}; Tracks", 100, -30, 30, 2000, -50, 50);

  h2["track_diff01_bary_abstheta"] = new TH2D("track_diff01_bary_abstheta", ";theta;x_{bary,0} - x_{bary,1}; Tracks", 100, -2, 30, 2000, -50, 50);
  h2["track_diff67_bary_abstheta"] = new TH2D("track_diff67_bary_abstheta", ";theta;x_{bary,6} - x_{bary,7}; Tracks", 100, -2, 30, 2000, -50, 50);
  h2["track_diff01_btpc_abstheta"] = new TH2D("track_diff01_btpc_abstheta", ";theta;x_{bary,0} - x_{bary,1}; Tracks", 100, -2, 30, 2000, -50, 50);
  h2["track_diff67_btpc_abstheta"] = new TH2D("track_diff67_btpc_abstheta", ";theta;x_{bary,6} - x_{bary,7}; Tracks", 100, -2, 30, 2000, -50, 50);
  h2["track_diff01_utpc_abstheta"] = new TH2D("track_diff01_utpc_abstheta", ";theta;x_{utpc,0} - x_{utpc,1}; Tracks", 100, -2, 30, 2000, -50, 50);
  h2["track_diff67_utpc_abstheta"] = new TH2D("track_diff67_utpc_abstheta", ";theta;x_{utpc,6} - x_{utpc,7}; Tracks", 100, -2, 30, 2000, -50, 50);
  h2["track_diff01_comb_abstheta"] = new TH2D("track_diff01_comb_abstheta", ";theta;x_{comb,0} - x_{comb,1}; Tracks", 100, -2, 30, 2000, -50, 50);
  h2["track_diff67_comb_abstheta"] = new TH2D("track_diff67_comb_abstheta", ";theta;x_{comb,6} - x_{comb,7}; Tracks", 100, -2, 30, 2000, -50, 50);

  h1["tpc_uncertainty"]   = new TH1D("tpc_uncertainty",   ";TPC unceratinty;TPC Clusters", 1000, -10, 40);
  h1["tpc_uncertainty_0"] = new TH1D("tpc_uncertainty_0", ";TPC unceratinty;TPC Clusters", 1000, -10, 40);
  h1["tpc_uncertainty_1"] = new TH1D("tpc_uncertainty_1", ";TPC unceratinty;TPC Clusters", 1000, -10, 40);
  h1["tpc_uncertainty_6"] = new TH1D("tpc_uncertainty_6", ";TPC unceratinty;TPC Clusters", 1000, -10, 40);
  h1["tpc_uncertainty_7"] = new TH1D("tpc_uncertainty_7", ";TPC unceratinty;TPC Clusters", 1000, -10, 40);
  h2["tpc_unc01"] = new TH2D("tpc_unc01", ";TPC uncertainty 0; TPC uncertainty 1;TPC Clusters", 1000, -1, 10, 1000, -1, 10);
  h2["tpc_unc67"] = new TH2D("tpc_unc67", ";TPC uncertainty 6; TPC uncertainty 7;TPC Clusters", 1000, -1, 10, 1000, -1, 10);
  h2["tpc_val01"] = new TH2D("tpc_val01", ";TPC val 0; TPC val 1;TPC Clusters", 100, -1, 10, 100, -1, 7);
  h2["tpc_val67"] = new TH2D("tpc_val67", ";TPC val 6; TPC val 7;TPC Clusters", 100, -1, 10, 100, -1, 7);

  TString name;
  for (ibo = 0; ibo < nboards; ibo++){
    name = Form("track_N1_theta_x_vs_residual_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x theta;x_{cluster} - x_{track, proj.};Tracks", 100, -35, 35, 200, -5.0, 5.0);

    name = Form("track_N1_theta_x_vs_residual_norm_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x_{cluster};x_{cluster} - x_{track, proj.};Tracks", 100, -20, 220, 200, -5.0, 5.0);

    name = Form("track_N1_x_vs_residual_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x_{cluster};x_{cluster} - x_{track, proj.};Tracks", 100, -20, 220, 200, -5.0, 5.0);

    name = Form("track_N1_x_vs_residual_10deg_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x_{cluster};x_{cluster} - x_{track, proj.};Tracks", 100, -20, 220, 200, -5.0, 5.0);

    name = Form("track_N1_theta_x_vs_utpc_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x theta;x_{utpc} - x_{track, proj.};Tracks", 100, -35, 35, 200, -5.0, 5.0);

    name = Form("track_N1_theta_x_vs_comp_%i", ibo);
    h1[name.Data()] = new TH1D(name+"1", "        ;x_{utpc} - x_{bary};Tracks",               200, -5.0, 5.0);
    h2[name.Data()] = new TH2D(name,     ";x theta;x_{utpc} - x_{bary};Tracks", 100, -35, 35, 200, -5.0, 5.0);

    name = Form("track_N1_theta_x_vs_comb_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x theta;x_{comb} - x_{track, proj.};Tracks", 100, -35, 35, 200, -5.0, 5.0);
  }

  h2["strip_position_vs_board"] = new TH2D("strip_position_vs_board", ";strip number;MMFE number;charge [fC]", 512, 0.5, 512.5, 8, -0.5, 7.5);

  for (ibo = 0; ibo < nboards; ibo++){
    h2[Form("strip_q_vs_ch_%i",    ibo)] = new TH2D(Form("strip_q_vs_ch_%i",    ibo), ";strip number;Charge [fC];strip",    512, 0.5, 512.5, 512,   0,  128);
    h2[Form("strip_pdo_vs_ch_%i",  ibo)] = new TH2D(Form("strip_pdo_vs_ch_%i",  ibo), ";strip number;PDO [counts];strip",   512, 0.5, 512.5, 512,   0, 2048);
    h2[Form("strip_tdo_vs_ch_%i",  ibo)] = new TH2D(Form("strip_tdo_vs_ch_%i",  ibo), ";strip number;TDO [counts];strip",   512, 0.5, 512.5, 256,   0,  256);
    h2[Form("strip_tdoc_vs_ch_%i", ibo)] = new TH2D(Form("strip_tdoc_vs_ch_%i", ibo), ";strip number;TDO corr. [ns];strip", 512, 0.5, 512.5, 110, -10,  100);
    h2[Form("strip_dbc_vs_ch_%i",  ibo)] = new TH2D(Form("strip_dbc_vs_ch_%i",  ibo), ";strip number;#Delta BC;strip",      512, 0.5, 512.5,  64,   0,   64);
    h2[Form("strip_time_vs_ch_%i", ibo)] = new TH2D(Form("strip_time_vs_ch_%i", ibo), ";strip number;Time [ns];strip",      512, 0.5, 512.5, 100, 300, 1000);
    h2[Form("strip_zpos_vs_ch_%i", ibo)] = new TH2D(Form("strip_zpos_vs_ch_%i", ibo), ";strip number;z_{drift} [mm];strip", 512, 0.5, 512.5, 150, -10,   20);
    h2[Form("strip_zdri_vs_ch_%i", ibo)] = new TH2D(Form("strip_zdri_vs_ch_%i", ibo), ";strip number;z_{drift} [mm];strip", 512, 0.5, 512.5, 150, -10,   20);
    h2[Form("strip_zres_vs_ch_%i", ibo)] = new TH2D(Form("strip_zres_vs_ch_%i", ibo), ";strip number;#Delta z [mm];strip",  512, 0.5, 512.5, 200, -15,   15);
  }
  h2["dups_vs_ch"] = new TH2D("dups_vs_ch", ";strip number;MMFE number;Duplicates", 512, 0.5, 512.5, 8, -0.5, 7.5);

  for (ibo = 0; ibo < nboards; ibo++){
    h2[Form("strip_zpos_vs_ztrack_%i", ibo)] = new TH2D(Form("strip_zpos_vs_ztrack_%i", ibo), ";z_{track};z_{drift} [mm];strip", 100, -5, 10, 100, -2, 8);
    h2[Form("strip_chi2_vs_tpcunc_%i", ibo)] = new TH2D(Form("strip_chi2_vs_tpcunc_%i", ibo), ";chi2/ndf;#sigma(uTPC x);strip",  100,  0, 1.5, 100,  0, 2);
    h2[Form("tpc_phi_vs_phi_%i",       ibo)] = new TH2D(Form("tpc_phi_vs_phi_%i",       ibo), ";#phi [deg];#phi [deg];TPC cluster", 200,  -2, 2, 200, -2, 2);
  }

  h2["clus_vs_board"]          = new TH2D("clus_vs_board",          ";MMFE number;clusters;Events",            8, -0.5, 7.5, 32, -0.5, 31.5);
  h2["hits_vs_board"]          = new TH2D("hits_vs_board",          ";MMFE number;strips;Events",              8, -0.5, 7.5, 32, -0.5, 31.5);
  h2["dups_vs_board"]          = new TH2D("dups_vs_board",          ";MMFE number;duplicate strips;Events",    8, -0.5, 7.5, 32, -0.5, 31.5);
  h2["hits_per_clus_vs_board"] = new TH2D("hits_per_clus_vs_board", ";MMFE number;hits in a cluster;Clusters", 8, -0.5, 7.5, 13, -0.5, 12.5);
  h2["timediff_vs_board"]      = new TH2D("timediff_vs_board",      ";MMFE number;#DeltaBCID;Events",          8, -0.5, 7.5, 100, -0.5, 99.5);
  h1["timediff_track"]         = new TH1D("timediff_track",         ";#DeltaBCID;Strips on-track",             100, -0.5, 99.5);
  h1["timediff_dupli"]         = new TH1D("timediff_dupli",         ";#DeltaBCID;Duplicate strips",            100, -0.5, 99.5);
  h2["timediff_vs_charge"]     = new TH2D("timediff_vs_charge",     ";charge [fC];#DeltaBCID;Events",          200, 0, 200, 100, -0.5, 99.5);

  h2["track_exp_hit"] = new TH2D("track_exp_hit", ";VMM;Board;Trigger", 10, -1.5, 8.5, 17, -0.75, 7.75);
  h2["track_obs_hit"] = new TH2D("track_obs_hit", ";VMM;Board;Trigger", 10, -1.5, 8.5, 17, -0.75, 7.75);

  // cataloging duplicates
  MMFE8Hits duplicates;

  // collecting clusters and the nominal fit
  std::vector<MMClusterList> clusters_perboard;
  MMClusterList clusters_all;
  MMClusterList clusters_road;
  MMClusterList clusters_x;
  MMTrack track;

  // N-1 fits and resolutions
  MMClusterList clusters_N1;
  MMTrack track_N1;

  // comparing two 4-hit tracks: different quads
  MMTrack track_0123;
  MMTrack track_4567;
  MMClusterList clus_0123;
  MMClusterList clus_4567;

  // comparing two 4-hit tracks: diff-quad Xs, diff-quad UV
  MMTrack track_0256;
  MMTrack track_1347;
  MMClusterList clus_0256;
  MMClusterList clus_1347;

  // comparing two 4-hit tracks: diff-quad Xs, same-quad UV
  MMTrack track_0236;
  MMTrack track_1457;
  MMClusterList clus_0236;
  MMClusterList clus_1457;

  // miscellaneous helpers
  double start_time = -1;
  double time_since_start = -1;
  double days_since_start = -1;
  double seconds_per_day  = 60*60*24;
  int EventNum4Hist = 0;
  double z_track = 0.0;
  double z_utpc  = 0.0;
  double z_utpc_check = 0.0;
  double offset = 0.0;
  double x_utpc  = 0.0;
  double x_clus  = 0.0;
  double x_comb  = 0.0;
  double sign   = 0.0;
  double vdrift = 1.0 / 20; // mm per ns
  double deltaT = 41.3 / vdrift;
  std::vector<double> xs;
  std::vector<double> zs_dbc;
  std::vector<double> zs;
  std::vector<double> xs_sus;
  std::vector<double> zs_sus;
  std::vector<double> zs_sus_dbc;
  std::vector<double> neighbor_xs;
  std::vector<double> neighbor_zs;
  double residual = 0.0;
  double quasires = 0.0;
  TCanvas* can;
  GeoPlane plane;
  double z_middle = 0.0;
  for (ibo = 0; ibo < nboards; ibo++)
    z_middle += GEOMETRY->Get(ibo).Origin().Z();
  z_middle /= (double)(nboards);

  // board positions for uTPC and trigger
  // mu-telescope.rc.fas.harvard.edu:/data/mm_2016/mm_ana/oct_ana.C
//   const std::vector<double> zboard = {0.4,
//                                       10.8,
//                                       32.4,
//                                       43.6,
//                                       113.6,
//                                       124.8,
//                                       146.7,
//                                       156.5};

//   const std::vector<double> zboard = {0.0,
//                                       11.2,
//                                       32.4,
//                                       43.6,
//                                       113.6,
//                                       124.8,
//                                       146.0,
//                                       157.2};
  //TMultiGraph* utpc_mg      = 0;
  //TGraph*      utpc_graph   = 0;
  //TGraph*      utpc_dbc     = 0;
  //TGraph*      utpc_sus     = 0;
  //TGraph*      utpc_sus_dbc = 0;
  //TGraph*      utpc_bary    = 0;
  //TGraph*      utpc_pred    = 0;
  //TGraph*      utpc_others  = 0;
  //TF1*         utpc_fit     = 0;
  //TLine*       utpc_ref     = 0;
  //TLine*       utpc_horiz1  = 0;
  //TLine*       utpc_horiz2  = 0;

  // output
  TFile* fout = new TFile(outputFileName, "RECREATE");
  fout->mkdir("event_displays");
  MMPlot();

  // progress
  std::chrono::time_point<std::chrono::system_clock> time_start;
  std::chrono::duration<double> elapsed_seconds;
  time_start = std::chrono::system_clock::now();

  // loop
  for(int evt = 0; evt < Nevent; evt++){

    // get the event
    DATA->GetEntry(evt);
    if(evt % 1000 == 0) {    
      elapsed_seconds = (std::chrono::system_clock::now() - time_start);
      progress(elapsed_seconds.count(), evt, Nevent);
    }
    if(GEOMETRY->RunNumber() < 0)
      GEOMETRY->SetRunNumber(DATA->RunNum);
    if (start_time < 0)
      start_time = DATA->mm_EventHits.Time();
    time_since_start = DATA->mm_EventHits.Time() - start_time;
    days_since_start = time_since_start / seconds_per_day;

    // scintillator diagnostics
    for (int i = 0; i < DATA->sc_EventHits.GetNHits(); i++){
      auto hit = DATA->sc_EventHits.Get(i);
      h2["scint_ch_vs_counts"]->Fill(hit.Channel(), hit.Count());
    }
    if(!DATA->sc_EventHits.IsGoodEvent())
      continue;
    if(DATA->RunNum == 3545 && DATA->sci_EventNum > 150000 && DATA->sci_EventNum < 154600) continue;
    if(DATA->RunNum == 3546 && DATA->sci_EventNum >  11400 && DATA->sci_EventNum <  12100) continue;
    if(DATA->RunNum == 3547 && DATA->sci_EventNum > 123000 && DATA->sci_EventNum < 124700) continue;
    if(DATA->RunNum == 3548 && DATA->sci_EventNum > 127300 && DATA->sci_EventNum < 129000) continue;
    if(DATA->RunNum == 3549 && DATA->sci_EventNum >  72000 && DATA->sci_EventNum <  75100) continue;

    // for combining adjacent runs
    EventNum4Hist = DATA->mm_EventNum;
    if (DATA->RunNum == 3519 && Nevent > 400000)
      EventNum4Hist = EventNum4Hist + 200000;
    else if (DATA->RunNum == 3520 && Nevent > 400000)
      EventNum4Hist = EventNum4Hist + 500000;

    // reset clusters and tracks
    track.Reset();
    track_N1.Reset();
    clusters_all.Reset();
    clusters_road.Reset();
    clusters_N1.Reset();
    for (auto clus_list: clusters_perboard)
      clus_list.Reset();
    clusters_perboard.clear();
    
    //if (evt > 50000)
    //  break;

    // calibrate
    PDOCalibrator->Calibrate(DATA->mm_EventHits);
    TDOCalibrator->Calibrate(DATA->mm_EventHits);
    PACMAN->SetEventPadTime(0);
    FILTERER->SetRunNumber(DATA->RunNum);
    
    // run pacman
    nboardshit = DATA->mm_EventHits.GetNBoards();
    for(i = 0; i < nboardshit; i++){
      if(DATA->mm_EventHits[i].GetNHits() == 0)
        continue;
      MMClusterList board_clusters = PACMAN->Cluster(DATA->mm_EventHits[i]);
      if(board_clusters.GetNCluster() > 0)
        clusters_perboard.push_back(board_clusters);

      // strips: PDO vs channel
      // its slow. only run if hella desired.
      for(ich = 0; ich < DATA->mm_EventHits[i].GetNHits(); ich++){
        auto hit = DATA->mm_EventHits[i][ich];
        ibo = GEOMETRY->Index(DATA->mm_EventHits[i].MMFE8());

        h2["timediff_vs_board"]->Fill(ibo, hit.DeltaBC());
        h2["timediff_vs_charge"]->Fill(hit.Charge(), hit.DeltaBC());

        if( !PACMAN->IsGoodHit(hit) )
          continue;

        h2[Form("strip_tdoc_vs_ch_%i", ibo)]->Fill(hit.Channel(), hit.Time());
        h2[Form("strip_dbc_vs_ch_%i",  ibo)]->Fill(hit.Channel(), hit.DeltaBC());
        h2[Form("strip_time_vs_ch_%i", ibo)]->Fill(hit.Channel(), hit.DriftTime(deltaT));
        h2[Form("strip_zpos_vs_ch_%i", ibo)]->Fill(hit.Channel(), hit.DriftTime(deltaT) * vdrift);

        h2[Form("strip_q_vs_ch_%i",   ibo)]->Fill(hit.Channel(), hit.Charge());
        h2[Form("strip_pdo_vs_ch_%i", ibo)]->Fill(hit.Channel(), hit.PDO());
        h2[Form("strip_tdo_vs_ch_%i", ibo)]->Fill(hit.Channel(), hit.TDO());
      }
    }

    // hits, duplicates, clusters per board
    // ------------------------------------
    for (ibo = 0; ibo < nboards; ibo++){
      test = -1;
      for (i = 0; i < nboardshit; i++){
        // hits on this board!
        if (GEOMETRY->Index(DATA->mm_EventHits[i].MMFE8()) == ibo){
          test = i;
          h2["clus_vs_board"]->Fill(ibo, clusters_perboard[i].size());
          h2["hits_vs_board"]->Fill(ibo, DATA->mm_EventHits[i].GetNHits());
          h2["dups_vs_board"]->Fill(ibo, DATA->mm_EventHits[i].GetNDuplicates());
          duplicates = DATA->mm_EventHits[i].GetDuplicates();
          for(ich = 0; ich < (int)(duplicates.GetNHits()); ich++)
            h2["dups_vs_ch"]->Fill(duplicates[ich].Channel(), ibo);
        }
      }
      // no hits on this board!
      if (test == -1){
        h2["clus_vs_board"]->Fill(ibo, 0);
        h2["hits_vs_board"]->Fill(ibo, 0);
        h2["dups_vs_board"]->Fill(ibo, 0);
      }
    }
    
    for (auto clus_list: clusters_perboard)
      for (auto clus: clus_list)
        h2["hits_per_clus_vs_board"]->Fill(GEOMETRY->Index(clus->MMFE8()), clus->GetNHits());
    
    // preselection quality to run tracking
    // require at least N boards hit
    if (clusters_perboard.size() < 5)
      continue;

    // flatten these lists
    for (auto clus_list: clusters_perboard)
      for (auto clus: clus_list)
        clusters_all.AddCluster(*clus);
    
    if (debug || DATA->mm_EventNum == 14036 || DATA->mm_EventNum == 46873 || DATA->mm_EventNum == 41380){
      track.Reset();
      can = Plot_Track2D(Form("track2D_%05d_all", DATA->mm_EventNum), track, *GEOMETRY, &clusters_all);
      fout->cd("event_displays");
      can->Write();
      delete can;
    }

    // fit it!
    for (auto botpair: DATA->sc_EventHits.GetBotPair()){
      clusters_road.Reset();
      clusters_road = FILTERER->FilterClustersScint(clusters_all, *GEOMETRY, botpair.first->Channel(), DATA->mm_EventNum, debug);
      track  = FITTER ->Fit(clusters_road, *GEOMETRY, DATA->mm_EventNum);
      for (auto clus: clusters_road)
        track.CountHit(GEOMETRY->Index(clus->MMFE8()));
    }

    // angles
    if (clusters_road.size() == 6) h2["track_angle_6"]->Fill(theta(track.SlopeX()), theta(track.SlopeY()));
    if (clusters_road.size() == 7) h2["track_angle_7"]->Fill(theta(track.SlopeX()), theta(track.SlopeY()));
    if (clusters_road.size() == 8) h2["track_angle_8"]->Fill(theta(track.SlopeX()), theta(track.SlopeY()));

    // clusters per track
    h2["track_hits_vs_time"]->Fill(days_since_start, clusters_road.size());
    if (GEOMETRY->IsFiducial(track))
      h2["track_hits_vs_time_fid"]->Fill(days_since_start, clusters_road.size());

    // require a good track for further analysis
    if (clusters_road.size() < 6)
      continue;

    // wtf
    if (clusters_road.size() == 8){
      double y = track.SlopeY()*z_middle + track.ConstY();
      h2["track_res2_vs_y"]  ->Fill(y,              sqrt(track.ResidualSq()));
      h2["track_res2_vs_my"] ->Fill(track.SlopeY(), sqrt(track.ResidualSq()));
      h2["track_res2_vs_cy"] ->Fill(track.ConstY(), sqrt(track.ResidualSq()));
    }

    // some sanity plots
    h2["track_param_x"]->Fill(track.SlopeX(), track.ConstX());
    h2["track_param_y"]->Fill(track.SlopeY(), track.ConstY());
    for (auto clus: clusters_road)
      h2["track_clusmult_vs_theta"]->Fill(atan(track.SlopeX())*180/3.14159, clus->size());

    // interlude: some scintillator plots
    for (auto toppair: DATA->sc_EventHits.GetTopPair()){
      for (auto botpair: DATA->sc_EventHits.GetBotPair()){
        h2["scint_bot_vs_top"]->Fill(botpair.first->Channel(), toppair.first->Channel());
        if (toppair.first->Channel() == 19) h2["track_angle_x_scint_19"]->Fill(theta(track.SlopeX()), botpair.first->Channel());
        if (toppair.first->Channel() == 20) h2["track_angle_x_scint_20"]->Fill(theta(track.SlopeX()), botpair.first->Channel());
        if (toppair.first->Channel() == 21) h2["track_angle_x_scint_21"]->Fill(theta(track.SlopeX()), botpair.first->Channel());

        if (botpair.first->Channel() == 0) h2["track_scint_0_vs_time"]->Fill(days_since_start, theta(track.SlopeX()));
        if (botpair.first->Channel() == 1) h2["track_scint_1_vs_time"]->Fill(days_since_start, theta(track.SlopeX()));
        if (botpair.first->Channel() == 2) h2["track_scint_2_vs_time"]->Fill(days_since_start, theta(track.SlopeX()));
        if (botpair.first->Channel() == 3) h2["track_scint_3_vs_time"]->Fill(days_since_start, theta(track.SlopeX()));
        if (botpair.first->Channel() == 4) h2["track_scint_4_vs_time"]->Fill(days_since_start, theta(track.SlopeX()));
        if (botpair.first->Channel() == 5) h2["track_scint_5_vs_time"]->Fill(days_since_start, theta(track.SlopeX()));
      }
    }

    // VMM-level efficiency
    // ------------------------------------------------------------------
    double xclus = 0.0;
    double xproj = 0.0;
    int vmm_min  = 0, vmm_max  = 0;
    int chan_min = 0, chan_max = 0;
    for (ibo = 0; ibo < 8; ibo++){
      
      // make track from other boards
      clusters_N1.Reset();
      track_N1.Reset();
      for (auto clus: clusters_road)
        if (clus->MMFE8Index() != ibo)
          clusters_N1.AddCluster(*clus);
      track_N1 = FITTER->Fit(clusters_N1, *GEOMETRY, DATA->mm_EventNum);
      if (!track_N1.IsFit())
        continue;
      if (track_N1.NX() + track_N1.NU() + track_N1.NV() < 6)
        continue;

      // find track position at board and corresponding channel
      // find it at both ends of the board, to be safe for UV
      xproj    = track_N1.SlopeX()*GEOMETRY->Get(ibo).Origin().Z() + track_N1.ConstX();
      chan_min = (int)(channel_from_x(xproj, ibo, GEOMETRY, true));
      chan_max = (int)(channel_from_x(xproj, ibo, GEOMETRY, false));
    
      if (chan_min <= 0 || chan_max <= 0)
        continue;

      vmm_min = (chan_min-1)/64;
      vmm_max = (chan_max-1)/64;
      if (vmm_min != vmm_max)
        continue;

      chan_min = (chan_min-1)%64 - 1;
      chan_max = (chan_max-1)%64 - 1;
    
      // avoid edge cases
      if (chan_min==1 || chan_min==2 || chan_min==63 || chan_min==64)
        continue;
      if (chan_max==1 || chan_max==2 || chan_max==63 || chan_max==64)
        continue;
      
      h2["track_exp_hit"]->Fill(vmm_min, ibo);
      for (auto clus: clusters_road)
        if (clus->MMFE8Index() == ibo)
          h2["track_obs_hit"]->Fill(vmm_min, ibo);

    }

    // N-1 histograms (residuals)
    // ------------------------------------------------------------------
    for (auto clus: clusters_road){
      
      clusters_N1.Reset();
      track_N1.Reset();

      for (auto clus_other: clusters_road){
        if (clus == clus_other)
          continue;
        clusters_N1.AddCluster(*clus_other);
      }
      
      track_N1 = FITTER->Fit(clusters_N1, *GEOMETRY, DATA->mm_EventNum);
      residual = GEOMETRY->GetResidualX(*clus, track_N1);
      quasires = GEOMETRY->GetResidualX(*clus, track);
      ibo      = GEOMETRY->Index(clus->MMFE8());
      plane    = GEOMETRY->Get(GEOMETRY->Index(clus->MMFE8()));
      xclus    = plane.Origin().X() + plane.LocalXatYbegin(clus->Channel());

      h2["track_N1_board_vs_residual"]                ->Fill(ibo,                   residual);
      h2["track_N1_board_vs_residual_norm"]           ->Fill(ibo,                   quasires / clus->ChannelUnc(track.SlopeX()));
      h2[Form("track_N1_theta_x_vs_residual_%i", ibo)]->Fill(theta(track.SlopeX()), residual);
      h2[Form("track_N1_x_vs_residual_%i",       ibo)]->Fill(xclus,                 residual);
      if (theta(track.SlopeX()) < -10)
        h2[Form("track_N1_x_vs_residual_10deg_%i", ibo)]->Fill(xclus, residual);
    }

    // barycenter board_i vs board_j
    int hit_0 = 0, hit_1 = 0, hit_6 = 0, hit_7 = 0;
    for (auto clus: clusters_road){
      if      (clus->MMFE8Index() == 0) hit_0 = 1;
      else if (clus->MMFE8Index() == 1) hit_1 = 1;
      else if (clus->MMFE8Index() == 6) hit_6 = 1;
      else if (clus->MMFE8Index() == 7) hit_7 = 1;
    }

    GeoPlane plane_i;
    GeoPlane plane_j;
    double x_i = 0;
    double x_j = 0;
    double dx = 0;

    if (hit_0 && hit_1){
      plane_i = GEOMETRY->Get(0);
      plane_j = GEOMETRY->Get(1);
      for (auto clus: clusters_road){
        if (clus->MMFE8Index() == 0) x_i = plane_i.Origin().X() + plane_i.LocalXatYbegin(clus->Channel());
        if (clus->MMFE8Index() == 1) x_j = plane_j.Origin().X() + plane_j.LocalXatYbegin(clus->Channel());
      }
      dx = x_i - x_j + AngleCorrection(0, 1, track.SlopeX(), GEOMETRY);
      h1["track_diff01_bary"]      ->Fill(dx);
      h2["track_diff01_bary_theta"]->Fill(theta(track.SlopeX()), dx);
      h2["track_diff01_bary_abstheta"]->Fill(fabs(theta(track.SlopeX())), dx);
    }
    if (hit_6 && hit_7){
      plane_i = GEOMETRY->Get(6);
      plane_j = GEOMETRY->Get(7);
      for (auto clus: clusters_road){
        if (clus->MMFE8Index() == 6) x_i = plane_i.Origin().X() + plane_i.LocalXatYbegin(clus->Channel());
        if (clus->MMFE8Index() == 7) x_j = plane_j.Origin().X() + plane_j.LocalXatYbegin(clus->Channel());
      }
      dx = x_i - x_j + AngleCorrection(6, 7, track.SlopeX(), GEOMETRY);
      h1["track_diff67_bary"]      ->Fill(dx);
      h2["track_diff67_bary_theta"]->Fill(theta(track.SlopeX()), dx);
      h2["track_diff67_bary_abstheta"]->Fill(fabs(theta(track.SlopeX())), dx);
    }

    // uTPC, finally
    // ------------------------------------------------------------------
    track = FITTER->Fit(clusters_road, *GEOMETRY, DATA->mm_EventNum);
    if (std::fabs(theta(track.SlopeX())) < 10)
      continue;
    if (clusters_road.size() < 6)
      continue;

    int ntpc = 0;
    double z_half = 0;
    double fiducial_x    =  1.8;
    double fiducial_z_hi =  6.6;
    double fiducial_z_lo = -0.8;

    int utpc_0 = 0, utpc_1 = 0, utpc_6 = 0, utpc_7 = 0;
    double x_utpc_0 = -999, x_utpc_1 = -999, x_utpc_6 = -999, x_utpc_7 = -999;
    double x_bary_0 = -999, x_bary_1 = -999, x_bary_6 = -999, x_bary_7 = -999;
    double x_comb_0 = -999, x_comb_1 = -999, x_comb_6 = -999, x_comb_7 = -999;
    double u_utpc_0 = -999, u_utpc_1 = -999, u_utpc_6 = -999, u_utpc_7 = -999;
    double u_bary_0 = -999, u_bary_1 = -999, u_bary_6 = -999, u_bary_7 = -999;
    double u_comb_0 = -999, u_comb_1 = -999, u_comb_6 = -999, u_comb_7 = -999;

    for (auto clus: clusters_road){
      
      // N-1 track
      clusters_N1.Reset();
      track_N1.Reset();
      for (auto clus_other: clusters_road){
        if (clus == clus_other)
          continue;
        clusters_N1.AddCluster(*clus_other);
      }
      track_N1 = FITTER->Fit(clusters_N1, *GEOMETRY, DATA->mm_EventNum);

      // sign of the drift
      ibo = GEOMETRY->Index(clus->MMFE8());
      sign = (ibo==0 || ibo==2 || ibo==4 || ibo==6) ? -1.0 : 1.0;

      plane  = GEOMETRY->Get(ibo);
      z_half = plane.Origin().Z();

      // utpc points
      ntpc = 0;
      xs.clear();
      zs_dbc.clear();
      zs.clear();
      xs_sus.clear();
      zs_sus.clear();
      zs_sus_dbc.clear();
      for (auto hit: *clus){

        x_clus  = plane.Origin().X() + plane.LocalXatYbegin(clus->Channel());
        x_utpc  = plane.Origin().X() + plane.LocalXatYbegin(hit->Channel());
        z_track = (x_utpc - track.ConstX()) / track.SlopeX();
        offset  = OffsetByBoard(ibo) + channel_offsets[ibo][hit->VMM()][(int)(hit->VMMChannel())];
        z_utpc  = zboard[ibo] + sign*(vdrift * hit->DriftTime(deltaT) + offset);

        // rose-colored lens for suspicious BCIDs: keep the best one
        if (hit->SuspiciousBCID()){
          hit->SetBCID(hit->BCID() - 1);
          z_utpc_check = zboard[ibo] + sign*(vdrift * hit->DriftTime(deltaT) + offset);
          if (fabs(z_track - z_utpc_check) > fabs(z_track - z_utpc))
            hit->SetBCID(hit->BCID() + 1);
        }
        z_utpc = zboard[ibo] + sign*(vdrift * hit->DriftTime(deltaT) + offset);

        if (fabs(x_clus - x_utpc) > fiducial_x)
          continue;
        if (vdrift*hit->DriftTime(deltaT) + offset > fiducial_z_hi)
          continue;
        if (vdrift*hit->DriftTime(deltaT) + offset < fiducial_z_lo)
          continue;

        // h2[Form("strip_zdri_vs_ch_%i",     ibo)]->Fill(hit->Channel(), sign*(z_utpc - zboard[ibo]));
        h2[Form("strip_zdri_vs_ch_%i",     ibo)]->Fill(hit->Channel(), vdrift*hit->DriftTime(deltaT) + offset);
        h2[Form("strip_zres_vs_ch_%i",     ibo)]->Fill(hit->Channel(), z_utpc - z_track);
        h2[Form("strip_zpos_vs_ztrack_%i", ibo)]->Fill(sign*(z_track - zboard[ibo]), vdrift*hit->DriftTime(deltaT) + offset);
        xs.push_back(x_utpc);
        zs.push_back(z_utpc);
        ntpc++;
      }

      double slope, offset, x_fit, x_track, x_unc;
      double chi2, ndf, prob, cov00, cov01, cov10, cov11;

      // success! 
      if (ntpc >= 3){

        // local fit
        int invert = 1;
        std::tie(slope, offset, chi2, ndf, prob, cov00, cov01, cov10, cov11) = fit_root(xs, zs, invert);
        x_fit    = (z_half - offset) / slope;
        x_track  = track.SlopeX()*z_half + track.ConstX();
        residual = x_fit - x_track;
        h2["track_N1_board_vs_prob"]->Fill(ibo, prob);
        h2["track_N1_board_vs_utpc"]->Fill(ibo, residual);
        h2[Form("track_N1_theta_x_vs_utpc_%i", ibo)]->Fill(theta(track.SlopeX()), residual);

        // uncertainty on x_fit
        x_unc = TPCUncertainty(z_half, slope, offset, cov00, cov11, cov01, invert);
        if (!invert){
          slope  = 1/slope;
          offset = -1*offset*slope;
        }
        h1["tpc_uncertainty"]->Fill(x_unc);
        if (ibo == 0) h1["tpc_uncertainty_0"]->Fill(x_unc);
        if (ibo == 1) h1["tpc_uncertainty_1"]->Fill(x_unc);
        if (ibo == 6) h1["tpc_uncertainty_6"]->Fill(x_unc);
        if (ibo == 7) h1["tpc_uncertainty_7"]->Fill(x_unc);

        h2[Form("strip_chi2_vs_tpcunc_%i", ibo)]->Fill(chi2/ndf, x_unc);
        h2[Form("tpc_phi_vs_phi_%i",       ibo)]->Fill(-track.SlopeX(), -1/slope);

        // comparison
        double unc_bary = clus->ChannelUnc(track.SlopeX());
        plane = GEOMETRY->Get(ibo);
        x_clus = plane.Origin().X() + plane.LocalXatYbegin(clus->Channel());
        h1[Form("track_N1_theta_x_vs_comp_%i", ibo)]->Fill(                       (x_fit - x_clus) / sqrt(pow(x_unc, 2.0) + pow(unc_bary, 2.0)));
        h2[Form("track_N1_theta_x_vs_comp_%i", ibo)]->Fill(theta(track.SlopeX()), (x_fit - x_clus) / sqrt(pow(x_unc, 2.0) + pow(unc_bary, 2.0)));

        // combination
        double fraction = pow(unc_bary, 2.0) / ( pow(unc_bary, 2.0) + pow(x_unc, 2.0) );
        double unc_comb = pow(fraction, 2.0)*pow(x_unc, 2.0) + pow(1-fraction, 2.0)*pow(unc_bary, 2.0);
        unc_comb = sqrt(unc_comb);
        x_comb   = fraction*x_fit + (1-fraction)*x_clus;

        // keeping track
        if (ibo == 0) { utpc_0 = 1; x_utpc_0 = x_fit; x_bary_0 = x_clus;   x_comb_0 = x_comb;
                                    u_utpc_0 = x_unc; u_bary_0 = unc_bary; u_comb_0 = unc_comb;
        }
        if (ibo == 1) { utpc_1 = 1; x_utpc_1 = x_fit; x_bary_1 = x_clus;   x_comb_1 = x_comb;
                                    u_utpc_1 = x_unc; u_bary_1 = unc_bary; u_comb_1 = unc_comb;
        }
        if (ibo == 6) { utpc_6 = 1; x_utpc_6 = x_fit; x_bary_6 = x_clus;   x_comb_6 = x_comb;
                                    u_utpc_6 = x_unc; u_bary_6 = unc_bary; u_comb_6 = unc_comb;
        }
        if (ibo == 7) { utpc_7 = 1; x_utpc_7 = x_fit; x_bary_7 = x_clus;   x_comb_7 = x_comb;
                                    u_utpc_7 = x_unc; u_bary_7 = unc_bary; u_comb_7 = unc_comb;
        }

        //if (ibo == 1) { utpc_1 = 1; x_utpc_1 = x_fit; x_bary_1 = x_clus; x_comb_1 = x_comb; }
        //if (ibo == 6) { utpc_6 = 1; x_utpc_6 = x_fit; x_bary_6 = x_clus; x_comb_6 = x_comb; }
        //if (ibo == 7) { utpc_7 = 1; x_utpc_7 = x_fit; x_bary_7 = x_clus; x_comb_7 = x_comb; }

        // board v board: utpc
        if ((ibo == 0 && utpc_1) || (ibo == 1 && utpc_0)){
          dx = x_utpc_0 - x_utpc_1 + AngleCorrection(0, 1, track.SlopeX(), GEOMETRY);
          h1["track_diff01_utpc"]->Fill(dx);
          h1["track_diff01_utpc_norm"]->Fill(dx / sqrt(pow(u_utpc_0, 2.0) + pow(u_utpc_1, 2.0)));
          h2["track_diff01_utpc_theta"]->Fill(theta(track.SlopeX()), dx);
          h2["track_diff01_utpc_abstheta"]->Fill(fabs(theta(track.SlopeX())), dx);
          h2["tpc_unc01"]->Fill(u_utpc_0, u_utpc_1);
          h2["tpc_val01"]->Fill(x_utpc_0, x_utpc_1);
        }
        if ((ibo == 6 && utpc_7) || (ibo == 7 && utpc_6)){
          dx = x_utpc_6 - x_utpc_7 + AngleCorrection(6, 7, track.SlopeX(), GEOMETRY);
          h1["track_diff67_utpc"]->Fill(dx);
          h1["track_diff67_utpc_norm"]->Fill(dx / sqrt(pow(u_utpc_6, 2.0) + pow(u_utpc_7, 2.0)));
          h2["track_diff67_utpc_theta"]->Fill(theta(track.SlopeX()), dx);
          h2["track_diff67_utpc_abstheta"]->Fill(fabs(theta(track.SlopeX())), dx);
          h2["tpc_unc67"]->Fill(u_utpc_6, u_utpc_7);
          h2["tpc_val67"]->Fill(x_utpc_6, x_utpc_7);
        }
        if ((ibo == 0 && utpc_6) || (ibo == 6 && utpc_0)){
          dx = x_utpc_0 - x_utpc_6 + AngleCorrection(0, 6, track.SlopeX(), GEOMETRY);
          h1["track_diff06_utpc"]->Fill(dx);
          h1["track_diff06_utpc_norm"]->Fill(dx / sqrt(pow(u_utpc_0, 2.0) + pow(u_utpc_6, 2.0)));
        }
        if ((ibo == 1 && utpc_7) || (ibo == 7 && utpc_1)){
          dx = x_utpc_1 - x_utpc_7 + AngleCorrection(1, 7, track.SlopeX(), GEOMETRY);
          h1["track_diff17_utpc"]->Fill(dx);
          h1["track_diff17_utpc_norm"]->Fill(dx / sqrt(pow(u_utpc_1, 2.0) + pow(u_utpc_7, 2.0)));
        }

        // board v board: bary
        if ((ibo == 0 && utpc_1) || (ibo == 1 && utpc_0)){
          dx = x_bary_0 - x_bary_1 + AngleCorrection(0, 1, track.SlopeX(), GEOMETRY);
          h1["track_diff01_btpc"]->Fill(dx);
          h1["track_diff01_btpc_norm"]->Fill(dx / sqrt(pow(u_bary_0, 2.0) + pow(u_bary_1, 2.0)));
          h2["track_diff01_btpc_theta"]->Fill(theta(track.SlopeX()), dx);
          h2["track_diff01_btpc_abstheta"]->Fill(fabs(theta(track.SlopeX())), dx);
        }
        if ((ibo == 6 && utpc_7) || (ibo == 7 && utpc_6)){
          dx = x_bary_6 - x_bary_7 + AngleCorrection(6, 7, track.SlopeX(), GEOMETRY);
          h1["track_diff67_btpc"]->Fill(dx);
          h1["track_diff67_btpc_norm"]->Fill(dx / sqrt(pow(u_bary_6, 2.0) + pow(u_bary_7, 2.0)));
          h2["track_diff67_btpc_theta"]->Fill(theta(track.SlopeX()), dx);
          h2["track_diff67_btpc_abstheta"]->Fill(fabs(theta(track.SlopeX())), dx);
        }
        if ((ibo == 0 && utpc_6) || (ibo == 6 && utpc_0)){
          dx = x_bary_0 - x_bary_6 + AngleCorrection(0, 6, track.SlopeX(), GEOMETRY);
          h1["track_diff06_btpc"]->Fill(dx);
          h1["track_diff06_btpc_norm"]->Fill(dx / sqrt(pow(u_bary_0, 2.0) + pow(u_bary_6, 2.0)));
        }
        if ((ibo == 1 && utpc_7) || (ibo == 7 && utpc_1)){
          dx = x_bary_1 - x_bary_7 + AngleCorrection(1, 7, track.SlopeX(), GEOMETRY);
          h1["track_diff17_btpc"]->Fill(dx);
          h1["track_diff17_btpc_norm"]->Fill(dx / sqrt(pow(u_bary_1, 2.0) + pow(u_bary_7, 2.0)));
        }

        // board v board: comb
        if ((ibo == 0 && utpc_1) || (ibo == 1 && utpc_0)){
          dx = x_comb_0 - x_comb_1 + AngleCorrection(0, 1, track.SlopeX(), GEOMETRY);
          h1["track_diff01_comb"]->Fill(dx);
          h2["track_diff01_comb_theta"]->Fill(theta(track.SlopeX()), dx);
          h2["track_diff01_comb_abstheta"]->Fill(fabs(theta(track.SlopeX())), dx);
        }
        if ((ibo == 6 && utpc_7) || (ibo == 7 && utpc_6)){
          dx = x_comb_6 - x_comb_7 + AngleCorrection(6, 7, track.SlopeX(), GEOMETRY);
          h1["track_diff67_comb"]->Fill(dx);
          h2["track_diff67_comb_theta"]->Fill(theta(track.SlopeX()), dx);
          h2["track_diff67_comb_abstheta"]->Fill(fabs(theta(track.SlopeX())), dx);
        }

      }

      // x
      // x_comb   = fraction*x_utpc + (1-fraction)*x_clus;
      // x_track  = track.SlopeX()*z_half + track.ConstX();
      // residual = x_comb - x_track;
      // if (false){
      //   std::cout << Form("Evt %3d, Bo %d :: unc(bary) = %6.2f, unc(utpc) = %6.2f => unc(comb) = %6.2f", evt, ibo, unc_bary, unc_utpc, unc_comb);
      //   std::cout << Form(" :: x(bary) = %6.2f, x(utpc) = %6.2f => x(comb) = %6.2f :: x(track) = %6.2f", x_clus, x_utpc, x_comb, x_track) << std::endl;
      // }
      // h2[Form("track_N1_theta_x_vs_comb_%i", ibo)]->Fill(theta(track.SlopeX()), residual);
      // if (ibo == 0) x_comb_0 = x_comb;
      // if (ibo == 1) x_comb_1 = x_comb;
      // if (ibo == 6) x_comb_6 = x_comb;
      // if (ibo == 7) x_comb_7 = x_comb;
    }
  
//     // board v board: tpc
//     if (x_utpc_0 > -100 && x_utpc_1 > -100){
//       dx = x_utpc_0 - x_utpc_1 + AngleCorrection(0, 1, track.SlopeX(), GEOMETRY);
//       h1["track_diff01_utpc"]      ->Fill(dx);
//       h2["track_diff01_utpc_theta"]->Fill(theta(track.SlopeX()), dx);
//     }
//     if (x_utpc_6 > -100 && x_utpc_7 > -100){
//       dx = x_utpc_6 - x_utpc_7 + AngleCorrection(6, 7, track.SlopeX(), GEOMETRY);
//       h1["track_diff67_utpc"]      ->Fill(dx);
//       h2["track_diff67_utpc_theta"]->Fill(theta(track.SlopeX()), dx);
//     }
    
//     // board v board: combo
//     if (hit_0 && hit_1){
//       dx = x_comb_0 - x_comb_1 + AngleCorrection(0, 1, track.SlopeX(), GEOMETRY);
//       h1["track_diff01_comb"]      ->Fill(dx);
//       h2["track_diff01_comb_theta"]->Fill(theta(track.SlopeX()), dx);
//     }
//     if (hit_6 && hit_7){
//       dx = x_comb_6 - x_comb_7 + AngleCorrection(6, 7, track.SlopeX(), GEOMETRY);
//       h1["track_diff67_comb"]      ->Fill(dx);
//       h2["track_diff67_comb_theta"]->Fill(theta(track.SlopeX()), dx);
//     }
    
  }

  // write to file
  fout->cd();
  fout->mkdir("histograms");
  fout->cd("histograms");

  for (auto kv: h1)
    kv.second->Write();
  for (auto kv: h2)
    kv.second->Write();

  fout->Close();    
  std::cout << std::endl;
}

void progress(double time_diff, int nprocessed, int ntotal) {
  double rate = (double)(nprocessed+1)/time_diff;
  std::cout.precision(1);
  std::cout << "\r > " << nprocessed << " / " << ntotal 
            << " | "   << std::fixed << 100*(double)(nprocessed)/(double)(ntotal) << "%"
            << " | "   << std::fixed << rate << "Hz"
            << " | "   << std::fixed << time_diff/60 << "m elapsed"
            << " | "   << std::fixed << (double)(ntotal-nprocessed)/(rate*60) << "m remaining"
            << std::flush;
  std::cout.precision(6);
}

double channel_from_x(double xpos, int board, GeoOctuplet* geo, int begin) {
  int SignChannel = geo->Get(board).SignChannel();
  double xorigin  = geo->Get(board).Origin().X();
  double alpha    = (begin) ? geo->Get(board).StripAlpha() : 0.0;
  double channel  = ((xpos - xorigin - tan(alpha)*200) / (SignChannel*0.4)) + 256.5;
  return channel;
}

double theta(double slope){
  return atan(slope)*180/3.14159;
}

int OkayForTPC(MMCluster* clus){
  if (clus->size() < 3)
    return 0;
  //for (auto hit: *clus)
  //  if (hit->TrigBCID() - hit->BCID() > 34 || hit->TrigBCID() - hit->BCID() < 26)
  //    return 0;
  //for (auto hit: *clus)
  //  if (hit->BCID() % 4 == 1)
  //    return 0;
  return 1;
}

std::tuple<double, double> fit(std::vector<double> xs, std::vector<double> zs){
  double slope  = 0.0;
  double offset = 0.0;
  double avg_x  = std::accumulate(xs.begin(), xs.end(), 0.0) / (double)(xs.size());
  double avg_z  = std::accumulate(zs.begin(), zs.end(), 0.0) / (double)(zs.size());
  double sum_sq_z = std::inner_product(zs.begin(), zs.end(), zs.begin(), 0.0);
  for (unsigned int i = 0; i < xs.size(); i++)
    slope += xs[i] * ( (zs[i] - avg_z) / (sum_sq_z - zs.size()*pow(avg_z, 2)));
  offset = avg_x - slope*avg_z;
  return std::make_pair(slope, offset);
}

std::tuple<double, double, double, double, double, double, double, double, double> fit_root(std::vector<double> xs, std::vector<double> zs, int invert){

  TGraphErrors* graph = 0;
  TF1* fit = 0;
  TFitResultPtr result = 0;

  double dummy = 0;
  double unc_z = 0.60;
  std::vector<double> exs = {};
  std::vector<double> ezs = {};
  for (auto z: zs){
    exs.push_back(0.0);
    ezs.push_back(unc_z);
    dummy = dummy + z;
  }

  if (invert) graph  = new TGraphErrors(int(xs.size()), &xs[0], &zs[0], &exs[0], &ezs[0]);
  else        graph  = new TGraphErrors(int(xs.size()), &zs[0], &xs[0], &ezs[0], &exs[0]);

  //result = graph->Fit("pol1", "QSF");
  //fit    = graph->GetFunction("pol1");
  fit = new TF1("the_fit", "[0] + x*[1]", -5.0, 215.0);
  fit->SetParameter(0, 0.0);
  fit->SetParameter(1, 0.0);
  // result = graph->Fit("the_fit", "QSF");
  result = graph->Fit("the_fit", "QS");

  double slope, offset, chi2, ndf, prob;
  slope  = fit->GetParameter(1);
  offset = fit->GetParameter(0);
  chi2   = fit->GetChisquare();
  ndf    = fit->GetNDF();
  prob   = fit->GetProb();

  TMatrixD cov = result->GetCovarianceMatrix();
  double cov00 = cov[0][0];
  double cov01 = cov[0][1];
  double cov10 = cov[1][0];
  double cov11 = cov[1][1];

  delete fit;
  delete graph;

  return std::make_tuple(slope, offset, chi2, ndf, prob, cov00, cov01, cov10, cov11);
}

double TPCUncertainty(double z, double m, double b, double sigma2_b, double sigma2_m, double sigma_bm, int invert){
  double unc2 = -1;
  if (invert){
    unc2 = sigma2_m*(z-b)*(z-b)/(m*m) + sigma2_b + 2*sigma_bm*(z-b)/m;
    //unc2 = unc2 / sigma2_m;
    unc2 = unc2 / (m*m);
  }
  else{
    unc2 = (z*z)*sigma2_m + sigma2_b + 2*z*sigma_bm;
  }
  return sqrt(unc2);
}

double AngleCorrection(int i, int j, double slope, GeoOctuplet* GEOMETRY){
  GeoPlane plane_i = GEOMETRY->Get(i);
  GeoPlane plane_j = GEOMETRY->Get(j);
  double dz = std::fabs(plane_i.Origin().Z() - plane_j.Origin().Z());
  double dx = dz * slope;
  return dx;
}

double OffsetByBoard(int board){
  // via oct_anaplus.C
  std::vector<double> offset = {0, 0, 0, 0, 
                                0, 0, 0, 0};
  offset[0]=-0.22 ; //0.53; //0.55+.13-.03;
  offset[1]=+0.23 ; // -0.36;  //-.257-.08+.08;
  offset[2]=-0.77 ;  //0.9+.25;
  offset[3]=+0.15 ;  //-.175-.07-.08;
  offset[4]=-0.02 ;  //0.44+.1;
  offset[5]=0.52  ;  //.02-.03;
  offset[6]=0.1   ; //.32+.07;
  offset[7]=0.70  ;  //.06+.02+0.1;
  if (board < 0 || board > 7)
    std::cout << "Get ready for a seg fault: " << board << std::endl;
  return offset[board];
}

double OffsetByChannel(int vmm, int ch, std::map< std::tuple<int,int>, double > offsets){
  auto key = make_tuple(vmm, ch);
  if (offsets.count(key))
    return offsets[key];
  else
    return 0.0;
}

std::vector<double> ChannelOffsets(std::string filename, int board, int this_vmm){
  if (board == 0 && this_vmm == 0)
    std::cout << "Channel offsets: " << filename << std::endl;
  int vmm_ch, vmm, ch;
  std::string line;
  std::ifstream file(filename);
  std::vector<double> corr = {0, 0, 0, 0,
                              0, 0, 0, 0};
  std::vector<double> corrections = {};

  // index channels from 1
  corrections.push_back(0);

  // std::map< std::tuple<int,int>, double > the_map = {};
  if(file.is_open()){
    while( getline(file, line) ){
      std::stringstream sline;
      sline << line;
      sline >> vmm_ch;
      sline >> corr[0];
      sline >> corr[1];
      sline >> corr[2];
      sline >> corr[3];
      sline >> corr[4];
      sline >> corr[5];
      sline >> corr[6];
      sline >> corr[7];
      vmm = (vmm_ch - 1) / 64;
      ch  = (vmm_ch - 1) % 64 + 1;
      if (vmm == this_vmm)
        corrections.push_back(corr[board]);
      // the_map[make_tuple(vmm, ch)] = corr[board];
      //for (int board = 0; board < 8; board++)
      //  the_map[make_tuple(board, vmm, ch)] = corr[board];
    }
  }
  file.close();
  return corrections;
  //return the_map;
}
