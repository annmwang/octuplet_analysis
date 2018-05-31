///
///  \file   TunaAnalysisTPSR.C
///
///  \author Tunaface
///          (tuna@cern.ch)
///
///  \date   2018 Feb
///

#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TF1.h"
#include <iostream>
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
double xpos(TPHit* hit, GeoOctuplet* geo);
double xpos(MMCluster* hit, GeoOctuplet* geo);
int triggerableX(MMClusterList clusters, int roadsize, int overlap);

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
  bool b_nate  = false;
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
    if (strncmp(argv[i+1],"-n",2)==0){
      b_nate = true;
    }
  }

  if(!b_input) std::cout << "Error at Input: please specify  input file (-i flag)" << std::endl;
  if(!b_out)   std::cout << "Error at Input: please specify output file (-o flag)" << std::endl;
  if(!b_input || !b_out) return 0;

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
  DATA->SetTP(1);
  int Nevent = DATA->GetNEntries();

  std::map< string, TH1D* > h1;
  std::map< string, TH2D* > h2;

  h2["track_res2_vs_y"]  = new TH2D("track_res2_vs_y",  " ;y;Sqrt(Sum of res. sq.); Tracks", 220, 0, 220,    100, 0.0,  5.0);
  h2["track_res2_vs_my"] = new TH2D("track_res2_vs_my", ";my;Sqrt(Sum of res. sq.); Tracks", 100, -4.0, 4.0, 100, 0.0,  5.0);
  h2["track_res2_vs_cy"] = new TH2D("track_res2_vs_cy", ";cy;Sqrt(Sum of res. sq.); Tracks", 100, -300, 500, 100, 0.0,  5.0);

  h2["scint_bot_vs_top"] = new TH2D("scint_bot_vs_top", ";SC bottom channel;SC top channel;Tracks", 8, -1.5, 6.5, 8, 14.5, 22.5);
  h2["scint_ch_vs_counts"] = new TH2D("scint_ch_vs_counts", ";SC CH;SC Counts;", 30, -1.5, 28.5, 400, 0, 400);

  h2["track_angle_x_scint_17"] = new TH2D("track_angle_x_scint_17", ";x angle;scint. bottom channel;Tracks, scint. top channel 17", 100, -30, 30, 8, -1.5, 6.5);
  h2["track_angle_x_scint_18"] = new TH2D("track_angle_x_scint_18", ";x angle;scint. bottom channel;Tracks, scint. top channel 18", 100, -30, 30, 8, -1.5, 6.5);
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
  h2["track_hits_vs_time"]      = new TH2D("track_hits_vs_time",     ";Time [days];N(clusters in track);", 1000, 0, 20, 10, -0.5, 9.5);
  h2["track_hits_vs_time_fid"]  = new TH2D("track_hits_vs_time_fid", ";Time [days];N(clusters in track);", 1000, 0, 20, 10, -0.5, 9.5);
  h2["track_hits_vs_time_art"]  = new TH2D("track_hits_vs_time_art", ";Time [days];N(clusters in track);", 1000, 0, 20, 10, -0.5, 9.5);
  h2["track_clusmult_vs_theta"] = new TH2D("track_clusmult_vs_theta", ";#theta#lower[0.5]{track} [degrees];hits in cluster", 100, -35, 35, 13, -0.5, 12.5);

  h2["track_param_x"] = new TH2D("track_param_x", ";x slope; x constant;Tracks",   100, -0.6, 0.6, 100,  -80, 280);
  h2["track_param_y"] = new TH2D("track_param_y", ";y slope; y constant;Tracks",   100,   -4,   4, 100, -300, 500);
  h2["track_angle_6"] = new TH2D("track_angle_6", ";#theta#lower[0.5]{track} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);
  h2["track_angle_7"] = new TH2D("track_angle_7", ";#theta#lower[0.5]{track} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);
  h2["track_angle_8"] = new TH2D("track_angle_8", ";#theta#lower[0.5]{track} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);

  h2["track_N1_board_vs_residual"] = new TH2D("track_N1_board_vs_residual", ";board;x_{cluster} - x_{track, proj.};Tracks", 8, -0.5, 7.5, 200, -5.0, 5.0);
  h2["track_N1_board_vs_art"]      = new TH2D("track_N1_board_vs_art",      ";board;x_{ART} - x_{track, proj.};Tracks",     8, -0.5, 7.5, 400, -20.0, 20.0);
  h2["track_N1_board_vs_utpc"]     = new TH2D("track_N1_board_vs_utpc",     ";board;x_{utpc} - x_{track, proj.};Tracks",    8, -0.5, 7.5, 200, -5.0, 5.0);

  TString name;
  for (ibo = 0; ibo < nboards; ibo++){
    name = Form("track_N1_theta_x_vs_residual_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x theta;x_{cluster} - x_{track, proj.};Tracks", 100, -35, 35, 200, -5.0, 5.0);

    name = Form("track_N1_x_vs_residual_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x_{cluster};x_{cluster} - x_{track, proj.};Tracks", 100, -20, 220, 200, -5.0, 5.0);

    name = Form("track_N1_x_vs_residual_10deg_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x_{cluster};x_{cluster} - x_{track, proj.};Tracks", 100, -20, 220, 200, -5.0, 5.0);

    name = Form("track_N1_theta_x_vs_art_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x theta;x_{ART} - x_{track, proj.};Tracks", 100, -35, 35, 200, -5.0, 5.0);

    name = Form("track_N1_theta_x_vs_art_vetodray_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x theta;x_{ART} - x_{track, proj.};Tracks", 100, -35, 35, 200, -5.0, 5.0);

    name = Form("track_N1_theta_x_vs_utpc_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x theta;x_{utpc} - x_{track, proj.};Tracks", 100, -35, 35, 200, -5.0, 5.0);
  }

  h2["strip_position_vs_board"] = new TH2D("strip_position_vs_board", ";strip number;MMFE number;charge [fC]", 512, 0.5, 512.5, 8, -0.5, 7.5);

  for (ibo = 0; ibo < nboards; ibo++){
    h2[Form("strip_q_vs_ch_%i",    ibo)] = new TH2D(Form("strip_q_vs_ch_%i",    ibo), ";strip number;Charge [fC];strip",    512, 0.5, 512.5, 512,   0,  128);
    h2[Form("strip_pdo_vs_ch_%i",  ibo)] = new TH2D(Form("strip_pdo_vs_ch_%i",  ibo), ";strip number;PDO [counts];strip",   512, 0.5, 512.5, 512,   0, 2048);
    h2[Form("strip_tdo_vs_ch_%i",  ibo)] = new TH2D(Form("strip_tdo_vs_ch_%i",  ibo), ";strip number;TDO [counts];strip",   512, 0.5, 512.5, 256,   0,  256);
    h2[Form("strip_tdoc_vs_ch_%i", ibo)] = new TH2D(Form("strip_tdoc_vs_ch_%i", ibo), ";strip number;TDO corr. [ns];strip", 512, 0.5, 512.5, 110, -10,  100);
    h2[Form("strip_dbc_vs_ch_%i",  ibo)] = new TH2D(Form("strip_dbc_vs_ch_%i",  ibo), ";strip number;#Delta BC;strip",      512, 0.5, 512.5,  64,   0,   64);
    h2[Form("strip_time_vs_ch_%i", ibo)] = new TH2D(Form("strip_time_vs_ch_%i", ibo), ";strip number;Time [ns];strip",      512, 0.5, 512.5, 100, 300, 1000);
    h2[Form("strip_zpos_vs_ch_%i", ibo)] = new TH2D(Form("strip_zpos_vs_ch_%i", ibo), ";strip number;z_{drift} [mm];strip", 512, 0.5, 512.5,  60, -10,   20);
    h2[Form("strip_zdri_vs_ch_%i", ibo)] = new TH2D(Form("strip_zdri_vs_ch_%i", ibo), ";strip number;z_{drift} [mm];strip", 512, 0.5, 512.5,  60, -10,   20);
    h2[Form("strip_zres_vs_ch_%i", ibo)] = new TH2D(Form("strip_zres_vs_ch_%i", ibo), ";strip number;#Delta z [mm];strip",  512, 0.5, 512.5, 200, -15,   15);
  }
  h2["dups_vs_ch"] = new TH2D("dups_vs_ch", ";strip number;MMFE number;Duplicates", 512, 0.5, 512.5, 8, -0.5, 7.5);

  for (ibo = 0; ibo < nboards; ibo++){
    h2[Form("strip_zpos_vs_ztrack_%i", ibo)] = new TH2D(Form("strip_zpos_vs_ztrack_%i", ibo), ";z_{track};z_{drift} [mm];strip", 100, -5, 10, 100, -2, 8);
  }

  h2["clus_vs_board"]          = new TH2D("clus_vs_board",          ";MMFE number;clusters;Events",            8, -0.5, 7.5, 32, -0.5, 31.5);
  h1["clus"]                   = new TH1D("clus",                   ";clusters;Events",                        50, -0.5, 49.5);
  h2["hits_vs_board"]          = new TH2D("hits_vs_board",          ";MMFE number;strips;Events",              8, -0.5, 7.5, 32, -0.5, 31.5);
  h2["dups_vs_board"]          = new TH2D("dups_vs_board",          ";MMFE number;duplicate strips;Events",    8, -0.5, 7.5, 32, -0.5, 31.5);
  h2["hits_per_clus_vs_board"] = new TH2D("hits_per_clus_vs_board", ";MMFE number;hits in a cluster;Clusters", 8, -0.5, 7.5, 13, -0.5, 12.5);
  h2["timediff_vs_board"]      = new TH2D("timediff_vs_board",      ";MMFE number;#DeltaBCID;Events",          8, -0.5, 7.5, 100, -0.5, 99.5);
  h1["timediff_track"]         = new TH1D("timediff_track",         ";#DeltaBCID;Strips on-track",             100, -0.5, 99.5);
  h1["timediff_dupli"]         = new TH1D("timediff_dupli",         ";#DeltaBCID;Duplicate strips",            100, -0.5, 99.5);
  h2["timediff_vs_charge"]     = new TH2D("timediff_vs_charge",     ";charge [fC];#DeltaBCID;Events",          200, 0, 200, 100, -0.5, 99.5);

  h2["track_exp_hit"] = new TH2D("track_exp_hit", ";VMM;Board;Trigger", 10, -1.5, 8.5, 17, -0.75, 7.75);
  h2["track_obs_hit"] = new TH2D("track_obs_hit", ";VMM;Board;Trigger", 10, -1.5, 8.5, 17, -0.75, 7.75);


  h1["trig_dbc_scint"]      = new TH1D("trig_dbc_scint",      ";#DeltaBC(trig, scint);", 8191, -4095.5, 4095.5);
  h2["trig_per_event"]      = new TH2D("trig_per_event",      ";N(triggers);N(boards in MM track)", 20, -0.5, 19.5, 7, 2.5, 9.5);
  h1["trig_ambig"]          = new TH1D("trig_ambig",          ";N(spurious ART, first) - N(spurious ART, other)", 13, -6.5, 6.5);
  h2["trig_dtheta_vsNX"]    = new TH2D("trig_dtheta_vsNX",    ";#theta(TP)-#theta(MM);N(ART, X);Events", 1000, -500, 500, 5, 0.5, 5.5);
  h2["trig_dtheta_vsNX_ok"] = new TH2D("trig_dtheta_vsNX_ok", ";#theta(TP)-#theta(MM);N(ART, X);Events", 1000, -500, 500, 5, 0.5, 5.5);
  h1["trig_mm"]             = new TH1D("trig_mm",             ";N(MM clusters);",                        10, -0.5, 9.5);
  h1["trig_art"]            = new TH1D("trig_art",            ";N(ART);",                                10, -0.5, 9.5);
  h1["trig_theta"]          = new TH1D("trig_theta",          ";#theta [deg.];",                         100, -35, 35);
  h1["trig_eventmm"]        = new TH1D("trig_eventmm",        ";MM Event Number;",                       1000, 0, 1000);
  h2["trig_dbc_vs_N"]       = new TH2D("trig_dbc_vs_N",       ";#DeltaBC;N(ART)",                        10, -0.5, 9.5, 7, 2.5, 9.5);
  h2["trig_dbc_vs_theta"]   = new TH2D("trig_dbc_vs_theta",   ";#DeltaBC;#theta [deg]",                  10, -0.5, 9.5, 500, -30, 20);
  h2["trig_dbc_vs_evt"]     = new TH2D("trig_dbc_vs_evt",     ";MM Event Number;#DeltaBC",               1000, -100, 1000, 10, -0.5, 9.5);
  h1["trig_dbc_pairs"]      = new TH1D("trig_dbc_pairs",      ";#DeltaBCID, pairs;",                     23, -11.5, 11.5);
  h2["trig_art_vs_evt"]     = new TH2D("trig_art_vs_evt",     ";MM Event Number;N(ART in trigger)",      1000, -100, 1000, 10, -0.5, 9.5);

  h2["track_angle_N_denom"] = new TH2D("track_angle_N_denom", ";#theta#lower[0.5]{xz} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 450, -30, 30, 200, -150, 150);
  h2["track_angle_N_numer"] = new TH2D("track_angle_N_numer", ";#theta#lower[0.5]{xz} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 450, -30, 30, 200, -150, 150);
  h2["track_angle_N_denom_vetodray"] = new TH2D("track_angle_N_denom_vetodray", ";#theta#lower[0.5]{xz} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 450, -30, 30, 200, -150, 150);
  h2["track_angle_N_numer_vetodray"] = new TH2D("track_angle_N_numer_vetodray", ";#theta#lower[0.5]{xz} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 450, -30, 30, 200, -150, 150);
  h2["track_pos_N_denom"]   = new TH2D("track_pos_N_denom",   ";xpos [mm];ypos [mm];Tracks", 200, -10, 220, 200, -100, 300);
  h2["track_pos_N_numer"]   = new TH2D("track_pos_N_numer",   ";xpos [mm];ypos [mm];Tracks", 200, -10, 220, 200, -100, 300);

  h2["trig_dx_vs_theta_0"]  = new TH2D("trig_dx_vs_theta_0", ";#theta(MM);#Deltax(ART, MMFE track) [mm]", 100, -35, 35, 1000, -16, 16);
  h2["trig_dx_vs_theta_1"]  = new TH2D("trig_dx_vs_theta_1", ";#theta(MM);#Deltax(ART, MMFE track) [mm]", 100, -35, 35, 1000, -16, 16);
  h2["trig_dx_vs_theta_6"]  = new TH2D("trig_dx_vs_theta_6", ";#theta(MM);#Deltax(ART, MMFE track) [mm]", 100, -35, 35, 1000, -16, 16);
  h2["trig_dx_vs_theta_7"]  = new TH2D("trig_dx_vs_theta_7", ";#theta(MM);#Deltax(ART, MMFE track) [mm]", 100, -35, 35, 1000, -16, 16);

  // cataloging duplicates
  MMFE8Hits duplicates;

  // collecting clusters and the nominal fit
  std::vector<MMClusterList> clusters_perboard;
  MMClusterList clusters_all;
  MMClusterList clusters_road;
  MMClusterList clusters_x;
  MMClusterList clusters_tp;
  MMTrack track;
  MMTrack track_tp;

  // N-1 fits and resolutions
  MMClusterList clusters_N1;
  MMTrack track_N1;

  // miscellaneous helpers
  double start_time = -1;
  double time_since_start = -1;
  double days_since_start = -1;
  double seconds_per_day  = 60*60*24;
  int EventNum4Hist = 0;
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
  TCanvas* can;
  GeoPlane plane;
  double z_middle = 0.0;
  for (ibo = 0; ibo < nboards; ibo++)
    z_middle += GEOMETRY->Get(ibo).Origin().Z();
  z_middle /= (double)(nboards);

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
      if (!b_nate)
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

        h2["timediff_vs_board"]->Fill(ibo, hit.TrigBCID() - hit.BCID());
        h2["timediff_vs_charge"]->Fill(hit.Charge(), hit.TrigBCID() - hit.BCID());

        if( !PACMAN->IsGoodHit(hit) )
          continue;

        h2[Form("strip_q_vs_ch_%i",   ibo)]->Fill(hit.Channel(), hit.Charge());
        h2[Form("strip_pdo_vs_ch_%i", ibo)]->Fill(hit.Channel(), hit.PDO());
        h2[Form("strip_tdo_vs_ch_%i", ibo)]->Fill(hit.Channel(), hit.TDO());

        h2[Form("strip_tdoc_vs_ch_%i", ibo)]->Fill(hit.Channel(), hit.Time());
        h2[Form("strip_dbc_vs_ch_%i",  ibo)]->Fill(hit.Channel(), hit.DeltaBC());
        //h2[Form("strip_time_vs_ch_%i", ibo)]->Fill(hit.Channel(), hit.DriftTime(deltaT));
        //h2[Form("strip_zpos_vs_ch_%i", ibo)]->Fill(hit.Channel(), hit.DriftTime(deltaT) * vdrift);
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

    h1["clus"]->Fill(clusters_all.size());

    // fit it!
    for (auto botpair: DATA->sc_EventHits.GetBotPair()){
      clusters_road.Reset();
      clusters_road = FILTERER->FilterClustersScint(clusters_all, *GEOMETRY, botpair.first->Channel(), DATA->mm_EventNum, debug);
      track  = FITTER ->Fit(clusters_road, *GEOMETRY, DATA->mm_EventNum);
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
        if (toppair.first->Channel() == 17) h2["track_angle_x_scint_17"]->Fill(theta(track.SlopeX()), botpair.first->Channel());
        if (toppair.first->Channel() == 18) h2["track_angle_x_scint_18"]->Fill(theta(track.SlopeX()), botpair.first->Channel());
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

      for (auto clus_other: clusters_road)
        if (clus != clus_other)
          clusters_N1.AddCluster(*clus_other);
      
      track_N1 = FITTER->Fit(clusters_N1, *GEOMETRY, DATA->mm_EventNum);
      residual = GEOMETRY->GetResidualX(*clus, track_N1);
      ibo      = GEOMETRY->Index(clus->MMFE8());
      plane    = GEOMETRY->Get(GEOMETRY->Index(clus->MMFE8()));
      xclus    = plane.Origin().X() + plane.LocalXatYbegin(clus->Channel());

      h2["track_N1_board_vs_residual"]                ->Fill(ibo,                   residual);
      h2[Form("track_N1_theta_x_vs_residual_%i", ibo)]->Fill(theta(track.SlopeX()), residual);
      h2[Form("track_N1_x_vs_residual_%i",       ibo)]->Fill(xclus,                 residual);
      if (theta(track.SlopeX()) < -10)
        h2[Form("track_N1_x_vs_residual_10deg_%i", ibo)]->Fill(xclus, residual);
    }

    // TP!
    // ------------------------------
    if (!track.IsTrigCand(3, 3))
      continue;
    if (track.PointZ(z_middle).X() > 170)
      continue;
    //if (!triggerableX(clusters_road, 8, 4))
    //  continue;

    int dB = DATA->tp_EventTracks.FetchSciOffset(DATA->RunNum);
    DATA->tp_EventTracks.SetSciBCID(DATA->sc_EventHits.TPBCID(), DATA->sc_EventHits.TPph());
    DATA->tp_EventTracks.SetSciOffset(dB);

    // count # clusters per board                                                                                                                 
    vector <int> clus_boards;
    bool poss_dray = false;
    for (auto clus : clusters_all){
      int m_ib = clus[0].MMFE8Index();
      if (std::find(clus_boards.begin(),clus_boards.end(),m_ib) == clus_boards.end()){
        clus_boards.push_back(m_ib);
      }
      else
        poss_dray = true;
    }

    TPTrack* neo = DATA->tp_EventTracks.Highlander(clusters_road, true, 10);
    //TPTrack* neo = DATA->tp_EventTracks.Highlander(clusters_road, true, 10000);
    if (neo)
      h1["trig_dbc_scint"]->Fill(DATA->tp_EventTracks.deltaBCID(neo->BCID()));

    // trigger-finding efficiency
    h2["track_hits_vs_time_art"]->Fill(days_since_start, clusters_road.size());
    h2["track_angle_N_denom"]->Fill(theta(track.SlopeX()), theta(track.SlopeY()));
    h2["track_pos_N_denom"]  ->Fill(track.PointZ(z_middle).X(), track.PointZ(z_middle).Y());

    if (!poss_dray) {
      h2["track_angle_N_denom_vetodray"]->Fill(theta(track.SlopeX()), theta(track.SlopeY()));
    }
    // tracks without friends
    if (false && !neo && std::fabs(theta(track.SlopeX())) < 1. && track.IsTrigCand(3, 3) && poss_dray){
      //std::cout << "event had no tp track, but had possible delta ray" << std::endl;
      //std::cout << Form("%7i :: %i %i", DATA->mm_EventNum, DATA->mm_EventHits.Time_S(), DATA->mm_EventHits.Time_NS()) << std::endl;
      fout->cd("event_displays");
      can = Plot_Track2D(Form("track2D_%05d_MMall", DATA->mm_EventNum), track, *GEOMETRY, &clusters_all);
      can->Write(0, TObject::kOverwrite);
      delete can;
      can = Plot_Track2D(Form("track2D_%05d_MMfit", DATA->mm_EventNum), track, *GEOMETRY, &clusters_road);
      can->Write(0, TObject::kOverwrite);
      delete can;
    }

    // onto the numerators
    if (!neo)
      continue;
    if (!neo->IsTrigCand())
      continue;
    h2["track_angle_N_numer"]->Fill(theta(track.SlopeX()), theta(track.SlopeY()));
    h2["track_pos_N_numer"]  ->Fill(track.PointZ(z_middle).X(), track.PointZ(z_middle).Y());

    if (!poss_dray){
      h2["track_angle_N_numer_vetodray"]->Fill(theta(track.SlopeX()), theta(track.SlopeY()));
    }

    // event dump for nathan!
    if (b_nate){
      if (neo->EventNumGBT() != -1)
        {
          std::cout << Form("ALL %8d ", neo->EventNumGBT());
          for (auto tr: DATA->tp_EventTracks)
            std::cout << Form("%8d ", tr->EventNum());
          std::cout << std::endl;
          for (auto tr: DATA->tp_EventTracks)
            std::cout << Form("FIT %8d ", tr->EventNum()) << std::endl;
          std::cout << Form("GBT %8d ", neo->EventNumGBT()) << std::endl;
          std::cout << std::endl;
        }
    }

    // residuals
    // ------------------------------
    double xtrack = 0;
    if (clusters_road.size() >= 7){
      for (auto art: *neo){
        
        clusters_N1.Reset();
        track_N1.Reset();
        
        for (auto clus_other: clusters_road)
          if (art->MMFE8() != clus_other->MMFE8())
            clusters_N1.AddCluster(*clus_other);
        
        track_N1 = FITTER->Fit(clusters_N1, *GEOMETRY, DATA->mm_EventNum);
        ibo      = art->MMFE8Index();
        xtrack   = track_N1.SlopeX()*zboard[ibo] + track_N1.ConstX();
        residual = xpos(art, GEOMETRY) - xtrack;
        
        h2["track_N1_board_vs_art"]                ->Fill(ibo,                   residual);
        h2[Form("track_N1_theta_x_vs_art_%i", ibo)]->Fill(theta(track.SlopeX()), residual);

        if (!poss_dray)
          h2[Form("track_N1_theta_x_vs_art_vetodray_%i", ibo)]->Fill(theta(track.SlopeX()), residual);

        if (fabs(residual) > 3. and !poss_dray and (ibo < 2 || ibo > 5) ){
          clusters_tp.Reset();
          track_tp.Reset();
          for (auto art: *neo)
            clusters_tp.AddCluster(MMCluster(MMHit(art->MMFE8(), art->VMM(), art->VMMChannel(), DATA->RunNum)));
          track_tp = FITTER->Fit(clusters_tp, *GEOMETRY, DATA->mm_EventNum);
          
          fout->cd("event_displays");
          can = Plot_Track2D(Form("track2D_%05d_MMall", DATA->mm_EventNum), track, *GEOMETRY, &clusters_all);
          can->Write(0, TObject::kOverwrite);
          delete can;
          can = Plot_Track2D(Form("track2D_%05d_MMfit", DATA->mm_EventNum), track, *GEOMETRY, &clusters_road);
          can->Write(0, TObject::kOverwrite);
          delete can;
          can = Plot_Track2D(Form("track2D_%05d_TPfit", DATA->mm_EventNum), track_tp, *GEOMETRY, &clusters_tp);
          can->Write(0, TObject::kOverwrite);
          delete can;
        }

      }
    }
    
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

double xpos(TPHit* hit, GeoOctuplet* geo){
  return geo->Get(hit->MMFE8Index()).LocalXatYbegin(hit->Channel())
    + geo->Get(hit->MMFE8Index()).Origin().X();
}

double xpos(MMCluster* hit, GeoOctuplet* geo){
  return geo->Get(hit->MMFE8Index()).LocalXatYbegin(hit->Channel())
    + geo->Get(hit->MMFE8Index()).Origin().X();
}

int triggerableX(MMClusterList clusters, int roadsize, int overlap){

  int tmp_x  = 0;
  int nroads = 0;
  int nx     = 0;
  std::vector<int> trigger_x = {};

  // collect the x hits
  tmp_x = 0;
  trigger_x.clear();
  for (auto clus: clusters){
    if (!clus->isX())
      continue;
    tmp_x = clus->Channel() - 1;
    if (clus->MMFE8Index() == 0 ||
        clus->MMFE8Index() == 3 ||
        clus->MMFE8Index() == 5 ||
        clus->MMFE8Index() == 6)
      tmp_x = 511 - tmp_x;
    tmp_x += 64;
    trigger_x.push_back( std::round(tmp_x) );
  }

  // evaluate the roads
  nroads = 512 / roadsize;
  for (int ro = 0; ro < nroads; ro++){
    nx = 0;
    for (auto x: trigger_x)
      if (x/roadsize == ro || (x-overlap)/roadsize == ro)
        nx++;
    if (nx >= 3)
      return 1;
  }
  return 0;
}
