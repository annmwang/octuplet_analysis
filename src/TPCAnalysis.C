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
std::tuple<double, double, double, double, double, double> fit_root(std::vector<double> xs, std::vector<double> zs);
int OkayForTPC(MMCluster* clus);

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
  h2["track_hits_vs_time"]      = new TH2D("track_hits_vs_time",     ";Time [days];N(clusters in track);", 1000, 0, 20, 10, -0.5, 9.5);
  h2["track_hits_vs_time_fid"]  = new TH2D("track_hits_vs_time_fid", ";Time [days];N(clusters in track);", 1000, 0, 20, 10, -0.5, 9.5);
  h2["track_clusmult_vs_theta"] = new TH2D("track_clusmult_vs_theta", ";#theta#lower[0.5]{track} [degrees];hits in cluster", 100, -35, 35, 13, -0.5, 12.5);

  h2["track_param_x"] = new TH2D("track_param_x", ";x slope; x constant;Tracks",   100, -0.6, 0.6, 100,  -80, 280);
  h2["track_param_y"] = new TH2D("track_param_y", ";y slope; y constant;Tracks",   100,   -4,   4, 100, -300, 500);
  h2["track_angle_6"] = new TH2D("track_angle_6", ";#theta#lower[0.5]{track} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);
  h2["track_angle_7"] = new TH2D("track_angle_7", ";#theta#lower[0.5]{track} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);
  h2["track_angle_8"] = new TH2D("track_angle_8", ";#theta#lower[0.5]{track} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);

  h2["track_N1_board_vs_residual"] = new TH2D("track_N1_board_vs_residual", ";board;x_{cluster} - x_{track, proj.};Tracks", 8, -0.5, 7.5, 200, -5.0, 5.0);
  h2["track_N1_board_vs_utpc"]     = new TH2D("track_N1_board_vs_utpc",     ";board;x_{utpc} - x_{track, proj.};   Tracks", 8, -0.5, 7.5, 200, -5.0, 5.0);

  TString name;
  for (ibo = 0; ibo < nboards; ibo++){
    name = Form("track_N1_theta_x_vs_residual_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x theta;x_{cluster} - x_{track, proj.};Tracks", 100, -35, 35, 200, -5.0, 5.0);

    name = Form("track_N1_x_vs_residual_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x_{cluster};x_{cluster} - x_{track, proj.};Tracks", 100, -20, 220, 200, -5.0, 5.0);

    name = Form("track_N1_x_vs_residual_10deg_%i", ibo);
    h2[name.Data()] = new TH2D(name, ";x_{cluster};x_{cluster} - x_{track, proj.};Tracks", 100, -20, 220, 200, -5.0, 5.0);

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
  double z_utpc_check  = 0.0;
  double t_utpc  = 0.0;
  double x_utpc  = 0.0;
  double x_clus  = 0.0;
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
    
    if (evt > 500)
      break;

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
        h2[Form("strip_time_vs_ch_%i", ibo)]->Fill(hit.Channel(), hit.DriftTime(deltaT));
        h2[Form("strip_zpos_vs_ch_%i", ibo)]->Fill(hit.Channel(), hit.DriftTime(deltaT) * vdrift);
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
    if (clusters_road.size() < 7)
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
      ibo      = GEOMETRY->Index(clus->MMFE8());
      plane    = GEOMETRY->Get(GEOMETRY->Index(clus->MMFE8()));
      xclus    = plane.Origin().X() + plane.LocalXatYbegin(clus->Channel());

      h2["track_N1_board_vs_residual"]                ->Fill(ibo,                   residual);
      h2[Form("track_N1_theta_x_vs_residual_%i", ibo)]->Fill(theta(track.SlopeX()), residual);
      h2[Form("track_N1_x_vs_residual_%i",       ibo)]->Fill(xclus,                 residual);
      if (theta(track.SlopeX()) < -10)
        h2[Form("track_N1_x_vs_residual_10deg_%i", ibo)]->Fill(xclus, residual);

      if (evt < 1000 && clusters_road.size() == 8 && ibo == 4 && residual > 4){
        can = Plot_Track2D(Form("track2D_%05d_wtf_all", DATA->mm_EventNum), track_N1, *GEOMETRY, &clusters_all);
        fout->cd("event_displays");
        can->Write();
        delete can;
        can = Plot_Track2D(Form("track2D_%05d_wtf_fit", DATA->mm_EventNum), track_N1, *GEOMETRY, &clusters_road);
        fout->cd("event_displays");
        can->Write();
        delete can;
      }


      if (DATA->mm_EventNum == 40 && false){
        plane = GEOMETRY->Get(GEOMETRY->Index(clus->MMFE8()));
        std::cout << Form("Event %2i | Board %2i | x = %7.4f z = %7.3f | residual: %7.4f | track w/out this cluster: mx = %7.4f cx = %7.4f my = %7.4f cy = %7.4f",
                          DATA->mm_EventNum, ibo,
                          plane.Origin().X() + plane.LocalXatYbegin(clus->Channel()), plane.Origin().Z(),
                          residual,
                          track_N1.SlopeX(), track_N1.ConstX(), track_N1.SlopeY(), track_N1.ConstY()) << std::endl;
      }
    }

    // uTPC, finally
    // ------------------------------------------------------------------
    track = FITTER->Fit(clusters_road, *GEOMETRY, DATA->mm_EventNum);
    if (std::fabs(theta(track.SlopeX())) < 10)
      continue;
    if (clusters_road.size() < 7)
      continue;

    int ntpc = 0;
    double fiducial_x    =  1.8;
    double fiducial_z_hi =  6.6;
    double fiducial_z_lo = -0.8;

    for (auto clus: clusters_road){
      
      if (!OkayForTPC(clus))
        continue;

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
      // if (ibo == 2 || ibo == 3 || ibo == 4 || ibo == 5) continue;
      sign = (ibo==0 || ibo==2 || ibo==4 || ibo==6) ? -1.0 : 1.0;

      plane = GEOMETRY->Get(ibo);

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
        t_utpc  = hit->DriftTime(deltaT);

        if (fabs(x_clus - x_utpc) > fiducial_x)
          continue;
        if (vdrift*t_utpc > fiducial_z_hi)
          continue;
        if (vdrift*t_utpc < fiducial_z_lo)
          continue;

        ntpc++;
        z_utpc  = zboard[ibo] + vdrift * t_utpc * sign;
        z_track = (x_utpc - track_N1.ConstX()) / track_N1.SlopeX();

        // rose-colored lens for suspicious BCIDs
        if (hit->SuspiciousBCID()){
          z_utpc_check = zboard[ibo] + vdrift * (t_utpc - 25) * sign;
          if (fabs(z_track - z_utpc_check) < fabs(z_track - z_utpc))
              z_utpc = z_utpc_check;
        }

        h2[Form("strip_zres_vs_ch_%i",     ibo)]->Fill(hit->Channel(), z_utpc - z_track);
        h2[Form("strip_zdri_vs_ch_%i",     ibo)]->Fill(hit->Channel(), sign*(z_utpc - zboard[ibo]));
        h2[Form("strip_zpos_vs_ztrack_%i", ibo)]->Fill(sign*(z_track - zboard[ibo]), vdrift*t_utpc);
        xs.push_back(x_utpc);
        zs.push_back(z_utpc);
      }

      if (ntpc < 3)
        continue;

      // local fit
      double slope, offset, x_fit, x_track, z_half;
      double cov00, cov01, cov10, cov11;
      std::tie(slope, offset, cov00, cov01, cov10, cov11) = fit_root(xs, zs);
      z_half   = plane.Origin().Z();
      x_fit    = slope*z_half + offset;
      x_track  = track_N1.SlopeX()*z_half + track_N1.ConstX();
      residual = x_fit - x_track;
      h2["track_N1_board_vs_utpc"]->Fill(ibo, residual);
      h2[Form("track_N1_theta_x_vs_utpc_%i", ibo)]->Fill(theta(track.SlopeX()), residual);

      // uncertainty on x_fit

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

std::tuple<double, double, double, double, double, double> fit_root(std::vector<double> xs, std::vector<double> zs){
  TGraph* graph = 0;
  TF1*    fit   = 0;
  graph = new TGraph(int(xs.size()), &xs[0], &zs[0]);
  TFitResultPtr result = graph->Fit("pol1", "QS");
  fit = graph->GetFunction("pol1");

  double slope, offset;
  slope  = 1.0/fit->GetParameter(1);
  offset = -1.0*fit->GetParameter(0)/fit->GetParameter(1);

  TMatrixD cov = result->GetCovarianceMatrix();
  double cov00 = cov[0][0];
  double cov01 = cov[0][1];
  double cov10 = cov[1][0];
  double cov11 = cov[1][1];
  std::cout << pow(fit->GetParError(0), 2) << " " << pow(fit->GetParError(1), 2) << std::endl;
  std::cout << cov00 << " " << cov01 << " " << cov10 << " " << cov11 << std::endl;
  std::cout << std::endl;

  delete fit;
  delete graph;

  return std::make_tuple(slope, offset, cov00, cov01, cov10, cov11);
}

