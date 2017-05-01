///
///  \file   ResolutionAnalysis.C
///
///  \author Tunaface
///          (tuna@cern.ch)
///
///  \date   2017 Jan
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

  if(!b_input) std::cout << "Error at Input: please specify  input file (-i flag)" << std::endl;
  if(!b_out)   std::cout << "Error at Input: please specify output file (-o flag)" << std::endl;
  if(!b_input || !b_out) return 0;

  const int nboards    = 8;
  int       nboardshit = 0;
  int i = 0, ibo = 0, ich = 0, test = 0;
  int jbo = 0;
  int debug = 0;

  // class defs
  PDOToCharge* PDOCalibrator;
  TDOToTime*   TDOCalibrator;
  MMPacmanAlgo*                PACMAN       = new MMPacmanAlgo(5,5.,0.5);
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
  int Nevent = DATA->GetNEntries();

  TH2D* scint_bot_vs_top;

  TH2D* track_hits_v_board;
  TH2D* track_param_x;
  TH2D* track_param_y;
  TH2D* track_angle_7;
  TH2D* track_angle_8;

  TH2D* track_mx_vs_cx_0123_4567;
  TH2D* track_my_vs_cy_0123_4567;

  TH2D* track_mx_vs_cx_0236_1457;
  TH2D* track_my_vs_cy_0236_1457;

  TH2D* track_mx_vs_cx_0256_1347;
  TH2D* track_my_vs_cy_0256_1347;

  TH2D* track_slope_x_scint_19;
  TH2D* track_slope_x_scint_20;
  TH2D* track_slope_x_scint_21;

  TH2D* track_angle_x_scint_19;
  TH2D* track_angle_x_scint_20;
  TH2D* track_angle_x_scint_21;

  TH2D* track_x_vs_y_0123_4567;
  TH2D* track_x_vs_y_0236_1457;
  TH2D* track_x_vs_y_0256_1347;

  TH2D* track_N1_board_vs_residual;
  std::vector<TH2D*> track_N1_slope_x_vs_residual;
  std::vector<TH2D*> track_N1_slope_y_vs_residual;

  TH2D* strip_position_vs_board;
  std::vector<TH2D*> strip_charge_vs_channel;
  std::vector<TH2D*> strip_pdo_vs_channel;
  std::vector<TH2D*> strip_tdo_vs_channel;
  std::vector<TH2D*> strip_tdoc_vs_channel;
  std::vector<TH2D*> strip_dbc_vs_channel;
  std::vector<TH2D*> strip_time_vs_channel;
  std::vector<TH2D*> strip_zpos_vs_channel;
  std::vector<TH1D*> utpc_zdiff;
  TH2D* clus_vs_board;
  TH2D* hits_vs_board;
  TH2D* dups_vs_board;
  TH2D* dups_vs_channel;
  TH2D* hits_per_clus_vs_board;
  TH2D* timediff_vs_board;
  TH1D* timediff_dupli;
  TH1D* timediff_track;
  TH2D* timediff_vs_charge;

  // utpc declarations, matched with Paolo
  double z_track = 0.0;
  double z_utpc  = 0.0;
  double t_utpc  = 0.0;
  double x_utpc  = 0.0;
  double sign    = 0.0;
  double vdrift = 1.0 / 20; // mm per ns
  double deltaT = 850.0;
  std::vector<double> utpc_xs;
  std::vector<double> utpc_zs;
  std::vector<double> utpc_neighbors_xs;
  std::vector<double> utpc_neighbors_zs;
  TMultiGraph* utpc_mg      = 0;
  TGraph*      utpc_graph   = 0;
  TGraph*      utpc_others  = 0;
  TF1*         utpc_fit     = 0;
  TLine*       utpc_line    = 0;
  TH1D*        utpc_dphi    = 0;
  TLatex*      dphi_tex     = 0;
  double       dphi         = 0.0;
  double       m1 = 0, m2 = 0;
  const std::vector<double> zboard = {  1.260,
                                       10.520,
                                       33.710,
                                       42.740,
                                      114.240,
                                      124.260,
                                      147.380,
                                      152.900};

  scint_bot_vs_top = new TH2D("scint_bot_vs_top", ";SC bottom channel;SC top channel;Tracks", 8, -1.5, 6.5, 8, 14.5, 22.5);

  track_angle_x_scint_19 = new TH2D("track_angle_x_scint_19", ";x angle;scint. bottom channel;Tracks, scint. top channel 19", 100, -30, 30, 8, -1.5, 6.5);
  track_angle_x_scint_20 = new TH2D("track_angle_x_scint_20", ";x angle;scint. bottom channel;Tracks, scint. top channel 20", 100, -30, 30, 8, -1.5, 6.5);
  track_angle_x_scint_21 = new TH2D("track_angle_x_scint_21", ";x angle;scint. bottom channel;Tracks, scint. top channel 21", 100, -30, 30, 8, -1.5, 6.5);

  track_slope_x_scint_19 = new TH2D("track_slope_x_scint_19", ";x slope;scint. bottom channel;Tracks, scint. top channel 19", 100, -0.6, 0.6, 8, -1.5, 6.5);
  track_slope_x_scint_20 = new TH2D("track_slope_x_scint_20", ";x slope;scint. bottom channel;Tracks, scint. top channel 20", 100, -0.6, 0.6, 8, -1.5, 6.5);
  track_slope_x_scint_21 = new TH2D("track_slope_x_scint_21", ";x slope;scint. bottom channel;Tracks, scint. top channel 21", 100, -0.6, 0.6, 8, -1.5, 6.5);

  track_hits_v_board = new TH2D("track_hits_v_board", ";Board;N(Clusters);Tracks", 8, -0.5, 7.5, 6, -0.5, 5.5);

  track_param_x = new TH2D("track_param_x", ";x slope; x constant;Tracks",   100, -0.6, 0.6, 100,  -80, 280);
  track_param_y = new TH2D("track_param_y", ";y slope; y constant;Tracks",   100,   -4,   4, 100, -300, 500);
  track_angle_7 = new TH2D("track_angle_7", ";#theta#lower[0.5]{xz} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);
  track_angle_8 = new TH2D("track_angle_8", ";#theta#lower[0.5]{xz} [degrees];#theta#lower[0.5]{yz} [degrees];Tracks", 100, -35, 35, 100, -220, 220);

  track_mx_vs_cx_0123_4567 = new TH2D("track_mx_vs_cx_0123_vs_4567", ";#Delta(x slope);#Delta(x constant);Events with 8 boards", 100, -0.6, 0.6, 400,   -40.0,   40.0);
  track_my_vs_cy_0123_4567 = new TH2D("track_my_vs_cy_0123_vs_4567", ";#Delta(y slope);#Delta(y constant);Events with 8 boards", 100,  -40,  40, 400, -2000.0, 2000.0);

  track_mx_vs_cx_0256_1347 = new TH2D("track_mx_vs_cx_0256_vs_1347", ";#Delta(x slope);#Delta(x constant);Events with 8 boards", 100, -0.08, 0.08, 400,    -8.0,    8.0);
  track_my_vs_cy_0256_1347 = new TH2D("track_my_vs_cy_0256_vs_1347", ";#Delta(y slope);#Delta(y constant);Events with 8 boards", 100,   -10,   10, 400,  -800.0,  800.0);

  track_mx_vs_cx_0236_1457 = new TH2D("track_mx_vs_cx_0236_vs_1457", ";#Delta(x slope);#Delta(x constant);Events with 8 boards", 100, -0.08, 0.08, 400,    -8.0,    8.0);
  track_my_vs_cy_0236_1457 = new TH2D("track_my_vs_cy_0236_vs_1457", ";#Delta(y slope);#Delta(y constant);Events with 8 boards", 100,   -20,   20, 400, -2000.0, 2000.0);

  // track_x_vs_y_0123_4567 = new TH2D("track_x_vs_y_0123_vs_4567", ";#Delta(x @ middle);#Delta(y @ middle);Events with 8 boards", 100, -40, 40, 100, -2000, 2000);
  // track_x_vs_y_0256_1347 = new TH2D("track_x_vs_y_0256_vs_1347", ";#Delta(x @ middle);#Delta(y @ middle);Events with 8 boards", 100,  -5,  5, 100,  -200,  200);
  // track_x_vs_y_0236_1457 = new TH2D("track_x_vs_y_0236_vs_1457", ";#Delta(x @ middle);#Delta(y @ middle);Events with 8 boards", 100,  -5,  5, 100,  -800,  800);

  track_x_vs_y_0123_4567 = new TH2D("track_x_vs_y_0123_vs_4567", ";#Delta(x @ middle);#Delta(y @ middle);Events with 8 boards", 400, -40, 40, 400, -2000, 2000);
  track_x_vs_y_0256_1347 = new TH2D("track_x_vs_y_0256_vs_1347", ";#Delta(x @ middle);#Delta(y @ middle);Events with 8 boards", 400, -40, 40, 400, -2000, 2000);
  track_x_vs_y_0236_1457 = new TH2D("track_x_vs_y_0236_vs_1457", ";#Delta(x @ middle);#Delta(y @ middle);Events with 8 boards", 400, -40, 40, 400, -2000, 2000);

  track_N1_board_vs_residual = new TH2D("track_N1_board_vs_residual",       ";board;x_{cluster} - x_{track, proj.};Tracks", 8, -0.5, 7.5, 200, -5.0, 5.0);
  for (ibo = 0; ibo < nboards; ibo++){
    track_N1_slope_x_vs_residual.push_back(new TH2D(Form("track_N1_slope_x_vs_residual_%i", ibo), ";x slope;x_{cluster} - x_{track, proj.};Tracks", 100, -0.6, 0.6, 200, -5.0, 5.0));
    track_N1_slope_y_vs_residual.push_back(new TH2D(Form("track_N1_slope_y_vs_residual_%i", ibo), ";y slope;x_{cluster} - x_{track, proj.};Tracks", 100, -2.5, 2.5, 200, -5.0, 5.0));
  }

  strip_position_vs_board = new TH2D("strip_position_vs_board", ";strip number;MMFE number;charge [fC]", 512, 0.5, 512.5, 8, -0.5, 7.5);

  for (ibo = 0; ibo < nboards; ibo++){
    strip_charge_vs_channel.push_back(new TH2D(Form("strip_charge_vs_channel_%i", ibo), ";strip number;Charge [fC];strip",    512, 0.5, 512.5, 512,   0,  128));
    strip_pdo_vs_channel.push_back   (new TH2D(Form("strip_pdo_vs_channel_%i",    ibo), ";strip number;PDO [counts];strip",   512, 0.5, 512.5, 512,   0, 2048));
    strip_tdo_vs_channel.push_back   (new TH2D(Form("strip_tdo_vs_channel_%i",    ibo), ";strip number;TDO [counts];strip",   512, 0.5, 512.5, 256,   0,  256));
    strip_tdoc_vs_channel.push_back  (new TH2D(Form("strip_tdoc_vs_channel_%i",   ibo), ";strip number;TDO corr. [ns];strip", 512, 0.5, 512.5, 110, -10,  100));
    strip_dbc_vs_channel.push_back   (new TH2D(Form("strip_dbc_vs_channel_%i",    ibo), ";strip number;#Delta BC;strip",      512, 0.5, 512.5,  64,   0,   64));
    strip_time_vs_channel.push_back  (new TH2D(Form("strip_time_vs_channel_%i",   ibo), ";strip number;Time [ns];strip",      512, 0.5, 512.5, 100, 300, 1000));
    strip_zpos_vs_channel.push_back  (new TH2D(Form("strip_zpos_vs_channel_%i",   ibo), ";strip number;z_{drift} [mm];strip", 512, 0.5, 512.5, 100,  10,   50));

    utpc_zdiff.push_back(new TH1D(Form("utpc_zdiff_%i", ibo), ";z_{#muTPC} - z_{track} [mm];tracks with 8 clusters", 200, -10, 10));

  }
  dups_vs_channel = new TH2D("dups_vs_channel", ";strip number;MMFE number;Duplicates", 512, 0.5, 512.5, 8, -0.5, 7.5);

  clus_vs_board          = new TH2D("clus_vs_board",          ";MMFE number;clusters;Events",            8, -0.5, 7.5, 32, -0.5, 31.5);
  hits_vs_board          = new TH2D("hits_vs_board",          ";MMFE number;strips;Events",              8, -0.5, 7.5, 32, -0.5, 31.5);
  dups_vs_board          = new TH2D("dups_vs_board",          ";MMFE number;duplicate strips;Events",    8, -0.5, 7.5, 32, -0.5, 31.5);
  hits_per_clus_vs_board = new TH2D("hits_per_clus_vs_board", ";MMFE number;hits in a cluster;Clusters", 8, -0.5, 7.5, 13, -0.5, 12.5);
  timediff_vs_board      = new TH2D("timediff_vs_board",      ";MMFE number;#DeltaBCID;Events",          8, -0.5, 7.5, 100, -0.5, 99.5);
  timediff_track         = new TH1D("timediff_track",         ";#DeltaBCID;Strips on-track",             100, -0.5, 99.5);
  timediff_dupli         = new TH1D("timediff_dupli",         ";#DeltaBCID;Duplicate strips",            100, -0.5, 99.5);
  timediff_vs_charge     = new TH2D("timediff_vs_charge",     ";charge [fC];#DeltaBCID;Events",          200, 0, 200, 100, -0.5, 99.5);

  utpc_dphi = new TH1D("utpc_dphi", ";#Delta#phi(clusters track, #muTPC tracklet);", 160, 0.0, 1.6);

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

  // loop
  for(int evt = 0; evt < Nevent; evt++){

    // get the event
    DATA->GetEntry(evt);
    if(evt % (Nevent/20) == 0)
      cout << "Processing event # " << evt << " | " << Nevent << endl;
    if(GEOMETRY->RunNumber() < 0)
      GEOMETRY->SetRunNumber(DATA->RunNum);
    if(!DATA->sc_EventHits.IsGoodEvent())
      continue;

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
    PACMAN->SetEventTrigBCID(DATA->mm_trig_BCID);
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
        ibo = GEOMETRY->Index(DATA->mm_EventHits[i].MMFE8());
        timediff_vs_board->Fill(ibo, DATA->mm_trig_BCID - DATA->mm_EventHits[i][ich].BCID());
        if( !PACMAN->IsGoodHit(DATA->mm_EventHits[i][ich]) )
          continue;
        timediff_vs_charge->Fill(DATA->mm_EventHits[i][ich].Charge(), DATA->mm_trig_BCID - DATA->mm_EventHits[i][ich].BCID());
        strip_charge_vs_channel[ibo]->Fill(DATA->mm_EventHits[i][ich].Channel(), DATA->mm_EventHits[i][ich].Charge());
        strip_pdo_vs_channel[ibo]   ->Fill(DATA->mm_EventHits[i][ich].Channel(), DATA->mm_EventHits[i][ich].PDO());
        strip_tdo_vs_channel[ibo]   ->Fill(DATA->mm_EventHits[i][ich].Channel(), DATA->mm_EventHits[i][ich].TDO());
        strip_dbc_vs_channel[ibo]   ->Fill(DATA->mm_EventHits[i][ich].Channel(), DATA->mm_trig_BCID - DATA->mm_EventHits[i][ich].BCID());
        strip_tdoc_vs_channel[ibo]  ->Fill(DATA->mm_EventHits[i][ich].Channel(), DATA->mm_EventHits[i][ich].Time());
        strip_time_vs_channel[ibo]  ->Fill(DATA->mm_EventHits[i][ich].Channel(),  (DATA->mm_trig_BCID - DATA->mm_EventHits[i][ich].BCID())*25 + DATA->mm_EventHits[i][ich].Time());
        strip_zpos_vs_channel[ibo]  ->Fill(DATA->mm_EventHits[i][ich].Channel(), ((DATA->mm_trig_BCID - DATA->mm_EventHits[i][ich].BCID())*25 + DATA->mm_EventHits[i][ich].Time()) * vdrift);
      }

      // pseudo-event display for clustering visualization
      if (DATA->mm_EventNum == 7821 && board_clusters.size() > 0){
        for(ich = 0; ich < DATA->mm_EventHits[i].GetNHits(); ich++){
          if( !PACMAN->IsGoodHit(DATA->mm_EventHits[i][ich]) )
            continue;
          strip_position_vs_board->Fill(DATA->mm_EventHits[i][ich].Channel(), GEOMETRY->Index(board_clusters[0].MMFE8()), DATA->mm_EventHits[i][ich].Charge());
        }
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
          clus_vs_board->Fill(ibo, clusters_perboard[i].size());
          hits_vs_board->Fill(ibo, DATA->mm_EventHits[i].GetNHits());
          dups_vs_board->Fill(ibo, DATA->mm_EventHits[i].GetNDuplicates());

          duplicates = DATA->mm_EventHits[i].GetDuplicates();
	  //	  cout << "Duplicates? " << duplicates.GetNHits() << endl;
          // for(ich = 0; ich < (int)(duplicates.GetNHits()); ich++){
	  //   cout << "i " << ich << endl;
	  //   cout << "channel " << duplicates[ich].Channel() << endl;
	  //   cout << "nhits " << duplicates[ich].GetNHits() << endl;
          //   // do you want to fill this with a weight of duplicates[ich].GetNHits()-1?
	  //   dups_vs_channel->Fill(duplicates[ich].Channel(), ibo, duplicates[ich].GetNHits()-1);
	  //   timediff_dupli->Fill(DATA->mm_trig_BCID - duplicates[ich].BCID());
          // }
        }
      }

      // no hits on this board!
      if (test == -1){
        clus_vs_board->Fill(ibo, 0);
        hits_vs_board->Fill(ibo, 0);
        dups_vs_board->Fill(ibo, 0);
      }
    }

    
    for (auto clus_list: clusters_perboard)
      for (auto clus: clus_list)
        hits_per_clus_vs_board->Fill(GEOMETRY->Index(clus->MMFE8()), clus->GetNHits());

    
    // preselection quality to run tracking
    // require at least 7 boards hit
    if (clusters_perboard.size() < 7)
      continue;
    //    cout << "fucked" << endl;
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
    if (clusters_road.size() == 7) track_angle_7->Fill(atan(track.SlopeX())*180/3.14159, atan(track.SlopeY())*180/3.14159);
    if (clusters_road.size() == 8) track_angle_8->Fill(atan(track.SlopeX())*180/3.14159, atan(track.SlopeY())*180/3.14159);

    // timing of on-track hits
    if (clusters_road.size() >= 7)
      for (auto clus: clusters_road)
        for (i = 0; i < (int)(clus->GetNHits()); i++)
          timediff_track->Fill(DATA->mm_trig_BCID - clus->Get(i).BCID());

    // require a very good track, i.e., 8 boards hit
    if (clusters_road.size() != 8)
      continue;

    // some sanity plots
    track_param_x->Fill(track.SlopeX(), track.ConstX());
    track_param_y->Fill(track.SlopeY(), track.ConstY());

    // interlude: some scintillator plots
    for (auto toppair: DATA->sc_EventHits.GetTopPair()){
      for (auto botpair: DATA->sc_EventHits.GetBotPair()){
        scint_bot_vs_top->Fill(botpair.first->Channel(), toppair.first->Channel());
        if (toppair.first->Channel() == 19){
          track_angle_x_scint_19->Fill(atan(track.SlopeX())*180/3.14159, botpair.first->Channel());
          track_slope_x_scint_19->Fill(     track.SlopeX(),              botpair.first->Channel());
        }
        else if (toppair.first->Channel() == 20){
          track_angle_x_scint_20->Fill(atan(track.SlopeX())*180/3.14159, botpair.first->Channel());
          track_slope_x_scint_20->Fill(     track.SlopeX(),              botpair.first->Channel());
        }
        else if (toppair.first->Channel() == 21){
          track_angle_x_scint_21->Fill(atan(track.SlopeX())*180/3.14159, botpair.first->Channel());
          track_slope_x_scint_21->Fill(     track.SlopeX(),              botpair.first->Channel());
        }
      }
    }

    if (debug || DATA->mm_EventNum == 14036 || DATA->mm_EventNum == 46873 || DATA->mm_EventNum == 41380){
      can = Plot_Track2D(Form("track2D_%05d_road", DATA->mm_EventNum), track, *GEOMETRY, &clusters_road);
      fout->cd("event_displays");
      can->Write();
      delete can;
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

      track_N1_board_vs_residual       ->Fill(ibo,            residual);
      track_N1_slope_x_vs_residual[ibo]->Fill(track.SlopeX(), residual);
      track_N1_slope_y_vs_residual[ibo]->Fill(track.SlopeY(), residual);

      if (DATA->mm_EventNum == 40 && false){
        plane = GEOMETRY->Get(GEOMETRY->Index(clus->MMFE8()));
        std::cout << Form("Event %2i | Board %2i | x = %7.4f z = %7.3f | residual: %7.4f | track w/out this cluster: mx = %7.4f cx = %7.4f my = %7.4f cy = %7.4f",
                          DATA->mm_EventNum, ibo,
                          plane.Origin().X() + plane.LocalXatYbegin(clus->Channel()), plane.Origin().Z(),
                          residual,
                          track_N1.SlopeX(), track_N1.ConstX(), track_N1.SlopeY(), track_N1.ConstY()) << std::endl;
      }
    }


    // 4-hit track vs. 4-hit track histograms
    // ------------------------------------------------------------------
    track_0123.Reset(); clus_0123.Reset();
    track_4567.Reset(); clus_4567.Reset();
    track_0256.Reset(); clus_0256.Reset();
    track_1347.Reset(); clus_1347.Reset();
    track_0236.Reset(); clus_0236.Reset();
    track_1457.Reset(); clus_1457.Reset();
    
    for (auto clus: clusters_road){
      ibo = GEOMETRY->Index(clus->MMFE8());
      if (ibo == 0 || ibo == 1 || ibo == 2 || ibo == 3) clus_0123.AddCluster(*clus);
      if (ibo == 4 || ibo == 5 || ibo == 6 || ibo == 7) clus_4567.AddCluster(*clus);
      if (ibo == 0 || ibo == 2 || ibo == 5 || ibo == 6) clus_0256.AddCluster(*clus);
      if (ibo == 1 || ibo == 3 || ibo == 4 || ibo == 7) clus_1347.AddCluster(*clus);
      if (ibo == 0 || ibo == 2 || ibo == 3 || ibo == 6) clus_0236.AddCluster(*clus);
      if (ibo == 1 || ibo == 4 || ibo == 5 || ibo == 7) clus_1457.AddCluster(*clus);
    }

    track_0123 = FITTER->Fit(clus_0123, *GEOMETRY, DATA->mm_EventNum);
    track_4567 = FITTER->Fit(clus_4567, *GEOMETRY, DATA->mm_EventNum);
    track_0256 = FITTER->Fit(clus_0256, *GEOMETRY, DATA->mm_EventNum);
    track_1347 = FITTER->Fit(clus_1347, *GEOMETRY, DATA->mm_EventNum);
    track_0236 = FITTER->Fit(clus_0236, *GEOMETRY, DATA->mm_EventNum);
    track_1457 = FITTER->Fit(clus_1457, *GEOMETRY, DATA->mm_EventNum);
    
    track_mx_vs_cx_0123_4567->Fill(track_0123.SlopeX() - track_4567.SlopeX(), track_0123.ConstX() - track_4567.ConstX());
    track_my_vs_cy_0123_4567->Fill(track_0123.SlopeY() - track_4567.SlopeY(), track_0123.ConstY() - track_4567.ConstY());

    track_mx_vs_cx_0256_1347->Fill(track_0256.SlopeX() - track_1347.SlopeX(), track_0256.ConstX() - track_1347.ConstX());
    track_my_vs_cy_0256_1347->Fill(track_0256.SlopeY() - track_1347.SlopeY(), track_0256.ConstY() - track_1347.ConstY());

    track_mx_vs_cx_0236_1457->Fill(track_0236.SlopeX() - track_1457.SlopeX(), track_0236.ConstX() - track_1457.ConstX());
    track_my_vs_cy_0236_1457->Fill(track_0236.SlopeY() - track_1457.SlopeY(), track_0236.ConstY() - track_1457.ConstY());    

    track_x_vs_y_0123_4567->Fill(track_0123.PointZ(z_middle).X() - track_4567.PointZ(z_middle).X(), 
                                 track_0123.PointZ(z_middle).Y() - track_4567.PointZ(z_middle).Y());
    track_x_vs_y_0256_1347->Fill(track_0256.PointZ(z_middle).X() - track_1347.PointZ(z_middle).X(), 
                                 track_0256.PointZ(z_middle).Y() - track_1347.PointZ(z_middle).Y());
    track_x_vs_y_0236_1457->Fill(track_0236.PointZ(z_middle).X() - track_1457.PointZ(z_middle).X(), 
                                 track_0236.PointZ(z_middle).Y() - track_1457.PointZ(z_middle).Y());


    // MICRO TPC BABY
    // -----------------------------------------------------------------------------
    track = FITTER->Fit(clusters_road, *GEOMETRY, DATA->mm_EventNum);

    // only angled tracks for now
    if (std::fabs(atan(track.SlopeX())*180/3.14159) < 15)
      continue;

    for (auto clus: clusters_road){

      ibo = GEOMETRY->Index(clus->MMFE8());
      
      if (ibo == 2 || ibo == 3 || ibo == 4 || ibo == 5) continue;
      sign = (ibo==0 || ibo==2 || ibo==4 || ibo==6) ? -1.0 : 1.0;
      
      utpc_xs.clear();
      utpc_zs.clear();

      plane = GEOMETRY->Get(ibo);

      for (ich = 0; ich < clus->GetNHits(); ich++){

        // only "good" hits for now
        if (abs(DATA->mm_trig_BCID - clus->Get(ich).BCID()) >= 35 ||
            abs(DATA->mm_trig_BCID - clus->Get(ich).BCID()) <= 26 ||
            clus->Get(ich).Charge() > 70 ||
            clus->Get(ich).Charge() < 2
            )
          continue;

        x_utpc  = plane.Origin().X() + plane.LocalXatYbegin(clus->Get(ich).Channel());
        t_utpc  = deltaT - ((DATA->mm_trig_BCID - clus->Get(ich).BCID())*25 + clus->Get(ich).Time());
        z_utpc  = zboard[ibo] + vdrift*t_utpc*sign;
        z_track = (x_utpc - track.ConstX()) / track.SlopeX();
        utpc_zdiff[ibo]->Fill(z_utpc - z_track);

        utpc_xs.push_back(x_utpc);
        utpc_zs.push_back(z_utpc);

      }

      if ((int)(utpc_xs.size()) >= 3){
        can = new TCanvas(Form("utpc_%06i_Board%i", DATA->mm_EventNum, ibo), Form("utpc_%06i_Board%i", DATA->mm_EventNum, ibo), 800, 800);
        utpc_graph = new TGraph(int(utpc_xs.size()), &utpc_xs[0], &utpc_zs[0]);
        utpc_graph->GetXaxis()->SetTitle("x_{#muTPC} [mm]");
        utpc_graph->GetYaxis()->SetTitle("z_{#muTPC} [mm]");
        utpc_graph->SetName(Form("utpc_gr_%06i_Board%i", DATA->mm_EventNum, ibo));
        utpc_graph->SetLineStyle(0);
        utpc_graph->GetHistogram()->SetLineStyle(0);
        fout->cd("event_displays");

        utpc_graph->Fit("pol1", "Q");
        utpc_fit = utpc_graph->GetFunction("pol1");
        utpc_fit->SetLineWidth(2);
        utpc_fit->SetLineColor(kRed);

        // draw other boards for reference
        utpc_neighbors_xs.clear();
        utpc_neighbors_zs.clear();
        for (auto clus_other: clusters_road){
          jbo = GEOMETRY->Index(clus_other->MMFE8());
          if (jbo == 2 || jbo == 3 || jbo == 4 || jbo == 5)
            continue;
          utpc_neighbors_xs.push_back(GEOMETRY->Get(jbo).Origin().X() + GEOMETRY->Get(jbo).LocalXatYend(clus_other->Channel()));
          utpc_neighbors_zs.push_back(GEOMETRY->Get(jbo).Origin().Z());
        }
        utpc_others = new TGraph(int(utpc_neighbors_xs.size()), &utpc_neighbors_xs[0], &utpc_neighbors_zs[0]);
        utpc_others->SetMarkerColor(210);

        can->Draw();
        can->SetFillColor(10);
        gStyle->SetOptFit(0000);

        utpc_mg = new TMultiGraph();
        utpc_mg->Add(utpc_graph);
        utpc_mg->Add(utpc_others);
        utpc_mg->Draw("ap");
        utpc_mg->GetXaxis()->SetTitle("x [mm]");
        utpc_mg->GetYaxis()->SetTitle("z [mm]");

        double zmin = *std::min_element(utpc_neighbors_zs.begin(), utpc_neighbors_zs.end());
        double zmax = *std::max_element(utpc_neighbors_zs.begin(), utpc_neighbors_zs.end());
        double xmin = track.ConstX() + zmin*track.SlopeX();
        double xmax = track.ConstX() + zmax*track.SlopeX();
        utpc_line = new TLine(xmin, zmin, xmax, zmax);
        utpc_line->SetLineColor(210);
        utpc_line->Draw();

        // calculating the angle from two slopes!!
        // this is also used in the trigger processor baby!!
        m1 = track.SlopeX();
        m2 = 1 / utpc_fit->GetParameter(1);
        dphi = (m1-m2) / (1 + m1*m2);
        dphi = atan(dphi);
        dphi = std::fabs(dphi);
        dphi_tex = new TLatex(0.2, 0.93, Form("#Delta#phi = %5.2f", dphi));
        dphi_tex->SetNDC();
        dphi_tex->Draw();
        utpc_dphi->Fill(dphi);

        // can->Write();
        delete can;
        delete utpc_mg;
        delete utpc_line;
        delete dphi_tex;
      }
    }
 
  }


  // write to file
  fout->cd();
  fout->mkdir("histograms");
  fout->cd("histograms");

  scint_bot_vs_top->Write();

  track_param_x->Write();
  track_param_y->Write();
  track_angle_7->Write();
  track_angle_8->Write();

  track_mx_vs_cx_0123_4567->Write();
  track_my_vs_cy_0123_4567->Write();

  track_mx_vs_cx_0256_1347->Write();
  track_my_vs_cy_0256_1347->Write();

  track_mx_vs_cx_0236_1457->Write();
  track_my_vs_cy_0236_1457->Write();

  track_x_vs_y_0123_4567->Write();
  track_x_vs_y_0256_1347->Write();
  track_x_vs_y_0236_1457->Write();

  track_slope_x_scint_19->Write();
  track_slope_x_scint_20->Write();
  track_slope_x_scint_21->Write();

  track_angle_x_scint_19->Write();
  track_angle_x_scint_20->Write();
  track_angle_x_scint_21->Write();

  for (auto hist: track_N1_slope_x_vs_residual)
    hist->Write();
  for (auto hist: track_N1_slope_y_vs_residual)
    hist->Write();
  track_N1_board_vs_residual->Write();

  strip_position_vs_board->Write();
  for (auto hist: strip_charge_vs_channel) hist->Write();
  for (auto hist: strip_pdo_vs_channel)    hist->Write();
  for (auto hist: strip_tdo_vs_channel)    hist->Write();
  for (auto hist: strip_tdoc_vs_channel)   hist->Write();
  for (auto hist: strip_dbc_vs_channel)    hist->Write();
  for (auto hist: strip_time_vs_channel)   hist->Write();
  for (auto hist: strip_zpos_vs_channel)   hist->Write();
  for (auto hist: utpc_zdiff)              hist->Write();
    
  clus_vs_board->Write();
  hits_vs_board->Write();
  dups_vs_board->Write();
  dups_vs_channel->Write();
  hits_per_clus_vs_board->Write();
  timediff_vs_board->Write();
  timediff_dupli->Write();
  timediff_track->Write();
  timediff_vs_charge->Write();

  utpc_dphi->Write();

  fout->Close();    
}
