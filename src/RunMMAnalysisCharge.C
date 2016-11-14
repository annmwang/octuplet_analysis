///
///  \file   RunMMAnalysis.C
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///

#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <stdlib.h>

#include "include/MMPlot.hh"
#include "include/PDOToCharge.hh"
#include "include/TDOToTime.hh"
#include "include/MMDataAnalysis.hh"
#include "include/MMPacmanAlgo.hh"
#include "include/GeoOctuplet.hh"
#include "include/SimpleTrackFitter.hh"
#include "include/ChargeTrackFitter.hh"
#include "include/MMEfficiency.hh"

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
  bool b_align   = false;

  // pacman options
  double pac_csize = 5.;
  double pac_sthresh = 15.;
  double pac_hthresh = 10.;
  double fit_Charge2Res = 1;
  double fit_Boards2Res = 1;

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

    if (strncmp(argv[i],"--PacmanCsize",13)==0) {
      pac_csize = stod(argv[i+1]);
    }
    if (strncmp(argv[i],"--PacmanSthresh",15)==0) {
      pac_sthresh = stod(argv[i+1]);
    }
    if (strncmp(argv[i],"--PacmanHthresh",15)==0) {
      pac_hthresh = stod(argv[i+1]);
    }
    if (strncmp(argv[i],"--Charge2Res",12)==0) {
      fit_Charge2Res = stod(argv[i+1]);
    }
    if (strncmp(argv[i],"--Boards2Res",12)==0) {
      fit_Boards2Res = stod(argv[i+1]);
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

  MMPacmanAlgo* PACMAN = new MMPacmanAlgo(pac_csize, pac_sthresh, pac_hthresh);

  GeoOctuplet* GEOMETRY = new GeoOctuplet();
  if(b_align)
    GEOMETRY->SetAlignment(AlignFileName);

  ChargeTrackFitter *FITTER = new ChargeTrackFitter(fit_Charge2Res, fit_Boards2Res);
  SimpleTrackFitter *SIMPLEFITTER = new SimpleTrackFitter();

  //MMEfficiency *EFFICIENCY = new MMEfficiency(GEOMETRY, FITTER);

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

  // event num vs. stuff
  vector<TH2D*> board_PDO_v_EVT;
  vector<TH2D*> board_Q_v_EVT;
  vector<TH2D*> board_NHit_v_EVT;
  vector<TH2D*> board_Ndup_v_EVT;
  vector<TH2D*> board_Nclus_v_EVT;
  vector<TH2D*> board_Nclusdup_v_EVT;

  // event object counting
  vector<TH1D*> board_NHit;
  vector<TH1D*> board_Ndup;
  vector<TH1D*> board_Nclus;
  vector<TH1D*> board_Nclusdup;

  // histograms for hit analysis
  vector<TH1D*> board_hit_PDO;
  vector<TH1D*> board_hit_CH;
  vector<TH1D*> board_hit_Q;
  vector<TH2D*> board_hit_PDO_v_CH;
  vector<TH2D*> board_hit_Q_v_CH;

  vector<TH1D*> board_duphit_PDO;
  vector<TH1D*> board_duphit_NCH;
  vector<TH1D*> board_duphit_CH;
  vector<TH1D*> board_duphit_Q;
  vector<TH2D*> board_duphit_PDO_v_CH;
  vector<TH2D*> board_duphit_Q_v_CH;

  vector<TH1D*> board_duphit_dPDO;
  vector<TH1D*> board_duphit_dTDO;
  vector<TH1D*> board_duphit_dBCID; 
  
  // histograms for cluster analysis
  vector<TH1D*> board_clus_CH;
  vector<TH1D*> board_clus_Q;
  vector<TH1D*> board_clusN_CH;
  vector<TH1D*> board_clusN_Q;
  vector<TH1D*> board_clus_NHit;
  vector<TH1D*> board_clus_Ndup;
  vector<TH1D*> board_clusN_Ndup;
  vector<TH2D*> board_clus_Q_v_CH;
  vector<TH2D*> board_clusN_Q_v_CH;

  vector<TH1D*> board_clusdup_CH;
  vector<TH1D*> board_clusdup_Q;
  vector<TH1D*> board_clusdup_NHit;
  vector<TH2D*> board_clusdup_Q_v_CH;

  // histograms for track analysis
  vector<TH1D*> board_itrack_resX;
  vector<TH1D*> board_otrack_resX;
  vector<TH2D*> board_itrack_resX_v_CH;
  vector<TH2D*> board_otrack_resX_v_CH;
  vector<TH2D*> board_itrack_resX_v_EVT;
  vector<TH2D*> board_otrack_resX_v_EVT;

  vector<TH1D*> track_sumresX2;
  TH2D*         track_sumresX2_v_Nclus;

  for(int i = 0; i < 8; i++){
    board_PDO_v_EVT.push_back(new TH2D(Form("b_PDOvEVT_%d",i),
				       Form("b_PDOvEVT_%d",i),
				       1024,0.,max_EventNum,
				       128,0.,1028.));
    board_Q_v_EVT.push_back(new TH2D(Form("b_QvEVT_%d",i),
				       Form("b_QvEVT_%d",i),
				       1024,0.,max_EventNum,
				       128,0.,150.));

    board_NHit_v_EVT.push_back(new TH2D(Form("b_NHitvEVT_%d",i),
				       Form("b_NHitvEVT_%d",i),
				       1024,0.,max_EventNum,
				       40,0.5,40.5));
    board_Ndup_v_EVT.push_back(new TH2D(Form("b_NdupvEVT_%d",i),
				       Form("b_NdupvEVT_%d",i),
				       1024,0.,max_EventNum,
				       20,-0.5,19.5));
    board_Nclus_v_EVT.push_back(new TH2D(Form("b_NclusvEVT_%d",i),
				       Form("b_NclusvEVT_%d",i),
				       1024,0.,max_EventNum,
				       7,-0.5,6.5));
    board_Nclusdup_v_EVT.push_back(new TH2D(Form("b_NclusdupvEVT_%d",i),
				       Form("b_NclusdupvEVT_%d",i),
				       1024,0.,max_EventNum,
				       5,-0.5,4.5));

    board_hit_PDO.push_back(new TH1D(Form("b_h_PDO_%d",i),
				     Form("b_h_PDO_%d",i),
				     1024,0.,1028.));
    board_hit_CH.push_back(new TH1D(Form("b_h_CH_%d",i),
				    Form("b_h_CH_%d",i),
				    512,0.5,512.5));
    board_hit_Q.push_back(new TH1D(Form("b_h_Q_%d",i),
				   Form("b_h_Q_%d",i),
				   128,0.0,128));
    board_hit_PDO_v_CH.push_back(new TH2D(Form("b_h_PDOvCH_%d",i),
					  Form("b_h_PDOvCH_%d",i),
					  512,0.5,512.5,
					  128,0.,1028.));
    board_hit_Q_v_CH.push_back(new TH2D(Form("b_h_QvCH_%d",i),
					Form("b_h_QvCH_%d",i),
					512,0.5,512.5,
					128,0.,128.));
    
    board_duphit_PDO.push_back(new TH1D(Form("b_dh_PDO_%d",i),
					Form("b_dh_PDO_%d",i),
					1024,.0,1028));;
    board_duphit_NCH.push_back(new TH1D(Form("b_dh_NCH_%d",i),
					Form("b_dh_NCH_%d",i),
					512,0.5,512.5));;
    board_duphit_CH.push_back(new TH1D(Form("b_dh_CH_%d",i),
				       Form("b_dh_CH_%d",i),
				       512,0.5,512.5));;
    board_duphit_Q.push_back(new TH1D(Form("b_dh_Q_%d",i),
				      Form("b_dh_Q_%d",i),
				      128,0.0,128));;
    board_duphit_PDO_v_CH.push_back(new TH2D(Form("b_dh_PDOvCH_%d",i),
					     Form("b_dh_PDOvCH_%d",i),
					     512,0.5,512.5,
					     128,0.,1028.));
    board_duphit_Q_v_CH.push_back(new TH2D(Form("b_dh_QvCH_%d",i),
					   Form("b_dh_QvCH_%d",i),
					   512,0.5,512.5,
					   128,0.,128.));

    board_NHit.push_back(new TH1D(Form("b_NHit_%d",i),
				  Form("b_NHit_%d",i),
				  40,0.5,40.5));
    board_Ndup.push_back(new TH1D(Form("b_Ndup_%d",i),
				  Form("b_Ndup_%d",i),
				  20,-0.5,19.5));
    board_Nclus.push_back(new TH1D(Form("b_Nclus_%d",i),
				   Form("b_Nclus_%d",i),
				   7,-0.5,6.5));

    board_Nclusdup.push_back(new TH1D(Form("b_Nclusdup_%d",i),
				      Form("b_Nclusdup_%d",i),
				      5,-0.5,4.5));
  
    board_clus_CH.push_back(new TH1D(Form("b_c_CH_%d",i),
				     Form("b_c_CH_%d",i),
				     512,0.5,512.5));
    board_clus_Q.push_back(new TH1D(Form("b_c_Q_%d",i),
				    Form("b_c_Q_%d",i),
				    128,0.0,250));
    board_clusN_CH.push_back(new TH1D(Form("b_cN_CH_%d",i),
				      Form("b_cN_CH_%d",i),
				      512,0.5,512.5));
    board_clusN_Q.push_back(new TH1D(Form("b_cN_Q_%d",i),
				     Form("b_cN_Q_%d",i),
				     128,0.0,250));
    board_clus_NHit.push_back(new TH1D(Form("b_c_NHit_%d",i),
				       Form("b_c_NHit_%d",i),
				       15,0.5,15.5));
    board_clus_Ndup.push_back(new TH1D(Form("b_c_Ndup_%d",i),
				       Form("b_c_Ndup_%d",i),
				       6,-0.5,5.5));
    board_clusN_Ndup.push_back(new TH1D(Form("b_cN_Ndup_%d",i),
					Form("b_cN_Ndup_%d",i),
					6,-0.5,5.5));
    board_clus_Q_v_CH.push_back(new TH2D(Form("b_c_QvCH_%d",i),
					 Form("b_c_QvCH_%d",i),
					 512,0.5,512.5,
					 128,0.,250.));
    board_clusN_Q_v_CH.push_back(new TH2D(Form("b_cN_QvCH_%d",i),
					  Form("b_cN_QvCH_%d",i),
					  512,0.5,512.5,
					  128,0.,250.));

    board_clusdup_CH.push_back(new TH1D(Form("b_cd_CH_%d",i),
					Form("b_cd_CH_%d",i),
					512,0.5,512.5));
    board_clusdup_Q.push_back(new TH1D(Form("b_cd_Q_%d",i),
				       Form("b_cd_Q_%d",i),
				       128,0.0,250));
    board_clusdup_NHit.push_back(new TH1D(Form("b_cd_NHit_%d",i),
					  Form("b_cd_NHit_%d",i),
					  15,0.5,15.5));
    board_clusdup_Q_v_CH.push_back(new TH2D(Form("b_cd_QvCH_%d",i),
					    Form("b_cd_QvCH_%d",i),
					    512,0.5,512.5,
					    128,0.,250.));

    board_itrack_resX.push_back(new TH1D(Form("b_it_resX_%d",i),
					 Form("b_it_resX_%d",i),
					 1001,-5.,5.));
    board_otrack_resX.push_back(new TH1D(Form("b_ot_resX_%d",i),
					 Form("b_ot_resX_%d",i),
					 1001,-5.,5.));
    board_itrack_resX_v_CH.push_back(new TH2D(Form("b_it_resXvCH_%d",i),
					      Form("b_it_resXvCH_%d",i),
					      512,0.5,512.5,
					      201,-5.,5.));
    board_otrack_resX_v_CH.push_back(new TH2D(Form("b_ot_resXvCH_%d",i),
					      Form("b_ot_resXvCH_%d",i),
					      512,0.5,512.5,
					      201,-5.,5.));
    board_itrack_resX_v_EVT.push_back(new TH2D(Form("b_it_resXvEVT_%d",i),
					       Form("b_it_resXvEVT_%d",i),
					       1024,0.,max_EventNum,
					       201,-5.,5.));
    board_otrack_resX_v_EVT.push_back(new TH2D(Form("b_ot_resXvEVT_%d",i),
					       Form("b_ot_resXvEVT_%d",i),
					       1024,0.,max_EventNum,
					       201,-5.,5.));
  }
  
  for(int i = 5; i <= 8; i++){
    track_sumresX2.push_back(new TH1D(Form("t_sumresX2_%d",i),
				      Form("t_sumresX2_%d",i),
				      1024,0.,5.));
  }
  track_sumresX2_v_Nclus = new TH2D("t_sumresX2_v_Nclus",
				    "t_sumresX2_v_Nclus",
				    4, 4.5, 8.5,
				    1024,0.,5.);

  TFile* fout = new TFile(outputFileName, "RECREATE");
  fout->mkdir("event_displays");
  fout->mkdir("cluster_displays");
  MMPlot();

  for(int evt = 0; evt < Nevent; evt++){
    DATA->GetEntry(evt);
    if(evt%(Nevent/10) == 0) 
      cerr << "Processing event # " << evt << " | " << Nevent << endl;

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

    // book histograms for MM hits
    int Nboard = DATA->mm_EventHits.GetNBoards();
    for(int i = 0; i < Nboard; i++){
      int b = ib[DATA->mm_EventHits[i].MMFE8()];
      int Nhit = DATA->mm_EventHits[i].GetNHits();
      board_NHit[b]->Fill(Nhit);
      board_Ndup[b]->Fill(DATA->mm_EventHits[i].GetNDuplicates());
      board_Ndup_v_EVT[b]->Fill(DATA->mm_EventNum, DATA->mm_EventHits[i].GetNDuplicates());
      board_NHit_v_EVT[b]->Fill(DATA->mm_EventNum, Nhit);
      for(int j = 0; j < Nhit; j++){
	const MMLinkedHit& hit = DATA->mm_EventHits[i][j];
	board_hit_PDO[b]->Fill(hit.PDO());
	board_hit_CH[b]->Fill(hit.Channel());
	board_hit_Q[b]->Fill(hit.Charge());
	board_hit_PDO_v_CH[b]->Fill(hit.Channel(),hit.PDO());
	board_hit_Q_v_CH[b]->Fill(hit.Channel(),hit.Charge());
	  
	board_PDO_v_EVT[b]->Fill(DATA->mm_EventNum, hit.PDO());
	board_Q_v_EVT[b]->Fill(DATA->mm_EventNum, hit.Charge());

	// hit has duplicate
	int Ndup = hit.GetNHits();
	if(Ndup > 1){
	  board_duphit_CH[b]->Fill(hit.Channel());
	  const MMLinkedHit* phit = &hit;
	  while(phit){
	    board_duphit_PDO[b]->Fill(phit->PDO());
	    board_duphit_NCH[b]->Fill(phit->Channel());
	    board_duphit_Q[b]->Fill(phit->Charge());
	    board_duphit_PDO_v_CH[b]->Fill(phit->Channel(),phit->PDO());
	    board_duphit_Q_v_CH[b]->Fill(phit->Channel(), phit->Charge());
	    phit = phit->GetNext();
	  }
	}
      }
    }
    
    vector<MMClusterList> all_clusters;
    for(int i = 0; i < Nboard; i++){
      if(DATA->mm_EventHits[i].GetNHits() == 0)
        continue;
      
      MMClusterList clusters = PACMAN->Cluster(DATA->mm_EventHits[i]);
      board_Nclus[ib[DATA->mm_EventHits[i].MMFE8()]]->Fill(clusters.GetNCluster());
      board_Nclusdup[ib[DATA->mm_EventHits[i].MMFE8()]]->Fill(clusters.GetNDuplicates());
      board_Nclus_v_EVT[ib[DATA->mm_EventHits[i].MMFE8()]]->Fill(DATA->mm_EventNum, clusters.GetNCluster());
      board_Nclusdup_v_EVT[ib[DATA->mm_EventHits[i].MMFE8()]]->Fill(DATA->mm_EventNum, clusters.GetNDuplicates());
      if(clusters.GetNCluster() > 0)
        all_clusters.push_back(clusters);
    }
   
    // track fitting: use all the tracks 
    int Ncl = all_clusters.size();
    MMClusterList fit_clusters;
    double total_charge = 0;
    for(int i = 0; i < Ncl; i++){
      int Nc = all_clusters[i].GetNCluster();
      for(int j = 0; j < Nc; j++){
        const MMCluster& clus = all_clusters[i][j];

        //if (clus.GetNHits() > 1) {
            total_charge += clus.Charge();
            fit_clusters.AddCluster( clus );
        //}

        int b = ib[clus.MMFE8()];
        board_clus_CH[b]->Fill(clus.Channel());
        board_clus_Q[b]->Fill(clus.Charge());
        board_clus_NHit[b]->Fill(clus.GetNHits());
        board_clus_Ndup[b]->Fill(clus.GetNDuplicates());
        board_clus_Q_v_CH[b]->Fill(clus.Channel(),clus.Charge());
        if(clus.GetNHits() > 1){
          board_clusN_CH[b]->Fill(clus.Channel());
          board_clusN_Q[b]->Fill(clus.Charge());
          board_clusN_Ndup[b]->Fill(clus.GetNDuplicates());
          board_clusN_Q_v_CH[b]->Fill(clus.Channel(),clus.Charge());
        }
        if(clus.GetNDuplicates() > 0){
          board_clusdup_CH[b]->Fill(clus.Channel());
          board_clusdup_Q[b]->Fill(clus.Charge());
          board_clusdup_NHit[b]->Fill(clus.GetNHits());
          board_clusdup_Q_v_CH[b]->Fill(clus.Channel(),clus.Charge());
        }
      }
    }

    // see if we have enough tracks
    bool isinteresting = false;
    int Nboards_hit = 0;
    for(int i = 0; i < Ncl; i++){
      // add highest charge cluster from each board;
      int Nc = all_clusters[i].GetNCluster();
      Nboards_hit += (Nc > 0);
      if (Nc > 1)
          isinteresting = true;
    }

    if(Nboards_hit < 6)
      continue;

    bool loud = DATA->mm_EventNum == 1422;
    MMTrack track_all = FITTER->Fit(fit_clusters, *GEOMETRY, total_charge, loud);
    //cout << track_all.SlopeX() << " "<< track_all.SlopeY() << " "<< track_all.ConstX() << " " << track_all.ConstY() << endl;
    const MMClusterList track_clusters = track_all.Clusters();
    int Nclus_fit = track_clusters.GetNCluster();

    double sumresX2 = 0.;
    for(int c = 0; c < Nclus_fit; c++){
      const MMCluster& clus = track_clusters.Get(c);
      int b = ib[clus.MMFE8()];
      // fill on-track residuals
      double resX = GEOMETRY->GetResidualX(clus, track_all);
      sumresX2 += resX*resX;
    }

    track_sumresX2[Nclus_fit-5]->Fill(sumresX2);
    track_sumresX2_v_Nclus->Fill(Nclus_fit, sumresX2);
    
    //if( sumresX2 > double(Nclus_fit-2)*0.1)
    //  continue;

    for(int c = 0; c < Nclus_fit; c++){
      const MMCluster& clus = track_clusters.Get(c);
      int b = ib[clus.MMFE8()];
      // fill on-track residuals
      double resX = GEOMETRY->GetResidualX(clus, track_all);
      board_itrack_resX[b]->Fill(resX);
      board_itrack_resX_v_CH[b]->Fill(clus.Channel(), resX);
      board_itrack_resX_v_EVT[b]->Fill(DATA->mm_EventNum, resX);

      // new cluster list without this cluster
      MMClusterList clus_list;
      for(int o = 0; o < Nclus_fit; o++)
	if(o != c)
	  clus_list.AddCluster(track_clusters.Get(o));
	
      MMTrack track = SIMPLEFITTER->Fit(clus_list, *GEOMETRY);

      resX = GEOMETRY->GetResidualX(clus, track);
      board_otrack_resX[b]->Fill(resX);
      board_otrack_resX_v_CH[b]->Fill(clus.Channel(), resX);
      board_otrack_resX_v_EVT[b]->Fill(DATA->mm_EventNum, resX);
    }

    // Get the Effenciency Values
    //EFFICIENCY->UpdateEfficiency(fit_clusters);

    if (evt > 2000)
        break;

    //if(Nclus_fit < 8)
    //  continue;

    if (isinteresting) {
      vector<MMFE8Hits> boardhits;
      int Nboards = DATA->mm_EventHits.GetNBoards();
      for (int i = 0; i < Nboards; i++)
          boardhits.push_back( DATA->mm_EventHits[i] );

      TCanvas* can = Plot_Cluster2D(Form("clusters2D_%d",DATA->mm_EventNum), *GEOMETRY, all_clusters, track_all, boardhits, &track_clusters); 
      fout->cd("cluster_displays");
      can->Write();
      delete can;
    }

    {
      TCanvas* can = Plot_Track2D(Form("track2D_%d",DATA->mm_EventNum), track_all, *GEOMETRY, &track_clusters); 
      fout->cd("event_displays");
      can->Write();
      delete can;
    }

    TCanvas* canY = Plot_Track2DY(Form("track2DY_%d",DATA->mm_EventNum), track_all, *GEOMETRY, &track_clusters); 
    fout->cd("event_displays");
    canY->Write();
    delete canY;
    
    TCanvas* can3D = Plot_Track3D(Form("track3D_%d",DATA->mm_EventNum), track_all, *GEOMETRY, &track_clusters); 
    fout->cd("event_displays");
    can3D->Write();
    delete can3D;


  }

  // print the Efficiency 
  /*
  for (int ib = 0; ib < GEOMETRY->GetNPlanes(); ib++) {
    cout << "Board " << GEOMETRY->MMFE8(ib) << ": " << EFFICIENCY->GetBoardEfficiency(ib) << endl;
  }
  */

  fout->cd();
  fout->mkdir("histograms");
  for(int i = 0; i < 8; i++){
    fout->cd("histograms");
    board_PDO_v_EVT[i]->Write();
    board_Q_v_EVT[i]->Write();
    board_NHit_v_EVT[i]->Write();
    board_Ndup_v_EVT[i]->Write();
    board_Nclus_v_EVT[i]->Write();
    board_Nclusdup_v_EVT[i]->Write();

    board_hit_PDO[i]->Write();
    board_hit_CH[i]->Write();
    board_hit_Q[i]->Write();
    board_hit_PDO_v_CH[i]->Write();
    board_hit_Q_v_CH[i]->Write();
    board_duphit_PDO[i]->Write();
    board_duphit_NCH[i]->Write();
    board_duphit_CH[i]->Write();
    board_duphit_Q[i]->Write();
    board_duphit_PDO_v_CH[i]->Write();
    board_duphit_Q_v_CH[i]->Write();
    board_NHit[i]->Write();
    board_Ndup[i]->Write();
    board_Nclus[i]->Write();
    board_Nclusdup[i]->Write();
    board_clus_CH[i]->Write();
    board_clus_Q[i]->Write();
    board_clusN_CH[i]->Write();
    board_clusN_Q[i]->Write();
    board_clus_NHit[i]->Write();
    board_clus_Ndup[i]->Write();
    board_clusN_Ndup[i]->Write();
    board_clus_Q_v_CH[i]->Write();
    board_clusN_Q_v_CH[i]->Write();
    board_clusdup_CH[i]->Write();
    board_clusdup_Q[i]->Write();
    board_clusdup_NHit[i]->Write();
    board_clusdup_Q_v_CH[i]->Write();
    board_itrack_resX[i]->Write();
    board_otrack_resX[i]->Write();
    board_itrack_resX_v_CH[i]->Write();
    board_otrack_resX_v_CH[i]->Write();
    board_itrack_resX_v_EVT[i]->Write();
    board_otrack_resX_v_EVT[i]->Write();
  }

  for(int i = 0; i < 4; i++){
    fout->cd("histograms");
    track_sumresX2[i]->Write();
  }
  track_sumresX2_v_Nclus->Write();

  fout->cd();
  fout->mkdir("plots");
  fout->cd("plots");
  
  string title = "               Run "+string(Form("%d",DATA->RunNum));
  TCanvas* can;

  can = Plot_Octuplet("c_board_PDO_v_EVT", board_PDO_v_EVT, "Event Number", "PDO [counts]", "Number of hits",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_Q_v_EVT", board_Q_v_EVT, "Event Number", "PDO [fC]", "Number of hits",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_NHit_v_EVT", board_NHit_v_EVT, "Event Number", "Number of hits", "Number of events",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_Ndup_v_EVT", board_Ndup_v_EVT, "Event Number", "Number of duplicate hits", "Number of events",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_Nclus_v_EVT", board_Nclus_v_EVT, "Event Number", "Number of clusters", "Number of events",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_Nclusdup_v_EVT", board_Nclusdup_v_EVT, "Event Number", "Number of clusters with duplicates", "Number of events",
		      iboards, title, true);
  can->Write();
  delete can;

  can = Plot_Octuplet("c_board_PDO", board_hit_PDO, "PDO [counts]", "Number of hits",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_CH", board_hit_CH, "Channel", "Number of hits",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_Q", board_hit_Q, "PDO [fC]", "Number of hits",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_PDO_v_CH", board_hit_PDO_v_CH, "Channel", "PDO [counts]", "Number of hits",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_Q_v_CH", board_hit_Q_v_CH, "Channel", "PDO [fC]", "Number of hits",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_duphit_PDO", board_duphit_PDO, "PDO [counts]", "Number of duplicate hits",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_duphit_NCH", board_duphit_NCH, "Channel", "Number of duplicate hits",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_duphit_CH", board_duphit_CH, "Channel", "Number of events with duplicate hit",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_duphit_Q", board_duphit_Q, "PDO [fC]", "Number of duplicate hits",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_duphit_PDO_v_CH", board_duphit_PDO_v_CH, "Channel", "PDO [counts]", "Number of duplicate hits",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_duphit_Q_v_CH", board_duphit_Q_v_CH, "Channel", "PDO [fC]", "Number of duplicate hits",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_NHit", board_NHit, "Number of Hits", "Number of events",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_Ndup", board_Ndup, "Number of channels with duplicate hits", "Number of events",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_Nclus", board_Nclus, "Number of clusters", "Number of events",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_Nclusdup", board_Nclusdup, "Number of clusters with duplicates", "Number of events",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clus_CH", board_clus_CH, "Channel", "Number of clusters",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clus_Q", board_clus_Q, "Charge [fC]", "Number of clusters",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clusN_CH", board_clusN_CH, "Channel", "Number of clusters",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clusN_Q", board_clusN_Q, "Charge [fC]", "Number of clusters",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clus_NHit", board_clus_NHit, "Number of hits", "Number of clusters",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clus_Ndup", board_clus_Ndup, "Number of duplicate hits", "Number of clusters",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clusN_Ndup", board_clusN_Ndup, "Number of duplicate hits", "Number of clusters",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clus_Q_v_CH", board_clus_Q_v_CH, "Channel", "Charge [fC]", "Number of clusters",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clusN_Q_v_CH", board_clusN_Q_v_CH, "Channel", "Charge [fC]", "Number of clusters",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clusdup_CH", board_clusdup_CH, "Channel", "Number of clusters",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clusdup_Q", board_clusdup_Q, "Charge [fC]", "Number of clusters",
		      iboards, title);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clusdup_NHit", board_clusdup_NHit, "Number of hits", "Number of clusters",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_clusdup_Q_v_CH", board_clusdup_Q_v_CH, "Channel", "Charge [fC]", "Number of clusters",
		      iboards, title);
  can->Write();
  delete can;
 
  can = Plot_Octuplet("c_board_itrack_resX", board_itrack_resX, "X residual [mm]", "Number of clusters",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_otrack_resX", board_otrack_resX, "X residual [mm]", "Number of clusters",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_itrack_resX_v_CH", board_itrack_resX_v_CH, "Channel", "X residual [mm]", "Number of clusters",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_otrack_resX_v_CH", board_otrack_resX_v_CH, "Channel", "X residual [mm]", "Number of clusters",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_itrack_resX_v_EVT", board_itrack_resX_v_EVT, "Event number", "X residual [mm]", "Number of clusters",
		      iboards, title, true);
  can->Write();
  delete can;
  can = Plot_Octuplet("c_board_otrack_resX_v_EVT", board_otrack_resX_v_EVT, "Event number", "X residual [mm]", "Number of clusters",
		      iboards, title, true);
  can->Write();
  delete can;

  can = Plot_2D("c_track_sumresX2_v_Nclus", track_sumresX2_v_Nclus, "Number of clusters", 
		"#Sigma ( X residual )^2 [mm^{2}]", "Number of tracks", title);
  can->Write();
  delete can;

  for(int i = 0; i < 4; i++){
    track_sumresX2[i]->Write();
    can = Plot_1D(Form("c_track_sumresX2_%d",i+5), track_sumresX2[i], 
		  "#Sigma ( X residual )^2 [mm^{2}]", 
		  "Number of tracks", title);
    can->Write();
    delete can;
  }

  fout->Close();
    
}
