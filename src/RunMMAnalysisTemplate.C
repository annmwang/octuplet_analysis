///
///  \file   RunMMAnalysisTemplate.C
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///

#include "TH1D.h"
#include "TH2D.h"
#include <iostream>

#include "include/MMPlot.hh"
#include "include/PDOToCharge.hh"
#include "include/TDOToTime.hh"
#include "include/MMDataAnalysis.hh"
#include "include/MMPacmanAlgo.hh"

using namespace std;

int main(int argc, char* argv[]){

  char inputFileName[400];
  char outputFileName[400];
  char PDOFileName[400];
  char TDOFileName[400];
  
  if ( argc < 5 ){
    cout << "Error at Input: please specify input/output .root files ";
    cout << " and (optional) PDO/TDO calibration files" << endl;
    cout << "Example:   ./RunMMAnalysisTemplate.x -i input.root -o output.root" << endl;
    cout << "Example:   ./RunMMAnalysisTemplate.x -i input.root -o output.root";
    cout << " -p PDOcalib.root -t TDOcalib.root" << endl;
    return 0;
  }

  bool b_input = false;
  bool b_out   = false;
  bool b_pdo   = false;
  bool b_tdo   = false;
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
  }

  if(!b_input){
    cout << "Error at Input: please specify input file (-i flag)" << endl;
    return 0;
  }

  if(!b_out){
    cout << "Error at Input: please specify output file (-o flag)" << endl;
    return 0;
  }

  // board ID's for run 3508
  map<int,int> ib;
  ib[111] = 0;
  ib[116] = 1;
  ib[101] = 2;
  ib[109] = 3;
  ib[112] = 4;
  ib[102] = 5;
  ib[107] = 6;
  ib[105] = 7;
  vector<int> iboards;
  iboards.push_back(111);
  iboards.push_back(116);
  iboards.push_back(101);
  iboards.push_back(109);
  iboards.push_back(112);
  iboards.push_back(102);
  iboards.push_back(107);
  iboards.push_back(105);
  
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

  MMPacmanAlgo* PACMAN = new MMPacmanAlgo();

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

  // histograms for analysis
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

  vector<TH1D*> board_NHit;
  vector<TH1D*> board_Ndup;
  vector<TH1D*> board_Nclus;
  vector<TH1D*> board_Nclusdup;
  
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

  for(int i = 0; i < 8; i++){
    board_hit_PDO.push_back(new TH1D(Form("b_h_PDO_%d",i),
				     Form("b_h_PDO_%d",i),
				     1028,0.,1028.));
    board_hit_CH.push_back(new TH1D(Form("b_h_CH_%d",i),
				    Form("b_h_CH_%d",i),
				    512,0.5,512.5));
    board_hit_Q.push_back(new TH1D(Form("b_h_Q_%d",i),
				   Form("b_h_Q_%d",i),
				   128,0.0,128));
    board_hit_PDO_v_CH.push_back(new TH2D(Form("b_h_PDOvCH_%d",i),
					  Form("b_h_PDOvCH_%d",i),
					  512,0.5,512.5,
					  1028,0.,1028.));
    board_hit_Q_v_CH.push_back(new TH2D(Form("b_h_QvCH_%d",i),
					Form("b_h_QvCH_%d",i),
					512,0.5,512.5,
					128,0.,128.));
    
    board_duphit_PDO.push_back(new TH1D(Form("b_dh_PDO_%d",i),
					Form("b_dh_PDO_%d",i),
					1028,0.0,1028));;
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
					 128,0.,128.));
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
					  128,0.0,128));
    board_clusdup_Q_v_CH.push_back(new TH2D(Form("b_cd_QvCH_%d",i),
					    Form("b_cd_QvCH_%d",i),
					    512,0.5,512.5,
					    128,0.,250.));
  }
  
  for(int evt = 0; evt < Nevent; evt++){
       DATA->GetEntry(evt);
    if(evt%1 == 10000)
      cout << "Processing event # " << evt << " | " << Nevent << endl;

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
      for(int j = 0; j < Nhit; j++){
	const MMLinkedHit& hit = DATA->mm_EventHits[i][j];
	board_hit_PDO[b]->Fill(hit.PDO());
	board_hit_CH[b]->Fill(hit.Channel());
	board_hit_Q[b]->Fill(hit.Charge());
	board_hit_PDO_v_CH[b]->Fill(hit.Channel(),hit.PDO());
	board_hit_Q_v_CH[b]->Fill(hit.Channel(),hit.Charge());
	  
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
      if(clusters.GetNCluster() > 0)
	all_clusters.push_back(clusters);
    }
    
    int Ncl = all_clusters.size();
    for(int i = 0; i < Ncl; i++){
      int Nc = all_clusters[i].GetNCluster();
      for(int j = 0; j < Nc; j++){
	const MMCluster& clus = all_clusters[i][j];
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
  }

  TFile* fout = new TFile(outputFileName, "RECREATE");
  fout->cd();

  setstyle();
  
  TCanvas* can;
  can = Plot_Octuplet1D("test", board_hit_PDO, "PDO", "Number of hits",
			iboards, "Title");
  can->Write();
  delete can;
  
  for(int i = 0; i < 8; i++){
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
  }
  fout->Close();
    
}
