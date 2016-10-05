#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include "include/MMPlot.hh"

void Test_Plot_Octuplet(string filename){
  MMPlot();


  vector<int> iboards;
  iboards.push_back(111);
  iboards.push_back(116);
  iboards.push_back(101);
  iboards.push_back(109);
  iboards.push_back(112);
  iboards.push_back(102);
  iboards.push_back(107);
  iboards.push_back(105);

  TFile* f = new TFile(filename.c_str(),"READ");

  vector<TH1D*> board_hit_PDO;
  for(int i = 0; i < 8; i++)
    board_hit_PDO.push_back((TH1D*)f->Get(Form("histograms/b_h_PDO_%d",i)));


  TCanvas* can = Plot_Octuplet("test", board_hit_PDO, "PDO", "Number of hits",
			       iboards, "Run 3508", true);


  vector<TH2D*> board_hit_PDOvCH;
  for(int i = 0; i < 8; i++)
    board_hit_PDOvCH.push_back((TH2D*)f->Get(Form("histograms/b_h_PDOvCH_%d",i)));


  TCanvas* can2D = Plot_Octuplet("test2D", board_hit_PDOvCH, "Channel", "PDO", "Number of hits",
				 iboards, "Run 3508", false);






}
