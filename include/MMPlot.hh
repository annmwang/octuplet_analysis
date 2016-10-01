#ifndef SETSTYLE_h
#define SETSTYLE_h

#include <vector>
#include <string>

#include <TStyle.h>
#include <TColor.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TList.h>
#include <TF1.h>

#include "include/GeoOctuplet.hh"
#include "include/MMClusterList.hh"

using namespace std;

static int g_2Dcount = 0;

void MMPlot(){

  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefX(0);
  gStyle->SetCanvasDefY(0);
    
  // For the Pad:
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
    
  // For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
    
  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.15);
    
  // use large Times-Roman fonts
  gStyle->SetTitleFont(132,"xyz");  // set the all 3 axes title font
  gStyle->SetTitleFont(132," ");    // set the pad title font
  gStyle->SetTitleSize(0.06,"xyz"); // set the 3 axes title size
  gStyle->SetTitleSize(0.06," ");   // set the pad title size
  gStyle->SetLabelFont(132,"xyz");

  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelColor(1,"xyz");
  gStyle->SetTextFont(132);
  gStyle->SetTextSize(0.08);
  gStyle->SetStatFont(132);
    
  // use bold lines and markers
  gStyle->SetMarkerStyle(8);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
	
  //..Get rid of X error bars
  gStyle->SetErrorX(1.);
    
  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
    
  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

 int NZPalette = 28;
 double zcolor_s[5] = { 0.00, 0.50, 0.70, 0.82, 1.00 };
 double zcolor_r[5] = { 0.00, 0.00, 0.74, 1.00, 1.00 };
 double zcolor_g[5] = { 0.00, 0.61, 0.82, 0.70, 1.00 };
 double zcolor_b[5] = { 0.31, 0.73, 0.08, 0.00, 1.00 };

 //TColor::CreateGradientColorTable(5, zcolor_s, zcolor_r,
 //				  zcolor_g, zcolor_b, NZPalette);

  gStyle->cd();
}

const TColor blue0(7000,0.749,0.78,0.933);
const TColor blue1(7001,0.424,0.467,0.651);
const TColor blue2(7002,0.255,0.302,0.522);
const TColor blue3(7003,0.114,0.165,0.396);
const TColor blue4(7004,0.024,0.063,0.251);
const TColor green0(7010,0.737,0.949,0.784);
const TColor green1(7011,0.435,0.722,0.498);
const TColor green2(7012,0.239,0.576,0.314);
const TColor green3(7013,0.082,0.439,0.161);
const TColor green4(7014,0,0.275,0.063);
const TColor red0(7020,1,0.796,0.776);
const TColor red1(7021,0.957,0.612,0.576);
const TColor red2(7022,0.765,0.361,0.318);
const TColor red3(7023,0.58,0.157,0.11);
const TColor red4(7024,0.365,0.035,0);
const TColor yellow0(7030,1,0.933,0.776);
const TColor yellow1(7031,0.957,0.843,0.576);
const TColor yellow2(7032,0.765,0.631,0.318);
const TColor yellow3(7033,0.58,0.443,0.11);
const TColor yellow4(7034,0.365,0.259,0);
const TColor purple0(7040,0.937,0.729,0.898);
const TColor purple1(7041,0.753,0.478,0.702);
const TColor purple2(7042,0.6,0.286,0.541);
const TColor purple3(7043,0.42,0.075,0.353);
const TColor purple4(7044,0.196,0,0.161);
const TColor cyan0(7050,0.714,0.898,0.918);
const TColor cyan1(7051,0.424,0.639,0.659);
const TColor cyan2(7052,0.247,0.49,0.51);
const TColor cyan3(7053,0.067,0.329,0.357);
const TColor cyan4(7054,0,0.153,0.169);
const TColor orange0(7060,1,0.882,0.776);
const TColor orange1(7061,1,0.808,0.639);
const TColor orange2(7062,0.839,0.608,0.4);
const TColor orange3(7063,0.584,0.329,0.106);
const TColor orange4(7064,0.275,0.129,0);
const TColor lime0(7070,0.941,0.992,0.769);
const TColor lime1(7071,0.882,0.961,0.612);
const TColor lime2(7072,0.706,0.8,0.38);
const TColor lime3(7073,0.455,0.557,0.098);
const TColor lime4(7074,0.204,0.263,0);

TCanvas* Plot_Graph(string scan, TGraph* graph, string X, string Y, string title){
  TCanvas *c1 = new TCanvas(scan.c_str(),scan.c_str(),600,500);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  c1->SetBottomMargin(0.15);
  c1->SetTopMargin(0.08);
  c1->Draw();
  c1->SetGridx();
  c1->SetGridy();

  double x[1];
  x[0] = 0.;
  double y[1];
  y[0] = 0.;
  TGraph* gr = new TGraph(1,x,y);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(graph);
  mg->Add(gr);

  mg->Draw("ap");
  mg->GetXaxis()->SetTitle(X.c_str());
  mg->GetXaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleOffset(1.08);
  mg->GetYaxis()->SetTitle(Y.c_str());
  mg->GetYaxis()->SetTitleOffset(1.3);
  mg->GetYaxis()->CenterTitle();
  graph->SetMarkerStyle(4);
  graph->SetMarkerColor(7013);
  graph->SetMarkerSize(3);
  
  // graph->Draw("ap");
  // graph->GetXaxis()->SetTitle(X.c_str());
  // graph->GetXaxis()->CenterTitle();
  // graph->GetXaxis()->SetTitleOffset(1.08);
  // graph->GetYaxis()->SetTitle(Y.c_str());
  // graph->GetYaxis()->SetTitleOffset(1.3);
  // graph->GetYaxis()->CenterTitle();
  // graph->SetMarkerStyle(4);
  // graph->SetMarkerColor(7013);
  // graph->SetMarkerSize(3);
  
  TLatex l;
  l.SetTextFont(132);	
  l.SetNDC();	
  l.SetTextSize(0.05);
  l.SetTextFont(132);
  l.DrawLatex(0.5,0.94,title.c_str());
  l.SetTextSize(0.045);
  l.SetTextFont(42);
  l.DrawLatex(0.02,0.94,"#bf{#it{ATLAS}} Internal - MMFE8+VMM2");

  TF1* func = (TF1*) graph->GetListOfFunctions()->First();
  if(func){
    func->SetLineColor(7022);
    func->SetLineWidth(3);
    // func->SetNpx(50000);
    // TH1F* hfunc = (TH1F*)func->GetHistogram();
    // hfunc->SetLineColor(7012);
    // hfunc->SetLineWidth(3);
    // hfunc->Draw("SAME");
  }
	
  return c1;
}

TCanvas* Plot_2D(string scan, TH2D* histo, string X, string Y, string Z, string title){
  TCanvas *c1 = new TCanvas(scan.c_str(),scan.c_str(),600,500);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.22);
  c1->SetBottomMargin(0.15);
  c1->SetTopMargin(0.08);
  c1->Draw();
  c1->SetGridx();
  c1->SetGridy();
  c1->SetLogz();
  
  histo->Draw("COLZ");
  histo->GetXaxis()->SetTitle(X.c_str());
  histo->GetXaxis()->SetTitleOffset(1.08);
  histo->GetXaxis()->CenterTitle();
  histo->GetYaxis()->SetTitle(Y.c_str());
  histo->GetYaxis()->SetTitleOffset(1.11);
  histo->GetYaxis()->CenterTitle();
  histo->GetZaxis()->SetTitle(Z.c_str());
  histo->GetZaxis()->SetTitleOffset(1.5);
  histo->GetZaxis()->CenterTitle();
  histo->GetZaxis()->SetRangeUser(0.9*histo->GetMinimum(),1.1*histo->GetMaximum());
  histo->Draw("COLZ");
  
  TLatex l;
  l.SetTextFont(132);	
  l.SetNDC();	
  l.SetTextSize(0.05);
  l.SetTextFont(132);
  l.DrawLatex(0.5,0.94,title.c_str());
  l.SetTextSize(0.045);
  l.SetTextFont(42);
  l.DrawLatex(0.02,0.94,"#bf{#it{ATLAS}} Internal - MMFE8+VMM2");
	
  return c1;
}

TCanvas* Plot_1D(string can, TH1D* hist, string X, string Y, string title = ""){
  TCanvas *c1 = new TCanvas(can.c_str(),can.c_str(),700,500);
  c1->SetRightMargin(0.05);
  c1->Draw();
  c1->SetGridx();
  c1->SetGridy();

  hist->Draw();
  hist->GetXaxis()->SetTitle(X.c_str());
  hist->GetXaxis()->SetTitleOffset(1.08);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->SetTitle(Y.c_str());
  hist->GetYaxis()->SetTitleOffset(1.13);
  hist->GetYaxis()->CenterTitle();

  hist->SetLineColor(7003);
  hist->SetLineWidth(3);
  hist->SetMarkerColor(7003);
  hist->SetMarkerSize(0);
  hist->SetFillColor(7000);
  hist->SetFillStyle(3002);

  TLatex l;
  l.SetTextFont(132);	
  l.SetNDC();	
  l.SetTextSize(0.05);
  l.SetTextFont(132);
  l.DrawLatex(0.5,0.94,title.c_str());
  l.SetTextSize(0.045);
  l.SetTextFont(42);
  l.DrawLatex(0.02,0.94,"#bf{#it{ATLAS}} Internal - MMFE8+VMM2");

  TF1* func = (TF1*) hist->GetListOfFunctions()->First();
  if(func){
    func->SetLineColor(7012);
    func->SetLineWidth(3);
    func->SetNpx(50000);
    TH1F* hfunc = (TH1F*)func->GetHistogram();
    hfunc->SetLineColor(7012);
    hfunc->SetLineWidth(3);
    hfunc->Draw("SAME");
  }

  return c1;
}

TCanvas* Plot_1D(string can, vector<TH1D*>& histo, string X, string Y, 
		 string title = "", const vector<string>& label = vector<string>()){
  TCanvas *c1 = new TCanvas(can.c_str(),can.c_str(),700,500);
  c1->SetRightMargin(0.05);
  c1->Draw();
  c1->SetGridx();
  c1->SetGridy();

  int Nh = histo.size();
  int imax = 0;
  int imin = 0;
  double max = 0;
  double min = -1;
  for(int i = 0; i < Nh; i++){
    if(histo[i]->GetMaximum() > max){
      imax = i;
      max = histo[i]->GetMaximum();
    }
    if(histo[i]->GetMinimum(0.) < min || min < 0){
      imin = i;
      min = histo[i]->GetMinimum(0.);
    }
  }

  histo[imax]->Draw();
  histo[imax]->GetXaxis()->SetTitle(X.c_str());
  histo[imax]->GetXaxis()->SetTitleOffset(1.08);
  histo[imax]->GetXaxis()->CenterTitle();
  histo[imax]->GetYaxis()->SetTitle(Y.c_str());
  histo[imax]->GetYaxis()->SetTitleOffset(1.13);
  histo[imax]->GetYaxis()->CenterTitle();
  histo[imax]->GetYaxis()->SetRangeUser(0.9*min,1.1*max);

  for(int i = 0; i < Nh; i++){
    histo[i]->SetLineColor(7003 + (i%8)*10);
    histo[i]->SetLineWidth(3);
    histo[i]->SetMarkerColor(7003 + (i%8)*10);
    histo[i]->SetMarkerSize(0);
    histo[i]->SetFillColor(7000 + (i%8)*10);
    histo[i]->SetFillStyle(3002);
    histo[i]->Draw("SAME");
    TF1* func = (TF1*) histo[i]->GetListOfFunctions()->First();
    if(func){
      //func->SetLineColorAlpha(kWhite, 0);
      func->SetLineWidth(0);
      func->SetNpx(50000);
      TH1F* hfunc = (TH1F*)func->GetHistogram();
      hfunc->SetLineColor(7002 + ((i+1)%8)*10);
      hfunc->SetLineWidth(3);
      hfunc->Draw("SAME");
    }
  }

  TLegend* leg = new TLegend(0.688,0.22,0.93,0.42);
  if(label.size() == histo.size()){
    leg->SetTextFont(132);
    leg->SetTextSize(0.045);
    leg->SetFillColor(kWhite);
    leg->SetLineColor(kWhite);
    leg->SetShadowColor(kWhite);
    for(int i = 0; i < Nh; i++)
      leg->AddEntry(histo[i],label[i].c_str());
    leg->SetLineColor(kWhite);
    leg->SetFillColor(kWhite);
    leg->SetShadowColor(kWhite);
    leg->Draw("SAME");
  }

  TLatex l;
  l.SetTextFont(132);	
  l.SetNDC();	
  l.SetTextSize(0.05);
  l.SetTextFont(132);
  l.DrawLatex(0.5,0.94,title.c_str());
  l.SetTextSize(0.045);
  l.SetTextFont(42);
  l.DrawLatex(0.02,0.94,"#bf{#it{ATLAS}} Internal - MMFE8+VMM2");

  return c1;
}

TCanvas* Plot_Track2D(string can, const MMTrack& track, const GeoOctuplet& geo, 
		      const MMClusterList* clusters = 0){
  TCanvas *c1 = new TCanvas(can.c_str(),can.c_str(),700,700);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.14);
  c1->SetLeftMargin(0.14);

  TMultiGraph* mg = new TMultiGraph();

  TGraph* gr_track = track.GetXZGraph(-5., 165.);
  mg->Add(gr_track);

  vector<TGraph*> gr_planes;
  int Nplane = geo.GetNPlanes();
  for(int i = 0; i < Nplane; i++){
    gr_planes.push_back(geo[i].GetXZGraph());
    mg->Add(gr_planes[i]);
  }

  c1->Draw();
  c1->cd();
  mg->Draw("ACP");
  mg->GetXaxis()->SetTitle("x position [mm]");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetTitle("z position [mm]");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.15);

  if(clusters){
    TGraphErrors* gr_clusters = geo.GetXZGraphErrors(*clusters);
    gr_clusters->Draw("P same");
    gr_clusters->Draw("e1 same");
  }
  
  return c1;
}

TCanvas* Plot_Track2DY(string can, const MMTrack& track, const GeoOctuplet& geo, 
		       const MMClusterList* clusters = 0){
  TCanvas *c1 = new TCanvas(can.c_str(),can.c_str(),700,700);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.14);
  c1->SetLeftMargin(0.14);

  TMultiGraph* mg = new TMultiGraph();
  
  TGraph* gr_track = geo.GetXZGraph(track);
  mg->Add(gr_track);

  vector<TGraph*> gr_planes;
  int Nplane = geo.GetNPlanes();
  for(int i = 0; i < Nplane; i++){
    gr_planes.push_back(geo[i].GetXZGraph());
    mg->Add(gr_planes[i]);
  }

  c1->Draw();
  c1->cd();
  mg->Draw("ALP");
  mg->GetXaxis()->SetTitle("x position [mm]");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetTitle("z position [mm]");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.15);

  if(clusters){
    TGraph* gr_clusters = geo.GetXZGraph(*clusters);
    gr_clusters->Draw("P same");
  }
  
  return c1;
}

TCanvas* Plot_Track3D(string can, const MMTrack& track, const GeoOctuplet& geo, 
		       const MMClusterList* clusters = 0){
  TCanvas *c1 = new TCanvas(can.c_str(),can.c_str(),700,700);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.14);
  c1->SetLeftMargin(0.14);
  c1->Draw();
  c1->cd();

  double x[3];
  double y[3];
  double z[3];
  x[0] = 215.;
  y[0] = -20.;
  z[0] = -5.;
  x[1] = -15.;
  y[1] = 210.;
  z[1] = -5.;
  x[2] = -15.;
  y[2] = -20.;
  z[2] = 165.;

  TGraph2D* gr_frame = new TGraph2D(3,x,y,z);
  gr_frame->SetName(Form("g2D_%d",g_2Dcount++));
  gr_frame->SetTitle("");
  gr_frame->SetMarkerSize(0);
  gr_frame->Draw("P");
  gr_frame->GetXaxis()->SetRangeUser(-15.,215.);
  gr_frame->GetYaxis()->SetRangeUser(-20.,210.);
  gr_frame->GetZaxis()->SetRangeUser(-5.,165.);
  gr_frame->GetXaxis()->SetTitle("x position [mm]");
  gr_frame->GetXaxis()->CenterTitle();
  gr_frame->GetZaxis()->SetTitle("z position [mm]");
  gr_frame->GetZaxis()->CenterTitle();
  gr_frame->GetZaxis()->SetTitleOffset(1.15);
  
  
  TGraph2D* gr_track = track.Get2DGraph(-5., 165.);
  gr_track->SetName(Form("g2D_%d",g_2Dcount++));
  gr_track->SetTitle("");
  gr_track->Draw("LINE same");

  TGraph2D* gr_trackhits = geo.Get2DGraph(track);
  gr_trackhits->SetName(Form("g2D_%d",g_2Dcount++));
  gr_trackhits->SetTitle("");
  gr_trackhits->Draw("P same");

  vector<TGraph2D*> gr_planes;
  int Nplane = geo.GetNPlanes();
  for(int i = 0; i < Nplane; i++){
    gr_planes.push_back(geo[i].Get2DGraph());
    gr_planes[i]->SetName(Form("g2D_%d",g_2Dcount++));
    gr_planes[i]->SetTitle("");
    gr_planes[i]->Draw("LINE same");
  }

  vector<TGraph2D*> gr_clus;
  if(clusters){
    int Nclus = clusters->GetNCluster();
    for(int i = 0; i < Nclus; i++){
      int index = geo.Index(clusters->Get(i).MMFE8());
      if(index < 0)
	continue;
      TGraph2D* gr = geo[index].GetChannelGraph(clusters->Get(i).Channel());
      gr_clus.push_back(gr);
      gr->SetName(Form("g2D_%d",g_2Dcount++));
      gr->SetTitle("");
      gr->Draw("LINE same");
    }
  }
  
  
  return c1;
}

TCanvas* Plot_Octuplet1D(string can, vector<TH1D*>& histo, string X, string Y,
			 vector<int>& iboard, string title = ""){
  TCanvas *c1 = new TCanvas(can.c_str(),can.c_str(),700,500);
  c1->SetRightMargin(0.05);
  c1->Draw();
  c1->Divide(1,8,0.,0.);
  c1->SetGridx();
  c1->SetGridy();

  int Nh = histo.size();
  for(int i = 0; i < Nh; i++){
    TVirtualPad* pad = c1->GetPad(i+1);
    pad->SetRightMargin(0.0);
    pad->SetLeftMargin(0.0);
    pad->SetTopMargin(0.0);
    pad->SetBottomMargin(0.0);
    c1->cd(i+1);
    histo[i]->Draw();
  }
  
 

  // histo[imax]->Draw();
  // histo[imax]->GetXaxis()->SetTitle(X.c_str());
  // histo[imax]->GetXaxis()->SetTitleOffset(1.08);
  // histo[imax]->GetXaxis()->CenterTitle();
  // histo[imax]->GetYaxis()->SetTitle(Y.c_str());
  // histo[imax]->GetYaxis()->SetTitleOffset(1.13);
  // histo[imax]->GetYaxis()->CenterTitle();
  // histo[imax]->GetYaxis()->SetRangeUser(0.9*min,1.1*max);

  // for(int i = 0; i < Nh; i++){
  //   histo[i]->SetLineColor(7003 + (i%8)*10);
  //   histo[i]->SetLineWidth(3);
  //   histo[i]->SetMarkerColor(7003 + (i%8)*10);
  //   histo[i]->SetMarkerSize(0);
  //   histo[i]->SetFillColor(7000 + (i%8)*10);
  //   histo[i]->SetFillStyle(3002);
  //   histo[i]->Draw("SAME");
  //   TF1* func = (TF1*) histo[i]->GetListOfFunctions()->First();
  //   if(func){
  //     //func->SetLineColorAlpha(kWhite, 0);
  //     func->SetLineWidth(0);
  //     func->SetNpx(50000);
  //     TH1F* hfunc = (TH1F*)func->GetHistogram();
  //     hfunc->SetLineColor(7002 + ((i+1)%8)*10);
  //     hfunc->SetLineWidth(3);
  //     hfunc->Draw("SAME");
  //   }
  // }

  // TLegend* leg = new TLegend(0.688,0.22,0.93,0.42);
  // if(label.size() == histo.size()){
  //   leg->SetTextFont(132);
  //   leg->SetTextSize(0.045);
  //   leg->SetFillColor(kWhite);
  //   leg->SetLineColor(kWhite);
  //   leg->SetShadowColor(kWhite);
  //   for(int i = 0; i < Nh; i++)
  //     leg->AddEntry(histo[i],label[i].c_str());
  //   leg->SetLineColor(kWhite);
  //   leg->SetFillColor(kWhite);
  //   leg->SetShadowColor(kWhite);
  //   leg->Draw("SAME");
  // }

  // TLatex l;
  // l.SetTextFont(132);	
  // l.SetNDC();	
  // l.SetTextSize(0.05);
  // l.SetTextFont(132);
  // l.DrawLatex(0.5,0.94,title.c_str());
  // l.SetTextSize(0.045);
  // l.SetTextFont(42);
  // l.DrawLatex(0.02,0.94,"#bf{#it{ATLAS}} Internal - MMFE8+VMM2");

  return c1;
}

#endif
