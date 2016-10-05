#include <iostream>

#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"

using namespace std;

void Hollywood(string filename, string canname){
  TFile* f = new TFile(filename.c_str(),"READ");
  if(!f->IsOpen())
    return;

  TCanvas* can = (TCanvas*) f->Get(canname.c_str());
  if(!can)
    return;

  gSystem->Exec("rm -rf hollywood");
  gSystem->Exec("mkdir hollywood");

  int gifcount = 0;
  
  int Nscene = 1000;
  for(int i = 0; i < Nscene; i++){
    double phi = -30. + (330.+30.)*float(i)/float(Nscene);
    double theta = 5. + 50.*(1.-pow(sin(phi/180.*acos(-1.)),2.));
    can->SetPhi(phi);
    can->SetTheta(theta);
    can->Draw();
    can->Modified();
    can->Update();
    can->SaveAs(Form("hollywood/hollywood_%.4d.gif",gifcount));
    gifcount++;
  }
  //can->Print("hollywood.gif++");

  gSystem->Exec("gifsicle --delay=10 --loop=5 --colors 256 hollywood/h*.gif > hollywood.gif");
  gSystem->Exec("rm -rf hollywood");
  
  delete can;
}
