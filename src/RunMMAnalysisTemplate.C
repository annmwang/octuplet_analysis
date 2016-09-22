///
///  \file   RunMMAnalysisTemplate.C
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///

#include <iostream>

#include "include/PDOToCharge.hh"
#include "include/TDOToTime.hh"
#include "include/MMDataAnalysis.hh"

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
    if (strncmp(argv[i],"-o",2)==0){
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

  for(int evt = 0; evt < Nevent; evt++){
    DATA->GetEntry(evt);

    PDOCalibrator->Calibrate(DATA->mm_EventHits);
    TDOCalibrator->Calibrate(DATA->mm_EventHits);
  
    int Nboard = DATA->mm_EventHits.GetNBoards();
    cout << "N boards with hits: " << Nboard << endl;

    for(int i = 0; i < Nboard; i++){
      int Nhit = DATA->mm_EventHits[i]->GetNHits();
      cout << "board " << i << "(" << DATA->mm_EventHits[i]->MMFE8() << ") ";
      cout << "has " << Nhit << " hits" << endl;
      for(int j = 0; j < Nhit; j++){
	cout << "(" << DATA->mm_EventHits[i]->Get(j)->Charge() << ", " << DATA->mm_EventHits[i]->Get(j)->Time() << ") ";
      }
      cout << endl;
    }
    cout << endl << endl;
  }
}
