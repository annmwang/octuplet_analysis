//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun May 21 23:41:20 2017 by ROOT version 6.06/00
// from TTree COMB_data/COMB_data
// found on file: data/Run3522_part/Run3522_combSCMMTPTIME.root
//////////////////////////////////////////////////////////

#ifndef MMDataBase_h
#define MMDataBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class MMDataBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           RunNum;
   Int_t           sci_EventNum;
   Int_t           sci_Time_sec;
   Int_t           sci_Time_nsec;
   Int_t           mm_EventNum;
   Int_t           mm_Time_sec;
   Int_t           mm_Time_nsec;
   Int_t           mm_trig_BCID;
   Int_t           N_sci;
   vector<int>     *sci_CH;
   vector<int>     *sci_count;
   Int_t           N_mm;
   vector<int>     *mm_VMM;
   vector<int>     *mm_CH;
   vector<int>     *mm_PDO;
   vector<int>     *mm_TDO;
   vector<int>     *mm_BCID;
   vector<int>     *mm_trigBCID;
   vector<int>     *mm_trigphase;
   vector<int>     *mm_MMFE8;
   vector<int>     *mm_FIFOcount;
   Int_t           tp_n;
   vector<int>     *tp_EventNum;
   vector<int>     *tp_cntr;
   vector<int>     *tp_Time_sec;
   vector<int>     *tp_Time_nsec;
   vector<int>     *tp_BCID;
   vector<int>     *tp_hit_n;
   vector<double>  *tp_mxlocal;
   vector<vector<int> > *tp_hit_MMFE8;
   vector<vector<int> > *tp_hit_VMM;
   vector<vector<int> > *tp_hit_CH;
   vector<vector<int> > *tp_hit_BCID;
   Int_t           tpsci_EventNum;
   Int_t           tpsci_BCID;
   Int_t           tpsci_ph;
   Int_t           tpsci_overflow;
   Int_t           tpsci_Time_sec;
   Int_t           tpsci_Time_nsec;

   // List of branches
   TBranch        *b_RunNum;   //!
   TBranch        *b_sci_EventNum;   //!
   TBranch        *b_sci_Time_sec;   //!
   TBranch        *b_sci_Time_nsec;   //!
   TBranch        *b_mm_EventNum;   //!
   TBranch        *b_mm_Time_sec;   //!
   TBranch        *b_mm_Time_nsec;   //!
   TBranch        *b_mm_trig_BCID;   //!
   TBranch        *b_N_sci;   //!
   TBranch        *b_sci_CH;   //!
   TBranch        *b_sci_count;   //!
   TBranch        *b_N_mm;   //!
   TBranch        *b_mm_VMM;   //!
   TBranch        *b_mm_CH;   //!
   TBranch        *b_mm_PDO;   //!
   TBranch        *b_mm_TDO;   //!
   TBranch        *b_mm_BCID;   //!
   TBranch        *b_mm_trigBCID;   //!
   TBranch        *b_mm_trigphase;   //!
   TBranch        *b_mm_MMFE8;   //!
   TBranch        *b_mm_FIFOcount;   //!
   TBranch        *b_tp_n;   //!
   TBranch        *b_tp_EventNum;   //!
   TBranch        *b_tp_cntr;   //!
   TBranch        *b_tp_Time_sec;   //!
   TBranch        *b_tp_Time_nsec;   //!
   TBranch        *b_tp_BCID;   //!
   TBranch        *b_tp_hit_n;   //!
   TBranch        *b_tp_mxlocal;   //!
   TBranch        *b_tp_hit_MMFE8;   //!
   TBranch        *b_tp_hit_VMM;   //!
   TBranch        *b_tp_hit_CH;   //!
   TBranch        *b_tp_hit_BCID;   //!
   TBranch        *b_tpsci_EventNum;   //!
   TBranch        *b_tpsci_BCID;   //!
   TBranch        *b_tpsci_ph;   //!
   TBranch        *b_tpsci_overflow;   //!
   TBranch        *b_tpsci_Time_sec;   //!
   TBranch        *b_tpsci_Time_nsec;   //!

   MMDataBase(TTree *tree=0);
   virtual ~MMDataBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

inline MMDataBase::MMDataBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data/Run3522_part/Run3522_combSCMMTPTIME.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data/Run3522_part/Run3522_combSCMMTPTIME.root");
      }
      f->GetObject("COMB_data",tree);

   }
   Init(tree);
}

inline MMDataBase::~MMDataBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

inline Int_t MMDataBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
inline Long64_t MMDataBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

inline void MMDataBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   sci_CH = 0;
   sci_count = 0;
   mm_VMM = 0;
   mm_CH = 0;
   mm_PDO = 0;
   mm_TDO = 0;
   mm_BCID = 0;
   mm_MMFE8 = 0;
   mm_FIFOcount = 0;
   tp_EventNum = 0;
   tp_cntr = 0;
   tp_Time_sec = 0;
   tp_Time_nsec = 0;
   tp_BCID = 0;
   tp_hit_n = 0;
   tp_mxlocal = 0;
   tp_hit_MMFE8 = 0;
   tp_hit_VMM = 0;
   tp_hit_CH = 0;
   tp_hit_BCID = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNum", &RunNum, &b_RunNum);
   fChain->SetBranchAddress("sci_EventNum", &sci_EventNum, &b_sci_EventNum);
   fChain->SetBranchAddress("sci_Time_sec", &sci_Time_sec, &b_sci_Time_sec);
   fChain->SetBranchAddress("sci_Time_nsec", &sci_Time_nsec, &b_sci_Time_nsec);
   fChain->SetBranchAddress("mm_EventNum", &mm_EventNum, &b_mm_EventNum);
   fChain->SetBranchAddress("mm_Time_sec", &mm_Time_sec, &b_mm_Time_sec);
   fChain->SetBranchAddress("mm_Time_nsec", &mm_Time_nsec, &b_mm_Time_nsec);
   fChain->SetBranchAddress("mm_trig_BCID", &mm_trig_BCID, &b_mm_trig_BCID);
   fChain->SetBranchAddress("N_sci", &N_sci, &b_N_sci);
   fChain->SetBranchAddress("sci_CH", &sci_CH, &b_sci_CH);
   fChain->SetBranchAddress("sci_count", &sci_count, &b_sci_count);
   fChain->SetBranchAddress("N_mm", &N_mm, &b_N_mm);
   fChain->SetBranchAddress("mm_VMM", &mm_VMM, &b_mm_VMM);
   fChain->SetBranchAddress("mm_CH", &mm_CH, &b_mm_CH);
   fChain->SetBranchAddress("mm_PDO", &mm_PDO, &b_mm_PDO);
   fChain->SetBranchAddress("mm_TDO", &mm_TDO, &b_mm_TDO);
   fChain->SetBranchAddress("mm_BCID", &mm_BCID, &b_mm_BCID);
   fChain->SetBranchAddress("mm_trigBCID", &mm_trigBCID, &b_mm_trigBCID);
   fChain->SetBranchAddress("mm_trigphase", &mm_trigphase, &b_mm_trigphase);
   fChain->SetBranchAddress("mm_MMFE8", &mm_MMFE8, &b_mm_MMFE8);
   fChain->SetBranchAddress("mm_FIFOcount", &mm_FIFOcount, &b_mm_FIFOcount);
   fChain->SetBranchAddress("tp_n", &tp_n, &b_tp_n);
   fChain->SetBranchAddress("tp_EventNum", &tp_EventNum, &b_tp_EventNum);
   fChain->SetBranchAddress("tp_cntr", &tp_cntr, &b_tp_cntr);
   fChain->SetBranchAddress("tp_Time_sec", &tp_Time_sec, &b_tp_Time_sec);
   fChain->SetBranchAddress("tp_Time_nsec", &tp_Time_nsec, &b_tp_Time_nsec);
   fChain->SetBranchAddress("tp_BCID", &tp_BCID, &b_tp_BCID);
   fChain->SetBranchAddress("tp_hit_n", &tp_hit_n, &b_tp_hit_n);
   fChain->SetBranchAddress("tp_mxlocal", &tp_mxlocal, &b_tp_mxlocal);
   fChain->SetBranchAddress("tp_hit_MMFE8", &tp_hit_MMFE8, &b_tp_hit_MMFE8);
   fChain->SetBranchAddress("tp_hit_VMM", &tp_hit_VMM, &b_tp_hit_VMM);
   fChain->SetBranchAddress("tp_hit_CH", &tp_hit_CH, &b_tp_hit_CH);
   fChain->SetBranchAddress("tp_hit_BCID", &tp_hit_BCID, &b_tp_hit_BCID);
   fChain->SetBranchAddress("tpsci_EventNum", &tpsci_EventNum, &b_tpsci_EventNum);
   fChain->SetBranchAddress("tpsci_BCID", &tpsci_BCID, &b_tpsci_BCID);
   fChain->SetBranchAddress("tpsci_ph", &tpsci_ph, &b_tpsci_ph);
   fChain->SetBranchAddress("tpsci_overflow", &tpsci_overflow, &b_tpsci_overflow);
   fChain->SetBranchAddress("tpsci_Time_sec", &tpsci_Time_sec, &b_tpsci_Time_sec);
   fChain->SetBranchAddress("tpsci_Time_nsec", &tpsci_Time_nsec, &b_tpsci_Time_nsec);
   Notify();
}

inline Bool_t MMDataBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

inline void MMDataBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
inline Int_t MMDataBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
