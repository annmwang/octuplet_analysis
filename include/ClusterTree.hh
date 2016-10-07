//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct  6 19:50:08 2016 by ROOT version 5.34/34
// from TTree ClusterTree/ClusterTree
// found on file: clusters.root
//////////////////////////////////////////////////////////

#ifndef ClusterTree_h
#define ClusterTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class ClusterTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           N_clus;
   vector<int>     *clus_MMFE8;
   vector<int>     *clus_Index;
   vector<double>  *clus_Charge;
   vector<double>  *clus_Time;
   vector<double>  *clus_Channel;

   // List of branches
   TBranch        *b_N_clus;   //!
   TBranch        *b_clus_MMFE8;   //!
   TBranch        *b_clus_Index;   //!
   TBranch        *b_clus_Charge;   //!
   TBranch        *b_clus_Time;   //!
   TBranch        *b_clus_Channel;   //!

   ClusterTree(TTree *tree=0);
   virtual ~ClusterTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

inline ClusterTree::ClusterTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("clusters.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("clusters.root");
      }
      f->GetObject("ClusterTree",tree);

   }
   Init(tree);
}

inline ClusterTree::~ClusterTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

inline Int_t ClusterTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
inline Long64_t ClusterTree::LoadTree(Long64_t entry)
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

inline void ClusterTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   clus_MMFE8 = 0;
   clus_Index = 0;
   clus_Charge = 0;
   clus_Time = 0;
   clus_Channel = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("N_clus", &N_clus, &b_N_clus);
   fChain->SetBranchAddress("clus_MMFE8", &clus_MMFE8, &b_clus_MMFE8);
   fChain->SetBranchAddress("clus_Index", &clus_Index, &b_clus_Index);
   fChain->SetBranchAddress("clus_Charge", &clus_Charge, &b_clus_Charge);
   fChain->SetBranchAddress("clus_Time", &clus_Time, &b_clus_Time);
   fChain->SetBranchAddress("clus_Channel", &clus_Channel, &b_clus_Channel);
   Notify();
}

inline Bool_t ClusterTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

inline void ClusterTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
inline Int_t ClusterTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
