//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct  7 12:24:46 2016 by ROOT version 5.34/25
// from TTree AlignmentTree/AlignmentTree
// found on file: align_3513_tranX.root
//////////////////////////////////////////////////////////

#ifndef AlignmentTree_h
#define AlignmentTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class AlignmentTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        tranX_1;
   Double_t        tranX_err_1;
   Double_t        tranX_2;
   Double_t        tranX_err_2;
   Double_t        tranX_3;
   Double_t        tranX_err_3;
   Double_t        tranX_4;
   Double_t        tranX_err_4;
   Double_t        tranX_5;
   Double_t        tranX_err_5;
   Double_t        tranX_6;
   Double_t        tranX_err_6;
   Double_t        tranX_7;
   Double_t        tranX_err_7;
   Double_t        tranY_1;
   Double_t        tranY_err_1;
   Double_t        tranY_2;
   Double_t        tranY_err_2;
   Double_t        tranY_3;
   Double_t        tranY_err_3;
   Double_t        tranY_4;
   Double_t        tranY_err_4;
   Double_t        tranY_5;
   Double_t        tranY_err_5;
   Double_t        tranY_6;
   Double_t        tranY_err_6;
   Double_t        tranY_7;
   Double_t        tranY_err_7;
   Double_t        tranZ_1;
   Double_t        tranZ_err_1;
   Double_t        tranZ_2;
   Double_t        tranZ_err_2;
   Double_t        tranZ_3;
   Double_t        tranZ_err_3;
   Double_t        tranZ_4;
   Double_t        tranZ_err_4;
   Double_t        tranZ_5;
   Double_t        tranZ_err_5;
   Double_t        tranZ_6;
   Double_t        tranZ_err_6;
   Double_t        tranZ_7;
   Double_t        tranZ_err_7;
   Double_t        rotX_1;
   Double_t        rotX_err_1;
   Double_t        rotX_2;
   Double_t        rotX_err_2;
   Double_t        rotX_3;
   Double_t        rotX_err_3;
   Double_t        rotX_4;
   Double_t        rotX_err_4;
   Double_t        rotX_5;
   Double_t        rotX_err_5;
   Double_t        rotX_6;
   Double_t        rotX_err_6;
   Double_t        rotX_7;
   Double_t        rotX_err_7;
   Double_t        rotY_1;
   Double_t        rotY_err_1;
   Double_t        rotY_2;
   Double_t        rotY_err_2;
   Double_t        rotY_3;
   Double_t        rotY_err_3;
   Double_t        rotY_4;
   Double_t        rotY_err_4;
   Double_t        rotY_5;
   Double_t        rotY_err_5;
   Double_t        rotY_6;
   Double_t        rotY_err_6;
   Double_t        rotY_7;
   Double_t        rotY_err_7;
   Double_t        rotZ_1;
   Double_t        rotZ_err_1;
   Double_t        rotZ_2;
   Double_t        rotZ_err_2;
   Double_t        rotZ_3;
   Double_t        rotZ_err_3;
   Double_t        rotZ_4;
   Double_t        rotZ_err_4;
   Double_t        rotZ_5;
   Double_t        rotZ_err_5;
   Double_t        rotZ_6;
   Double_t        rotZ_err_6;
   Double_t        rotZ_7;
   Double_t        rotZ_err_7;

   // List of branches
   TBranch        *b_tranX_1;   //!
   TBranch        *b_tranX_err_1;   //!
   TBranch        *b_tranX_2;   //!
   TBranch        *b_tranX_err_2;   //!
   TBranch        *b_tranX_3;   //!
   TBranch        *b_tranX_err_3;   //!
   TBranch        *b_tranX_4;   //!
   TBranch        *b_tranX_err_4;   //!
   TBranch        *b_tranX_5;   //!
   TBranch        *b_tranX_err_5;   //!
   TBranch        *b_tranX_6;   //!
   TBranch        *b_tranX_err_6;   //!
   TBranch        *b_tranX_7;   //!
   TBranch        *b_tranX_err_7;   //!
   TBranch        *b_tranY_1;   //!
   TBranch        *b_tranY_err_1;   //!
   TBranch        *b_tranY_2;   //!
   TBranch        *b_tranY_err_2;   //!
   TBranch        *b_tranY_3;   //!
   TBranch        *b_tranY_err_3;   //!
   TBranch        *b_tranY_4;   //!
   TBranch        *b_tranY_err_4;   //!
   TBranch        *b_tranY_5;   //!
   TBranch        *b_tranY_err_5;   //!
   TBranch        *b_tranY_6;   //!
   TBranch        *b_tranY_err_6;   //!
   TBranch        *b_tranY_7;   //!
   TBranch        *b_tranY_err_7;   //!
   TBranch        *b_tranZ_1;   //!
   TBranch        *b_tranZ_err_1;   //!
   TBranch        *b_tranZ_2;   //!
   TBranch        *b_tranZ_err_2;   //!
   TBranch        *b_tranZ_3;   //!
   TBranch        *b_tranZ_err_3;   //!
   TBranch        *b_tranZ_4;   //!
   TBranch        *b_tranZ_err_4;   //!
   TBranch        *b_tranZ_5;   //!
   TBranch        *b_tranZ_err_5;   //!
   TBranch        *b_tranZ_6;   //!
   TBranch        *b_tranZ_err_6;   //!
   TBranch        *b_tranZ_7;   //!
   TBranch        *b_tranZ_err_7;   //!
   TBranch        *b_rotX_1;   //!
   TBranch        *b_rotX_err_1;   //!
   TBranch        *b_rotX_2;   //!
   TBranch        *b_rotX_err_2;   //!
   TBranch        *b_rotX_3;   //!
   TBranch        *b_rotX_err_3;   //!
   TBranch        *b_rotX_4;   //!
   TBranch        *b_rotX_err_4;   //!
   TBranch        *b_rotX_5;   //!
   TBranch        *b_rotX_err_5;   //!
   TBranch        *b_rotX_6;   //!
   TBranch        *b_rotX_err_6;   //!
   TBranch        *b_rotX_7;   //!
   TBranch        *b_rotX_err_7;   //!
   TBranch        *b_rotY_1;   //!
   TBranch        *b_rotY_err_1;   //!
   TBranch        *b_rotY_2;   //!
   TBranch        *b_rotY_err_2;   //!
   TBranch        *b_rotY_3;   //!
   TBranch        *b_rotY_err_3;   //!
   TBranch        *b_rotY_4;   //!
   TBranch        *b_rotY_err_4;   //!
   TBranch        *b_rotY_5;   //!
   TBranch        *b_rotY_err_5;   //!
   TBranch        *b_rotY_6;   //!
   TBranch        *b_rotY_err_6;   //!
   TBranch        *b_rotY_7;   //!
   TBranch        *b_rotY_err_7;   //!
   TBranch        *b_rotZ_1;   //!
   TBranch        *b_rotZ_err_1;   //!
   TBranch        *b_rotZ_2;   //!
   TBranch        *b_rotZ_err_2;   //!
   TBranch        *b_rotZ_3;   //!
   TBranch        *b_rotZ_err_3;   //!
   TBranch        *b_rotZ_4;   //!
   TBranch        *b_rotZ_err_4;   //!
   TBranch        *b_rotZ_5;   //!
   TBranch        *b_rotZ_err_5;   //!
   TBranch        *b_rotZ_6;   //!
   TBranch        *b_rotZ_err_6;   //!
   TBranch        *b_rotZ_7;   //!
   TBranch        *b_rotZ_err_7;   //!

   AlignmentTree(TTree *tree=0);
   virtual ~AlignmentTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

inline AlignmentTree::AlignmentTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("align_3513_tranX.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("align_3513_tranX.root");
      }
      f->GetObject("AlignmentTree",tree);

   }
   Init(tree);
}

inline AlignmentTree::~AlignmentTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

inline Int_t AlignmentTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
inline Long64_t AlignmentTree::LoadTree(Long64_t entry)
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

inline void AlignmentTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("tranX_1", &tranX_1, &b_tranX_1);
   fChain->SetBranchAddress("tranX_err_1", &tranX_err_1, &b_tranX_err_1);
   fChain->SetBranchAddress("tranX_2", &tranX_2, &b_tranX_2);
   fChain->SetBranchAddress("tranX_err_2", &tranX_err_2, &b_tranX_err_2);
   fChain->SetBranchAddress("tranX_3", &tranX_3, &b_tranX_3);
   fChain->SetBranchAddress("tranX_err_3", &tranX_err_3, &b_tranX_err_3);
   fChain->SetBranchAddress("tranX_4", &tranX_4, &b_tranX_4);
   fChain->SetBranchAddress("tranX_err_4", &tranX_err_4, &b_tranX_err_4);
   fChain->SetBranchAddress("tranX_5", &tranX_5, &b_tranX_5);
   fChain->SetBranchAddress("tranX_err_5", &tranX_err_5, &b_tranX_err_5);
   fChain->SetBranchAddress("tranX_6", &tranX_6, &b_tranX_6);
   fChain->SetBranchAddress("tranX_err_6", &tranX_err_6, &b_tranX_err_6);
   fChain->SetBranchAddress("tranX_7", &tranX_7, &b_tranX_7);
   fChain->SetBranchAddress("tranX_err_7", &tranX_err_7, &b_tranX_err_7);
   fChain->SetBranchAddress("tranY_1", &tranY_1, &b_tranY_1);
   fChain->SetBranchAddress("tranY_err_1", &tranY_err_1, &b_tranY_err_1);
   fChain->SetBranchAddress("tranY_2", &tranY_2, &b_tranY_2);
   fChain->SetBranchAddress("tranY_err_2", &tranY_err_2, &b_tranY_err_2);
   fChain->SetBranchAddress("tranY_3", &tranY_3, &b_tranY_3);
   fChain->SetBranchAddress("tranY_err_3", &tranY_err_3, &b_tranY_err_3);
   fChain->SetBranchAddress("tranY_4", &tranY_4, &b_tranY_4);
   fChain->SetBranchAddress("tranY_err_4", &tranY_err_4, &b_tranY_err_4);
   fChain->SetBranchAddress("tranY_5", &tranY_5, &b_tranY_5);
   fChain->SetBranchAddress("tranY_err_5", &tranY_err_5, &b_tranY_err_5);
   fChain->SetBranchAddress("tranY_6", &tranY_6, &b_tranY_6);
   fChain->SetBranchAddress("tranY_err_6", &tranY_err_6, &b_tranY_err_6);
   fChain->SetBranchAddress("tranY_7", &tranY_7, &b_tranY_7);
   fChain->SetBranchAddress("tranY_err_7", &tranY_err_7, &b_tranY_err_7);
   fChain->SetBranchAddress("tranZ_1", &tranZ_1, &b_tranZ_1);
   fChain->SetBranchAddress("tranZ_err_1", &tranZ_err_1, &b_tranZ_err_1);
   fChain->SetBranchAddress("tranZ_2", &tranZ_2, &b_tranZ_2);
   fChain->SetBranchAddress("tranZ_err_2", &tranZ_err_2, &b_tranZ_err_2);
   fChain->SetBranchAddress("tranZ_3", &tranZ_3, &b_tranZ_3);
   fChain->SetBranchAddress("tranZ_err_3", &tranZ_err_3, &b_tranZ_err_3);
   fChain->SetBranchAddress("tranZ_4", &tranZ_4, &b_tranZ_4);
   fChain->SetBranchAddress("tranZ_err_4", &tranZ_err_4, &b_tranZ_err_4);
   fChain->SetBranchAddress("tranZ_5", &tranZ_5, &b_tranZ_5);
   fChain->SetBranchAddress("tranZ_err_5", &tranZ_err_5, &b_tranZ_err_5);
   fChain->SetBranchAddress("tranZ_6", &tranZ_6, &b_tranZ_6);
   fChain->SetBranchAddress("tranZ_err_6", &tranZ_err_6, &b_tranZ_err_6);
   fChain->SetBranchAddress("tranZ_7", &tranZ_7, &b_tranZ_7);
   fChain->SetBranchAddress("tranZ_err_7", &tranZ_err_7, &b_tranZ_err_7);
   fChain->SetBranchAddress("rotX_1", &rotX_1, &b_rotX_1);
   fChain->SetBranchAddress("rotX_err_1", &rotX_err_1, &b_rotX_err_1);
   fChain->SetBranchAddress("rotX_2", &rotX_2, &b_rotX_2);
   fChain->SetBranchAddress("rotX_err_2", &rotX_err_2, &b_rotX_err_2);
   fChain->SetBranchAddress("rotX_3", &rotX_3, &b_rotX_3);
   fChain->SetBranchAddress("rotX_err_3", &rotX_err_3, &b_rotX_err_3);
   fChain->SetBranchAddress("rotX_4", &rotX_4, &b_rotX_4);
   fChain->SetBranchAddress("rotX_err_4", &rotX_err_4, &b_rotX_err_4);
   fChain->SetBranchAddress("rotX_5", &rotX_5, &b_rotX_5);
   fChain->SetBranchAddress("rotX_err_5", &rotX_err_5, &b_rotX_err_5);
   fChain->SetBranchAddress("rotX_6", &rotX_6, &b_rotX_6);
   fChain->SetBranchAddress("rotX_err_6", &rotX_err_6, &b_rotX_err_6);
   fChain->SetBranchAddress("rotX_7", &rotX_7, &b_rotX_7);
   fChain->SetBranchAddress("rotX_err_7", &rotX_err_7, &b_rotX_err_7);
   fChain->SetBranchAddress("rotY_1", &rotY_1, &b_rotY_1);
   fChain->SetBranchAddress("rotY_err_1", &rotY_err_1, &b_rotY_err_1);
   fChain->SetBranchAddress("rotY_2", &rotY_2, &b_rotY_2);
   fChain->SetBranchAddress("rotY_err_2", &rotY_err_2, &b_rotY_err_2);
   fChain->SetBranchAddress("rotY_3", &rotY_3, &b_rotY_3);
   fChain->SetBranchAddress("rotY_err_3", &rotY_err_3, &b_rotY_err_3);
   fChain->SetBranchAddress("rotY_4", &rotY_4, &b_rotY_4);
   fChain->SetBranchAddress("rotY_err_4", &rotY_err_4, &b_rotY_err_4);
   fChain->SetBranchAddress("rotY_5", &rotY_5, &b_rotY_5);
   fChain->SetBranchAddress("rotY_err_5", &rotY_err_5, &b_rotY_err_5);
   fChain->SetBranchAddress("rotY_6", &rotY_6, &b_rotY_6);
   fChain->SetBranchAddress("rotY_err_6", &rotY_err_6, &b_rotY_err_6);
   fChain->SetBranchAddress("rotY_7", &rotY_7, &b_rotY_7);
   fChain->SetBranchAddress("rotY_err_7", &rotY_err_7, &b_rotY_err_7);
   fChain->SetBranchAddress("rotZ_1", &rotZ_1, &b_rotZ_1);
   fChain->SetBranchAddress("rotZ_err_1", &rotZ_err_1, &b_rotZ_err_1);
   fChain->SetBranchAddress("rotZ_2", &rotZ_2, &b_rotZ_2);
   fChain->SetBranchAddress("rotZ_err_2", &rotZ_err_2, &b_rotZ_err_2);
   fChain->SetBranchAddress("rotZ_3", &rotZ_3, &b_rotZ_3);
   fChain->SetBranchAddress("rotZ_err_3", &rotZ_err_3, &b_rotZ_err_3);
   fChain->SetBranchAddress("rotZ_4", &rotZ_4, &b_rotZ_4);
   fChain->SetBranchAddress("rotZ_err_4", &rotZ_err_4, &b_rotZ_err_4);
   fChain->SetBranchAddress("rotZ_5", &rotZ_5, &b_rotZ_5);
   fChain->SetBranchAddress("rotZ_err_5", &rotZ_err_5, &b_rotZ_err_5);
   fChain->SetBranchAddress("rotZ_6", &rotZ_6, &b_rotZ_6);
   fChain->SetBranchAddress("rotZ_err_6", &rotZ_err_6, &b_rotZ_err_6);
   fChain->SetBranchAddress("rotZ_7", &rotZ_7, &b_rotZ_7);
   fChain->SetBranchAddress("rotZ_err_7", &rotZ_err_7, &b_rotZ_err_7);
   Notify();
}

inline Bool_t AlignmentTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

inline void AlignmentTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

inline Int_t AlignmentTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

