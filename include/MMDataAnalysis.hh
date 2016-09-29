///
///  \file   MMDataAnalysis.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMDataAnalysis_HH
#define MMDataAnalysis_HH

#include "include/MMDataBase.hh"
#include "include/SCEventHits.hh"
#include "include/MMEventHits.hh"

class MMDataAnalysis : public MMDataBase {

public:
  MMDataAnalysis(TTree *tree=0);
  ~MMDataAnalysis();

  virtual Int_t GetNEntries();
  virtual Int_t GetEntry(Long64_t entry);

  MMEventHits mm_EventHits;
  SCEventHits sc_EventHits;

private:
  int m_Nentry;
  
};

#endif

inline MMDataAnalysis::MMDataAnalysis(TTree *tree)
  : MMDataBase(tree)
{
  if(tree)
    m_Nentry = tree->GetEntries();
  else
    m_Nentry = 0;
}
  
inline MMDataAnalysis::~MMDataAnalysis() {}

inline Int_t MMDataAnalysis::GetNEntries(){
  return m_Nentry;
}

inline Int_t MMDataAnalysis::GetEntry(Long64_t entry){
  if(entry < 0 || entry >= m_Nentry)
    return false;

  int ret = MMDataBase::GetEntry(entry);

  // clear previous event micromega hits;
  mm_EventHits = MMEventHits();

  // fill event hits
  for(int i = 0; i < N_mm; i++){
    MMHit hit(mm_MMFE8->at(i),
	      mm_VMM->at(i),
	      mm_CH->at(i));
    hit.SetPDO(mm_PDO->at(i));
    hit.SetTDO(mm_TDO->at(i));
    hit.SetBCID(mm_BCID->at(i));
    hit.SetFIFOcount(mm_FIFOcount->at(i));

    mm_EventHits += hit;
  }

  // clear previous event scintillator hits
  sc_EventHits = SCEventHits();

  // fill event hits
  for(int i = 0; i < N_sci; i++){
    SCHit hit(sci_CH->at(i),
	      sci_count->at(i));

    sc_EventHits += hit;
  }

  return ret;
}
