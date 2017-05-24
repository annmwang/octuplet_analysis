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
#include "include/TPEventTracks.hh"

class MMDataAnalysis : public MMDataBase {

public:
  MMDataAnalysis(TTree *tree=0);
  ~MMDataAnalysis();

  virtual Int_t GetNEntries();
  virtual Int_t GetEntry(Long64_t entry);

  MMEventHits mm_EventHits;
  SCEventHits sc_EventHits;
  TPEventTracks tp_EventTracks;

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
	      mm_CH->at(i),
	      RunNum);
    hit.SetPDO(mm_PDO->at(i));
    hit.SetTDO(mm_TDO->at(i));
    hit.SetBCID(mm_BCID->at(i));
    hit.SetFIFOcount(mm_FIFOcount->at(i));

    mm_EventHits += hit;
  }

  // clear previous event scintillator hits
  sc_EventHits = SCEventHits(RunNum);

  // fill event hits
  for(int i = 0; i < N_sci-1; i++){
   
    SCHit hit(sci_CH->at(i),
	      sci_count->at(i));
    sc_EventHits += hit;
    
  }
  sc_EventHits.SetTPBCID(tpsci_BCID);
  if (RunNum >= 3522)
    sc_EventHits.SetTPph(tpsci_ph);
  if (RunNum >= 3522){
  tp_EventTracks = TPEventTracks();
  // fill event tracks
  for (int i = 0; i < tp_hit_n->size(); i++){
    TPTrack tp_track = TPTrack();
    for(int j = 0; j < tp_hit_n->at(i); j++){
      TPHit thit(tp_hit_MMFE8->at(i)[j],
                 tp_hit_VMM->at(i)[j],
                 tp_hit_CH->at(i)[j],
                 tp_hit_BCID->at(i)[j],
                 RunNum);
      tp_track += thit;
    }
    tp_track.SetMxLocal(tp_mxlocal->at(i));
    tp_track.SetBCID(tp_BCID->at(i));
    
    tp_EventTracks.AddTrack(tp_track);
  }
  }
  // }
  return ret;
}
