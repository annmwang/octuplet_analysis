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
#include "include/MMEventHits.hh"

class MMDataAnalysis : public MMDataBase {

public:
  MMDataAnalysis(TTree *tree=0);
  ~MMDataAnalysis();

  virtual Int_t GetEntry(Long64_t entry);

  MMEventHits mm_EventHits;

private:
  int m_Nentry;
  
};

#endif

inline MMDataAnalysis::MMDataAnalysis(TTree *tree=0)
  : MMDataBase(tree)
{
  if(tree)
    m_Nentry = tree->GetEntries();
  else
    m_Nentry = 0;
}
  
inline MMDataAnalysis::~MMDataAnalysis() {}

inline Int_t MMDataAnalysis::GetEntry(Long64_t entry){
  if(entry < 0 || entry >= m_Nentry)
    return false;

  MMDataBase::GetEntry(entry);

  // clear previous event micromega hits;
  mm_EventHits = MMEventHits();

  
}

inline bool MMEventHits::AddHit(const MMHit& hit){
  int Nboard = GetNBoards();
  for(int i = 0; i < Nboard; i++)
    if(m_boards.AddHit(hit))
      return true;

  m_boards.push_back(MMFE8Hits(hit));
  return true;
}

inline bool MMEventHits::AddLinkedHit(const MMLinkedHit& hit){
  int Nboard = GetNBoards();
  for(int i = 0; i < Nboard; i++)
    if(m_boards.AddLinkedHit(hit))
      return true;
  
  m_boards.push_back(MMFE8Hits(hit));
  return true;
}

inline bool MMEventHits::AddHits(const MMFE8Hits& hits){
  int Nboard = GetNBoards();
  for(int i = 0; i < Nboard; i++)
    if(m_boards[i].IsSameMMFE8(hits))
      if(!m_boards[i].AddHits(hits))
	 return true;
	 
  m_boards.push_back(hits);
  return true;
}

inline bool MMEventHits::operator += (const MMHit& hit){
  return AddHit(hit);
}
inline bool MMEventHits::operator += (const MMFE8Hits& hits){
  return AddHits(hits);
}
  
inline int MMEventHits::MMFE8(int iboard) const {
  if(iboard < 0 || iboard >= GetNBoards())
    return -1;

  return m_boards[0].MMFE8();
}

inline int MMEventHits::GetNBoards() const {
  return int(m_boards.size());
}

inline const MMFE8Hits* MMEventHits::Get(int iboard) const {
  if(iboard < 0 || iboard >= GetNBoards())
    return nullptr;

  return &m_boards[iboard];
}

inline const MMFE8Hits* MMEventHits::operator [] (int iboard) const {
  return Get(ihit);
}

inline int MMEventHits::GetNDuplicates() const {
  int Ndup = 0;
  int N = GetNBoards();
  for(int i = 0; i < N; i++)
    Ndup += m_boards[i].GetNDuplicates();
  
  return Ndup;
}

inline MMEventHits MMEventHits::GetDuplicates() const {
  MMEventHits dups;
  int N = GetNBoards();
  for(int i = 0; i < N; i++)
    dups += m_boards[i].GetDuplicates();
  
  return dups;
}

inline MMEventHits operator + (const MMEventHits& evt_hits, 
			       const MMFE8Hits& hits){
  MMEventHits ret(evt_hits);
  ret += hits;
  return ret;
}

inline MMEventHits operator + (const MMEventHits& evt_hits, 
			       const MMHit& hit){
  MMEventHits ret(evt_hits);
  ret += hit;
  return ret;
}

inline MMEventHits operator + (const MMFE8Hits& hits, 
			       const MMEventHits& evt_hits){
  MMEventHits ret(evt_hits);
  ret += hits;
  return ret;
}

inline MMEventHits operator + (const MMHit& hit, 
			       const MMEventHits& evt_hits){
  MMEventHits ret(evt_hits);
  ret += hit;
  return ret;
}
