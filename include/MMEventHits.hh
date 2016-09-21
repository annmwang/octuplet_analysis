///
///  \file   MMEventHits.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMEventHits_HH
#define MMEventHits_HH

#include "include/MMFE8Hits.hh"

class MMEventHits {

public:
  MMEventHits();
  MMEventHits(const MMHit& hit);
  MMEventHits(const MMFE8Hits& hits);
  MMEventHits(const MMEventHits& event_hits);
  
  ~MMEventHits();

  bool AddHit(const MMHit& hit);
  bool AddHits(const MMFE8Hits& hits);
  bool AddLinkedHit(const MMLinkedHit& hit);

  bool operator += (const MMHit& hit);
  bool operator += (const MMFE8Hits& hits);
  
  int GetNBoards() const;
  int MMFE8(int iboard) const;
  const MMFE8Hits* Get(int iboard) const;
  const MMFE8Hits* operator [] (int iboard) const;

  int GetNDuplicates() const;

  MMEventHits GetDuplicates() const;
  
private:
  vector<MMFE8Hits> m_boards;

};

MMFE8Hits operator + (const MMEventHits& evt_hits, const MMFE8Hits& hits);
MMFE8Hits operator + (const MMEventHits& evt_hits, const MMHit& hit);
MMFE8Hits operator + (const MMFE8Hits& hits, const MMEventHits& evt_hits);
MMFE8Hits operator + (const MMHit& hit, const MMEventHits& evt_hits);

#endif

inline MMEventHits::MMEventHits() {}

inline MMEventHits::MMEventHits(const MMHit& hit){
  AddHit(hit);
}
inline MMEventHits::MMEventHits(const MMFE8Hits& hits){
  AddHits(hits);
}

inline MMEventHits::MMEventHits(const MMEventHits& evt_hits){
  int Nboard = evt_hits.GetNBoards();
  for(int i = 0; i < Nboard; i++)
    AddHits(*evt_hits[i]);
}
  
inline MMEventHits::~MMEventHits() {}

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
