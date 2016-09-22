///
///  \file   MMFE8Hits.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMFE8Hits_HH
#define MMFE8Hits_HH

#include "include/MMLinkedHit.hh"

class MMFE8Hits {

public:
  MMFE8Hits();
  MMFE8Hits(const MMHit& hit);
  MMFE8Hits(const MMFE8Hits& hits);
  MMFE8Hits(const MMLinkedHit& hit);
  
  ~MMFE8Hits();

  bool IsSameMMFE8(const MMHit& hit) const;
  bool IsSameMMFE8(const MMFE8Hits& hits) const;

  bool AddHit(const MMHit& hit);
  bool AddHits(const MMFE8Hits& hits);
  bool AddLinkedHit(const MMLinkedHit& hit);

  bool operator += (const MMHit& hit);
  bool operator += (const MMFE8Hits& hits);
  
  int MMFE8() const;

  int GetNHits() const;
  const MMLinkedHit* Get(int ihit) const;
  const MMLinkedHit* operator [] (int ihit) const;

  int GetNDuplicates() const;

  MMFE8Hits GetDuplicates() const;
  
private:
  vector<MMLinkedHit> m_hits;

  friend class PDOToCharge;
  friend class TDOToTime;
};

MMFE8Hits operator + (const MMFE8Hits& hits_a, const MMFE8Hits& hits_b);
MMFE8Hits operator + (const MMFE8Hits& hits, const MMHit& hit);
MMFE8Hits operator + (const MMHit& hit, const MMFE8Hits& hits);
MMFE8Hits operator + (const MMHit& hit_a, const MMHit& hit_b);

#endif

inline MMFE8Hits::MMFE8Hits(){ }

inline MMFE8Hits::MMFE8Hits(const MMHit& hit){
  AddHit(hit);
}

inline MMFE8Hits::MMFE8Hits(const MMFE8Hits& hits){
  AddHits(hits);
}

inline MMFE8Hits::MMFE8Hits(const MMLinkedHit& hit){
  AddLinkedHit(hit);
}
  
inline MMFE8Hits::~MMFE8Hits(){ }

inline bool MMFE8Hits::IsSameMMFE8(const MMHit& hit) const {
  if(GetNHits() == 0)
    return true;

  if(MMFE8() == hit.MMFE8() && hit.MMFE8() > 0)
    return true;
  return false;
}

inline bool MMFE8Hits::IsSameMMFE8(const MMFE8Hits& hits) const {
  if(GetNHits() == 0)
    return true;
  
  if(MMFE8() == hits.MMFE8() && hits.MMFE8() > 0)
    return true;
  return false;
}

inline bool MMFE8Hits::AddHit(const MMHit& hit){
  if(IsSameMMFE8(hit)){
    int N = GetNHits();
    for(int i = 0; i < N; i++){
      if(hit.Channel() == m_hits[i].Channel()){
	m_hits[i].AddHit(hit);
	return true;
      }	
    }
    m_hits.push_back(MMLinkedHit(hit));
    return true;
  }
  return false;
}

inline bool MMFE8Hits::AddLinkedHit(const MMLinkedHit& hit){
  if(IsSameMMFE8(hit)){
    int N = GetNHits();
    for(int i = 0; i < N; i++){
      if(hit.Channel() == m_hits[i].Channel()){
	m_hits[i].AddLinkedHit(hit);
	return true;
      }	
    }
    m_hits.push_back(MMLinkedHit(hit));
    return true;
  }
  return false;
}

inline bool MMFE8Hits::AddHits(const MMFE8Hits& hits){
  int N = hits.GetNHits();
  double ret = true;
  for(int i = 0; i < N; i++)
    ret = AddLinkedHit(*hits[i]) && ret;
  return ret;
}

inline bool MMFE8Hits::operator += (const MMHit& hit){
  return AddHit(hit);
}
inline bool MMFE8Hits::operator += (const MMFE8Hits& hits){
  return AddHits(hits);
}
  
inline int MMFE8Hits::MMFE8() const {
  if(GetNHits() == 0)
    return -1;

  return m_hits[0].MMFE8();
}

inline int MMFE8Hits::GetNHits() const {
  return int(m_hits.size());
}

inline const MMLinkedHit* MMFE8Hits::Get(int ihit) const {
  if(ihit < 0 || ihit >= GetNHits())
    return nullptr;

  return &m_hits[ihit];
}

inline const MMLinkedHit* MMFE8Hits::operator [] (int ihit) const {
  return Get(ihit);
}

inline int MMFE8Hits::GetNDuplicates() const {
  int Ndup = 0;
  int N = GetNHits();
  for(int i = 0; i < N; i++)
    if(m_hits[i].GetNHits() > 1)
      Ndup++;
  return Ndup;
}
inline MMFE8Hits MMFE8Hits::GetDuplicates() const {
  MMFE8Hits dups;
  int N = GetNHits();
  for(int i = 0; i < N; i++)
    if(m_hits[i].GetNHits() > 1)
      dups.AddLinkedHit(m_hits[i]);
  return dups;
}

inline MMFE8Hits operator + (const MMFE8Hits& hits_a,
			     const MMFE8Hits& hits_b){
  MMFE8Hits ret(hits_a);
  ret += hits_b;
  return ret;
}

inline MMFE8Hits operator + (const MMFE8Hits& hits,
			     const MMHit& hit){
  MMFE8Hits ret(hits);
  ret += hit;
  return ret;
}

inline MMFE8Hits operator + (const MMHit& hit,
			     const MMFE8Hits& hits){
  MMFE8Hits ret(hit);
  ret += hits;
  return ret;
}

inline MMFE8Hits operator + (const MMHit& hit_a,
			     const MMHit& hit_b){
  MMFE8Hits ret(hit_a);
  ret += hit_b;
  return ret;
}
