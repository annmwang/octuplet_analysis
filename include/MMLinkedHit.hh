///
///  \file   MMLinkedHit.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMLinkedHit_HH
#define MMLinkedHit_HH

#include "include/MMHit.hh"

class MMLinkedHit : public MMHit {

public:
  MMLinkedHit(const MMHit& hit);
  MMLinkedHit(const MMLinkedHit& hit);
  
  ~MMLinkedHit();

  int GetNHits() const;

  const MMLinkedHit* GetNext() const;

  void AddHit(const MMHit& hit);
  void AddLinkedHit(const MMHit& hit);
  
private:
  MMLinkedHit* m_next;

};

#endif

inline MMLinkedHit::MMLinkedHit(const MMHit& hit)
  : MMHit(hit)
{
  m_next = nullptr;
}

inline MMLinkedHit::MMLinkedHit(const MMLinkedHit& hit)
  : MMHit(hit)
{
  if(hit.GetNHits() > 1)
    m_next = new MMLinkedHit(*hit.GetNext());
  else
    m_next = nullptr;
}
  
inline MMLinkedHit::~MMLinkedHit(){
  if(m_next)
    delete m_next;
}

inline int MMLinkedHit::GetNHits() const {
  if(m_next)
    return m_next->GetNHits()+1;
  return 1;
}

inline const MMLinkedHit* MMLinkedHit::GetNext() const {
  return m_next;
}

inline void MMLinkedHit::AddHit(const MMHit& hit) {
  if(m_next)
    m_next->AddHit(hit);
  else
    m_next = new MMLinkedHit(hit);
}

inline void MMLinkedHit::AddLinkedHit(const MMLinkedHit& hit) {
  if(m_next)
    m_next->AddLinkedHit(hit);
  else
    m_next = new MMLinkedHit(hit);
}
