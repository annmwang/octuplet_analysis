///
///  \file   TPLinkedTrack.hh
///
///  \author Ann Miao Wang
///
///  \date   2017 May
///


#ifndef TPLinkedTrack_HH
#define TPLinkedTrack_HH

#include "include/TPTrack.hh"

class TPLinkedTrack : public TPTrack {

public:
  TPLinkedTrack(const TPTrack& track);
  TPLinkedTrack(const TPLinkedTrack& track);
  
  ~TPLinkedTrack();

  int GetNTracks() const;

  void AddTrack(const TPTrack& track);
  void AddLinkedTrack(const TPLinkedTrack& track);

  const TPLinkedTrack* GetNext() const;
  
private:
  TPLinkedTrack* m_next;

};

#endif

inline TPLinkedTrack::TPLinkedTrack(const TPTrack& track)
  : TPTrack(track)
{
  m_next = nullptr;
}

inline TPLinkedTrack::TPLinkedTrack(const TPLinkedTrack& track)
  : TPTrack(track)
{
  if(track.GetNTracks() > 1)
    m_next = new TPLinkedTrack(*track.GetNext());
  else
    m_next = nullptr;
}
  
inline TPLinkedTrack::~TPLinkedTrack(){
  if(m_next){
    delete m_next;
  }
}

inline int TPLinkedTrack::GetNTracks() const {
  if(m_next)
    return m_next->GetNTracks()+1;
  return 1;
}

inline const TPLinkedTrack* TPLinkedTrack::GetNext() const {
  return m_next;
}

inline void TPLinkedTrack::AddTrack(const TPTrack& track) {
  if (m_next)
    m_next->AddTrack(track);
  else
    m_next = new TPLinkedTrack(track);
}
inline void TPLinkedTrack::AddLinkedTrack(const TPLinkedTrack& track) {
  if(m_next)
    m_next->AddLinkedTrack(track);
  else
    m_next = new TPLinkedTrack(track);
}

