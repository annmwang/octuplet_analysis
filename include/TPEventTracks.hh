///
///  \file   TPEventTracks.hh
///
///  \author Ann Miao Wang
///
///  \date   2017 May
///


#ifndef TPEventTracks_HH
#define TPEventTracks_HH

#include "include/TPTrack.hh"

class TPEventTracks : public TPTrack {

public:
  TPEventTracks();
  TPEventTracks(const TPTrack& track);
  TPEventTracks(const TPEventTracks& track);
  
  ~TPEventTracks();

  int GetNTracks() const;

  void AddTrack(const TPTrack& track);
  void AddLinkedTrack(const TPEventTracks& track);

  const TPEventTracks* GetNext() const;
  
private:
  TPEventTracks* m_next;

};

#endif

inline TPEventTracks::TPEventTracks()
{
  m_next = nullptr;
}

inline TPEventTracks::TPEventTracks(const TPTrack& track)
  : TPTrack(track)
{
  m_next = nullptr;
}

inline TPEventTracks::TPEventTracks(const TPEventTracks& track)
  : TPTrack(track)
{
  if(track.GetNTracks() > 1)
    m_next = new TPEventTracks(*track.GetNext());
  else
    m_next = nullptr;
}
  
inline TPEventTracks::~TPEventTracks(){
  if(m_next){
    delete m_next;
  }
}

inline int TPEventTracks::GetNTracks() const {
  if(m_next)
    return m_next->GetNTracks()+1;
  return 1;
}

inline const TPEventTracks* TPEventTracks::GetNext() const {
  return m_next;
}

inline void TPEventTracks::AddTrack(const TPTrack& track) {
  if (m_next)
    m_next->AddTrack(track);
  else
    m_next = new TPEventTracks(track);
}
inline void TPEventTracks::AddLinkedTrack(const TPEventTracks& track) {
  if(m_next)
    m_next->AddLinkedTrack(track);
  else
    m_next = new TPEventTracks(track);
}

