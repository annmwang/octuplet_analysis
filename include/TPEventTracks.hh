///
///  \file   TPEventTracks.hh
///
///  \author AW
///
///  \date   2017 May
///


#ifndef TPEventTracks_HH
#define TPEventTracks_HH

#include "include/TPLinkedTrack.hh"

class TPEventTracks {

public:
  TPEventTracks();
  TPEventTracks(const TPTrack& track);
  
  ~TPEventTracks();

  bool AddTrack(const TPTrack& track);

  bool operator += (const TPTrack& track);
  
  int GetNTrack() const;
  TPLinkedTrack Get() const;
  TPLinkedTrack operator [] (int itrack) const;

  
private:
  std::vector<TPTrack*> m_track;

  friend class PDOToCharge;
  friend class TDOToTime;
};

inline TPEventTracks::TPEventTracks() {}

inline TPEventTracks::TPEventTracks(const TPTrack& track){
  AddTrack(track);
}

inline TPEventTracks::~TPEventTracks(){
  int N = GetNTrack();
  for(int i = 0; i < N; i++)
    delete m_track[i];
}

inline bool TPEventTracks::AddTrack(const TPTrack& track){
  int N = GetNTrack();
  for(int i = 0; i < N; i++){
    if(track.GetNHits() > m_track[i]->GetNHits()){
      m_track.insert(m_track.begin()+i, new TPTrack(track));
      return true;
    }
    if(track.GetNHits() == m_track[i]->GetNHits()){
      m_track.insert(m_track.begin()+i+1, new TPTrack(track));
      return true;
    }
  }
  m_track.push_back(new TPLinkedTrack(track));
  return true;
}


inline bool TPEventTracks::operator += (const TPTrack& track){
  return AddTrack(track);
}

  
inline int TPEventTracks::GetNTrack() const {
  return int(m_track.size());
}

inline TPLinkedTrack TPEventTracks::Get() const {
  TPLinkedTrack tr = TPLinkedTrack(*m_track[0]);
  for (int i = 1; i < GetNTrack(); i++){
    tr.AddTrack(*m_track[i]);
  }
  return tr;
}

inline TPLinkedTrack TPEventTracks::operator [] (int itrack) const {
  return Get();
}

#endif
