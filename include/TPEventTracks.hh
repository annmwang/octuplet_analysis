///
///  \file   TPEventTracks.hh
///
///  \author AW
///
///  \date   2017 May
///


#ifndef TPEventTracks_HH
#define TPEventTracks_HH

#include "include/TPTrack.hh"
#include "include/MMClusterList.hh"

class TPEventTracks {

public:
  TPEventTracks();
  TPEventTracks(const TPTrack& track);
  
  ~TPEventTracks();

  bool AddTrack(const TPTrack& track);

  bool operator += (const TPTrack& track);
  
  int GetNTrack() const;
  size_t size() const;

  void SetSciBCID(int sciBCID, int sciph);
  void SetSciOffset(int offset);

  double SciBCID();
  int SciOffset();

  TPTrack const& Get(int itrack) const;
  TPTrack const& operator [] (int itrack) const;

  std::vector<TPTrack*>::iterator begin();
  std::vector<TPTrack*>::iterator end();

  TPTrack* Highlander(const MMClusterList& clusters, bool scimatch=false, int cut=10);
  double deltaBCID(int BCID);
  
private:
  std::vector<TPTrack*> m_track;

  double m_sciBCID; 
  int m_offset;

  friend class PDOToCharge;
  friend class TDOToTime;
};

inline TPEventTracks::TPEventTracks() {
  m_sciBCID = -1;
  m_offset = -1;
}

inline TPEventTracks::TPEventTracks(const TPTrack& track){
  m_sciBCID = -1;
  m_offset = -1;
  AddTrack(track);
}

inline TPEventTracks::~TPEventTracks(){
  int N = GetNTrack();
  for(int i = 0; i < N; i++)
    delete m_track[i];
}

inline std::vector<TPTrack*>::iterator TPEventTracks::begin(){
  return m_track.begin();
}

inline std::vector<TPTrack*>::iterator TPEventTracks::end(){
  return m_track.end();
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
  m_track.push_back(new TPTrack(track));
  return true;
}


inline bool TPEventTracks::operator += (const TPTrack& track){
  return AddTrack(track);
}

  
inline int TPEventTracks::GetNTrack() const {
  return int(m_track.size());
}

inline size_t TPEventTracks::size() const {
  return m_track.size();
}

inline void TPEventTracks::SetSciBCID(int sciBCID, int sciph) {
  m_sciBCID = sciBCID+sciph/16.;
}

inline void TPEventTracks::SetSciOffset(int offset) {
  m_offset = offset;
}

inline double TPEventTracks::SciBCID() {
  return m_sciBCID;
}

inline int TPEventTracks::SciOffset() {
  return m_offset;
}

inline TPTrack const& TPEventTracks::Get(int itrack) const {
  return *m_track[itrack];
}

inline TPTrack const& TPEventTracks::operator [] (int itrack) const {
  return Get(itrack);
}

inline TPTrack* TPEventTracks::Highlander(const MMClusterList& clusters, bool scimatch, int cut){
  //
  // Choose the TPTrack with most hits in common with MMClusterList provided.
  // If >1 tracks satisfy this, pick the first (theyre sorted by N(hits)).
  //
  if (m_track.size() == 0)
    return nullptr;

  int mostmatch = -1, nmatch = -1;
  int nscimatch = -1;
  if (scimatch && m_sciBCID > 0 && m_offset != -1)
    nscimatch = 0;
  TPTrack* the_one = nullptr;

  for (size_t it = 0; it < m_track.size(); ++it){
    nmatch = 0;
    if (scimatch && m_sciBCID > 0 && m_offset != -1){
      if (fabs(deltaBCID(m_track[it]->BCID())) < cut)
        nscimatch++;
      else
        continue;
    }
    for (size_t ia = 0; ia < m_track[it]->size(); ++ia){
      if (clusters.ContainsTP(m_track[it]->Get(ia))){
        nmatch++;
      }
    }
    m_track[it]->SetNMatch(nmatch);
    if (nmatch > mostmatch){
      mostmatch = nmatch;
      the_one = m_track[it];
    }
  }

  if (nscimatch == 0){
    return nullptr;
  }
  return the_one;
}

inline double TPEventTracks::deltaBCID(int BCID) {
  if ((m_sciBCID - BCID) < 0)
    return m_sciBCID - BCID + m_offset;
  else
    return m_sciBCID - BCID + (m_offset - 4096);
}

#endif
