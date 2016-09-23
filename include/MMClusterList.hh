///
///  \file   MMClusterList.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMClusterList_HH
#define MMClusterList_HH

#include "include/MMCluster.hh"

class MMClusterList {

public:
  MMClusterList();
  MMClusterList(const MMCluster& clus);
  
  ~MMClusterList();

  void AddCluster(const MMCluster& clus);
  void AddHit(const MMHit& hit, int iclus);
  void AddLinkedHit(const MMLinkedHit& hit, int iclus);
  
  int GetNCluster() const;
  MMCluster const& Get(int i) const;
  MMCluster const& operator [] (int i) const;

  bool Contains(const MMHit& hit) const;

private:
  std::vector<MMCluster> m_clusters;
  
};

inline MMClusterList::MMClusterList() {}

inline MMClusterList::MMClusterList(const MMCluster& clus){
  AddCluster(clus);
}
  
inline MMClusterList::~MMClusterList() {}

inline void MMClusterList::AddCluster(const MMCluster& clus){
  int N = GetNCluster();
  for(int i = 0; i < N; i++){
    if(clus.Charge() < m_clusters[i].Charge()){
      m_clusters.insert(m_clusters.begin()+i, clus);
      return;
    }
  }
  m_clusters.push_back(clus);
}

inline void MMClusterList::AddHit(const MMHit& hit, int iclus){
  if(iclus < 0 || iclus >= GetNCluster())
    return;

  m_clusters[iclus].AddHit(hit);
}

inline void MMClusterList::AddLinkedHit(const MMLinkedHit& hit, int iclus){
  if(iclus < 0 || iclus >= GetNCluster())
    return;

  m_clusters[iclus].AddLinkedHit(hit);
}
  
inline int MMClusterList::GetNCluster() const {
  return m_clusters.size();
}

inline MMCluster const& MMClusterList::Get(int i) const {
  return m_clusters[i];
}

inline MMCluster const& MMClusterList::operator [] (int i) const {
  return Get(i);
}

inline bool MMClusterList::Contains(const MMHit& hit) const {
  int Nclus = GetNCluster();
  for(int i = 0; i < Nclus; i++)
    if(Get(i).Contains(hit))
      return true;
  return false;
}

#endif




