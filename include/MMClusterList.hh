///
///  \file   MMFE8Hits.hh
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
  
  int GetNCluster() const;
  const MMCluster* Get(int i) const;
  const MMCluster* operator [] (int i) const;

private:
  std::vector<MMCluster> m_clusters;
  
};

#endif

inline MMClusterList::MMClusterList() {}

inline MMClusterList::MMClusterList(const MMCluster& clus){
  AddCluster(clus);
}
  
inline MMClusterList::~MMClusterList() {}

inline double MMCluster::Channel() const {
  double ch = 0;
  int Nhit = GetNHits();
  for(int i = 0; i < Nhit; i++)
    ch += double(Get(i)->Channel())*Get(i)->Charge();
  ch /= Charge();
  return ch;
}

inline double MMCluster::Charge() const {
  double Q = 0;
  int Nhit = GetNHits();
  for(int i = 0; i < Nhit; i++)
    Q += double(Get(i)->Charge());
 
  return Q;
}

inline double MMCluster::Time() const {
  // to do
  return 0;
}
  
inline int MMCluster::NHoles() const {
  int Nhit = GetNHits();
  int clus_size = Get(Nhit-1)->Channel()-Get(0)->Channel()+1;

  return clus_size - Nhit;
}


