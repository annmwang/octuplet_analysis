///
///  \file   MMCluster.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMCluster_HH
#define MMCluster_HH

#include "include/MMFE8Hits.hh"

class MMCluster : public MMFE8Hits {

public:
  MMCluster();
  MMCluster(const MMHit& hit);
  MMCluster(const MMFE8Hits& hits);
  MMCluster(const MMLinkedHit& hit);
  MMCluster(const MMCluster& clus);
  
  ~MMCluster();
  
  double Channel() const;
  double Charge() const;
  double Time() const;
  
  int NHoles() const;
};

#endif

inline MMCluster::MMCluster()
  : MMFE8Hits() {}

inline MMCluster::MMCluster(const MMHit& hit)
  : MMFE8Hits(hit) {}

inline MMCluster::MMCluster(const MMFE8Hits& hits)
  : MMFE8Hits(hits) {}

inline MMCluster::MMCluster(const MMLinkedHit& hit)
  : MMFE8Hits(hit) {}

inline MMCluster::MMCluster(const MMCluster& clus)
  : MMFE8Hits(clus) {}
  
inline MMCluster::~MMCluster() {}

inline double MMCluster::Channel() const {
  double ch = 0;
  int Nhit = GetNHits();
  for(int i = 0; i < Nhit; i++)
    ch += double(Get(i).Channel())*Get(i).Charge();
  ch /= Charge();
  return ch;
}

inline double MMCluster::Charge() const {
  double Q = 0;
  int Nhit = GetNHits();
  for(int i = 0; i < Nhit; i++)
    Q += double(Get(i).Charge());
 
  return Q;
}

inline double MMCluster::Time() const {
  // to do
  return 0;
}
  
inline int MMCluster::NHoles() const {
  int Nhit = GetNHits();
  int clus_size = Get(Nhit-1).Channel()-Get(0).Channel()+1;

  return clus_size - Nhit;
}


