///
///  \file   CombinatoricTracktFitter.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2017 Jan
///

#ifndef CombinatoricTrackFitter_HH
#define CombinatoricTrackFitter_HH

#include "include/SimpleTrackFitter.hh"

class CombinatoricTrackFitter : public SimpleTrackFitter {

public:
  CombinatoricTrackFitter();
  ~CombinatoricTrackFitter();

  virtual MMTrack Fit(const MMClusterList& clusters, 
		      const GeoOctuplet& geometry,
		      const int evt = -1);
  
private:
  MMClusterList* m_clusters;
  const GeoOctuplet* m_geometry;
};

inline CombinatoricTrackFitter::CombinatoricTrackFitter(){
  m_clusters = nullptr;
  m_geometry = nullptr;
}

inline CombinatoricTrackFitter::~CombinatoricTrackFitter(){
  
}

inline MMTrack CombinatoricTrackFitter::Fit(const MMClusterList& clusters, 
                                      const GeoOctuplet&   geometry,
                                      const int evt){
  // get geometry pointer
  m_geometry = &geometry;

  // get cluster list
  int N = clusters.GetNCluster();
  
  int Nclus[8];
  std::vector<int> iclus[8]; 
  for(int i = 0; i < 8; i++)
    Nclus[i] = 0;

  // loop over each board to find 
  // number of clusters on each board
  for(int i = 0; i < N; i++){
    int index = m_geometry->Index(clusters[i].MMFE8());
    if(index >= 0 && index < 8){
      Nclus[index]++;
      iclus[index].push_back(i);
    }
  }
     
  // number of different combinatorics to check
  int Ncomb = 1;
  for(int i = 0; i < 8; i++)
    if(Nclus[i] > 0)
      Ncomb *= Nclus[i];

  // loop over all combinatorics
  MMTrack min_track;

  for(int c = 0; c < Ncomb; c++){
    m_clusters = new MMClusterList();
    int key = c;
    for(int b = 0; b < 8; b++){
      if(Nclus[b] == 0)
	continue;
      if(Nclus[b] == 1){
	m_clusters->AddCluster(clusters[iclus[b][0]]);
	continue;
      }
      
      m_clusters->AddCluster(clusters[iclus[b][key%Nclus[b]]]);
      key /= Nclus[b];
    }

    MMTrack track = SimpleTrackFitter::Fit(*m_clusters, geometry, evt);
    if(track.ResidualSq() < min_track.ResidualSq() ||
       min_track.IsFit() == false)
      min_track = track;

    delete m_clusters;
    m_clusters = nullptr;
  }

  m_geometry = nullptr;
  m_clusters = nullptr;

  return min_track;
}

#endif
