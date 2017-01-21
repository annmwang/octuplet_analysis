///
///  \file   CombinatoricChi2TracktFitter.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2017 Jan
///

#ifndef CombinatoricChi2TrackFitter_HH
#define CombinatoricChi2TrackFitter_HH

#include "include/SimpleTrackFitter.hh"

class CombinatoricChi2TrackFitter : public SimpleTrackFitter {

public:
  CombinatoricChi2TrackFitter();
  ~CombinatoricChi2TrackFitter();

  virtual MMTrack Fit(const MMClusterList& clusters, 
		      const GeoOctuplet& geometry,
		      const int evt = -1);
  
private:
  MMClusterList* m_clusters;
  const GeoOctuplet* m_geometry;

  MMTrack RecursiveFit(const MMClusterList& clusters, 
		       const GeoOctuplet& geometry,
		       const int evt = -1);
};

inline CombinatoricChi2TrackFitter::CombinatoricChi2TrackFitter(){
  m_clusters = nullptr;
  m_geometry = nullptr;
}

inline CombinatoricChi2TrackFitter::~CombinatoricChi2TrackFitter(){
  
}

inline MMTrack CombinatoricChi2TrackFitter::RecursiveFit(const MMClusterList& clusters, 
							 const GeoOctuplet&   geometry,
							 const int evt){

  // get cluster list
  int N = clusters.GetNCluster();

  MMTrack all_track = SimpleTrackFitter::Fit(clusters, geometry, evt);
 
  if(N <= 5)
    return all_track;

  MMTrack min_track;
  
  for(int i = 0; i < N; i++){
    m_clusters = new MMClusterList();
    // fill cluster list with all hits 
    // except this one
    for(int j = 0; j < N; j++)
      if(j != i)
	m_clusters->AddCluster(clusters[j]);
    
    MMTrack track = RecursiveFit(*m_clusters, geometry, evt);
    if((track.ResidualSq() < min_track.ResidualSq()) ||
       (track.IsFit() == true && min_track.IsFit() == false))
      min_track = track;
       
    delete m_clusters;
    m_clusters = nullptr;
  }

  if(!min_track.IsFit())
    return all_track;

  if(min_track.ResidualSq() / (min_track.NX()+min_track.NU()+min_track.NV()-4.) <
     all_track.ResidualSq() / (all_track.NX()+all_track.NU()+all_track.NV()-4.))
    return min_track;
  else
    return all_track;

}

inline MMTrack CombinatoricChi2TrackFitter::Fit(const MMClusterList& clusters, 
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

  // for(int i = 0; i < 8; i++)
  //   cout << "N " << i << " " << Nclus[i] << endl;

  // cout << "Ncomb = " << Ncomb << endl;

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

    MMTrack track = RecursiveFit(*m_clusters, geometry, evt);

    if(track.IsFit() &&
       ((min_track.ResidualSq() / (min_track.NX()+min_track.NU()+min_track.NV()-4.) >
	 track.ResidualSq() / (track.NX()+track.NU()+track.NV()-4.)) || 
	min_track.IsFit() == false))
       min_track = track;

    delete m_clusters;
    m_clusters = nullptr;
  }

  m_geometry = nullptr;
  m_clusters = nullptr;

  return min_track;
}

#endif
