///
///  \file   HighQTracktFitter.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2017 Jan
///

#ifndef HighQTrackFitter_HH
#define HighQTrackFitter_HH

#include "include/SimpleTrackFitter.hh"

class HighQTrackFitter : public SimpleTrackFitter {

public:
  HighQTrackFitter();
  ~HighQTrackFitter();

  virtual MMTrack Fit(const MMClusterList& clusters, 
		      const GeoOctuplet& geometry,
		      const int evt = -1);
  
private:
  MMClusterList* m_clusters;
  const GeoOctuplet* m_geometry;
};

inline HighQTrackFitter::HighQTrackFitter(){
  m_clusters = nullptr;
  m_geometry = nullptr;
}

inline HighQTrackFitter::~HighQTrackFitter(){
  
}

inline MMTrack HighQTrackFitter::Fit(const MMClusterList& clusters, 
                                      const GeoOctuplet&   geometry,
                                      const int evt){
  // get geometry pointer
  m_geometry = &geometry;

  // get cluster list
  int N = clusters.GetNCluster();
  m_clusters = new MMClusterList();
  
  // loop over each board to find 
  // highest charge cluster
  for(int ib = 0; ib < 8; ib++){
    double Qmax = -1;
    int imax = -1;
    for(int i = 0; i < N; i++){
      if(m_geometry->Index(clusters[i].MMFE8()) == ib){
	double Q = clusters[i].Charge();
	if(Q > Qmax){
	  Qmax = Q;
	  imax = i;
	}
      }
    }
    if(imax >= 0)
      m_clusters->AddCluster(clusters[imax]);
  }

  MMTrack track = SimpleTrackFitter::Fit(*m_clusters, geometry, evt);

  m_geometry = nullptr;
  delete m_clusters;
  m_clusters = nullptr;

  return track;
}

#endif
