///
///  \file   MMEfficiency.hh
///
///  \author Gray Putnam
///          (graylcputnam@gmail.com)
///
///  \date   2016 Oct
///


#ifndef MMEfficiency_HH
#define MMEfficiency_HH

#include "include/MMClusterList.hh"
#include "include/GeoOctuplet.hh"
#include "include/SimpleTrackFitter.hh"

class MMEfficiency {
public:
  MMEfficiency(GeoOctuplet *geo, SimpleTrackFitter *fitter, double tol);
  ~MMEfficiency();

  void UpdateEfficiency(const MMClusterList &clusters); 
  double GetBoardEfficiency(int i) const;
  void GetEfficiency(double *ret) const;
  double GetMMFE8Efficiency(int mmfe8) const;

private:
  GeoOctuplet *m_geo;
  SimpleTrackFitter *m_fitter;
  int *m_plane_hits;
  int *m_plane_misses;
  double m_tol;
};

#endif

inline MMEfficiency::MMEfficiency(GeoOctuplet *geo, SimpleTrackFitter *fitter, double tol = 2.5) {
  m_geo = geo;
  m_fitter = fitter;
  m_tol = tol;

  int n = m_geo->GetNPlanes();
  m_plane_hits = new int[n];
  m_plane_misses = new int[n];
  for (int i = 0; i < n; i ++) {
    m_plane_hits[i] = 0;
    m_plane_misses[i] = 0;
  }
}

inline MMEfficiency::~MMEfficiency() {
  delete[] m_plane_hits;
  delete[] m_plane_misses;
}

void MMEfficiency::UpdateEfficiency(const MMClusterList &clusters) {
  int Nplanes = m_geo->GetNPlanes();
  int Nclusters = clusters.GetNCluster();

  MMClusterList temp_clusters;
  for (int ip = 0; ip < Nplanes; ip ++) {
    int mmfe8 = m_geo->MMFE8(ip);
    for (int j = 0; j < Nclusters; j ++) {
      if (clusters[j].MMFE8() != mmfe8) {
        temp_clusters.AddCluster( clusters[j] );
      }
    }   
    MMTrack track = m_fitter->Fit(temp_clusters, *m_geo);

    bool ishit = false;
    for (int j = 0; j < Nclusters; j ++) {
      if (clusters[j].MMFE8() == mmfe8) {
        if ( fabs(m_geo->GetResidualX( clusters[j], track)) < m_tol) {
          ishit = true;
          break;
        }
      }
    }
    ishit ? m_plane_hits[ip] ++ : m_plane_misses[ip] ++;
    temp_clusters.ClearClusters();
  }
}

inline double MMEfficiency::GetMMFE8Efficiency(int mmfe8) const {
    return GetBoardEfficiency(m_geo->Index(mmfe8));
}

inline double MMEfficiency::GetBoardEfficiency(int i) const {
    if (m_plane_hits[i] == 0 && m_plane_misses[i] == 0) 
        return 0;
    return double(m_plane_hits[i]) / double( m_plane_hits[i] + m_plane_misses[i]);
}

void MMEfficiency::GetEfficiency(double *ret) const {
  int NPlanes = m_geo->GetNPlanes();
  for (int ip = 0; ip < NPlanes; ip++) {
    ret[ip] = GetBoardEfficiency(ip);
  }
}

