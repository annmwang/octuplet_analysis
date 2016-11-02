///
///  \file   ChargeTrackFitter.hh
///
///  \author Gray Putnam
///          (crogan@cern.ch)
///
///  \date   2016 Nov
///

#ifndef ChargeTrackFitter_HH
#define ChargeTrackFitter_HH

#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

#include "include/MMClusterList.hh"
#include "include/GeoOctuplet.hh"

class ChargeTrackFitter {

public:
  ChargeTrackFitter();
  ~ChargeTrackFitter();

  MMTrack Fit(const MMClusterList clusters, const GeoOctuplet& geometry, double total_charge);
  bool IsGoodCluster(const MMClusterList& clusters) const;

private:
  ROOT::Math::Minimizer* m_minimizer;
  ROOT::Math::Functor* m_functor;
  double m_total_charge;

  MMClusterList* m_clusters;
  const GeoOctuplet*   m_geometry;
  double EvaluateMetric(const double* param);
  MMTrack _DoFit(const MMClusterList& fit_clusters);

};

inline ChargeTrackFitter::ChargeTrackFitter(){
  m_minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
  m_minimizer->SetMaxFunctionCalls(10000000);
  m_minimizer->SetMaxIterations(100000);
  m_minimizer->SetTolerance(0.001);
  m_minimizer->SetPrintLevel(0);
  
  m_functor = new ROOT::Math::Functor(this, &ChargeTrackFitter::EvaluateMetric, 4);
  m_minimizer->SetFunction(*m_functor);

  m_minimizer->SetVariable(0, "c_x", 0., 0.001);
  m_minimizer->SetVariable(1, "s_x", 0., 0.001);
  m_minimizer->SetVariable(2, "c_y", 0., 0.001);
  m_minimizer->SetVariable(3, "s_y", 0., 0.001);

  m_clusters = nullptr;
  m_geometry = nullptr;
}

inline ChargeTrackFitter::~ChargeTrackFitter(){
  delete m_minimizer;
  delete m_functor;
}

inline MMTrack ChargeTrackFitter::Fit(const MMClusterList clusters, const GeoOctuplet& geometry, double total_charge){
  // return track
  MMTrack track;
  MMClusterList *track_clusters = new MMClusterList();
  bool track_initialized = false;

  // get geometry pointer
  m_geometry = &geometry;
  
  m_total_charge = total_charge;

  // get cluster list
  int N = clusters.GetNCluster();
  m_clusters = new MMClusterList();
  for(int i = 0; i < N; i++)
    if(m_geometry->Index(clusters[i].MMFE8()) >= 0)
      m_clusters->AddCluster(clusters[i]);

  // do fits
  // try 1 track
  // try all combinatorics

  // each cluster can be in one of two states in the track or not
  // so in total we have 2**Nclusters tracks here
  // we really should have less than 64 clusters, so I'm not going to worry about integer overflow
  unsigned NClusters = (unsigned)m_clusters->GetNCluster();

  for (unsigned i = 0; i < 1 << NClusters; i++) {
    // the jth cluster is included if the jth bit is 1
    // this will capture all the possibilities
    MMClusterList *fit_clusters = new MMClusterList();

    for (unsigned j = 1; j < NClusters; j++)
      if (i & (1 << j)) {
        fit_clusters->AddCluster( m_clusters->Get(j) );
      }
    
    if (IsGoodCluster(*fit_clusters)) {
      MMTrack temp_track = _DoFit(*fit_clusters);
      if (!track_initialized || temp_track.GetFitScore() < track.GetFitScore()) {
        track = temp_track;
        //cout << track.SlopeX() << track.SlopeY() << track.ConstX() << track.ConstY() << endl;
        delete track_clusters;
        track_clusters = fit_clusters;
        track_initialized = true;
      }
      else
        delete fit_clusters;
    }
    else
      delete fit_clusters;
  }
  m_geometry = nullptr;
  delete m_clusters;
  m_clusters = nullptr;
  total_charge = 0;

  track.SetClusters( *track_clusters );
  delete track_clusters;
  //cout << track.SlopeX() << " "<< track.SlopeY() << " "<< track.ConstX() << " " << track.ConstY() << endl;
  return track;
}

inline bool ChargeTrackFitter::IsGoodCluster(const MMClusterList& clusters) const {
  vector<int> planes;
  vector<int> uv_planes;

  int nClusters = clusters.GetNCluster();
  for (int i = 0; i < nClusters; i++) {
    MMCluster clus = clusters.Get(i);
    int MMFE8 = clus.MMFE8();
    GeoPlane plane = m_geometry->Get( m_geometry->Index(MMFE8) );
    bool is_uv = fabs(plane.LocalXatYend(32) - plane.LocalXatYbegin(32)) > 1e-2;
    if (is_uv) {
        int size = uv_planes.size();
        int j;
        for (j = 0; j < size; j ++)
            if (uv_planes[j] == MMFE8)
                return false;
        uv_planes.push_back(MMFE8);
    }
    else {
        int size = planes.size();
        int j;
        for (j = 0; j < size; j ++)
            if (planes[j] == MMFE8)
                return false;
        planes.push_back(MMFE8);
    }
  }

  return (uv_planes.size() > 1 && planes.size() > 2);
}

inline MMTrack ChargeTrackFitter::_DoFit(const MMClusterList& fit_clusters) {
  MMTrack track;
  m_minimizer->SetVariableValue(0, 0.);
  m_minimizer->SetVariableValue(1, 0.);
  m_minimizer->SetVariableValue(2, 0.);
  m_minimizer->SetVariableValue(3, 0.);

  m_minimizer->Minimize();

  // get the line score
  double cluster_charge = 0;
  int Nclus = m_clusters->GetNCluster();
  for (int i = 0; i < Nclus; i ++)
    cluster_charge += m_clusters->Get(i).Charge();
  double score = m_minimizer->MinValue() + pow(m_total_charge - cluster_charge, 2); 
  track.SetFitScore(score);
 
  const double* param = m_minimizer->X();
  track.SetConstX(param[0]);
  track.SetConstY(param[2]);
  track.SetSlopeX(param[1]);
  track.SetSlopeY(param[3]);

  return track;
}

inline double ChargeTrackFitter::EvaluateMetric(const double* param){
  MMTrack track;
  track.SetConstX(param[0]);
  track.SetConstY(param[2]);
  track.SetSlopeX(param[1]);
  track.SetSlopeY(param[3]);

  double chi2 = 0;
  int Nclus = m_clusters->GetNCluster();
  for(int i = 0; i < Nclus; i++){
    const GeoPlane& plane = 
      m_geometry->Get(m_geometry->Index(m_clusters->Get(i).MMFE8()));
    double x_diff = 
      plane.LocalXatYend(track) - 
      plane.LocalXatYend(m_clusters->Get(i).Channel());

    chi2 += x_diff*x_diff;
  }
  return chi2;
}

#endif
