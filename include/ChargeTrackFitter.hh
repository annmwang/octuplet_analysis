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
  ChargeTrackFitter(double Charge2Res, double Boards2Res);
  ~ChargeTrackFitter();

  MMTrack Fit(const MMClusterList clusters, const GeoOctuplet& geometry, double total_charge, bool loud);
  bool IsGoodCluster(const MMClusterList& clusters) const;

private:
  ROOT::Math::Minimizer* m_minimizer;
  ROOT::Math::Functor* m_functor;
  double m_total_charge;
  double m_charge2res;
  double m_boards2res;

  MMClusterList* m_clusters;
  MMClusterList* m_fit_clusters;
  const GeoOctuplet*   m_geometry;
  double EvaluateMetric(const double* param);
  MMTrack _DoFit(MMClusterList *fit_clusters);
  bool m_loud;

};

inline ChargeTrackFitter::ChargeTrackFitter(double Charge2Res, double Boards2Res){
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
  m_charge2res = Charge2Res;
  m_boards2res = Boards2Res; 
}

inline ChargeTrackFitter::~ChargeTrackFitter(){
  delete m_minimizer;
  delete m_functor;
}

inline MMTrack ChargeTrackFitter::Fit(const MMClusterList clusters, const GeoOctuplet& geometry, double total_charge, bool loud=false){
  // return track
  MMTrack track;
  MMClusterList *track_clusters = new MMClusterList();
  bool track_initialized = false;
  m_loud = loud;

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

  for (unsigned i = 1; i < 1 << NClusters; i++) {
    // the jth cluster is included if the jth bit is 1
    // this will capture all the possibilities
    MMClusterList *fit_clusters = new MMClusterList();

    for (unsigned j = 0; j < NClusters; j++)
      if (i & (1 << j)) {
        fit_clusters->AddCluster( m_clusters->Get(j) );
      }
    
    if (IsGoodCluster(*fit_clusters)) {
      if (m_loud) { 
          cout << "BEGIN" << endl;
    for (unsigned j = 0; j < NClusters; j++) {
      if (i & (1 << j)) {
        cout << "CLUSTER: " << j << " CHARGE: " << m_clusters->Get(j).Charge() << " VMM: " << m_geometry->Index(m_clusters->Get(j).MMFE8()) << " CH: " << m_clusters->Get(j).Channel() << endl;
      }
    }
      }
      MMTrack temp_track = _DoFit(fit_clusters);
      if (m_loud) {
        cout << "SCORE: " << temp_track.GetFitScore() << endl << endl;
      }
      if (temp_track.GetFitScore() > 0 && (!track_initialized || temp_track.GetFitScore() < track.GetFitScore())) {
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

  //cout << "FINAL SCORE: " << track.GetFitScore() << endl << endl << endl;
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

  return (uv_planes.size() > 1 && uv_planes.size() + planes.size() > 4);
}

inline MMTrack ChargeTrackFitter::_DoFit(MMClusterList *fit_clusters) {
  MMTrack track;
  // setup: get total charge and averagre x
  double init_x = 0;
  double cluster_charge = 0;
  int Nclus = fit_clusters->GetNCluster();
  for (int i = 0; i < Nclus; i ++) { 
    cluster_charge += fit_clusters->Get(i).Charge();
    init_x += m_geometry->Get(m_geometry->Index(fit_clusters->Get(i).MMFE8())).
        LocalXatYend( fit_clusters->Get(i).Channel() );
  }

  //init_x = init_x / double(Nclus);
  init_x = 0.;
  m_minimizer->SetVariableValue(0, init_x);
  m_minimizer->SetVariableValue(1, 0.);
  m_minimizer->SetVariableValue(2, 0.);
  m_minimizer->SetVariableValue(3, 0.);

  m_fit_clusters = fit_clusters;
  m_minimizer->Minimize();
  m_fit_clusters = NULL;

  if (m_loud)
    cout << "RES: " << m_minimizer->MinValue() << " Q: " << pow((m_total_charge - cluster_charge), 2) << " B: " << pow((Nclus - m_geometry->GetNPlanes()), 2) << endl;

  double score = m_minimizer->MinValue() + pow(m_charge2res*(m_total_charge - cluster_charge), 2)
     + pow(m_boards2res*(Nclus - m_geometry->GetNPlanes()), 2); 
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
  int Nclus = m_fit_clusters->GetNCluster();
  for(int i = 0; i < Nclus; i++){
    const GeoPlane& plane = 
      m_geometry->Get(m_geometry->Index(m_fit_clusters->Get(i).MMFE8()));

    double x_diff = plane.GetResidualX(m_fit_clusters->Get(i).Channel(), track);
    double y_diff = plane.GetResidualY(m_fit_clusters->Get(i).Channel(), track);
    if (m_loud) 
        cout << "PlaneX: " << plane.LocalXatYend(m_fit_clusters->Get(i).Channel()) << " TrackY: " << plane.Intersection(track).Y() << " TrackX: " << plane.LocalXatYend(track) << " YDiff: " << y_diff << " XDiff: "  << x_diff << " CLUS: " << i << " VMM: " << m_geometry->Index(m_fit_clusters->Get(i).MMFE8()) << endl;
    //chi2 += x_diff*x_diff + y_diff*y_diff;
    chi2 += x_diff*x_diff;
  }
  return chi2;
}

#endif
