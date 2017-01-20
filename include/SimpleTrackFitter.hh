///
///  \file   SimpleTrackFitter.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///

#ifndef SimpleTrackFitter_HH
#define SimpleTrackFitter_HH

#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

#include "include/MMClusterList.hh"
#include "include/GeoOctuplet.hh"

class SimpleTrackFitter {

public:
  SimpleTrackFitter();
  ~SimpleTrackFitter();

  virtual MMTrack Fit(const MMClusterList& clusters, 
		      const GeoOctuplet& geometry,
		      const int evt = -1);

private:
  ROOT::Math::Minimizer* m_minimizer;
  ROOT::Math::Functor* m_functor;

  MMClusterList* m_clusters;
  const GeoOctuplet*   m_geometry;
  double EvaluateMetric(const double* param);

};

inline SimpleTrackFitter::SimpleTrackFitter(){
  m_minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
  m_minimizer->SetMaxFunctionCalls(10000000);
  m_minimizer->SetMaxIterations(100000);
  m_minimizer->SetTolerance(0.001);
  m_minimizer->SetPrintLevel(0);
  
  m_functor = new ROOT::Math::Functor(this, &SimpleTrackFitter::EvaluateMetric, 4);
  m_minimizer->SetFunction(*m_functor);

  m_minimizer->SetVariable(0, "c_x", 100., 0.01);
  m_minimizer->SetVariable(1, "s_x", 0.,   0.001);
  m_minimizer->SetVariable(2, "c_y", 100., 0.01);
  m_minimizer->SetVariable(3, "s_y", 0.,   0.001);

  m_clusters = nullptr;
  m_geometry = nullptr;
}

inline SimpleTrackFitter::~SimpleTrackFitter(){
  delete m_minimizer;
  delete m_functor;
}

inline MMTrack SimpleTrackFitter::Fit(const MMClusterList& clusters, 
                                      const GeoOctuplet&   geometry,
                                      const int evt){
  // return track
  MMTrack track;

  // get geometry pointer
  m_geometry = &geometry;
  
  // get cluster list
  int N = clusters.GetNCluster();
  m_clusters = new MMClusterList();
  for(int i = 0; i < N; i++)
    if(m_geometry->Index(clusters[i].MMFE8()) >= 0){
      m_clusters->AddCluster(clusters[i]);
      track.CountHit(clusters[i].MMFE8Index());
    }

  // need at least two X and two U/V planes to fit
  if(track.NX() < 2 || track.NU()+track.NV() < 2)
    return track;

  // do fit
  if(m_clusters->GetNCluster() > 0){
    m_minimizer->SetVariableValue(0, 0.);
    m_minimizer->SetVariableValue(1, 0.);
    m_minimizer->SetVariableValue(2, 0.);
    m_minimizer->SetVariableValue(3, 0.);

    bool ok = m_minimizer->Minimize();
    if (!ok && evt != -1)
      std::cout << " Fit failed on Event " << evt
                << " | N(clus) = "  << (int)(m_clusters->size())
                << " | N(X) = "     << track.NX()
                << " | N(U) = "     << track.NU()
                << " | N(V) = "     << track.NV()
                << std::endl;

    
    const double* param = m_minimizer->X();
    track.SetRes2(m_minimizer->MinValue());
    track.SetConstX(param[0]);
    track.SetConstY(param[2]);
    track.SetSlopeX(param[1]);
    track.SetSlopeY(param[3]);
    track.SetIsFit();
  }

  m_geometry = nullptr;
  delete m_clusters;
  m_clusters = nullptr;

  return track;
}

inline double SimpleTrackFitter::EvaluateMetric(const double* param){
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
    double diff = 
      plane.LocalXatYend(track) - 
      plane.LocalXatYend(m_clusters->Get(i).Channel());

    chi2 += diff*diff;
  }
  return chi2;
}

#endif
