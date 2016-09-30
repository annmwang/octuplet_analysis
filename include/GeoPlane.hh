///
///  \file   GeoPlane.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef GeoPlane_HH
#define GeoPlane_HH

#include "include/MMTrack.hh"

class GeoPlane {

public:
  GeoPlane();
  GeoPlane(const GeoPlane& plane);
  ~GeoPlane();

  TVector3 const& Origin() const;
  TVector3 const& nX() const;
  TVector3 const& nY() const;
  TVector3 const& nZ() const;
  double StripAlpha() const;

  int SignChannel() const;
  
  // local x-coordinate at chamber edges
  double LocalXatYbegin(const MMTrack& track,
			double zlocal = 0.) const;
  double LocalXatYend(const MMTrack& track,
		      double zlocal = 0.) const;
  double LocalXatYbegin(double channel) const;
  double LocalXatYend(double channel) const;
  

  TVector3 Intersection(const MMTrack& track,
			double zlocal = 0.) const;

  void SetOrigin(const TVector3& p);
  void SetStripAlpha(double alpha);
  void SetSignChannel(int sign);

  TGraph* GetXZGraph() const;
  
private:
  TVector3 m_Origin;

  TVector3 m_nX;
  TVector3 m_nY;
  TVector3 m_nZ;

  double m_Alpha;

  int m_SignChannel;

};

inline GeoPlane::GeoPlane(){
  m_Origin.SetXYZ(0.,0.,0.);
  m_nX.SetXYZ(1.,0.,0.);
  m_nY.SetXYZ(0.,1.,0.);
  m_nZ.SetXYZ(0.,0.,1.);
  m_Alpha = 0.;
  m_SignChannel = 1;
}

inline GeoPlane::GeoPlane(const GeoPlane& plane){
  m_Origin = plane.Origin();
  m_nX = plane.nX();
  m_nY = plane.nY();
  m_nZ = plane.nZ();
  m_Alpha = plane.StripAlpha();
  SetSignChannel(plane.SignChannel());
}

inline GeoPlane::~GeoPlane() {}

inline TVector3 const& GeoPlane::Origin() const {
  return m_Origin;
}

inline TVector3 const& GeoPlane::nX() const {
  return m_nX;
}

inline TVector3 const& GeoPlane::nY() const {
  return m_nY;
}

inline TVector3 const& GeoPlane::nZ() const {
  return m_nZ;
}

inline double GeoPlane::StripAlpha() const {
  return m_Alpha;
}

inline int GeoPlane::SignChannel() const {
  return m_SignChannel;
}

double GeoPlane::LocalXatYbegin(const MMTrack& track,
				double zlocal) const {
  TVector3 p = Intersection(track, zlocal) - 
    (m_Origin+zlocal*m_nZ.Unit());
  double x = p.Dot(m_nX.Unit());
  double y = p.Dot(m_nY.Unit());

  double dx = tan(m_Alpha)*(y+100.);

  return x+dx;
}

double GeoPlane::LocalXatYbegin(double channel) const {
  return m_SignChannel*(channel-256.5)*0.4 + tan(m_Alpha)*200.;
}

double GeoPlane::LocalXatYend(const MMTrack& track,
			      double zlocal) const {
  TVector3 p = Intersection(track, zlocal) - 
    (m_Origin+zlocal*m_nZ.Unit());
  double x = p.Dot(m_nX.Unit());
  double y = p.Dot(m_nY.Unit());

  double dx = tan(m_Alpha)*(y-100.);

  return x+dx;
}

double GeoPlane::LocalXatYend(double channel) const {
  return m_SignChannel*(channel-256.5)*0.4;
}

inline TVector3 GeoPlane::Intersection(const MMTrack& track, 
				       double zlocal) const {
  TVector3 L0(0.,0.,0.);
  L0.SetX(track.ConstX());
  L0.SetY(track.ConstY());

  TVector3 L(0.,0.,1.);
  L.SetX(track.SlopeX());
  L.SetY(track.SlopeY());

  double z = (m_Origin+zlocal*m_nZ.Unit()-L0).Dot(m_nZ.Unit())/L.Dot(m_nZ.Unit());

  return track.PointZ(z);
}

inline void GeoPlane::SetOrigin(const TVector3& p){
  m_Origin = p;
}

inline void GeoPlane::SetStripAlpha(double alpha){
  m_Alpha = alpha;
}

inline void GeoPlane::SetSignChannel(int sign){
  if(sign > 0)
    m_SignChannel = 1;
  if(sign < 0)
    m_SignChannel = -1;
}

inline TGraph* GeoPlane::GetXZGraph() const {
  double x[2], z[2];
  x[0] = (m_Origin + LocalXatYend(0)*m_nX).X();
  x[1] = (m_Origin + LocalXatYend(513)*m_nX).X();
  z[0] = (m_Origin + LocalXatYend(0)*m_nX).Z();
  z[1] = (m_Origin + LocalXatYend(513)*m_nX).Z();

  TGraph* gr = new TGraph(2,x,z);
  gr->SetLineColor(kBlue+2);
  gr->SetMarkerSize(0.);
  gr->SetMarkerColor(kBlue+2);
  return gr;
}

#endif



