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

  double GetResidualX(double channel,
		      const MMTrack& track);

  TGraph*   GetXZGraph() const;
  TGraph2D* Get2DGraph() const;
  TGraph2D* GetChannelGraph(double channel) const;

  void TranslateX(double x);
  void TranslateY(double y);
  void TranslateZ(double z);

  void RotateX(double phix);
  void RotateY(double phiy);
  void RotateZ(double phiz);
  
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

inline double GeoPlane::LocalXatYbegin(const MMTrack& track,
				       double zlocal) const {
  TVector3 p = Intersection(track, zlocal) - 
    (m_Origin+zlocal*m_nZ.Unit());
  double x = p.Dot(m_nX.Unit());
  double y = p.Dot(m_nY.Unit());
  
  double dx = tan(m_Alpha)*(y+100.);

  return x+dx;
}

inline double GeoPlane::LocalXatYbegin(double channel) const {
  return m_SignChannel*(channel-256.5)*0.4 + tan(m_Alpha)*200.;
}

inline double GeoPlane::LocalXatYend(const MMTrack& track,
			      double zlocal) const {
  TVector3 p = Intersection(track, zlocal) - 
    (m_Origin+zlocal*m_nZ.Unit());
  double x = p.Dot(m_nX.Unit());
  double y = p.Dot(m_nY.Unit());

  double dx = tan(m_Alpha)*(y-100.);

  return x+dx;
}

inline double GeoPlane::LocalXatYend(double channel) const {
  return m_SignChannel*(channel-256.5)*0.4;
}

inline double GeoPlane::GetResidualX(double channel,
			      const MMTrack& track){
  return LocalXatYend(track) - LocalXatYend(channel);
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

inline TGraph2D* GeoPlane::Get2DGraph() const {
  double x[5], y[5], z[5];
  TVector3 p0 = m_Origin +
    LocalXatYend(0)*m_nX +
    100.*m_nY;
  TVector3 p1 = m_Origin +
    LocalXatYend(513)*m_nX +
    100.*m_nY;
  TVector3 p2 = m_Origin +
    LocalXatYbegin(513)*m_nX +
    -100.*m_nY;
  TVector3 p3 = m_Origin +
    LocalXatYbegin(0)*m_nX +
    -100.*m_nY;
  
  x[0] = p0.X();
  y[0] = p0.Y();
  z[0] = p0.Z();
  x[1] = p1.X();
  y[1] = p1.Y();
  z[1] = p1.Z();
  x[2] = p2.X();
  y[2] = p2.Y();
  z[2] = p2.Z();
  x[3] = p3.X();
  y[3] = p3.Y();
  z[3] = p3.Z();
  x[4] = p0.X();
  y[4] = p0.Y();
  z[4] = p0.Z();

  TGraph2D* gr = new TGraph2D(5,x,y,z);
  gr->SetLineColor(kBlue+2);
  gr->SetMarkerSize(0.);
  gr->SetMarkerColor(kBlue+2);
  return gr;
}


/*
inline TGraph2D* GeoPlane::Get2DGraph() const {
  TVector3 p0 = m_Origin +
    LocalXatYend(0)*m_nX +
    100.*m_nY;
  TVector3 p1 = m_Origin +
    LocalXatYend(513)*m_nX +
    100.*m_nY;
  TVector3 p2 = m_Origin +
    LocalXatYbegin(513)*m_nX +
    -100.*m_nY;
  TVector3 p3 = m_Origin +
    LocalXatYbegin(0)*m_nX +
    -100.*m_nY;
  
  double x[29], y[29], z[29];
  x[0] = p0.X();
  y[0] = p0.Y();
  z[0] = p0.Z();
  x[1] = p1.X();
  y[1] = p1.Y();
  z[1] = p1.Z();
  x[2] = p2.X();
  y[2] = p2.Y();
  z[2] = p2.Z();
  x[3] = p3.X();
  y[3] = p3.Y();
  z[3] = p3.Z();
  x[4] = p0.X();
  y[4] = p0.Y();
  z[4] = p0.Z();
  
  for(int i = 0; i < 8; i++){
    TVector3 s0 = m_Origin +
      LocalXatYbegin(i*64+1)*m_nX +
      -100.*m_nY;
    TVector3 s1 = m_Origin +
      LocalXatYend(i*64+1)*m_nX +
      100.*m_nY;

    if(i%2 == 0){
      x[i*3+5] = s0.X();
      y[i*3+5] = s0.Y();
      z[i*3+5] = s0.Z();
      x[i*3+1+5] = s1.X();
      y[i*3+1+5] = s1.Y();
      z[i*3+1+5] = s1.Z();
    } else {
      x[i*3+5] = s1.X();
      y[i*3+5] = s1.Y();
      z[i*3+5] = s1.Z();
      x[i*3+1+5] = s0.X();
      y[i*3+1+5] = s0.Y();
      z[i*3+1+5] = s0.Z();
    }
    if(i > 0){
      x[i*3+2+5] = x[i*3-3+5];
      y[i*3+2+5] = y[i*3-3+5];
      z[i*3+2+5] = z[i*3-3+5];
    } else {
      x[i*3+2+5] = x[i*3+1+5];
      y[i*3+2+5] = y[i*3+1+5];
      z[i*3+2+5] = z[i*3+1+5];
    }
  }

  TGraph2D* gr = new TGraph2D(29,x,y,z);
  gr->SetLineColor(kBlue-9);
  gr->SetMarkerSize(0.);
  gr->SetMarkerColor(kBlue-9);
  return gr;
}
*/

inline TGraph2D* GeoPlane::GetChannelGraph(double channel) const {
  double x[2], y[2], z[2];
  TVector3 p0 = m_Origin +
    LocalXatYend(channel)*m_nX +
    100.*m_nY;
  TVector3 p1 = m_Origin +
    LocalXatYbegin(channel)*m_nX +
    -100.*m_nY;
  
  x[0] = p0.X();
  y[0] = p0.Y();
  z[0] = p0.Z();
  x[1] = p1.X();
  y[1] = p1.Y();
  z[1] = p1.Z();

  TGraph2D* gr = new TGraph2D(2,x,y,z);
  gr->SetLineColor(kBlack);
  gr->SetMarkerSize(0.);
  gr->SetMarkerColor(kBlack);
  gr->SetLineWidth(2);
  return gr;
}

inline void GeoPlane::TranslateX(double x){
  m_Origin.SetX(m_Origin.X()+x);
}

inline void GeoPlane::TranslateY(double y){
  m_Origin.SetY(m_Origin.Y()+y);
}

inline void GeoPlane::TranslateZ(double z){
  m_Origin.SetZ(m_Origin.Z()+z);
}

inline void GeoPlane::RotateX(double phix){
  m_nX.RotateX(phix);
  m_nY.RotateX(phix);
  m_nZ.RotateX(phix);
  m_nX = m_nX.Unit();
  m_nY = m_nY.Unit();
  m_nZ = m_nZ.Unit();
}

inline void GeoPlane::RotateY(double phiy){
  m_nX.RotateY(phiy);
  m_nY.RotateY(phiy);
  m_nZ.RotateY(phiy);
  m_nX = m_nX.Unit();
  m_nY = m_nY.Unit();
  m_nZ = m_nZ.Unit();
}

inline void GeoPlane::RotateZ(double phiz){
  m_nX.RotateZ(phiz);
  m_nY.RotateZ(phiz);
  m_nZ.RotateZ(phiz);
  m_nX = m_nX.Unit();
  m_nY = m_nY.Unit();
  m_nZ = m_nZ.Unit();
}

#endif



