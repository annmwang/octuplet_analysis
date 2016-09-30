///
///  \file   MMTrack.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMTrack_HH
#define MMTrack_HH

#include "TVector3.h"
#include "TGraph.h"

class MMTrack {

public:
  MMTrack();
  MMTrack(const MMTrack& track);
  ~MMTrack();

  double ConstX() const;
  double ConstY() const;
  double SlopeX() const;
  double SlopeY() const;

  TVector3 PointX(double x) const;
  TVector3 PointY(double y) const;
  TVector3 PointZ(double z) const;

  void SetConstX(double CX);
  void SetConstY(double CY);
  void SetSlopeX(double SX);
  void SetSlopeY(double SY);
  
  TGraph* GetXZGraph(double zmin, double zmax) const;

private:
  double m_CX;
  double m_CY;
  double m_SX;
  double m_SY;

};

inline MMTrack::MMTrack(){
  m_CX = 0.;
  m_CY = 0.;
  m_SX = 0.;
  m_SY = 0.;
}

inline MMTrack::MMTrack(const MMTrack& track){
  m_CX = track.ConstX();
  m_CY = track.ConstY();
  m_SX = track.SlopeX();
  m_SY = track.SlopeX();
}

inline MMTrack::~MMTrack() {}

inline double MMTrack::ConstX() const {
  return m_CX;
}

inline double MMTrack::ConstY() const {
  return m_CY;
}

inline double MMTrack::SlopeX() const {
  return m_SX;
}

inline double MMTrack::SlopeY() const {
  return m_SY;
}

inline TVector3 MMTrack::PointX(double x) const {
  TVector3 p;
  p.SetX(x);
  if(m_SX != 0){
    p.SetY(m_CY + m_SY*(x-m_CX)/m_SX);
    p.SetZ((x-m_CX)/m_SX);
  } else {
    p.SetY(0.);
    p.SetZ(0.);
  }
  
  return p;
}

inline TVector3 MMTrack::PointY(double y) const {
  TVector3 p;
  p.SetY(y);
  if(m_SY != 0){
    p.SetX(m_CX + m_SX*(y-m_CY)/m_SY);
    p.SetZ((y-m_CY)/m_SY);
  } else {
    p.SetX(0.);
    p.SetZ(0.);
  }
  
  return p;
}

inline TVector3 MMTrack::PointZ(double z) const {
  TVector3 p;
  p.SetX(m_CX + m_SX*z);
  p.SetY(m_CY + m_SY*z);
  p.SetZ(z);

  return p;
}

inline void MMTrack::SetConstX(double CX){
  m_CX = CX;
}

inline void MMTrack::SetConstY(double CY){
  m_CY = CY;
}

inline void MMTrack::SetSlopeX(double SX){
  m_SX = SX;
}

inline void MMTrack::SetSlopeY(double SY){
  m_SY = SY;
}

inline TGraph* MMTrack::GetXZGraph(double zmin, double zmax) const {
  double x[2], z[2];
  x[0] = m_CX + m_SX*zmin;
  x[1] = m_CX + m_SX*zmax;
  z[0] = zmin;
  z[1] = zmax;

  TGraph* gr = new TGraph(2,x,z);
  gr->SetLineColor(kRed+2);
  gr->SetMarkerSize(0.);
  gr->SetMarkerColor(kRed+2);
  gr->SetLineWidth(2);
  return gr;
}

#endif



