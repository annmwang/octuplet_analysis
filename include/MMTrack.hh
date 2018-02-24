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
#include "TGraph2D.h"

class MMTrack {

public:
  MMTrack();
  MMTrack(const MMTrack& track);
  ~MMTrack();

  bool IsFit() const;
  bool IsTrigCand(int m_nx = 2, int m_nuv = 2) const;
  void SetIsFit(bool isfit = true);

  double ConstX() const;
  double ConstY() const;
  double SlopeX() const;
  double SlopeY() const;
  int NX() const;
  int NU() const;
  int NV() const;
  int NBoard(int b) const;

  double ResidualSq() const;

  TVector3 PointX(double x) const;
  TVector3 PointY(double y) const;
  TVector3 PointZ(double z) const;

  void SetRes2(double Res2);
  void SetConstX(double CX);
  void SetConstY(double CY);
  void SetSlopeX(double SX);
  void SetSlopeY(double SY);
  void CountHit(int immfe8);
  
  TGraph*   GetXZGraph(double zmin, double zmax) const;
  TGraph2D* Get2DGraph(double zmin, double zmax) const;

  void Reset();

private:
  // track parameters
  double m_CX;
  double m_CY;
  double m_SX;
  double m_SY;
  
  // number of hits counting for all boards
  int m_NX;
  int m_NU;
  int m_NV;
  int m_N0; int m_N1; int m_N2; int m_N3;
  int m_N4; int m_N5; int m_N6; int m_N7;

  // residual^2 sum from fit
  double m_Res2;

  // flag indicating whether track has been fit or not
  bool m_IsFit;

};

inline MMTrack::MMTrack(){
  m_CX = 0.;
  m_CY = 0.;
  m_SX = 0.;
  m_SY = 0.;
  m_NX = 0;
  m_NU = 0;
  m_NV = 0;
  m_N0 = 0; m_N1 = 0; m_N2 = 0; m_N3 = 0;
  m_N4 = 0; m_N5 = 0; m_N6 = 0; m_N7 = 0;
  m_Res2 = 0.;
  m_IsFit = false;
}

inline MMTrack::MMTrack(const MMTrack& track){
  m_CX = track.ConstX();
  m_CY = track.ConstY();
  m_SX = track.SlopeX();
  m_SY = track.SlopeX();
  m_NX = track.NX();
  m_NU = track.NU();
  m_NV = track.NV();
  m_N0 = track.NBoard(0); m_N1 = track.NBoard(1); m_N2 = track.NBoard(2); m_N3 = track.NBoard(3);
  m_N4 = track.NBoard(4); m_N5 = track.NBoard(5); m_N6 = track.NBoard(6); m_N7 = track.NBoard(7);
  m_Res2 = track.ResidualSq();
  m_IsFit = track.IsFit();
}

inline MMTrack::~MMTrack() {}

inline void MMTrack::Reset(){
  m_CX = 0.;
  m_CY = 0.;
  m_SX = 0.;
  m_SY = 0.;
  m_NX = 0;
  m_NU = 0;
  m_NV = 0;
  m_N0 = 0; m_N1 = 0; m_N2 = 0; m_N3 = 0;
  m_N4 = 0; m_N5 = 0; m_N6 = 0; m_N7 = 0;
}

inline bool MMTrack::IsFit() const {
  return m_IsFit;
}

inline bool MMTrack::IsTrigCand(int m_nx, int m_nuv) const {
  if (m_N0 + m_N1 < 1 ||
      m_N6 + m_N7 < 1 ||
      m_NU + m_NV < m_nuv ||
      m_NX < m_nx ) 
    return false;
  else
    return true;
}

inline void MMTrack::SetIsFit(bool isfit){
  m_IsFit = isfit;
}

inline double MMTrack::ResidualSq() const { return m_Res2; }
inline double MMTrack::ConstX() const { return m_CX; }
inline double MMTrack::ConstY() const { return m_CY; }
inline double MMTrack::SlopeX() const { return m_SX; }
inline double MMTrack::SlopeY() const { return m_SY; }
inline    int MMTrack::NX()     const { return m_NX; }
inline    int MMTrack::NU()     const { return m_NU; }
inline    int MMTrack::NV()     const { return m_NV; }

inline int MMTrack::NBoard(int b) const { 
    if      (b==0) return m_N0;
    else if (b==1) return m_N1;
    else if (b==2) return m_N2;
    else if (b==3) return m_N3;
    else if (b==4) return m_N4;
    else if (b==5) return m_N5;
    else if (b==6) return m_N6;
    else if (b==7) return m_N7;
    return -1;
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

inline void MMTrack::SetRes2(double Res2){ m_Res2 = Res2; }
inline void MMTrack::SetConstX(double CX){ m_CX = CX; }
inline void MMTrack::SetConstY(double CY){ m_CY = CY; }
inline void MMTrack::SetSlopeX(double SX){ m_SX = SX; }
inline void MMTrack::SetSlopeY(double SY){ m_SY = SY; }

inline void MMTrack::CountHit(int immfe8) {
    if      (immfe8 == 0 || immfe8 == 1 ||
             immfe8 == 6 || immfe8 == 7)   m_NX++;
    else if (immfe8 == 2 || immfe8 == 4)   m_NU++;
    else if (immfe8 == 3 || immfe8 == 5)   m_NV++;

    if      (immfe8 == 0) m_N0++;
    else if (immfe8 == 1) m_N1++;
    else if (immfe8 == 2) m_N2++;
    else if (immfe8 == 3) m_N3++;
    else if (immfe8 == 4) m_N4++;
    else if (immfe8 == 5) m_N5++;
    else if (immfe8 == 6) m_N6++;
    else if (immfe8 == 7) m_N7++;
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

inline TGraph2D* MMTrack::Get2DGraph(double zmin, double zmax) const {
  double x[2], y[2], z[2];
  x[0] = m_CX + m_SX*zmin;
  x[1] = m_CX + m_SX*zmax;
  y[0] = m_CY + m_SY*zmin;
  y[1] = m_CY + m_SY*zmax;
  z[0] = zmin;
  z[1] = zmax;

  TGraph2D* gr = new TGraph2D(2,x,y,z);
  gr->SetLineColor(kRed+2);
  gr->SetMarkerSize(0.);
  gr->SetMarkerColor(kRed+2);
  gr->SetLineWidth(2);
  return gr;
}

#endif



