///
///  \file   GeoOctuplet.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef GeoOctuplet_HH
#define GeoOctuplet_HH

#include <map>

#include "TGraphErrors.h"

#include "include/GeoPlane.hh"
#include "include/MMClusterList.hh"

class GeoOctuplet {

public:
  GeoOctuplet();
  GeoOctuplet(int RunNumber);
  GeoOctuplet(const GeoOctuplet& oct);
  ~GeoOctuplet();

  int GetNPlanes() const;
  GeoPlane const& Get(int index) const;
  GeoPlane const& operator [] (int index) const;

  int Index(int mmfe8) const;
  int MMFE8(int index) const;

  std::vector<int> MMFE8list() const;
  
  int RunNumber() const;
  void SetRunNumber(int RunNumber);

  TGraph* GetXZGraph(const MMClusterList& clusters) const;
  TGraphErrors* GetXZGraphErrors(const MMClusterList& clusters) const;
  TGraph*   GetXZGraph(const MMTrack& track) const;
  TGraph2D* Get2DGraph(const MMTrack& track) const;

private:
  void Init();

  std::vector<GeoPlane*> m_planes;

  int m_RunNum;
  mutable std::map<int,int> m_MMFE82Index;
  mutable std::map<int,int> m_Index2MMFE8;

};

inline GeoOctuplet::GeoOctuplet(){
  m_RunNum = -1;
  Init();
}

inline GeoOctuplet::GeoOctuplet(int RunNumber){
  m_RunNum = -1;
  SetRunNumber(RunNumber);
  Init();
}

inline GeoOctuplet::GeoOctuplet(const GeoOctuplet& oct){
  SetRunNumber(oct.RunNumber());
  int N = oct.GetNPlanes();
  for(int i = 0; i < N; i++)
    m_planes.push_back(new GeoPlane(oct[i]));
}

inline GeoOctuplet::~GeoOctuplet(){
  int N = GetNPlanes();
  for(int i = 0; i < N; i++)
    delete m_planes[i];
}

inline void GeoOctuplet::Init(){
  int i = 0;
  TVector3 origin;

  // plane 0
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100., 0.);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(0.);
  m_planes[i]->SetSignChannel(-1);

  i++;
  // plane 1
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100., 11.2);
  m_planes[i]->SetOrigin(origin);
  ;  m_planes[i]->SetStripAlpha(0.);
  m_planes[i]->SetSignChannel(1);

  i++;
  // plane 2
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100.+17.9, 32.4);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(-0.0261799);
  m_planes[i]->SetSignChannel(1);

  i++;
  // plane 3
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100.+17.9, 43.6);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(0.0261799);
  m_planes[i]->SetSignChannel(-1);

  i++;
  // plane 4
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100.+17.9, 113.6);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(-0.0261799);
  m_planes[i]->SetSignChannel(1);

  i++;
  // plane 5
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100.+17.9, 124.8);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(0.0261799);
  m_planes[i]->SetSignChannel(-1);

  i++;
  // plane 6
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100., 146.0);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(0.);
  m_planes[i]->SetSignChannel(-1);

  i++;
  // plane 7
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100., 157.2);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(0.);
  m_planes[i]->SetSignChannel(1);
}

inline int GeoOctuplet::GetNPlanes() const {
  return int(m_planes.size());
}

inline GeoPlane const& GeoOctuplet::Get(int index) const {
  return *m_planes[index];
}

inline GeoPlane const& GeoOctuplet::operator [] (int index) const {
  return Get(index);
}

inline int GeoOctuplet::Index(int mmfe8) const {
  if(m_MMFE82Index.count(mmfe8) <= 0)
    return -1;
  return m_MMFE82Index[mmfe8];
}

inline int GeoOctuplet::MMFE8(int index) const {
  if(m_Index2MMFE8.count(index) <= 0)
    return -1;
  return m_Index2MMFE8[index];
}

inline std::vector<int> GeoOctuplet::MMFE8list() const {
  int N = m_Index2MMFE8.size();
  std::vector<int> ret;
  for(int i = 0; i < N; i++)
    ret.push_back(m_Index2MMFE8[i]);
  return ret;
}

inline int GeoOctuplet::RunNumber() const {
  return m_RunNum;
}

inline void GeoOctuplet::SetRunNumber(int RunNumber) {

  if(RunNumber >= 3504 && RunNumber <= 3506){
    m_MMFE82Index.clear(); 
    m_MMFE82Index[110] = 0;
    m_MMFE82Index[109] = 1;
    m_MMFE82Index[101] = 2;
    m_MMFE82Index[103] = 3;
    m_MMFE82Index[106] = 4;
    m_MMFE82Index[102] = 5;
    m_MMFE82Index[107] = 6;
    m_MMFE82Index[105] = 7;

    m_Index2MMFE8.clear();
    m_Index2MMFE8[0] = 110;
    m_Index2MMFE8[1] = 109;
    m_Index2MMFE8[2] = 101;
    m_Index2MMFE8[3] = 103;
    m_Index2MMFE8[4] = 106;
    m_Index2MMFE8[5] = 102;
    m_Index2MMFE8[6] = 107;
    m_Index2MMFE8[7] = 105;
    
    m_RunNum = RunNumber;
  }

  if(RunNumber >= 3507 && RunNumber <= 3512){
    m_MMFE82Index.clear(); 
    m_MMFE82Index[111] = 0;
    m_MMFE82Index[116] = 1;
    m_MMFE82Index[101] = 2;
    m_MMFE82Index[109] = 3;
    m_MMFE82Index[112] = 4;
    m_MMFE82Index[102] = 5;
    m_MMFE82Index[107] = 6;
    m_MMFE82Index[105] = 7;

    m_Index2MMFE8.clear();
    m_Index2MMFE8[0] = 111;
    m_Index2MMFE8[1] = 116;
    m_Index2MMFE8[2] = 101;
    m_Index2MMFE8[3] = 109;
    m_Index2MMFE8[4] = 112;
    m_Index2MMFE8[5] = 102;
    m_Index2MMFE8[6] = 107;
    m_Index2MMFE8[7] = 105;

    m_RunNum = RunNumber;
  }

  if(RunNumber >= 3513){
    m_MMFE82Index.clear(); 
    m_MMFE82Index[111] = 0;
    m_MMFE82Index[116] = 1;
    m_MMFE82Index[101] = 2;
    m_MMFE82Index[109] = 3;
    m_MMFE82Index[117] = 4;
    m_MMFE82Index[102] = 5;
    m_MMFE82Index[107] = 6;
    m_MMFE82Index[105] = 7;

    m_Index2MMFE8.clear();
    m_Index2MMFE8[0] = 111;
    m_Index2MMFE8[1] = 116;
    m_Index2MMFE8[2] = 101;
    m_Index2MMFE8[3] = 109;
    m_Index2MMFE8[4] = 117;
    m_Index2MMFE8[5] = 102;
    m_Index2MMFE8[6] = 107;
    m_Index2MMFE8[7] = 105;

    m_RunNum = RunNumber;
  }
}

inline TGraph* GeoOctuplet::GetXZGraph(const MMClusterList& clusters) const {
  std::vector<double> vx;
  std::vector<double> vz;
  int Nclus = clusters.GetNCluster();
  for(int i = 0; i < Nclus; i++){
    int index = Index(clusters[i].MMFE8());
    if(index < 0)
      continue;
    double ch = clusters[i].Channel();
    GeoPlane& plane = *m_planes[index];
    TVector3 p = plane.Origin() + plane.LocalXatYbegin(ch)*plane.nX();
    vx.push_back(p.X());
    vz.push_back(p.Z());
  }

  int N = vx.size();
  double x[N];
  double z[N];
  for(int i = 0; i < N; i++){
    x[i] = vx[i];
    z[i] = vz[i];
  }
  TGraph* gr = new TGraph(N,x,z);
  gr->SetMarkerSize(1);
  gr->SetMarkerColor(kBlack);

  return gr;
}

inline TGraphErrors* GeoOctuplet::GetXZGraphErrors(const MMClusterList& clusters) const {
  std::vector<double> vx;
  std::vector<double> vdx;
  std::vector<double> vz;
  int Nclus = clusters.GetNCluster();
  for(int i = 0; i < Nclus; i++){
    int index = Index(clusters[i].MMFE8());
    if(index < 0)
      continue;
    double ch = clusters[i].Channel();
    GeoPlane& plane = *m_planes[index];
    TVector3 p1 = plane.Origin() + plane.LocalXatYend(ch)*plane.nX();
    TVector3 p2 = plane.Origin() + plane.LocalXatYbegin(ch)*plane.nX();
    vx.push_back((p1.X()+p2.X())/2.);
    vz.push_back((p1.Z()+p2.Z())/2.);
    vdx.push_back(fabs(p1.X()-p2.X())/2.);
  }

  int N = vx.size();
  double x[N];
  double dx[N];
  double z[N];
  for(int i = 0; i < N; i++){
    x[i]  = vx[i];
    z[i]  = vz[i];
    dx[i] = vdx[i];
  }
  TGraphErrors* gr = new TGraphErrors(N,x,z,dx,0);
  gr->SetMarkerSize(1);
  gr->SetMarkerColor(kBlack);
  gr->SetLineWidth(2);


  return gr;
}

inline TGraph* GeoOctuplet::GetXZGraph(const MMTrack& track) const {
  std::vector<double> vx;
  std::vector<double> vz;

  TVector3 p, ip;
  p = track.PointZ(-5.);
  vx.push_back(p.X());
  vz.push_back(p.Z());

  int Nplane = GetNPlanes();
  for(int i = 0; i < Nplane; i++){
    GeoPlane& plane = *m_planes[i];
    ip = plane.Intersection(track, -5.6);
    p = track.PointZ(ip.Z());
    vx.push_back(p.X());
    vz.push_back(p.Z());
    
    p = plane.Origin() - 5.6*plane.nZ() + 
      plane.LocalXatYbegin(track,-5.6)*plane.nX();
    vx.push_back(p.X());
    vz.push_back(p.Z());

    p = plane.Origin() + 5.6*plane.nZ() + 
      plane.LocalXatYbegin(track, 5.6)*plane.nX();
    vx.push_back(p.X());
    vz.push_back(p.Z());

    ip = plane.Intersection(track, 5.6);
    p = track.PointZ(ip.Z());
    vx.push_back(p.X());
    vz.push_back(p.Z());
  }

  p = track.PointZ(165.);
  vx.push_back(p.X());
  vz.push_back(p.Z());

  int N = vx.size();
  double x[N];
  double z[N];
  for(int i = 0; i < N; i++){
    x[i] = vx[i];
    z[i] = vz[i];
  }
  TGraph* gr = new TGraph(N,x,z);
  gr->SetLineColor(kRed+2);
  gr->SetMarkerSize(0.);
  gr->SetMarkerColor(kRed+2);
  gr->SetLineWidth(2);
  return gr;
}

inline TGraph2D* GeoOctuplet::Get2DGraph(const MMTrack& track) const {
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> vz;

  TVector3 p;
  
  int Nplane = GetNPlanes();
  for(int i = 0; i < Nplane; i++){
    GeoPlane& plane = *m_planes[i];
    p = plane.Intersection(track);  
    vx.push_back(p.X());
    vy.push_back(p.Y());
    vz.push_back(p.Z());
  }

  int N = vx.size();
  double x[N];
  double y[N];
  double z[N];
  for(int i = 0; i < N; i++){
    x[i] = vx[i];
    y[i] = vy[i];
    z[i] = vz[i];
  }
  TGraph2D* gr = new TGraph2D(N,x,y,z);
  gr->SetLineColor(kRed+2);
  gr->SetMarkerSize(1);
  gr->SetMarkerColor(kRed+2);
  gr->SetLineWidth(2);
  return gr;
}

#endif



