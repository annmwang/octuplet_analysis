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
#include "include/AlignmentTree.hh"

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

  double GetResidualX(const MMCluster& clus,
                      const MMTrack& track);

  double GetQuadraticSumOfResidualsX(MMClusterList clusters,
                                     const MMTrack& track,
                                     bool normalize=true);

  TGraph* GetXZGraph(const MMClusterList& clusters) const;
  TGraphErrors* GetXZGraphErrors(const MMClusterList& clusters) const;
  TGraph*   GetXZGraph(const MMTrack& track) const;
  TGraph2D* Get2DGraph(const MMTrack& track) const;

  void TranslateX(double x, int iplane);
  void TranslateY(double y, int iplane);
  void TranslateZ(double z, int iplane);

  void RotateX(double phix, int iplane);
  void RotateY(double phiy, int iplane);
  void RotateZ(double phiz, int iplane);

  void SetAlignment(const std::string& align_filename);

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
  //origin.SetXYZ(102.3, 100., 0.);
  origin.SetXYZ(102.3, 100., -2.7);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(0.);
  m_planes[i]->SetSignChannel(-1);

  i++;
  // plane 1
  m_planes.push_back(new GeoPlane());
  //origin.SetXYZ(102.3, 100., 11.2);
  origin.SetXYZ(102.3, 100., 11.2+2.7);
  m_planes[i]->SetOrigin(origin);
  ;  m_planes[i]->SetStripAlpha(0.);
  m_planes[i]->SetSignChannel(1);

  i++;
  // plane 2
  m_planes.push_back(new GeoPlane());
  //origin.SetXYZ(102.3, 100.+17.9, 32.4);
  origin.SetXYZ(102.3, 100.+17.9, 32.4-2.7);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(-0.0261799);
  m_planes[i]->SetSignChannel(1);

  i++;
  // plane 3
  m_planes.push_back(new GeoPlane());
  //origin.SetXYZ(102.3, 100.+17.9, 43.6);
  origin.SetXYZ(102.3, 100.+17.9, 43.6+2.7);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(0.0261799);
  m_planes[i]->SetSignChannel(-1);

  i++;
  // plane 4
  m_planes.push_back(new GeoPlane());
  //origin.SetXYZ(102.3, 100.+17.9, 113.6);
  origin.SetXYZ(102.3, 100.+17.9, 113.6-2.7);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(-0.0261799);
  m_planes[i]->SetSignChannel(1);

  i++;
  // plane 5
  m_planes.push_back(new GeoPlane());
  //origin.SetXYZ(102.3, 100.+17.9, 124.8);
  origin.SetXYZ(102.3, 100.+17.9, 124.8+2.7);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(0.0261799);
  m_planes[i]->SetSignChannel(-1);

  i++;
  // plane 6
  m_planes.push_back(new GeoPlane());
  //origin.SetXYZ(102.3, 100., 146.0);
  origin.SetXYZ(102.3, 100., 146.0-2.7);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(0.);
  m_planes[i]->SetSignChannel(-1);

  i++;
  // plane 7
  m_planes.push_back(new GeoPlane());
  //origin.SetXYZ(102.3, 100., 157.2);
  origin.SetXYZ(102.3, 100., 157.2+2.7);
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

inline double GeoOctuplet::GetResidualX(const MMCluster& clus,
                                        const MMTrack& track){
  int iplane = Index(clus.MMFE8());
  if(iplane < 0)
    return 0.;

  return m_planes[iplane]->GetResidualX(clus.Channel(), track);
}

inline double GeoOctuplet::GetQuadraticSumOfResidualsX(MMClusterList clusters,
                                                       const MMTrack& track,
                                                       bool normalize){
  if (clusters.size() == 0)
    return 0.0;

  double sum = 0.0;
  for (auto clus: clusters){
    sum += pow(GetResidualX(*clus, track), 2.0);
  }
  sum = pow(sum, 0.5);

  if (normalize) return sum / (double)(clusters.size());
  else           return sum;
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

inline void GeoOctuplet::TranslateX(double x, int iplane){
  if(iplane < 0 || iplane >= GetNPlanes())
    return;
  m_planes[iplane]->TranslateX(x);
}

inline void GeoOctuplet::TranslateY(double y, int iplane){
  if(iplane < 0 || iplane >= GetNPlanes())
    return;
  m_planes[iplane]->TranslateY(y);
}

inline void GeoOctuplet::TranslateZ(double z, int iplane){
  if(iplane < 0 || iplane >= GetNPlanes())
    return;
  m_planes[iplane]->TranslateZ(z);
}

inline void GeoOctuplet::RotateX(double phix, int iplane){
  if(iplane < 0 || iplane >= GetNPlanes())
    return;
  m_planes[iplane]->RotateX(phix);
}

inline void GeoOctuplet::RotateY(double phiy, int iplane){
  if(iplane < 0 || iplane >= GetNPlanes())
    return;
  m_planes[iplane]->RotateY(phiy);
}

inline void GeoOctuplet::RotateZ(double phiz, int iplane){
  if(iplane < 0 || iplane >= GetNPlanes())
    return;
  m_planes[iplane]->RotateZ(phiz);
}

inline void GeoOctuplet::SetAlignment(const std::string& align_filename){
  TFile f(align_filename.c_str(), "READ");
  if(!f.IsOpen())
    return;
  TTree* T = (TTree*) f.Get("AlignmentTree");
  if(!T)
    return;
  
  AlignmentTree align(T);
  align.GetEntry(0);

  TranslateX(align.tranX_1, 1);
  TranslateX(align.tranX_2, 2);
  TranslateX(align.tranX_3, 3);
  TranslateX(align.tranX_4, 4);
  TranslateX(align.tranX_5, 5);
  TranslateX(align.tranX_6, 6);
  TranslateX(align.tranX_7, 7);
  TranslateY(align.tranY_1, 1);
  TranslateY(align.tranY_2, 2);
  TranslateY(align.tranY_3, 3);
  TranslateY(align.tranY_4, 4);
  TranslateY(align.tranY_5, 5);
  TranslateY(align.tranY_6, 6);
  TranslateY(align.tranY_7, 7);
  TranslateZ(align.tranZ_1, 1);
  TranslateZ(align.tranZ_2, 2);
  TranslateZ(align.tranZ_3, 3);
  TranslateZ(align.tranZ_4, 4);
  TranslateZ(align.tranZ_5, 5);
  TranslateZ(align.tranZ_6, 6);
  TranslateZ(align.tranZ_7, 7);

  RotateX(align.rotX_1, 1);
  RotateX(align.rotX_2, 2);
  RotateX(align.rotX_3, 3);
  RotateX(align.rotX_4, 4);
  RotateX(align.rotX_5, 5);
  RotateX(align.rotX_6, 6);
  RotateX(align.rotX_7, 7);
  RotateY(align.rotY_1, 1);
  RotateY(align.rotY_2, 2);
  RotateY(align.rotY_3, 3);
  RotateY(align.rotY_4, 4);
  RotateY(align.rotY_5, 5);
  RotateY(align.rotY_6, 6);
  RotateY(align.rotY_7, 7);
  RotateZ(align.rotZ_1, 1);
  RotateZ(align.rotZ_2, 2);
  RotateZ(align.rotZ_3, 3);
  RotateZ(align.rotZ_4, 4);
  RotateZ(align.rotZ_5, 5);
  RotateZ(align.rotZ_6, 6);
  RotateZ(align.rotZ_7, 7);

  f.Close();
}

#endif



