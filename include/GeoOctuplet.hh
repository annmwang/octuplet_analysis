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

#include "include/GeoPlane.hh"

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
  
  int RunNumber() const;
  void SetRunNumber(int RunNumber);

private:
  void Init();

  vector<GeoPlane*> m_planes;

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
  m_planes[i]->SetStripAlpha(0.);
  m_planes[i]->SetSignChannel(1);

  i++;
  // plane 2
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100.-17.9, 32.4);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(-0.0261799);
  m_planes[i]->SetSignChannel(1);

  i++;
  // plane 3
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100.-17.9, 43.6);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(0.0261799);
  m_planes[i]->SetSignChannel(1);

  i++;
  // plane 4
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100.-17.9, 113.6);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(-0.0261799);
  m_planes[i]->SetSignChannel(1);

  i++;
  // plane 5
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100.-17.9, 124.8);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(0.0261799);
  m_planes[i]->SetSignChannel(1);

  i++;
  // plane 6
  m_planes.push_back(new GeoPlane());
  origin.SetXYZ(102.3, 100., 146.0);
  m_planes[i]->SetOrigin(origin);
  m_planes[i]->SetStripAlpha(0.);
  m_planes[i]->SetSignChannel(1);

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

inline int GeoOctuplet::RunNumber() const {
  return m_RunNum;
}

inline void GeoOctuplet::SetRunNumber(int RunNumber) {

  if(RunNumber == 3505){
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

  if(RunNumber == 3508){
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
}


#endif



