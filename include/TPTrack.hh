///
///  \file   TPTrack.hh
///
///  \author Ann Miao Wang
///        
///  \date   2017 May
///


#ifndef TPTrack_HH
#define TPTrack_HH

#include "include/TPHit.hh"
#include <numeric>

class TPTrack {

public:
  TPTrack();
  TPTrack(const TPHit& hit);
  TPTrack(const TPTrack& track);
  
  ~TPTrack();

  bool AddHit(const TPHit& hit);
  bool AddTrack(const TPTrack& track);

  bool operator += (const TPHit& hit);
  bool operator += (const TPTrack& track);
  
  bool Contains(const TPHit& hit) const;
  int GetIndex(const TPHit& hit) const;

  void SetMxLocal(double mxl);
  void SetBCID(int bcid);
  void SetTime(int time, int time_ns);    
  void SetEventNum(int evt);    
  void SetEventNumGBT(int evt_gbt);
  void SetNMatch(int nmatch);

  double MxLocal() const;
  int BCID() const;
  double Time() const;
  int EventNum() const;
  int EventNumGBT() const;
  int NMatch() const;

  int NX();
  int NU();
  int NV();
  bool IsTrigCand(int m_nx = 2, int m_nuv = 2);
  bool WrapsAround();
  int BCIDWindow();
  int BCIDEarliest();
  int BCIDEarliestA();
  int BCIDEarliestB();
  int BCIDLatest();
  double BCIDAverage();
  double BCIDAverageA();
  double BCIDAverageB();

  size_t size() const;
  int GetNHits() const;
  TPHit const& Get(int ihit) const;
  TPHit const& operator [] (int ihit) const;

  std::vector<TPHit*>::iterator begin();
  std::vector<TPHit*>::iterator end();

private:
  std::vector<TPHit*> m_hits;

  double m_mxlocal;
  int m_BCID;
  double m_time;
  int m_evt;
  int m_evt_gbt;
  int m_nmatch;
  std::vector<int> m_bcids;
  std::vector<int> m_bcids_a;
  std::vector<int> m_bcids_b;

};

TPTrack operator + (const TPTrack& hits_a, const TPTrack& hits_b);
TPTrack operator + (const TPTrack& hits, const TPHit& hit);
TPTrack operator + (const TPHit& hit, const TPTrack& hits);
TPTrack operator + (const TPHit& hit_a, const TPHit& hit_b);

inline TPTrack::TPTrack() {
  m_mxlocal = -999.;
  m_BCID = -1;
  m_time = -1;
  m_evt = -1;
  m_evt_gbt = -1;
  m_nmatch = -1;
}

inline TPTrack::TPTrack(const TPHit& hit){
  m_mxlocal = -999.;
  m_BCID = -1;
  m_time = -1;
  m_evt = -1;
  m_evt_gbt = -1;
  m_nmatch = -1;
  AddHit(hit);
}

inline TPTrack::TPTrack(const TPTrack& track){
  m_mxlocal = track.MxLocal();
  m_BCID = track.BCID();
  m_time = track.Time();
  m_evt = track.EventNum();
  m_evt_gbt = track.EventNumGBT();
  m_nmatch = track.NMatch();
  AddTrack(track);
}

inline TPTrack::~TPTrack(){
  int N = GetNHits();
  for(int i = 0; i < N; i++)
    delete m_hits[i];
}

inline std::vector<TPHit*>::iterator TPTrack::begin(){
  return m_hits.begin();
}

inline std::vector<TPHit*>::iterator TPTrack::end(){
  return m_hits.end();
}

inline bool TPTrack::AddHit(const TPHit& hit){
  TPHit* newhit = nullptr;
  newhit = new TPHit(hit);
  m_hits.push_back(newhit);
  return true;
}

inline bool TPTrack::AddTrack(const TPTrack& track){
  int N = track.GetNHits();
  double ret = true;
  for(int i = 0; i < N; i++)
    ret = AddHit(track[i]) && ret;
  return ret;
}

inline bool TPTrack::Contains(const TPHit& hit) const {
  int Nhit = GetNHits();
  for(int i = 0; i < Nhit; i++){
    if(Get(i).Channel() == hit.Channel()
       && Get(i).MMFE8() == hit.MMFE8()
       && Get(i).BCID() == hit.BCID())
      return true;
  }
  return false;
}

inline int TPTrack::GetIndex(const TPHit& hit) const {
  int Nhit = GetNHits();
  for(int i = 0; i < Nhit; i++){
    if(Get(i).Channel() == hit.Channel()
       && Get(i).MMFE8() == hit.MMFE8()
       && Get(i).BCID() == hit.BCID())
      return i;
  }
  return -1;
}

inline void TPTrack::SetMxLocal(double mxl){
  m_mxlocal = mxl;
}

inline void TPTrack::SetBCID(int bcid){
  m_BCID = bcid;
}

inline void TPTrack::SetTime(int time, int time_ns){
  m_time = time+time_ns/pow(10,9.);
}

inline void TPTrack::SetEventNum(int evt){
  m_evt = evt;
}

inline void TPTrack::SetEventNumGBT(int evt_gbt){
  m_evt_gbt = evt_gbt;
}

inline void TPTrack::SetNMatch(int nmatch){
  m_nmatch = nmatch;
}

inline bool TPTrack::operator += (const TPHit& hit){
  return AddHit(hit);
}

inline bool TPTrack::operator += (const TPTrack& track){
  return AddTrack(track);
}

inline double TPTrack::MxLocal() const{
  return m_mxlocal;
}

inline int TPTrack::BCID() const{
  return m_BCID;
}

inline double TPTrack::Time() const{
  return m_time;
}

inline int TPTrack::EventNum() const{
  return m_evt;
}

inline int TPTrack::EventNumGBT() const{
  return m_evt_gbt;
}

inline int TPTrack::NMatch() const{
  return m_nmatch;
}

inline int TPTrack::NX() {
  int Nhit = GetNHits();
  int nb = 0;
  for (int i = 0; i < Nhit; i++){
    if (Get(i).MMFE8Index() < 2 || Get(i).MMFE8Index() > 5)
      nb++;
  }
  return nb;
}

inline int TPTrack::NU() {
  int Nhit = GetNHits();
  int nb = 0;
  for (int i = 0; i < Nhit; i++){
    if (Get(i).MMFE8Index() == 2 || Get(i).MMFE8Index() == 4)
      nb++;
  }
  return nb;
}

inline int TPTrack::NV() {
  int Nhit = GetNHits();
  int nb = 0;
  for (int i = 0; i < Nhit; i++){
    if (Get(i).MMFE8Index() == 3 || Get(i).MMFE8Index() == 5)
      nb++;
  }
  return nb;
}

inline int TPTrack::GetNHits() const {
  return int(m_hits.size());
}

inline size_t TPTrack::size() const {
  return m_hits.size();
}

inline TPHit const& TPTrack::Get(int ihit) const {
  return *m_hits[ihit];
}

inline TPHit const& TPTrack::operator [] (int ihit) const {
  return Get(ihit);
}

inline TPTrack operator + (const TPTrack& hits_a,
			     const TPTrack& hits_b){
  TPTrack ret(hits_a);
  ret += hits_b;
  return ret;
}

inline TPTrack operator + (const TPTrack& hits,
			     const TPHit& hit){
  TPTrack ret(hits);
  ret += hit;
  return ret;
}

inline TPTrack operator + (const TPHit& hit,
			     const TPTrack& hits){
  TPTrack ret(hit);
  ret += hits;
  return ret;
}

inline TPTrack operator + (const TPHit& hit_a,
			     const TPHit& hit_b){
  TPTrack ret(hit_a);
  ret += hit_b;
  return ret;
}

inline bool TPTrack::IsTrigCand(int m_nx, int m_nuv){
  bool x01 = false, x67 = false;
  for (int i = 0; i < m_hits.size(); i++){
    if      (Get(i).MMFE8Index() == 0 || Get(i).MMFE8Index() == 1) x01 = true;
    else if (Get(i).MMFE8Index() == 6 || Get(i).MMFE8Index() == 7) x67 = true;
  }
  return (x01 && x67 && NU()+NV() >= m_nuv && NX() > m_nx);
}

inline bool TPTrack::WrapsAround(){
  if (m_bcids.size() == 0){
    std::cout << "WARNING: Using TPTrack::WrapsAround with no BCIDs." << std::endl;
    return false;
  }
  int bcid_min = *std::min_element(m_bcids.begin(), m_bcids.end());
  int bcid_max = *std::max_element(m_bcids.begin(), m_bcids.end());
  if (bcid_min < 100 && bcid_max > 4000)
    return true;
  else
    return false;
}

inline int TPTrack::BCIDWindow() {

  m_bcids.clear();
  for (int i = 0; i < m_hits.size(); i++)
    m_bcids.push_back(Get(i).BCID());

  if (std::find(m_bcids.begin(), m_bcids.end(), -1) != m_bcids.end())
    return -1;

  if (WrapsAround()){
    m_bcids.clear();
    for (int i = 0; i < m_hits.size(); i++)
      m_bcids.push_back((Get(i).BCID() < 100) ? Get(i).BCID()+4096 : Get(i).BCID());
  }

  int bcid_min = *std::min_element(m_bcids.begin(), m_bcids.end());
  int bcid_max = *std::max_element(m_bcids.begin(), m_bcids.end());
  int bcid_win = bcid_max - bcid_min + 1;

  m_bcids.clear();
  return bcid_win;
}

inline int TPTrack::BCIDEarliest() {
  m_bcids.clear();
  for (int i = 0; i < m_hits.size(); i++)
    m_bcids.push_back(Get(i).BCID());

  if (std::find(m_bcids.begin(), m_bcids.end(), -1) != m_bcids.end())
    return -1;

  if (WrapsAround()){
    m_bcids.clear();
    for (int i = 0; i < m_hits.size(); i++)
      m_bcids.push_back((Get(i).BCID() < 100) ? Get(i).BCID()+4096 : Get(i).BCID());
  }

  int bcid_min = *std::min_element(m_bcids.begin(), m_bcids.end());

  m_bcids.clear();
  return bcid_min;
}

inline int TPTrack::BCIDEarliestA() {
  m_bcids_a.clear();
  m_bcids.clear();
  for (int i = 0; i < m_hits.size(); i++){
    m_bcids.push_back(Get(i).BCID());
    if (Get(i).MMFE8Index() < 4)
      m_bcids_a.push_back(Get(i).BCID());
  }

  if (std::find(m_bcids_a.begin(), m_bcids_a.end(), -1) != m_bcids_a.end())
    return -1;

  if (WrapsAround()){
    m_bcids_a.clear();
    for (int i = 0; i < m_hits.size(); i++)
      m_bcids_a.push_back((Get(i).BCID() < 100) ? Get(i).BCID()+4096 : Get(i).BCID());
  }

  int bcid_min = *std::min_element(m_bcids_a.begin(), m_bcids_a.end());
  
  m_bcids.clear();
  m_bcids_a.clear();
  return bcid_min;
}

inline int TPTrack::BCIDEarliestB() {
  m_bcids_b.clear();
  m_bcids.clear();
  for (int i = 0; i < m_hits.size(); i++){
    m_bcids.push_back(Get(i).BCID());
    if (Get(i).MMFE8Index() > 3)
      m_bcids_b.push_back(Get(i).BCID());
  }
  if (std::find(m_bcids_b.begin(), m_bcids_b.end(), -1) != m_bcids_b.end())
    return -1;

  if (WrapsAround()){
    m_bcids_b.clear();
    for (int i = 0; i < m_hits.size(); i++)
      m_bcids_b.push_back((Get(i).BCID() < 100) ? Get(i).BCID()+4096 : Get(i).BCID());
  }

  int bcid_min = *std::min_element(m_bcids_b.begin(), m_bcids_b.end());

  m_bcids.clear();
  m_bcids_b.clear();
  return bcid_min;
}

inline double TPTrack::BCIDAverage() {
  m_bcids.clear();
  for (int i = 0; i < m_hits.size(); i++)
    m_bcids.push_back(Get(i).BCID());

  if (std::find(m_bcids.begin(), m_bcids.end(), -1) != m_bcids.end())
    return -1;

  if (WrapsAround()){
    m_bcids.clear();
    for (int i = 0; i < m_hits.size(); i++)
      m_bcids.push_back((Get(i).BCID() < 100) ? Get(i).BCID()+4096 : Get(i).BCID());
  }

  double sum = std::accumulate(m_bcids.begin(), m_bcids.end(), 0.0);
  double avg = sum / (double)(m_bcids.size());

  m_bcids.clear();
  return avg;
}

inline double TPTrack::BCIDAverageA() {
  m_bcids.clear();
  m_bcids_a.clear();
  for (int i = 0; i < m_hits.size(); i++){
    m_bcids.push_back(Get(i).BCID());
    if (Get(i).MMFE8Index() < 4)
      m_bcids_a.push_back(Get(i).BCID());
  }

  if (std::find(m_bcids_a.begin(), m_bcids_a.end(), -1) != m_bcids_a.end())
    return -1;

  if (WrapsAround()){
    m_bcids_a.clear();
    for (int i = 0; i < m_hits.size(); i++)
      m_bcids_a.push_back((Get(i).BCID() < 100) ? Get(i).BCID()+4096 : Get(i).BCID());
  }

  double sum = std::accumulate(m_bcids_a.begin(), m_bcids_a.end(), 0.0);
  double avg = sum / (double)(m_bcids_a.size());

  m_bcids_a.clear();
  m_bcids.clear();
  return avg;
}

inline double TPTrack::BCIDAverageB() {
  m_bcids.clear();
  m_bcids_b.clear();
  for (int i = 0; i < m_hits.size(); i++){
    m_bcids.push_back(Get(i).BCID());
    if (Get(i).MMFE8Index() > 3)
      m_bcids_b.push_back(Get(i).BCID());
  }

  if (std::find(m_bcids_b.begin(), m_bcids_b.end(), -1) != m_bcids_b.end())
    return -1;

  if (WrapsAround()){
    m_bcids_b.clear();
    for (int i = 0; i < m_hits.size(); i++)
      m_bcids_b.push_back((Get(i).BCID() < 100) ? Get(i).BCID()+4096 : Get(i).BCID());
  }

  double sum = std::accumulate(m_bcids_b.begin(), m_bcids_b.end(), 0.0);
  double avg = sum / (double)(m_bcids_b.size());

  m_bcids_b.clear();
  m_bcids.clear();
  return avg;
}

inline int TPTrack::BCIDLatest() {
  m_bcids.clear();
  for (int i = 0; i < m_hits.size(); i++)
    m_bcids.push_back(Get(i).BCID());

  if (std::find(m_bcids.begin(), m_bcids.end(), -1) != m_bcids.end())
    return -1;

  if (WrapsAround()){
    m_bcids.clear();
    for (int i = 0; i < m_hits.size(); i++)
      m_bcids.push_back((Get(i).BCID() < 100) ? Get(i).BCID()+4096 : Get(i).BCID());
  }

  int bcid_max = *std::max_element(m_bcids.begin(), m_bcids.end());

  m_bcids.clear();
  return bcid_max;
}

#endif
