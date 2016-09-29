#ifndef TDOToTime_HH
#define TDOToTime_HH

#include "include/TDOcalibBase.hh"

using namespace std;

///////////////////////////////////////////////
// TDOToTime class
//
// Class takes TDOcalibBase input
// and provide methods for calibrated
// TDO to time conversion
///////////////////////////////////////////////

class TDOToTime {

public:
  TDOToTime();
  
  TDOToTime(const string& TDOcalib_filename);

  ~TDOToTime();

  // returns charge in fC
  double GetTime(double TDO, int MMFE8, int VMM, int CH) const;

  double GetTimeDefault(double TDO) const;

  // returns chi2 from PDO v charge fit
  double GetFitChi2(int MMFE8, int VMM, int CH) const;

  // returns probability from PDO v charge fit
  double GetFitProb(int MMFE8, int VMM, int CH) const;

  void Calibrate(MMEventHits& evt_hits) const;
  void Calibrate(MMFE8Hits& hits) const;
  void Calibrate(MMHit& hit) const;

private:
  mutable map<pair<int,int>, int> m_MMFE8VMM_to_index;
  mutable vector<map<int,int> > m_CH_to_index;

  double m_Cdef;
  double m_Sdef;
  
  vector<double> m_C;
  vector<double> m_S;
  vector<double> m_chi2;
  vector<double> m_prob;

  void PrintError(int MMFE8, int VMM, int CH) const {
    cout << "TDOToTime ERROR: ";
    cout << "No parameters for requested MMFE8 = " << MMFE8;
    cout << " VMM = " << VMM << endl;
    cout << " CH = " << CH << endl;
  }

};

#endif

inline TDOToTime::TDOToTime(){
  m_Cdef = 12.;
  m_Sdef = 1.3;
}
  
inline TDOToTime::TDOToTime(const string& TDOcalib_filename){

  m_Cdef = 12.;
  m_Sdef = 1.3;

  TChain tree("TDO_calib");
  tree.AddFile(TDOcalib_filename.c_str());
  TDOcalibBase base(&tree);
    
  int Nentry = base.fChain->GetEntries();

  for(int i = 0; i < Nentry; i++){
    base.GetEntry(i);

    pair<int,int> key(base.MMFE8, base.VMM);

    if(m_MMFE8VMM_to_index.count(key) == 0){
      m_MMFE8VMM_to_index[key] = m_CH_to_index.size();
      m_CH_to_index.push_back(map<int,int>());
    }

    int index = m_MMFE8VMM_to_index[key];

    if(m_CH_to_index[index].count(base.CH) == 0){
      m_CH_to_index[index][base.CH] = m_prob.size();
      m_C.push_back(base.C);
      m_S.push_back(base.S);
      m_chi2.push_back(base.chi2);
      m_prob.push_back(base.prob);
    } else {
      int c = m_CH_to_index[index][base.CH];
      m_C[c] = base.C;
      m_S[c] = base.S;
      m_chi2[c] = base.chi2;
      m_prob[c] = base.prob;
    }
  }
}

inline TDOToTime::~TDOToTime(){}

// returns charge in fC
inline double TDOToTime::GetTime(double TDO, int MMFE8, int VMM, int CH) const {
  pair<int,int> key(MMFE8,VMM);
  if(m_MMFE8VMM_to_index.count(key) == 0){
    //PrintError(MMFE8,VMM,CH);
    return GetTimeDefault(TDO);
  }
  int i = m_MMFE8VMM_to_index[key];

  if(m_CH_to_index[i].count(CH) == 0){
    //PrintError(MMFE8,VMM,CH);
    return GetTimeDefault(TDO);
  }
  int c = m_CH_to_index[i][CH];

  return (TDO-m_C[c])/m_S[c];
}

inline double TDOToTime::GetTimeDefault(double TDO) const {
  return (TDO-m_Cdef)/m_Sdef;
}

// returns chi2 from PDO v charge fit
inline double TDOToTime::GetFitChi2(int MMFE8, int VMM, int CH) const {
  pair<int,int> key(MMFE8,VMM);
  if(m_MMFE8VMM_to_index.count(key) == 0){
    PrintError(MMFE8,VMM,CH);
    return 0.;
  }
  int i = m_MMFE8VMM_to_index[key];

  if(m_CH_to_index[i].count(CH) == 0){
    PrintError(MMFE8,VMM,CH);
    return 0.;
  }
  int c = m_CH_to_index[i][CH];
  return m_chi2[c];
}

// returns probability from PDO v charge fit
inline double TDOToTime::GetFitProb(int MMFE8, int VMM, int CH) const {
  pair<int,int> key(MMFE8,VMM);
  if(m_MMFE8VMM_to_index.count(key) == 0){
    PrintError(MMFE8,VMM,CH);
    return 0.;
  }
  int i = m_MMFE8VMM_to_index[key];

  if(m_CH_to_index[i].count(CH) == 0){
    PrintError(MMFE8,VMM,CH);
    return 0.;
  }
  int c = m_CH_to_index[i][CH];
  return m_prob[c];
}

inline void TDOToTime::Calibrate(MMEventHits& evt_hits) const {
  int NBoards = evt_hits.GetNBoards();
  for(int i = 0; i < NBoards; i++)
    Calibrate(*evt_hits.m_boards[i]);
}

inline void TDOToTime::Calibrate(MMFE8Hits& hits) const {
  int NHits = hits.GetNHits();
  for(int i = 0; i < NHits; i++){
    MMLinkedHit* hit_ptr = hits.m_hits[i];
    int Ndup = hit_ptr->GetNHits();
    while(Ndup > 0){
      Calibrate(*hit_ptr);
      hit_ptr = hit_ptr->m_next;
      Ndup--;
    }
  }
}

inline void TDOToTime::Calibrate(MMHit& hit) const {
  hit.SetTime(GetTime(hit.TDO(), hit.MMFE8(), hit.VMM(), hit.VMMChannel()));
}
