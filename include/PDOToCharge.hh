#ifndef PDOToCharge_HH
#define PDOToCharge_HH

#include "include/PDOcalibBase.hh"
#include "include/MMEventHits.hh"

using namespace std;

///////////////////////////////////////////////
// PDOToCharge class
//
// Class takes PDOcalibBase input
// and provide methods for calibrated
// PDO to charge conversion
///////////////////////////////////////////////

class PDOToCharge {

public:
  PDOToCharge();
  PDOToCharge(const string& PDOcalib_filename);
  ~PDOToCharge();

  // returns charge in fC
  double GetCharge(double PDO, int MMFE8, int VMM, int CH) const;

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
  
  vector<double> m_c0;
  vector<double> m_A2;
  vector<double> m_t02;
  vector<double> m_d21;
  vector<double> m_chi2;
  vector<double> m_prob;

  void PrintError(int MMFE8, int VMM, int CH) const {
    cout << "PDOToCharge ERROR: ";
    cout << "No parameters for requested MMFE8 = " << MMFE8;
    cout << " VMM = " << VMM << endl;
    cout << " CH = " << CH << endl;
  }

};

#endif

inline PDOToCharge::PDOToCharge() {}

inline PDOToCharge::PDOToCharge(const string& PDOcalib_filename){

  TChain tree("PDO_calib");
  tree.AddFile(PDOcalib_filename.c_str());
  PDOcalibBase base(&tree);
  
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
      m_c0.push_back(base.c0);
      m_A2.push_back(base.A2);
      m_t02.push_back(base.t02);
      m_d21.push_back(base.d21);
      m_chi2.push_back(base.chi2);
      m_prob.push_back(base.prob);
    } else {
      int c = m_CH_to_index[index][base.CH];
      m_c0[c] = base.c0;
      m_A2[c] = base.A2;
      m_t02[c] = base.t02;
      m_d21[c] = base.d21;
      m_chi2[c] = base.chi2;
      m_prob[c] = base.prob;
    }
  }
}

inline PDOToCharge::~PDOToCharge(){}

// returns charge in fC
inline double PDOToCharge::GetCharge(double PDO, int MMFE8, int VMM, int CH) const {
    pair<int,int> key(MMFE8,VMM);
    if(m_MMFE8VMM_to_index.count(key) == 0){
      //PrintError(MMFE8,VMM,CH);
      return -1.;
    }
    int i = m_MMFE8VMM_to_index[key];

    if(m_CH_to_index[i].count(CH) == 0){
      //PrintError(MMFE8,VMM,CH);
      return -1.;
    }
    int c = m_CH_to_index[i][CH];

    // Pedestal unrealistic value                                                                                                                          
    if (fabs(m_c0[c]-m_A2[c]*m_d21[c]*(m_d21[c]+2.*m_t02[c])) > 200.){
      //PrintError(MMFE8,VMM,CH);
      return -1.;
    }

    // Gain unrealistic value                                                                                                                              
    if ((2.*m_A2[c]*m_d21[c] > 20.) | (2.*m_A2[c]*m_d21[c] < 5.)){
      //PrintError(MMFE8,VMM,CH);
      return -1.;
    }

    // PDO above fit saturation
    if(PDO >= m_c0[c])
      return m_t02[c];

    // PDO in quadratic part
    if(PDO > m_c0[c] + m_A2[c]*m_d21[c]*m_d21[c])
      return -1.*sqrt(max(0., (PDO-m_c0[c])/m_A2[c])) + m_t02[c];

    // linear part
    return 0.5*( (PDO-m_c0[c])/m_A2[c]/m_d21[c] + m_d21[c] + 2.*m_t02[c] );
  }

// returns chi2 from PDO v charge fit
inline double PDOToCharge::GetFitChi2(int MMFE8, int VMM, int CH) const {
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
inline double PDOToCharge::GetFitProb(int MMFE8, int VMM, int CH) const {
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

inline void PDOToCharge::Calibrate(MMEventHits& evt_hits) const {
  int NBoards = evt_hits.GetNBoards();
  for(int i = 0; i < NBoards; i++)
    Calibrate(*evt_hits.m_boards[i]);
}

inline void PDOToCharge::Calibrate(MMFE8Hits& hits) const {
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

inline void PDOToCharge::Calibrate(MMHit& hit) const {
  hit.SetCharge(GetCharge(hit.PDO(), hit.MMFE8(), hit.VMM(), hit.VMMChannel()));
}
