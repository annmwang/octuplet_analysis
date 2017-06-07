///
///  \file   SCEventHits.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef SCEventHits_HH
#define SCEventHits_HH

#include "include/SCHit.hh"

class SCEventHits {

public:
  SCEventHits();
  SCEventHits(int RunNumber);
  SCEventHits(const SCHit& hit);
  SCEventHits(const SCEventHits& evt_hits);
  
  ~SCEventHits();

  bool AddHit(const SCHit& hit);

  bool operator += (const SCHit& hit);
  
  int GetNHits() const;
  
  SCHit const& Get(int ihit) const;
  SCHit const& operator [] (int ihit) const;

  bool IsGoodEvent();
  
  int NBotPair();
  int NTopPair();
  int TPBCID();
  int TPph();

  void BotXY(double& x, double& y);
  void TopXY(double& x, double& y);
  void SetTPBCID(int BCID);
  void SetTPph(int ph);

  double BotTDC();
  double TopTDC();

  std::vector<std::pair<SCHit*,SCHit*> > GetTopPair();
  std::vector<std::pair<SCHit*,SCHit*> > GetBotPair();

  
private:
  std::vector<SCHit*> m_hits;

  SCHit* m_padw;
  SCHit* m_pade;

  // unpaired hits
  std::vector<SCHit*> m_bothits;
  std::vector<SCHit*> m_tophits;

  // matched pairs
  std::vector<std::pair<SCHit*,SCHit*> > m_bot_pairs;
  std::vector<std::pair<SCHit*,SCHit*> > m_top_pairs;

  static double sc_bot_atdc_corr[6];
  static double sc_top_atdc_corr[6];

  int m_RunNumber;
  int m_tpBCID;
  int m_tpPh;
};

inline SCEventHits::SCEventHits(){
  m_padw = nullptr;
  m_pade = nullptr;
}

inline SCEventHits::SCEventHits(int RunNumber){
  m_padw = nullptr;
  m_pade = nullptr;
  m_RunNumber = RunNumber;
}

inline SCEventHits::SCEventHits(const SCHit& hit){
  m_padw = nullptr;
  m_pade = nullptr;
  AddHit(hit);
}

inline SCEventHits::SCEventHits(const SCEventHits& evt_hits){
  m_padw = nullptr;
  m_pade = nullptr;
  int N = evt_hits.GetNHits();
  for(int i = 0; i < N; i++)
    AddHit(evt_hits[i]);
}
  
inline SCEventHits::~SCEventHits(){
  int N = GetNHits();
  for(int i = 0; i < N; i++)
    delete m_hits[i];
}

inline bool SCEventHits::AddHit(const SCHit& hit){
  SCHit* newhit = nullptr;
  int N = GetNHits();
  for(int i = 0; i < N; i++){
    if(hit.Channel() < m_hits[i]->Channel()){
      newhit = new SCHit(hit);
      m_hits.insert(m_hits.begin()+i, newhit);
      break;
    }
  }
  if(!newhit){
    newhit = new SCHit(hit);
    m_hits.push_back(newhit);
  }

  if(newhit->Channel() == 12){
    m_padw = newhit;
  }

  if(newhit->Channel() == 13){
    m_pade = newhit;
  }

  if(newhit->Channel() < 12){
    newhit->CorrectCount(m_RunNumber);
    if (newhit->PassCountReqs(m_RunNumber)){
      int Nbot = m_bothits.size();
      for(int i = 0; i < Nbot; i++){
	if(abs(newhit->Channel()-m_bothits[i]->Channel()) == 6){
	  if(newhit->Channel() < m_bothits[i]->Channel()) { 
	    m_bot_pairs.push_back(std::pair<SCHit*,SCHit*>(newhit, m_bothits[i]));
	  }
	  else { 
	    m_bot_pairs.push_back(std::pair<SCHit*,SCHit*>(m_bothits[i], newhit));
	  }
	  return true;
	}
      }
      m_bothits.push_back(newhit);
    }
  }
  if(newhit->Channel() >= 16 && newhit->Channel() < 28){
    newhit->CorrectCount(m_RunNumber);
    if (newhit->PassCountReqs(m_RunNumber)){
      int Ntop = m_tophits.size();
      for(int i = 0; i < Ntop; i++){
	if(abs(newhit->Channel()-m_tophits[i]->Channel()) == 6){
	  if(newhit->Channel() < m_tophits[i]->Channel()){
	    m_top_pairs.push_back(std::pair<SCHit*,SCHit*>(newhit, m_tophits[i]));
	  }
	  else {
	    m_top_pairs.push_back(std::pair<SCHit*,SCHit*>(m_tophits[i], newhit));
	  }
	  return true;
	}
      }
      m_tophits.push_back(newhit);
    }
  }  
  return true;
}

inline bool SCEventHits::IsGoodEvent(){
  if(!m_padw || !m_pade)
    return false;
  if (!m_padw->PassCountReqs(m_RunNumber) || !m_pade->PassCountReqs(m_RunNumber))
    return false;
  if(NBotPair() != 1 || NTopPair() != 1)
    return false;

  double topx, topy, botx, boty;
  BotXY(botx,boty);
  TopXY(topx,topy);

  if(fabs(topx) > 100. ||
     fabs(botx) > 100.)
    return false;

  double length = sqrt(274.3*274.3 +
		       (botx-topx)*(botx-topx)+
		       (boty-topy)*(boty-topy))/100.;
  // 0.5 ns per count
  double timediff = TopTDC() - BotTDC() - length*2.*3.3;
  if(fabs(timediff) > 20.)
    return false;
  return true;
}

inline int SCEventHits::NBotPair(){
  return m_bot_pairs.size();
}

inline int SCEventHits::NTopPair(){
  return m_top_pairs.size();
}

inline int SCEventHits::TPBCID(){
  return m_tpBCID;
}

inline int SCEventHits::TPph(){
  return m_tpPh;
}

inline void SCEventHits::BotXY(double& x, double& y){
  if(m_bot_pairs.size() <= 0){
    x = 0; y = 0;
    return;
  }
  double tdcdiff = m_bot_pairs[0].first->Count()-
    m_bot_pairs[0].second->Count();

  x = -tdcdiff*3.144;
  y = 50. - m_bot_pairs[0].first->Channel()*20.;
}

inline void SCEventHits::TopXY(double& x, double& y){
  if(m_top_pairs.size() <= 0){
    x = 0; y = 0;
    return;
  }
  double tdcdiff = m_top_pairs[0].first->Count()-
    m_top_pairs[0].second->Count();
  
  x = -tdcdiff*3.144;
  y = 50. - (m_top_pairs[0].first->Channel()-16)*20.;
}

inline void SCEventHits::SetTPBCID(int bcid){
  m_tpBCID = bcid;
}

inline void SCEventHits::SetTPph(int ph){
  m_tpPh = ph;
}

inline double SCEventHits::BotTDC(){
  if(m_bot_pairs.size() <= 0){
    return 0.;
  }
  double tdcav = (m_bot_pairs[0].first->Count()+
		  m_bot_pairs[0].second->Count())/2.;
  tdcav += sc_bot_atdc_corr[m_bot_pairs[0].first->Channel()];
  return tdcav;
}

inline double SCEventHits::TopTDC(){
  if(m_top_pairs.size() <= 0){
    return 0.;
  }
  return sc_top_atdc_corr[m_top_pairs[0].first->Channel()] +
    (m_top_pairs[0].first->Count()+m_top_pairs[0].second->Count())/2.;
}

inline std::vector<std::pair<SCHit*,SCHit*> > SCEventHits::GetTopPair(){
  return m_top_pairs;
}

inline std::vector<std::pair<SCHit*,SCHit*> > SCEventHits::GetBotPair(){
  return m_bot_pairs;
}

inline bool SCEventHits::operator += (const SCHit& hit){
  return AddHit(hit);
}

inline int SCEventHits::GetNHits() const {
  return int(m_hits.size());
}

inline SCHit const& SCEventHits::Get(int ihit) const {
  return *m_hits[ihit];
}

inline SCHit const& SCEventHits::operator [] (int ihit) const {
  return Get(ihit);
}

double SCEventHits::sc_bot_atdc_corr[6] = { -2.8, -2.7, -2.6,
					    -2.6, -2.6, -2.9 };
double SCEventHits::sc_top_atdc_corr[6] = { -1.4, -1.4, -1.4,
					    -1.4, -1.4, -1.4 };
    
#endif
