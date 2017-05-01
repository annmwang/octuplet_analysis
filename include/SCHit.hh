///
///  \file   SCHit.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef SCHit_HH
#define SCHit_HH

class SCHit {

public:
  SCHit();
  SCHit(int ch, double count = -1);
  SCHit(const SCHit& hit);
  ~SCHit();
  
  int Channel() const;
  double Count() const;

  void SetChannel(int ch);
  void SetCount(double count);
  void CorrectCount();
  
  int PassCountReqs(int RunNumber);

private:
  int m_CH;
  double m_Count;
  
};

#endif

inline SCHit::SCHit(){
  m_CH = -1;
  m_Count = -1;
}

inline SCHit::SCHit(int ch, double count){
  m_CH = ch;
  m_Count = count;
}

inline SCHit::SCHit(const SCHit& hit){
  m_CH = hit.Channel();
  m_Count = hit.Count();
}
  
inline SCHit::~SCHit() {}

inline int SCHit::Channel() const {
  return m_CH;
}

inline double SCHit::Count() const {
  return m_Count;
}

inline void SCHit::SetChannel(int ch){
  m_CH = ch;
}

inline void SCHit::SetCount(double count){
  m_Count = count;
}

inline void SCHit::CorrectCount(){
    if      (m_CH ==  0) m_Count += -0.6;
    else if (m_CH ==  1) m_Count += -1.5;
    else if (m_CH ==  2) m_Count +=  0.8;
    else if (m_CH ==  3) m_Count += -1.5;
    else if (m_CH ==  4) m_Count +=  1.4;
    else if (m_CH ==  5) m_Count +=  0.8;
    else if (m_CH ==  6) m_Count += -0.8;
    else if (m_CH ==  7) m_Count +=  0.1;
    else if (m_CH ==  8) m_Count += -0.0;
    else if (m_CH ==  9) m_Count += -0.8;
    else if (m_CH == 10) m_Count +=  1.4;
    else if (m_CH == 11) m_Count +=  0.5;
    // ----------------------------------
    // no correction to middle scint.
    // ----------------------------------
    else if (m_CH == 16) m_Count += -1.8;
    else if (m_CH == 17) m_Count +=  2.4;
    else if (m_CH == 18) m_Count +=  1.2;
    else if (m_CH == 19) m_Count +=  0.6;
    else if (m_CH == 20) m_Count +=  0.5;
    else if (m_CH == 21) m_Count += -1.6;
    else if (m_CH == 22) m_Count +=  1.1;
    else if (m_CH == 23) m_Count +=  2.0;
    else if (m_CH == 24) m_Count += -0.1;
    else if (m_CH == 25) m_Count += -2.0;
    else if (m_CH == 26) m_Count +=  0.3;
    else if (m_CH == 27) m_Count += -1.4;
}

inline int SCHit::PassCountReqs(int RunNumber) {
  if (RunNumber == 3513) {
    if      (m_CH < 12) return (int)(m_Count >  160.0 && m_Count  < 220.0);
    else if (m_CH > 15) return (int)(m_Count >  180.0 && m_Count  < 240.0);
    else                return (int)(m_Count >= 200.0 && m_Count <= 250.0);
  }
  else if (RunNumber >= 3516) {
    if      (m_CH < 12) return (int)(m_Count >  110.0 && m_Count  < 170.0);
    else if (m_CH > 15) return (int)(m_Count >  130.0 && m_Count  < 190.0);
    else                return (int)(m_Count >= 150.0 && m_Count <= 200.0);
  }
  else{
    std::cout << "Need to add RunNumber settings to include/SCHit.hh! Error!" << std::endl;
    return -1;
  }
  // probably shouldnt compare Counts() with float boundaries,
  // since they want to be ints, but who cares for now
}

