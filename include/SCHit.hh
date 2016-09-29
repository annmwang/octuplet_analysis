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

