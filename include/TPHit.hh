///
///  \file   TPHit.hh
///
///  \author Ann Miao Wang
///
///  \date   2017 May
///


#ifndef TPHit_HH
#define TPHit_HH

class TPHit {

public:
  TPHit();
  TPHit(int mmfe8, int vmm = -1, double ch = -1, int bcid = -1, int RunNumber = -1);
  TPHit(const TPHit& hit);
  ~TPHit();

  int MMFE8() const;
  int MMFE8Index() const;
  int VMM() const;
  double Channel() const;
  double VMMChannel() const;
  int BCID() const;
  int RunNumber() const;
  int isX() const;
  int isU() const;
  int isV() const;

  void SetMMFE8(int mmfe8);
  void SetMMFE8Index(int RunNumber);
  void SetVMM(int vmm);
  void SetChannel(double ch);
  void SetBCID(int bcid);
  void SetRunNumber(int RunNumber);
  
private:
  int m_MMFE8;
  int m_MMFE8index;
  int m_VMM;
  double m_CH;
  int m_BCID;
  int m_RunNumber;
};

inline TPHit::TPHit(){
  m_MMFE8 = -1;
  m_VMM = -1;
  m_CH = -1;
  m_BCID = -1;
  m_RunNumber = -1;

  SetMMFE8Index(m_RunNumber);
}

inline TPHit::TPHit(int mmfe8, int vmm, double ch, int bcid, int RunNumber){
  m_MMFE8 = mmfe8;
  m_VMM = vmm;
  m_CH = ch;
  m_BCID = bcid;
  m_RunNumber = RunNumber;

  SetMMFE8Index(RunNumber);
}

inline TPHit::TPHit(const TPHit& hit){
  m_MMFE8 = hit.MMFE8();
  m_MMFE8index = hit.MMFE8Index();
  m_VMM = hit.VMM();
  m_CH = hit.VMMChannel();
  m_BCID = hit.BCID();
  m_RunNumber = hit.RunNumber();
}
  
inline TPHit::~TPHit(){

}

inline int TPHit::MMFE8() const {
  return m_MMFE8;
}

inline int TPHit::MMFE8Index() const {
  return m_MMFE8index;
}

inline int TPHit::VMM() const {
  return m_VMM;
}

inline double TPHit::Channel() const {
  return std::max(-1.,64.*m_VMM + m_CH);
}

inline double TPHit::VMMChannel() const {
  return m_CH;
}

inline int TPHit::BCID() const {
  return m_BCID;
}

inline int TPHit::RunNumber() const {
  return m_RunNumber;
}

inline int TPHit::isX() const {
  return (MMFE8Index() == 0 ||
          MMFE8Index() == 1 ||
          MMFE8Index() == 6 ||
          MMFE8Index() == 7);
}

inline int TPHit::isU() const {
  return (MMFE8Index() == 2 ||
          MMFE8Index() == 4);
}

inline int TPHit::isV() const {
  return (MMFE8Index() == 3 ||
          MMFE8Index() == 5);
}

inline void TPHit::SetMMFE8(int mmfe8){
  m_MMFE8 = mmfe8;
}

inline void TPHit::SetVMM(int vmm){
  m_VMM = vmm;
}

inline void TPHit::SetChannel(double ch){
  m_CH = ch;
}

inline void TPHit::SetBCID(int bcid){
  m_BCID = bcid;
}

inline void TPHit::SetRunNumber(int RunNumber){
  m_RunNumber = RunNumber;
}

inline void TPHit::SetMMFE8Index(int RunNumber) {
  // only valid for Run 3513.
  // lets put this info into the ntuple?
  m_MMFE8index = -1;
  if (RunNumber == 3513) {
    if      (m_MMFE8 == 111) m_MMFE8index = 0;
    else if (m_MMFE8 == 116) m_MMFE8index = 1;
    else if (m_MMFE8 == 101) m_MMFE8index = 2;
    else if (m_MMFE8 == 109) m_MMFE8index = 3;
    else if (m_MMFE8 == 117) m_MMFE8index = 4;
    else if (m_MMFE8 == 102) m_MMFE8index = 5;
    else if (m_MMFE8 == 107) m_MMFE8index = 6;
    else if (m_MMFE8 == 105) m_MMFE8index = 7;
  }
  else if (RunNumber == 3516) {
    if      (m_MMFE8 == 111) m_MMFE8index = 0;
    else if (m_MMFE8 == 116) m_MMFE8index = 1;
    else if (m_MMFE8 == 117) m_MMFE8index = 2;
    else if (m_MMFE8 == 119) m_MMFE8index = 3;
    else if (m_MMFE8 == 106) m_MMFE8index = 4;
    else if (m_MMFE8 == 107) m_MMFE8index = 5;
    else if (m_MMFE8 == 118) m_MMFE8index = 6;
    else if (m_MMFE8 == 105) m_MMFE8index = 7;
  }
  else if (RunNumber >= 3518 && RunNumber < 3524) {
    if      (m_MMFE8 == 118) m_MMFE8index = 0;
    else if (m_MMFE8 == 116) m_MMFE8index = 1;
    else if (m_MMFE8 == 102) m_MMFE8index = 2;
    else if (m_MMFE8 == 119) m_MMFE8index = 3;
    else if (m_MMFE8 == 106) m_MMFE8index = 4;
    else if (m_MMFE8 == 107) m_MMFE8index = 5;
    else if (m_MMFE8 == 117) m_MMFE8index = 6;
    else if (m_MMFE8 == 105) m_MMFE8index = 7;
  }
  else if (RunNumber == 3524) {
    if      (m_MMFE8 == 118) m_MMFE8index = 0;
    else if (m_MMFE8 == 116) m_MMFE8index = 1;
    else if (m_MMFE8 == 102) m_MMFE8index = 2;
    else if (m_MMFE8 == 119) m_MMFE8index = 3;
    else if (m_MMFE8 == 106) m_MMFE8index = 4;
    else if (m_MMFE8 == 107) m_MMFE8index = 5;
    else if (m_MMFE8 == 101) m_MMFE8index = 6;
    else if (m_MMFE8 == 105) m_MMFE8index = 7;
  }
  else if (RunNumber >= 3525) {
    if      (m_MMFE8 == 118) m_MMFE8index = 0;
    else if (m_MMFE8 == 111) m_MMFE8index = 1;
    else if (m_MMFE8 == 120) m_MMFE8index = 2;
    else if (m_MMFE8 == 119) m_MMFE8index = 3;
    else if (m_MMFE8 == 106) m_MMFE8index = 4;
    else if (m_MMFE8 == 107) m_MMFE8index = 5;
    else if (m_MMFE8 == 101) m_MMFE8index = 6;
    else if (m_MMFE8 == 105) m_MMFE8index = 7;
  }
  else {
    std::cout << "Need to add RunNumber settings to include/TPHit.hh! Error!" << std::endl;
  }
}

#endif

