///
///  \file   MMHit.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMHit_HH
#define MMHit_HH

class MMHit {

public:
  MMHit();
  MMHit(int mmfe8, int vmm = -1, double ch = -1, int RunNumber = -1);
  MMHit(const MMHit& hit);
  ~MMHit();

  int MMFE8() const;
  int MMFE8Index() const;
  int VMM() const;
  double Channel() const;
  double VMMChannel() const;
  int PDO() const;
  int TDO() const;
  int BCID() const;
  int TrigBCID() const;
  int TrigPhase() const;
  int FIFOcount() const;
  int RunNumber() const;

  double Charge() const;
  double Time() const;

  bool IsChargeCalib() const;
  bool IsTimeCalib() const;
  
  void SetMMFE8(int mmfe8);
  void SetMMFE8Index(int RunNumber);
  void SetVMM(int vmm);
  void SetChannel(double ch);
  void SetPDO(int pdo);
  void SetTDO(int tdo);
  void SetBCID(int bcid);
  void SetTrigBCID(int bcid);
  void SetTrigPhase(int phase);
  void SetRunNumber(int RunNumber);
  void SetFIFOcount(int fifo);
  void SetCharge(double q);
  void SetTime(double t);
  
private:
  int m_MMFE8;
  int m_MMFE8index;
  int m_VMM;
  double m_CH;
  int m_PDO;
  int m_TDO;
  int m_BCID;
  int m_trigBCID;
  int m_trigphase;
  int m_FIFOcount;
  int m_RunNumber;

  double m_charge;
  double m_time;

  bool m_charge_calib;
  bool m_time_calib;

  friend class PDOToCharge;
  friend class TDOToTime;
};

#endif

inline MMHit::MMHit(){
  m_MMFE8 = -1;
  m_VMM = -1;
  m_CH = -1;
  m_PDO = -1;
  m_TDO = -1;
  m_BCID = -1;
  m_trigBCID = -1;
  m_trigphase = -1;
  m_FIFOcount = -1;
  m_RunNumber = -1;
  m_charge = -1;
  m_time = -1;

  m_charge_calib = false;
  m_time_calib = false;

  SetMMFE8Index(m_RunNumber);
}

inline MMHit::MMHit(int mmfe8, int vmm, double ch, int RunNumber){
  m_MMFE8 = mmfe8;
  m_VMM = vmm;
  m_CH = ch;
  m_PDO = -1;
  m_TDO = -1;
  m_BCID = -1;
  m_trigBCID = -1;
  m_trigphase = -1;
  m_FIFOcount = -1;
  m_RunNumber = RunNumber;
  m_charge = -1;
  m_time = -1;

  m_charge_calib = false;
  m_time_calib = false;

  SetMMFE8Index(RunNumber);
}

inline MMHit::MMHit(const MMHit& hit){
  m_MMFE8 = hit.MMFE8();
  m_MMFE8index = hit.MMFE8Index();
  m_VMM = hit.VMM();
  m_CH = hit.VMMChannel();
  m_PDO = hit.PDO();
  m_TDO = hit.TDO();
  m_BCID = hit.BCID();
  m_trigBCID = hit.TrigBCID();
  m_trigphase = hit.TrigPhase();
  m_FIFOcount = hit.FIFOcount();
  m_RunNumber = hit.RunNumber();

  m_charge = hit.Charge();
  m_time = hit.Time();

  m_charge_calib = hit.IsChargeCalib();
  m_time_calib = hit.IsTimeCalib();
}
  
inline MMHit::~MMHit(){

}

inline int MMHit::MMFE8() const {
  return m_MMFE8;
}

inline int MMHit::MMFE8Index() const {
  return m_MMFE8index;
}

inline int MMHit::VMM() const {
  return m_VMM;
}

inline double MMHit::Channel() const {
  return std::max(-1.,64.*m_VMM + m_CH);
}

inline double MMHit::VMMChannel() const {
  return m_CH;
}

inline int MMHit::PDO() const {
  return m_PDO;
}

inline int MMHit::TDO() const {
  return m_TDO;
}

inline int MMHit::BCID() const {
  return m_BCID;
}

inline int MMHit::TrigBCID() const {
  return m_trigBCID;
}

inline int MMHit::TrigPhase() const {
  return m_trigphase;
}

inline int MMHit::FIFOcount() const {
  return m_FIFOcount;
}

inline int MMHit::RunNumber() const {
  return m_RunNumber;
}

inline double MMHit::Charge() const {
  return m_charge;
}

inline double MMHit::Time() const {
  return m_time;
}

inline bool MMHit::IsChargeCalib() const {
  return m_charge_calib;
}

inline bool MMHit::IsTimeCalib() const {
  return m_time_calib;
}

inline void MMHit::SetMMFE8(int mmfe8){
  m_MMFE8 = mmfe8;
}

inline void MMHit::SetVMM(int vmm){
  m_VMM = vmm;
}

inline void MMHit::SetChannel(double ch){
  m_CH = ch;
}

inline void MMHit::SetPDO(int pdo){
  m_PDO = pdo;
  m_charge = -1;
  m_charge_calib = false;
}

inline void MMHit::SetTDO(int tdo){
  m_TDO = tdo;
  m_time = -1;
  m_time_calib = false;
}

inline void MMHit::SetBCID(int bcid){
  m_BCID = bcid;
}

inline void MMHit::SetTrigBCID(int bcid){
  m_trigBCID = bcid;
}

inline void MMHit::SetTrigPhase(int phase){
  m_trigphase = phase;
}

inline void MMHit::SetRunNumber(int RunNumber){
  m_RunNumber = RunNumber;
}

inline void MMHit::SetFIFOcount(int fifo){
  m_FIFOcount = fifo;
}

inline void MMHit::SetCharge(double q){
  m_charge = q;
  m_charge_calib = true;
}

inline void MMHit::SetTime(double t){
  m_time = t;
  m_time_calib = true;
}

inline void MMHit::SetMMFE8Index(int RunNumber) {
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
  else if (RunNumber >= 3525 && RunNumber < 3540) {
    if      (m_MMFE8 == 118) m_MMFE8index = 0;
    else if (m_MMFE8 == 111) m_MMFE8index = 1;
    else if (m_MMFE8 == 120) m_MMFE8index = 2;
    else if (m_MMFE8 == 119) m_MMFE8index = 3;
    else if (m_MMFE8 == 106) m_MMFE8index = 4;
    else if (m_MMFE8 == 107) m_MMFE8index = 5;
    else if (m_MMFE8 == 101) m_MMFE8index = 6;
    else if (m_MMFE8 == 105) m_MMFE8index = 7;
  }
  else if (RunNumber >= 3540) {
    if      (m_MMFE8 == 119) m_MMFE8index = 0;
    else if (m_MMFE8 == 124) m_MMFE8index = 1;
    else if (m_MMFE8 == 122) m_MMFE8index = 2;
    else if (m_MMFE8 == 126) m_MMFE8index = 3;
    else if (m_MMFE8 == 106) m_MMFE8index = 4;
    else if (m_MMFE8 == 109) m_MMFE8index = 5;
    else if (m_MMFE8 == 125) m_MMFE8index = 6;
    else if (m_MMFE8 == 123) m_MMFE8index = 7;
  }
  else {
    std::cout << "Need to add RunNumber settings to include/MMHit.hh! Error! You gave: " << RunNumber << std::endl;
  }
}

