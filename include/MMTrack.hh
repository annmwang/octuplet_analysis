///
///  \file   MMTrack.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMTrack_HH
#define MMTrack_HH

class MMTrack {

public:
  MMHit();
  MMHit(int mmfe8, int vmm = -1, int ch = -1);
  MMHit(const MMHit& hit);
  ~MMHit();

  int MMFE8() const;
  int VMM() const;
  int Channel() const;
  int VMMChannel() const;
  int PDO() const;
  int TDO() const;
  int BCID() const;
  int FIFOcount() const;

  double Charge() const;
  double Time() const;

  bool IsChargeCalib() const;
  bool IsTimeCalib() const;
  
  void SetMMFE8(int mmfe8);
  void SetVMM(int vmm);
  void SetChannel(int ch);
  void SetPDO(int pdo);
  void SetTDO(int tdo);
  void SetBCID(int bcid);
  void SetFIFOcount(int fifo);
  void SetCharge(double q);
  void SetTime(double t);
  
private:
  int m_MMFE8;
  int m_VMM;
  int m_CH;
  int m_PDO;
  int m_TDO;
  int m_BCID;
  int m_FIFOcount;

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
  m_FIFOcount = -1;
  m_charge = -1;
  m_time = -1;

  m_charge_calib = false;
  m_time_calib = false;
}

inline MMHit::MMHit(int mmfe8, int vmm, int ch){
  m_MMFE8 = mmfe8;
  m_VMM = vmm;
  m_CH = ch;
  m_PDO = -1;
  m_TDO = -1;
  m_BCID = -1;
  m_FIFOcount = -1;
  m_charge = -1;
  m_time = -1;

  m_charge_calib = false;
  m_time_calib = false;
}

inline MMHit::MMHit(const MMHit& hit){
  m_MMFE8 = hit.MMFE8();
  m_VMM = hit.VMM();
  m_CH = hit.VMMChannel();
  m_PDO = hit.PDO();
  m_TDO = hit.TDO();
  m_BCID = hit.BCID();
  m_FIFOcount = hit.FIFOcount();

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

inline int MMHit::VMM() const {
  return m_VMM;
}

inline int MMHit::Channel() const {
  return std::max(-1,64*m_VMM + m_CH);
}

inline int MMHit::VMMChannel() const {
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

inline int MMHit::FIFOcount() const {
  return m_FIFOcount;
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

inline void MMHit::SetChannel(int ch){
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

