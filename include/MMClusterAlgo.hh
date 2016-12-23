///
///  \file   MMClusterAlgo.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMClusterAlgo_HH
#define MMClusterAlgo_HH

#include "include/MMClusterList.hh"

class MMClusterAlgo {

public:
  MMClusterAlgo();
  
  ~MMClusterAlgo() {}

  // virtual function must be implemented in derived classes
  virtual MMClusterList Cluster(const MMFE8Hits& hits) = 0;
  
  void SetMaxBCIDDiff(int diff);
  void SetMinBCIDDiff(int diff);
  void SetEventTrigBCID(double trig_bcid);

  void SetMaxTimeDiff(double max_tdo);
  void SetEventPadTime(double pad_time);

  bool IsGoodHit(const MMHit& hit);

private:
  int m_max_BCID_diff;
  int m_min_BCID_diff;
  int m_trig_BCID;

  double m_max_time_diff;
  double m_pad_time;
  
};

inline MMClusterAlgo::MMClusterAlgo(){
  // default max/min difference between trig and event BCID
  m_max_BCID_diff = 40;
  m_min_BCID_diff = 15;
  m_trig_BCID = 0;

  // default max TDO/time upper cut
  m_max_time_diff = 440;
  m_pad_time = 0;
}

inline void MMClusterAlgo::SetMaxBCIDDiff(int diff){
  m_max_BCID_diff = diff;
}

inline void MMClusterAlgo::SetMinBCIDDiff(int diff){
  m_min_BCID_diff = diff;
}
inline void MMClusterAlgo::SetEventTrigBCID(double trig_bcid){
  m_trig_BCID = trig_bcid;
}

inline void MMClusterAlgo::SetMaxTimeDiff(double max_tdo){
  m_max_time_diff = max_tdo;
}

inline void MMClusterAlgo::SetEventPadTime(double pad_time){
  m_pad_time = pad_time;
}

inline bool MMClusterAlgo::IsGoodHit(const MMHit& hit){
  if(!hit.IsChargeCalib())
    return false;
  if(!hit.IsTimeCalib())
    return false;
  if(hit.MMFE8() <= 0)
    return false;
  if(hit.PDO() < 0)
    return false;
  if(hit.TDO() < 0)
    return false;
  if(hit.Charge() < 0)
    return false;

  // BCID
  if(fabs(hit.BCID() - m_trig_BCID) > m_max_BCID_diff)
    return false;
  if(fabs(hit.BCID() - m_trig_BCID) < m_min_BCID_diff)
    return false;

  // TDO - add back
  // if(hit.Time() - 0.5*m_pad_time + 120. > m_max_time_diff)
  //   return false;

  return true;
}

#endif




