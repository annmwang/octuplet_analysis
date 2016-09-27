///
///  \file   MMPacmanAlgo.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMPacmanAlgo_HH
#define MMPacmanAlgo_HH

#include "include/MMClusterAlgo.hh"

class MMPacmanAlgo : public MMClusterAlgo {

public:
  MMPacmanAlgo(int clus_size = 5,
	       double seed_thresh = 15.,
	       double hit_thresh  = 10.);
  
  ~MMPacmanAlgo() {}

  MMClusterList Cluster(const MMFE8Hits& hits);

  void SetClusterSize(int clus_size);
  void SetSeedThreshold(double thresh);
  void SetHitThreshold(double thresh);

private:
  int m_clus_size;
  double m_seed_thresh;
  double m_hit_thresh;
  
};

#endif
  
inline MMPacmanAlgo::MMPacmanAlgo(int clus_size,
				  double seed_thresh,
				  double hit_thresh){
  m_clus_size = clus_size;
  m_seed_thresh = seed_thresh;
  m_hit_thresh = hit_thresh;
}

inline MMClusterList MMPacmanAlgo::Cluster(const MMFE8Hits& hits){
  MMClusterList cluster_list;
  
  // forward step
  int Nhit = hits.GetNHits();
  for(int i = 0; i < Nhit; i++){
    if(!IsGoodHit(hits[i]))
      continue;
   
    // new cluster if seed above thresh
    if(hits[i].Charge() > m_seed_thresh){
      MMCluster cluster(hits[i]);
      int last_channel = hits[i].Channel();
      // look for additional hits forward
      for(int j = i+1; j < Nhit; j++){
	if(!IsGoodHit(hits[j])){
	  continue;
	}
	if(hits[j].Channel() <= last_channel+m_clus_size){
	  i = j; // move index so we don't look for seeds in this channel
	  if(hits[j].Charge() > m_hit_thresh){
	    cluster.AddLinkedHit(hits[j]);
	    last_channel = hits[j].Channel();
	  }
	} else {
	  break;
	}
      }
      cluster_list.AddCluster(cluster);
    }
  }
  
  // backward step
  int Nclus = cluster_list.GetNCluster();

  for(int c = 0; c < Nclus; c++){
    // get index of first hit in cluster
    int i = hits.GetIndex(cluster_list[c][0]);
    int first_channel = cluster_list[c][0].Channel();
    for(int j = i-1; j >= 0; j--){
      if(!IsGoodHit(hits[j]))
	continue;
      if(cluster_list.Contains(hits[j]))
	break; // already in another cluster
      if(hits[j].Channel() >= first_channel-m_clus_size){
	if(hits[j].Charge() > m_hit_thresh){
	  cluster_list.AddLinkedHit(hits[j], c);
	  first_channel = hits[j].Channel();
	}
      } else {
	break;
      }
    }
  }
  
  return cluster_list;
}

inline void MMPacmanAlgo::SetClusterSize(int clus_size){
  if(clus_size >= 1)
    m_clus_size = clus_size;
}

inline void MMPacmanAlgo::SetSeedThreshold(double thresh){
  m_seed_thresh = thresh;
}

inline void MMPacmanAlgo::SetHitThreshold(double thresh){
  m_hit_thresh = thresh;
}




