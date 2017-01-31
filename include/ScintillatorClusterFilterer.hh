///
///  \file   ScintillatorClusterFilterer.hh
///
///  \author Tuna
///          (tuna@cern.ch)
///
///  \date   2017 Jan
///

#ifndef ScintillatorClusterFilterer_HH
#define ScintillatorClusterFilterer_HH

#include <numeric>
#include "include/GeoPlane.hh"
#include "include/SimpleTrackFitter.hh"

class ScintillatorClusterFilterer : public SimpleTrackFitter {

public:
  ScintillatorClusterFilterer();
  ~ScintillatorClusterFilterer();

  MMClusterList FilterClustersScint(MMClusterList& clusters, 
                                    GeoOctuplet& geometry,
                                    int scint_bot,
                                    int evt = -1,
                                    int debug = 0
                                    );

private:

  MMClusterList FilterClusters(MMClusterList& clusters,
                               double slope,
                               double road_xx,
                               double road_uv);
  
  MMClusterList ChooseClosestPerBoard(MMClusterList& clusters,
                                      MMTrack track);
  
  int LargestBoardSeparation(MMClusterList& clusters);

  double CalculateAverageSlope(MMClusterList& clusters);
  double xpos(MMCluster* clus);
  double zpos(MMCluster* clus);
  double xdiff(MMCluster* clus, MMTrack track);

  MMClusterList* m_clusters;
  const GeoOctuplet* m_geometry;
  GeoPlane m_plane;
  MMTrack m_track_input;
  std::vector<double> m_scint_slope_mean;
  std::vector<double> m_diff_smallest;

  double m_slope_scint;
  double m_slope_avgd;

  MMClusterList m_clusters_within_road_scint;
  MMClusterList m_clusters_within_road_avgd;

  double m_slope;
  std::vector<double> m_slopes;
  std::vector<double> m_weights;

  int m_nchilds;
  int m_nboards;
  int m_mmfe;
  double m_sum_resid;
  double m_sum_resid_test;

  int m_ibo;
  double m_road;
  double m_xproj;

  int m_event;
  int m_debug;
  int iter;
};

inline ScintillatorClusterFilterer::ScintillatorClusterFilterer(){
  m_clusters = nullptr;
  m_geometry = nullptr;

  // xz-slopes from paolo
  m_scint_slope_mean.push_back(-0.287);
  m_scint_slope_mean.push_back(-0.194);
  m_scint_slope_mean.push_back(-0.105);
  m_scint_slope_mean.push_back(-0.017);
  m_scint_slope_mean.push_back( 0.070);
  m_scint_slope_mean.push_back( 0.158);
}

inline ScintillatorClusterFilterer::~ScintillatorClusterFilterer(){
}

inline MMClusterList ScintillatorClusterFilterer::FilterClustersScint(MMClusterList& clusters, 
                                                                           GeoOctuplet&   geometry,
                                                                           int scint_bot, int evt, int debug){
  
  // setup
  m_event = evt;
  m_debug = debug;
  m_clusters_within_road_scint.Reset();
  m_clusters_within_road_avgd.Reset();
  m_slope_scint = m_scint_slope_mean[scint_bot];
  m_geometry = &geometry;

  // filter all clusters based on scintillator slope and road size
  for (auto clus: FilterClusters(clusters, m_slope_scint, 10.0, 15.0))
    m_clusters_within_road_scint.AddCluster(*clus);

  if (m_debug)
    std::cout << Form("%40s | %6i | %2i clusters within scint. road, scint slope %.4f",
                      "Scint::FitScint", m_event, (int)(m_clusters_within_road_scint.size()), m_slope_scint) << std::endl;

  // calculate average slope of filtered clusters
  m_slope_avgd = CalculateAverageSlope(m_clusters_within_road_scint);

  // start over, with filtering based on averaged slope and tighter road
  for (auto clus: FilterClusters(clusters, m_slope_avgd, 5.0, 5.0))
    m_clusters_within_road_avgd.AddCluster(*clus);

  if (m_debug)
    std::cout << Form("%40s | %6i | Summary: %2i clusters before, %2i clusters after, avgd slope = %.4f",
                      "Scint::FitScint", m_event, (int)(clusters.size()), (int)(m_clusters_within_road_avgd.size()), m_slope_avgd) << std::endl;

  // done!
  return m_clusters_within_road_avgd;

}

inline MMClusterList ScintillatorClusterFilterer::FilterClusters(MMClusterList& clusters,
                                                                      double slope,
                                                                      double road_xx,
                                                                      double road_uv){
  // setup
  MMTrack track_road;
  MMClusterList clusters_road;
  MMClusterList clusters_test;
  MMClusterList clusters_test_1perboard;
  m_nchilds = -1;
  m_nboards = -1;
  m_mmfe    = 99;
  m_sum_resid = 1000;

  if (m_debug){
    std::cout << Form("%40s | %6i | slope: %7.4f | road_xx: %7.3f | road_uv: %7.2f",
                      "Scint::FilterClusters", m_event, slope, road_xx, road_uv) << std::endl;
  }

  // loop over clusters to find a seed
  for (auto clus: clusters){

    // reset stuff
    track_road.Reset();
    clusters_test.Reset();
    clusters_test_1perboard.Reset();

    // configure the reference track to land on the seed candidate
    track_road.SetSlopeX(slope);
    track_road.SetConstX(xpos(clus) - slope*zpos(clus));

    // check for children
    for (auto testclus: clusters){
      m_road = (testclus->isX()) ? road_xx : road_uv;
      m_xproj = track_road.SlopeX()*zpos(testclus) + track_road.ConstX();
      if (fabs(m_xproj - xpos(testclus)) < m_road){
        clusters_test.AddCluster(*testclus);
      }
      if (m_debug){
        std::cout << Form("%40s | %6i | cand @ %7.2f %7.2f | test @ %7.2f %7.2f | track m = %7.4f b = %7.2f | proj @ %7.2f | diff %7.2f",
                          "Scint::FilterClusters", m_event, xpos(clus), zpos(clus), xpos(testclus), zpos(testclus), 
                          track_road.SlopeX(), track_road.ConstX(), m_xproj, fabs(m_xproj - xpos(testclus))) << std::endl;
      }
    }

    if (m_debug)
      std::cout << Form("%40s | %6i | %i clusters for seed candidate (unlimited per board)",
                        "Scint::FilterClusters", m_event, (int)(clusters_test.size())) << std::endl;

    // keep at most 1 cluster per board
    for (auto clus_close: ChooseClosestPerBoard(clusters_test, track_road))
      clusters_test_1perboard.AddCluster(*clus_close);

    if (m_debug)
      std::cout << Form("%40s | %6i | %i clusters for seed candidate (max 1 per board)",
                        "Scint::FilterClusters", m_event, (int)(clusters_test_1perboard.size())) << std::endl;

    // check the sum of residuals in case there are many seed candidates
    // NB: dont need to normalize by N(children), since we only compare residuals
    //     when the seeds have the same number of children
    m_sum_resid_test = 0.0;
    for (auto testclus: clusters_test_1perboard){
      m_xproj = track_road.SlopeX()*zpos(testclus) + track_road.ConstX();
      m_sum_resid_test += fabs(m_xproj - xpos(testclus));
    }

    // does this qualify as a new seed?
    if (((int)(clusters_test_1perboard.size())  > m_nchilds) ||
        ((int)(clusters_test_1perboard.size()) == m_nchilds && LargestBoardSeparation(clusters_test_1perboard)  > m_nboards) ||
        ((int)(clusters_test_1perboard.size()) == m_nchilds && LargestBoardSeparation(clusters_test_1perboard) == m_nboards && m_sum_resid_test < m_sum_resid)
        ){
      clusters_road.Reset();
      for (auto child: clusters_test_1perboard)
        clusters_road.AddCluster(*child);
      m_nboards   = LargestBoardSeparation(clusters_road);
      m_nchilds   = (int)(clusters_road.size());
      m_mmfe      = m_geometry->Index(clus->MMFE8());
      m_sum_resid = m_sum_resid_test;
      if (m_debug)
        std::cout << Form("%40s | %6i | !!!! NEW SEED !!!! board = %i, children = %i, sum(resid.) = %.4f",
                          "Scint::FilterClusters", m_event, m_geometry->Index(clus->MMFE8()), m_nchilds, m_sum_resid) << std::endl;
    }
    else {
      if (m_debug)
        std::cout << Form("%40s | %6i | Seed skipped on board %i. %i vs. %i children, %i vs. %i z-separation of boards, %.4f vs. %.4f sum(resid)",
                          "Scint::FilterClusters", m_event, m_geometry->Index(clus->MMFE8()),
                          (int)(clusters_test_1perboard.size()), m_nchilds,
                          LargestBoardSeparation(clusters_test_1perboard), m_nboards,
                          m_sum_resid_test, m_sum_resid) << std::endl;
    }
  }

  if (m_debug)
    std::cout << Form("%40s | %6i | Seed, ultimately: board %i, %i clusters",
                      "Scint::FilterClusters", m_event, m_mmfe, (int)(clusters_road.size())) << std::endl;

  return clusters_road;
}

inline MMClusterList ScintillatorClusterFilterer::ChooseClosestPerBoard(MMClusterList& clusters,
                                                                             MMTrack track){
  MMClusterList closest; 
  closest.Reset();

  // to start: everything is far away
  m_diff_smallest.clear();
  for (iter = 0; iter < m_geometry->GetNPlanes(); iter++)
    m_diff_smallest.push_back(9e2);

  // tally the closest distance between cluster and track per board
  for (auto clus: clusters){
    m_ibo = m_geometry->Index(clus->MMFE8());
    m_diff_smallest[m_ibo] = std::min(m_diff_smallest[m_ibo], xdiff(clus, track));
  }

  if (m_debug){
    for (iter = 0; iter < m_geometry->GetNPlanes(); iter++)
      std::cout << Form("%40s | %6i | %2i %7.3f", 
                        "Scint::ChooseClosestPerBoard", m_event, iter, m_diff_smallest[iter]) << std::endl;
  }

  // save the closest (beware floating point comparisons)
  for (auto clus: clusters){
    m_ibo = m_geometry->Index(clus->MMFE8());
    if (     xdiff(clus, track) < m_diff_smallest[m_ibo] ||
        fabs(xdiff(clus, track) - m_diff_smallest[m_ibo]) < 1e-6
       ){
      closest.AddCluster(*clus);
    }
  }

  if (m_debug){
    double sum_of_resid = 0.0;
    for (auto resid: m_diff_smallest)
      if (resid < 5e2)
        sum_of_resid += resid;
    std::cout << Form("%40s | %6i | Sum of residuals: %7.3f",
                      "Scint::ChooseClosestPerBoard", m_event, sum_of_resid) << std::endl;
  }

  return closest;
}

inline double ScintillatorClusterFilterer::xpos(MMCluster* clus){
  m_plane = m_geometry->Get(m_geometry->Index(clus->MMFE8()));
  return (m_plane.Origin().X() + m_plane.LocalXatYbegin(clus->Channel()) +
          m_plane.Origin().X() + m_plane.LocalXatYend(  clus->Channel())) / 2.0;
}

inline double ScintillatorClusterFilterer::zpos(MMCluster* clus){
  m_plane = m_geometry->Get(m_geometry->Index(clus->MMFE8()));
  return m_plane.Origin().Z();
}

inline double ScintillatorClusterFilterer::xdiff(MMCluster* clus, MMTrack track){
  return fabs(xpos(clus) - (track.SlopeX()*zpos(clus) + track.ConstX()));
}

inline int ScintillatorClusterFilterer::LargestBoardSeparation(MMClusterList& clusters){
  int sep = 0;
  for (auto clus1: clusters){
    for (auto clus2: clusters){
      sep = std::max(sep,
                     abs(m_geometry->Index(clus1->MMFE8()) -
                         m_geometry->Index(clus2->MMFE8())));
    }
  }
  return sep;
}

inline double ScintillatorClusterFilterer::CalculateAverageSlope(MMClusterList& clusters){
  m_slopes.clear();
  m_weights.clear();

  if (m_debug){
    std::cout << Form("%40s | %6i | %i input clusters", 
                      "Scint::CalculateAverageSlope",  m_event, (int)(clusters.size())) << std::endl;
  }

  for (auto clus1: clusters){
    for (auto clus2: clusters){

      // to avoid double-counting, require the pairs be z-ordered
      if (m_geometry->Index(clus1->MMFE8()) > m_geometry->Index(clus2->MMFE8()))
        continue;

      // require clusters to be >20mm away
      if (fabs(zpos(clus1) - zpos(clus2)) < 20)
        continue;

      // weight the average by z-distance
      m_slopes.push_back( (xpos(clus1)-xpos(clus2)) / (zpos(clus1)-zpos(clus2)) );
      m_weights.push_back( fabs(zpos(clus1)-zpos(clus2)) );

      // debug
      if (m_debug)
        std::cout << Form("%40s | %6i | %2i @ %7.3f %7.3f | %2i @ %7.3f %7.3f | %7.4f %7.4f", 
                          "Scint::CalculateAverageSlope", m_event,
                          m_geometry->Index(clus1->MMFE8()), xpos(clus1), zpos(clus1),
                          m_geometry->Index(clus2->MMFE8()), xpos(clus2), zpos(clus2),
                          m_slopes.back(), m_weights.back()) << std::endl;
    }
  }

  // calculate the average
  m_slope = 0.0;
  for (iter = 0; iter < (int)(m_slopes.size()); iter++)
    m_slope += m_weights[iter]*m_slopes[iter];
  m_slope /= std::accumulate(m_weights.begin(), m_weights.end(), 0.0);

  if (m_debug)
    std::cout << Form("%40s | %6i | %7.4f", "Scint::CalculateAverageSlope", m_event, m_slope) << std::endl;

  return m_slope;
}

#endif
