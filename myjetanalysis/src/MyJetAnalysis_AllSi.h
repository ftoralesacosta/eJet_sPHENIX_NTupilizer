#ifndef MYJETANALYSIS_ALLSI_H
#define MYJETANALYSIS_ALLSI_H

#include <fun4all/SubsysReco.h>

#include <g4eval/JetTruthEval.h>

#include <set>

#include <memory>
#include <string>
#include <utility>  // std::pair, std::make_pair

#if ! defined(__CINT__) || defined(__CLING__)
#include <array>
#endif  // #ifndef __CINT__

class PHCompositeNode;
class JetEvalStack;
class TTree;
class TH1;

class Jet;
class JetMap;
class SvtxTrackMap;
class PHG4TruthInfoContainer;

/// \class MyJetAnalysis_AllSi
class MyJetAnalysis_AllSi : public SubsysReco
{
 public:
  MyJetAnalysis_AllSi(
      const std::string &recojetname = "AntiKt_Track_r10",
      const std::string &truthjetname = "AntiKt_Truth_r10",
      const std::string &outputfilename = "myjetanalysis.root");

  virtual ~MyJetAnalysis_AllSi();

  //! set eta range
  void
  setEtaRange(double low, double high)
  {
    m_etaRange.first = low;
    m_etaRange.second = high;
  }
  //! set eta range
  void
  setPtRange(double low, double high)
  {
    m_ptRange.first = low;
    m_ptRange.second = high;
  }

  void
  setelectronEmin(double Emin)
  {
    m_eEmin = Emin;  
  }
  
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:

#if ! defined(__CINT__) || defined(__CLING__)

  //! cache the jet evaluation modules
  std::shared_ptr<JetEvalStack> m_jetEvalStack;

  Jet* max_truth_jet_by_track_fastsim(Jet*);
  int get_track_fastsim_contribution(Jet* recojet, Jet* truthjet);
  std::set<PHG4Particle*> all_truth_particles(Jet* truthjet);

  JetMap* _truthjets;
  SvtxTrackMap* _trackmap;
  PHG4TruthInfoContainer* _truthinfo;

  std::string m_recoJetName;
  std::string m_truthJetName;
  std::string m_outputFileName;

  //! eta range
  std::pair<double, double> m_etaRange;

  //! pT range
  std::pair<double, double> m_ptRange;

  //electron Energy min
  double m_eEmin;
  
  //! max track-jet matching radius
  double m_trackJetMatchingRadius;

  //! Output histograms
  TH1 *m_hInclusiveE;
  TH1 *m_hInclusiveEta;
  TH1 *m_hInclusivePhi;
  TH1 *m_hInclusiveNJets;
  
  //! Output Tree variables
  TTree *m_T;

  int m_event;

  //Electron Truth Variables
  float m_etruthEta;
  float m_etruthPhi;
  float m_etruthE;
  float m_etruthPt;
  float m_etruthpX;
  float m_etruthpY;
  float m_etruthpZ;
  int m_etruthPID;

  int m_njets;
  int m_ntruthjets;
  
  enum {MaxNumJets = 10};
  std::array<int,MaxNumJets> m_id;
  std::array<int,MaxNumJets> m_nComponent;
  std::array<float,MaxNumJets> m_eta;
  std::array<float,MaxNumJets> m_phi;
  std::array<float,MaxNumJets> m_e;
  std::array<float,MaxNumJets> m_pt;
  
  std::array<int,MaxNumJets> m_truthID;
  std::array<int,MaxNumJets> m_truthNComponent;
  std::array<float,MaxNumJets> m_truthEta;
  std::array<float,MaxNumJets> m_truthPhi;
  std::array<float,MaxNumJets> m_truthE;
  std::array<float,MaxNumJets> m_truthPt;  

  // // Tracks
  // // ! number of matched tracks
  // std::array<int,MaxNumJets> m_nMatchedTrack;

  // enum
  // {
  //   //! max number of tracks
  //   kMaxMatchedTrack = 1000
  // };
  
  // std::array<std::array<float, kMaxMatchedTrack>, MaxNumJets > m_trackdR;
  // std::array<std::array<float, kMaxMatchedTrack>, MaxNumJets > m_trackpT;

#endif  // #ifndef __CINT__
};

#endif  // MYJETANALYSIS_ALLSI_H
