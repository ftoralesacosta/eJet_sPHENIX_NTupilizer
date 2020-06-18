#ifndef MYJETANALYSIS_H
#define MYJETANALYSIS_H

#include <fun4all/SubsysReco.h>

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

/// \class MyJetAnalysis
class MyJetAnalysis : public SubsysReco
{
 public:
  MyJetAnalysis(
      const std::string &recojetname = "AntiKt_Track_r08",
      const std::string &truthjetname = "AntiKt_Truth_r08",
      const std::string &outputfilename = "myjetanalysis.root");

  virtual ~MyJetAnalysis();

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


  std::string m_recoJetName;
  std::string m_truthJetName;
  std::string m_outputFileName;

  //! eta range
  std::pair<double, double> m_etaRange;

  //! pT range
  std::pair<double, double> m_ptRange;

  //electron Energy min
  double m_eEmin;

  //electron-jet matching radius for veto
  double m_electronJetMatchingRadius;
  
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
  float m_electron_truthEta;
  float m_electron_truthPhi;
  float m_electron_truthE;
  float m_electron_truthPt;
  float m_electron_truthpX;
  float m_electron_truthpY;
  float m_electron_truthpZ;
  int m_electron_truthPID;

  int m_njets;
  int m_ntruthjets;
  int m_nAlltruthjets;
  
  enum {MaxNumJets = 20};
  std::array<int,MaxNumJets> m_id;
  std::array<int,MaxNumJets> m_nComponent;
  std::array<float,MaxNumJets> m_eta;
  std::array<float,MaxNumJets> m_phi;
  std::array<float,MaxNumJets> m_e;
  std::array<float,MaxNumJets> m_pt;
  
  std::array<int,MaxNumJets> m_matched_truthID;
  std::array<int,MaxNumJets> m_matched_truthNComponent;
  std::array<float,MaxNumJets> m_matched_truthEta;
  std::array<float,MaxNumJets> m_matched_truthPhi;
  std::array<float,MaxNumJets> m_matched_truthE;
  std::array<float,MaxNumJets> m_matched_truthPt;  

  std::array<int,MaxNumJets> m_all_truthID;
  std::array<int,MaxNumJets> m_all_truthNComponent;
  std::array<float,MaxNumJets> m_all_truthEta;
  std::array<float,MaxNumJets> m_all_truthPhi;
  std::array<float,MaxNumJets> m_all_truthE;
  std::array<float,MaxNumJets> m_all_truthPt;  
  
  //Matched Constituens
  // enum {kMaxConstituents = 100};
  // std::array<std::array<int, kMaxConstituents >, MaxNumJets > m_Constituent_truthPID;
  // std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_Constituent_truthEta;
  // std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_Constituent_truthPhi;
  // std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_Constituent_truthPt;
  // std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_Constituent_truthE;


  
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

#endif  // MYJETANALYSIS_H
