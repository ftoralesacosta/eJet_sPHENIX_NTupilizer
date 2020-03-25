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
  int m_id;
  int m_nComponent;
  float m_eta;
  float m_phi;
  float m_e;
  float m_pt;

  int m_subleading_event;
  int m_subleading_id;
  int m_subleading_nComponent;
  float m_subleading_eta;
  float m_subleading_phi;
  float m_subleading_e;
  float m_subleading_pt;
  
  int m_truthID;
  int m_truthNComponent;
  float m_truthEta;
  float m_truthPhi;
  float m_truthE;
  float m_truthPt;

  int m_subleading_truthID;
  int m_subleading_truthNComponent;
  float m_subleading_truthEta;
  float m_subleading_truthPhi;
  float m_subleading_truthE;
  float m_subleading_truthPt;
  
  //! number of matched tracks
  int m_nMatchedTrack;
  int m_subleading_nMatchedTrack;

  //Electron Truth Variables
  float m_etruthEta;
  float m_etruthPhi;
  float m_etruthE;
  float m_etruthPt;

  float m_etruthpX;
  float m_etruthpY;
  float m_etruthpZ;
  int m_etruthPID;
  int m_etruthParentID;

  enum
  {
    //! max number of tracks
    kMaxMatchedTrack = 1000
  };
  std::array<float, kMaxMatchedTrack> m_trackdR;
  std::array<float, kMaxMatchedTrack> m_trackpT;
  std::array<float, kMaxMatchedTrack> m_subleading_trackdR;
  std::array<float, kMaxMatchedTrack> m_subleading_trackpT;

#endif  // #ifndef __CINT__
};

#endif  // MYJETANALYSIS_H
