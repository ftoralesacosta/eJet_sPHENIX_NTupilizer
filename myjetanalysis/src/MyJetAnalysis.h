#ifndef MYJETANALYSIS_H
#define MYJETANALYSIS_H

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>
#include <utility>  // std::pair, std::make_pair
#include <limits>

#if ! defined(__CINT__) || defined(__CLING__)
#include <array>
#endif  // #ifndef __CINT__

#define NaN numeric_limits<float>::signaling_NaN()

class PHCompositeNode;
class JetEvalStack;
class CaloEvalStack;
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

  
  double
  get_jet_radius_from_string(std::string jetname){
    //Assumes string compatible with G4_Jets.C
    std::string substring = jetname.substr(jetname.find("_r")+2);
    return .1*stof(substring);
  }

  void use_initial_vertex(const bool b = true) {initial_vertex = b;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:

#if ! defined(__CINT__) || defined(__CLING__)

  //! cache the jet evaluation modules
  std::shared_ptr<JetEvalStack> m_jetEvalStack;
//  std::shared_ptr<CaloEvalStack> m_eemcEvalStack;
//  std::shared_ptr<CaloEvalStack> m_cemcEvalStack;
  std::string m_recoJetName;
  std::string m_truthJetName;
  std::string m_outputFileName;

  //! eta range
  std::pair<double, double> m_etaRange;
  //! pT range
  std::pair<double, double> m_ptRange;

  //Jet Resolution Parameter
  double m_jet_R;

  //Electrons
  double m_electronJetMatchingRadius;
  double m_eEmin;
  
  //! flag to use initial vertex in track evaluator
  bool initial_vertex = false;

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
  float m_electron_truthP;
  int m_electron_truthPID;

  //Electron Reco Variables
  float m_electron_recoEta;
  float m_electron_recoPhi;
  float m_electron_recoE;
  float m_electron_recoPt;
  float m_electron_recoP;

  int m_njets;
  int m_nAlltruthjets;

  enum {MaxNumJets = 20, kMaxConstituents = 100};  
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

  //PID float in order to fill branches with NaN correctly.
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_Constituent_recoP;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_Constituent_recoEta;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_Constituent_recoPhi;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_Constituent_recoPt;
  
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_matched_Constituent_truthPID;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_matched_Constituent_truthCharge;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_matched_Constituent_truthEta;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_matched_Constituent_truthPhi;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_matched_Constituent_truthPt;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_matched_Constituent_truthE;

  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_all_Constituent_truthPID;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_all_Constituent_truthCharge;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_all_Constituent_truthEta;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_all_Constituent_truthPhi;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_all_Constituent_truthPt;
  std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_all_Constituent_truthE;
#endif  // #ifndef __CINT__
};

#endif  // MYJETANALYSIS_H
