#ifndef __PHOTONJET_H__
#define __PHOTONJET_H__

#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TTree.h>
#include <math.h>
class MyClass {

  //Physics Variables
  Int_t           m_event;
  Int_t           id;
  Int_t           nComponent;
  Float_t         eta;
  Float_t         phi;
  Float_t         e;
  Float_t         pt;
  Int_t           truthID;
  Int_t           truthNComponent;
  Float_t         truthEta;
  Float_t         truthPhi;
  Float_t         truthE;
  Float_t         truthPt;
  Int_t           nMatchedTrack;

  //Binning
  std::vector<float> pT_bins;
  std::vector<float> pT_Centers;
  std::vector<float> eta_bins;
  std::vector<float> eta_Centers;
  float Fit_pT_Min;
  float Fit_pT_Max;

  //JES Histos
  std::vector<std::vector<TH1F*> > vectors_TH1F_Calib;
  std::vector<TH2F*> vector_TH2F_Calib;
 
public:
  void get_binning(std::vector<float> &edges, std::vector<float> &centers,
		   float low, float high, float width);

  void get_pT_binning(float low, float high, float width);
  void get_eta_binning(float low, float high, float width);
  void set_pT_FitRange(float low, float high);
  void make_TH1F_Calib(TString reco_or_corr);
  void make_TH2F(TString reco_or_corr);
};

#endif  // __PHOTONJET_H__
