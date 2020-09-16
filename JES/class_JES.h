#ifndef __PHOTONJET_H__
#define __PHOTONJET_H__

#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TTree.h>
#include <math.h>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

class MyClass {

  //Binning
  std::vector<float> pT_bins;
  std::vector<float> pT_Centers;
  std::vector<float> eta_bins;
  std::vector<float> eta_Centers;
  float Fit_pT_Min;
  float Fit_pT_Max;

  //Reco Jet Cuts
  int min_N;
  float min_E;
  float min_Eta;//FIXME: Maybe just set them here? Do the same for binning?
  std::vector<std::vector<int> > jet_list;
  
  //JES Histos
  std::vector<TH1F*> pT_Differences;
  std::vector<TH1F*> pT_Ratios;
  std::vector<TH1F*> pT_Slices;
  std::vector<TH1F*> Phi_Deltas;
  std::vector<TH1F*> Eta_Deltas;

  std::vector<TH2F*> Ratio_TH2F_v;
  std::vector<TH2F*> Diff_TH2F_v;
  std::vector<TH2F*> Slice_TH2F_v;

  TString reco_or_corr;
  
public:
      
  void set_binning(std::vector<float> &edges, std::vector<float> &centers,
		   float low, float high, float width);

  void set_pT_binning(float low, float high, float width);
  void set_eta_binning(float low, float high, float width);
  void set_pT_FitRange(float low, float high);
  TString set_reco_or_corr(TString input_string);

  std::vector<TH1F*> create_TH1F(
     TString root_name, TString title,float *binning);

  std::vector<TH2F*> create_TH2F(
     TString root_name, TString title,float *binning);

  void initialize_histograms();

  void apply_cuts(TTreeReader Tree);
  
  void reco_loop(TTreeReader Tree);

};

#endif
