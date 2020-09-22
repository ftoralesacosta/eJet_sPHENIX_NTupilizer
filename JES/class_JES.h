#ifndef __PHOTONJET_H__
#define __PHOTONJET_H__

//#include <vector>
#include <math.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

class MyClass {

  //Binning
  std::vector<float> E_Bin_Edges; //coarse Energy slices for analysis (gaus)
  std::vector<float> E_Bin_Centers;
  std::vector<float> eta_bins;
  std::vector<float> eta_Centers;

  //Fitting Parameters
  float stdv_range = 1.5;
  float Fit_E_Min;
  float Fit_E_Max;

  float th1_binning[3] = {16,-2,2};; //n,min,max. use for Truth-Reco & Reco/Truth TH1
  int n_slice_bins = 10; //# bins in recoE TH1, within slices of truthE

  /* Jet Cuts */
  int min_N;// Jet constituents
  float min_E;
  float max_Eta;//FIXME: Maybe just set them here? Do the same for binning?
  
  //JES Histos
  std::vector<TH1F*> E_Differences;
  std::vector<TH1F*> E_Ratios;
  std::vector<TH1F*> E_Slices;
  std::vector<TH1F*> Phi_Deltas;
  std::vector<TH1F*> Eta_Deltas;

  //TH2s: TH1s vs Jet E
  std::vector<TH2F*> Ratio_TH2F_v;
  std::vector<TH2F*> Diff_TH2F_v;
  std::vector<TH2F*> Slice_TH2F_v;
  //FIXME: Can these all be put into a struct
  //Don't use vectors. It makes vector[2]->Fill(E_diff) hard to read

  
  TString reco_or_corr;

public:
      
  void set_binning(std::vector<float> &edges, std::vector<float> &centers,
		   float low, float high, float width);

  void set_E_binning(float low, float high, float width);
  void set_eta_binning(float low, float high, float width);
  void set_E_FitRange(float low, float high);
  TString set_reco_or_corr(TString input_string);

  std::vector<TH1F*> create_TH1F(
     TString root_name, TString title);

  std::vector<TH2F*> create_TH2F(
     TString root_name, TString title);

  void initialize_histograms();

  void reconstructed_loop(TTreeReader *Tree);

  void fill_histograms(TLorentzVector *truth, TLorentzVector *reco, int Ebin);

  void fit_histograms();
  
  void write_histograms(TFile *out_file);
  
  //void reco_loop(TTreeReader Tree);

};
#endif
