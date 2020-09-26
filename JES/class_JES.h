#ifndef __PHOTONJET_H__
#define __PHOTONJET_H__

//#include <vector>
#include <math.h>

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TMath.h>
#include <TTree.h>
#include <TLorentzVector.h>
//#include <LorentzVector.>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

class MyClass {

  //Binning
  int num_E_Bins;
  float E_Max = 80.; //GeV
  float E_Min = 2.; 
  float E_Bin_Width = 1.0;
  std::vector<float> E_Bin_Edges; //coarse Energy slices for analysis (gaus)
  std::vector<float> E_Bin_Centers;

  int num_eta_Bins;
  float eta_Max = 10;
  float eta_Min = -10;
  float eta_Bin_Width = 1.0;
  std::vector<float> eta_Bin_Edges;
  std::vector<float> eta_Bin_Centers;

  //Fitting Parameters
  float stdv_range = 1.5;
  float Fit_E_Min;
  float Fit_E_Max;

  float th1_binning[3] = {16,-2,2}; //n,min,max. use diff$ratio
  int n_slice_bins = 10; //# bins in recoE TH1, within slices of truthE
  
  /* Jet Cuts */
  int min_N = 3;// Jet constituents
  
  //JES Histos (flattened "2D" vectors for E and eta)
  std::vector<TH1F*> E_Differences;
  std::vector<TH1F*> E_Ratios;
  std::vector<TH1F*> E_Slices;
  std::vector<TH1F*> Phi_Deltas;
  std::vector<TH1F*> Eta_Deltas;
  
  std::vector<std::vector<TH1F*>*> all_TH1F_vecs = 
  /* all_TH1F_vecs.push_back(E_Differences); */
  /* all_TH1F_vecs.push_back(E_Ratios); */
  /* all_TH1F_vecs.push_back(E_Slices); */
  /* all_TH1F_vecs.push_back(Phi_Deltas); */
  /* all_TH1F_vecs.push_back(Eta_Deltas); */
  {
    &E_Differences,
    &E_Ratios,
    &E_Slices,
    &Phi_Deltas,
    &Eta_Deltas
  };
  //store the means,sigma of the histograms for plotting
  std::vector<std::pair<float,float> > mean_E_Differences;
  std::vector<std::pair<float,float> > mean_E_Ratios;
  std::vector<std::pair<float,float> > mean_E_Slices;
  std::vector<std::pair<float,float> > mean_Phi_Deltas;
  std::vector<std::pair<float,float> > mean_Eta_Deltas;

  std::vector<std::vector<std::pair<float,float> > > all_mean_vecs =
    {
      mean_E_Differences,
      mean_E_Ratios,
      mean_E_Slices,
      mean_Phi_Deltas,
      mean_Eta_Deltas,
    };


  //TH2s: TH1s vs Jet E
  std::vector<TH2F*> Ratio_TH2F_v;
  std::vector<TH2F*> Diff_TH2F_v;
  std::vector<TH2F*> Slice_TH2F_v;

  std::vector<TH2F*> Ratio_vs_Truth;
  std::vector<TH2F*> Diff_vs_Truth;
  std::vector<TH2F*> Reco_vs_Truth;
  
  std::map<TString ,std::vector<TH2F*> > TH2F_map =
    {
      {"Ratio", Ratio_vs_Truth},
      {"Diff", Diff_vs_Truth},
      {"Reco", Reco_vs_Truth}
    };
/* FIXME: don't need TH2f vectors. Just TH2fs in eta and E */
  
  TString reco_or_corr;

public:
      
  void set_binning(std::vector<float> &edges, std::vector<float> &centers,
		   float low, float high, float width);

  void set_E_binning();
  void set_eta_binning();
  void set_E_FitRange(float low, float high);
  TString set_reco_or_corr(TString input_string);

  std::vector<TH1F*> create_TH1F(
     TString root_name, TString title);

  std::vector<TH2F*> create_TH2F(
     TString root_name, TString title);

  void initialize_histograms();

  void reconstructed_loop(TTreeReader *Tree);

  void fill_histograms(TLorentzVector truth, TLorentzVector reco, int etabin, int Ebin);

  void fit_histograms();
  
  void write_histograms(TFile *out_file);

  void resolution_vs_eta();

  void resolution_vs_E();
  //void reco_loop(TTreeReader Tree);

};
#endif
