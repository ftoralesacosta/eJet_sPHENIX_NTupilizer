#include "class_JES.h"
#include <iostream>

void MyClass::get_binning(std::vector <float> &edges, std::vector <float> &centers,
			  float low, float high, float width)
{
  if (fmod((high-low),width) != 0)
    fprintf(stderr,"WARNING: Choice of binning leaves remainder \n");

  int n_bins = int((high-low)/width);
  for (int i; i < n_bins+1; ++i)
    edges.push_back(low+i*width);
  for (int i; i < n_bins; ++i)
    centers.push_back(low + width*(0.5+i));

  //myobj.get_binning(myobj.pT_bins,myobj.pT_Centers,0,10,2);
}

void MyClass::get_pT_binning(float low, float high, float width)
{
  get_binning(pT_bins,pT_Centers,low,high,width);
}

void MyClass::get_eta_binning(float low, float high, float width)
{
  get_binning(eta_bins,eta_Centers,low,high,width);
}

void MyClass::set_pT_FitRange(float low, float high)
{
  Fit_pT_Min = low;
  Fit_pT_Max = high;
}

void MyClass::make_TH1F_Calib(TString reco_or_corr = "Reconstructed")
{ 
  std::vector<TH1F*> pT_Differences;
  std::vector<TH1F*> pT_Ratios;
  std::vector<TH1F*> pT_Slices;
  std::vector<TH1F*> Phi_Deltas;
  std::vector<TH1F*> Eta_Deltas;
  
  TString reco_or_corr_pT = TString("p_{T}^{") + reco_or_corr + TString("}");
  
  for (auto it = pT_bins.begin(); std::next(it,1) != pT_bins.end(); ++it){
    TString root_name = Form("pT%1.0f_%1.0f",*it*10,*std::next(it,1)*10);
    TString title_range = TString(Form("(%1.2f < p_{T}^{Truth} < %1.2f)",*it,*std::next(it,1)));

    pT_Differences.push_back(new TH1F(reco_or_corr + (TString("_pT_Difference_") + root_name).Data(),
				      (reco_or_corr_pT + TString(" - p_{T}^{Truth} ")
				       + title_range).Data(),16,-2,2));

    pT_Ratios.push_back(new TH1F((reco_or_corr + TString("_pT_Ratio_") + root_name).Data(),
				 (reco_or_corr_pT + TString(" / p_{T}^{Truth} ")
				  + title_range).Data(),16,-2,2));
    
    pT_Slices.push_back(new TH1F((reco_or_corr + TString("_pT_Slices") + root_name).Data(),
				 (reco_or_corr_pT + title_range).Data(),16,-2,2));

    Phi_Deltas.push_back(new TH1F((reco_or_corr + TString("_delta_phi_") + root_name).Data(),
				  (TString("#varphi^{") + reco_or_corr+ TString("} - #varphi^{Truth}")
				   + title_range).Data(),16,-2,2));

    Eta_Deltas.push_back(new TH1F((reco_or_corr + TString("_delta_eta_") + root_name).Data(),
				  (TString("#eta^{") + reco_or_corr+ TString("} - #eta^{Truth}")
				   + title_range).Data(),16,-2,2));

  }
  vectors_TH1F_Calib.push_back(pT_Differences);
  vectors_TH1F_Calib.push_back(pT_Ratios);
  vectors_TH1F_Calib.push_back(pT_Slices);
  vectors_TH1F_Calib.push_back(Phi_Deltas);
  vectors_TH1F_Calib.push_back(Eta_Deltas);
}

void MyClass::make_TH2F(TString reco_or_corr = "Reconstructed")
{

  TString reco_or_corr_pT = TString("p_{T}^{") + reco_or_corr + TString("}");
  TH2F* H_Diff_pT = new TH2F((TString("2D") + reco_or_corr + TString(pT_Residuals)).Data(),)

    TH2F* H_Diff_Tpt = new TH2F(Form("%s_H_2D_%s_pT_Residuals",r_c_string.Data(),r_cut_string.Data()),
                                Form("p_{T}^{%s} - p_{T}^{Truth} vs. p_{T}^{Truth}, %s",r_c_string.Data(),cut_string.Data()),38,0.5,10,50,-10,10);
    TH2F* H_Ratio_Tpt = new TH2F(Form("%s_H_2D_%s_pT_Ratios",r_c_string.Data(),r_cut_string.Data()),
                                 Form("p_{T}^{%s} / p_{T}^{Truth} vs. p_{T}^{Truth}, %s",r_c_string.Data(),cut_string.Data()),38,0.5,10,60,0,2);
    TH2F* H_Reco_Tpt = new TH2F(Form("H_2D_%s_%s_Tpt",r_c_string.Data(),r_cut_string.Data()),
                                Form("p_{T}^{%s} vs. p_{T}^{Truth}, %s",r_c_string.Data(),cut_string.Data()),38,0.5,10,38,0.5,10);

    H_Diff_Tpt->GetXaxis()->SetTitle("p_{T}^{Truth}");
    H_Diff_Tpt->GetYaxis()->SetTitle(Form("p_{T}^{%s} - p_{T}^{Truth}",r_c_string.Data()));

    H_Ratio_Tpt->GetXaxis()->SetTitle("p_{T}^{Truth}");
    H_Ratio_Tpt->GetYaxis()->SetTitle(Form("#frac{p_{T}^{%s}}{p_{T}^{Truth}}",r_c_string.Data()));

    H_Reco_Tpt->GetXaxis()->SetTitle("p_{T}^{Truth}");
    H_Reco_Tpt->GetYaxis()->SetTitle(Form("p_{T}^{%s}",r_c_string.Data()));

    vector_TH2F_Calib.push_back(H_Diff_Tpt);
    vector_TH2F_Calib.push_back(H_Ratio_Tpt);
    vector_TH2F_Calib.push_back(H_Reco_Tpt);
}



int main()
{
  MyClass myobj;
  myobj.get_pT_binning(0,10,2);
  myobj.get_eta_binning(0,10,1);
  myobj.set_pT_FitRange(2,20);
  myobj.make_TH1F_Calib();
  return 0;
}
