#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"
#include  <TLegend.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLatex.h>
#include <TColor.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
//#include "/gpfs/mnt/gpfs02/sphenix/user/fernando/jet_analysis/macros-rec/macros/sPHENIXStyle/sPhenixStyle.h"
//#include "/gpfs/mnt/gpfs02/sphenix/user/fernando/jet_analysis/macros-rec/macros/sPHENIXStyle/sPhenixStyle.C"
#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

const int MAX_INPUT_LENGTH = 200;


std::vector<std::vector<TH1F*> > make_TH1F_Calib(int N_pT_Bins, float *pT_Edges, TString cut_string, TString r_cut_string, bool is_corrected = false){
  //N-Dimensionality of calibration is mostly contained here. The naming scheme is currently pT and Eta, and can be expanded. This function should be called in a loop of eta bins and can be added to a nested loop for more dependancies.


  std::vector<std::vector<TH1F*> > v_TH1F_Calib;

  std::vector<TH1F*> pT_Differences;
  std::vector<TH1F*> pT_Ratios;
  std::vector<TH1F*> Reco_pT_Slices; //reco pT in "slices" of truth pT
  std::vector<TH1F*> Phi_Deltas;
  std::vector<TH1F*> Eta_Deltas;

  TString r_c_string = "Reco";

  if (is_corrected)
    r_c_string = "Corrected";

    for (int ipt = 0; ipt < N_pT_Bins; ipt++){    

      float pT_Truth = pT_Edges[ipt];
      float pT_Truth_plus = pT_Edges[ipt+1];


      pT_Differences.push_back(new TH1F(Form("%s_pT_Difference_pT%1.0f_%1.0f_%s",r_c_string.Data(),
					     pT_Truth*10,pT_Truth_plus*10,r_cut_string.Data()),
					Form("p_{T}^{%s} - p_{T}^{Truth} (%1.2f < p_{T}^{Truth} < %1.2f, %s )",
					     r_c_string.Data(),pT_Truth,pT_Truth_plus,cut_string.Data()), 150, -75,75));

      pT_Ratios.push_back(new TH1F(Form("%s_pT_Ratio_pT%1.0f_%1.0f_%s",r_c_string.Data(),
					pT_Truth*10,pT_Truth*10,pT_Truth_plus*10,r_cut_string.Data()),
				   Form("p_{T}^{%s} / p_{T}^{Truth} (%1.2f < p_{T}^{Truth} < %1.2f, %s )",
					r_c_string.Data(),pT_Truth,pT_Truth_plus,cut_string.Data()), 60,0,2));

      Reco_pT_Slices.push_back(new TH1F(Form("%s_pT_for_Truth_pT%1.0f_%1.0f_%s",r_c_string.Data(),
					     pT_Truth*10,pT_Truth*10,pT_Truth_plus*10,r_cut_string.Data()),
					Form("p_{T}^{%s}  (%1.2f < p_{T}^{Truth} < %1.2f, %s )",
					     r_c_string.Data(),pT_Truth,pT_Truth_plus,cut_string.Data()), 150,0,75));

      Phi_Deltas.push_back(new TH1F(Form("%s_dPhi_for_Truth_pT%1.0f_%1.0f_%s",r_c_string.Data(),
					 pT_Truth*10,pT_Truth*10,pT_Truth_plus*10,r_cut_string.Data()),
				    Form("#Delta#varphi^{%s}  (%1.2f < p_{T}^{Truth} < %1.2f, %s )",
					 r_c_string.Data(),pT_Truth,pT_Truth_plus,cut_string.Data()), 60,-0.3,0.3));

      Eta_Deltas.push_back(new TH1F(Form("%s_dEta_for_Truth_pT%1.0f_%1.0f_%s",r_c_string.Data(),
					 pT_Truth*10,pT_Truth*10,pT_Truth_plus*10,r_cut_string.Data()),
				    Form("#Delta|#eta|^{%s}  (%1.2f < p_{T}^{Truth} < %1.2f, %s )",
					 r_c_string.Data(),pT_Truth,pT_Truth_plus,cut_string.Data()), 60,-0.3,0.3));
    }

    v_TH1F_Calib.push_back(pT_Differences);
    v_TH1F_Calib.push_back(pT_Ratios);
    v_TH1F_Calib.push_back(Reco_pT_Slices);
    v_TH1F_Calib.push_back(Phi_Deltas);
    v_TH1F_Calib.push_back(Eta_Deltas);

    return v_TH1F_Calib;

  }

std::vector<TH2F*>  make_TH2F(TString cut_string, TString r_cut_string,bool is_corrected = false){
    
    TString r_c_string = "Reco";

    if (is_corrected)
      r_c_string = "Corrected";

    std::vector<TH2F*> v_TH2F_Calib;
    
    TH2F* H_Diff_Tpt = new TH2F(Form("%s_H_2D_%s_pT_Residuals",r_c_string.Data(),r_cut_string.Data()),
				Form("p_{T}^{%s} - p_{T}^{Truth} vs. p_{T}^{Truth}, %s",r_c_string.Data(),cut_string.Data()),130,10,75,150,-75,75);
    TH2F* H_Ratio_Tpt = new TH2F(Form("%s_H_2D_%s_pT_Ratios",r_c_string.Data(),r_cut_string.Data()),
				 Form("p_{T}^{%s} / p_{T}^{Truth} vs. p_{T}^{Truth}, %s",r_c_string.Data(),cut_string.Data()),130,10,75,60,0,2);
    TH2F* H_Reco_Tpt = new TH2F(Form("H_2D_%s_%s_Tpt",r_c_string.Data(),r_cut_string.Data()),
				Form("p_{T}^{%s} vs. p_{T}^{Truth}, %s",r_c_string.Data(),cut_string.Data()),130,10,75,130,10,75);

    H_Diff_Tpt->GetXaxis()->SetTitle("p_{T}^{Truth}");
    H_Diff_Tpt->GetYaxis()->SetTitle(Form("p_{T}^{%s} - p_{T}^{Truth}",r_c_string.Data()));

    H_Ratio_Tpt->GetXaxis()->SetTitle("p_{T}^{Truth}");
    H_Ratio_Tpt->GetYaxis()->SetTitle(Form("#frac{p_{T}^{%s}}{p_{T}^{Truth}}",r_c_string.Data()));

    H_Reco_Tpt->GetXaxis()->SetTitle("p_{T}^{Truth}");
    H_Reco_Tpt->GetYaxis()->SetTitle(Form("p_{T}^{%s}",r_c_string.Data()));

    v_TH2F_Calib.push_back(H_Diff_Tpt);
    v_TH2F_Calib.push_back(H_Ratio_Tpt);
    v_TH2F_Calib.push_back(H_Reco_Tpt);
    
    return v_TH2F_Calib;

  }



std::vector< std::vector<float> >  Fit_Gauss(int N_pT_Bins, std::vector<TH1F*> Histo) {

  //returns std::vectors in following order:
    std::vector<float> Gauss_Mean;
    std::vector<float> Numerical_Mean;
    std::vector<float> Gauss_Sigma;
    std::vector<float> Numerical_Sigma;
    std::vector<float> Gauss_Mean_Error;
    std::vector<float> Gauss_Sigma_Error;

    for (int ipt = 0; ipt < N_pT_Bins; ipt++){
      
      //obtain fit range for gauss
      float c = 1.5;
      float mean = Histo[ipt]->GetMean();
      float stdev = Histo[ipt]->GetStdDev();
      Numerical_Mean.push_back(mean);
      Numerical_Sigma.push_back(stdev);
 
      float num_min = mean-(c*stdev);
      float num_max = mean+(c*stdev);
     
      gStyle->SetStatY(0.85);
      TF1* num_fit = new TF1("nfit","gaus",num_min,num_max);

      Histo[ipt]->Fit(num_fit,"QR");

      float gaus_mean = num_fit->GetParameters()[1];
      float gaus_sigma = num_fit->GetParameters()[2];
      float mean_error = num_fit->GetParError(1);
      float sigma_error = num_fit->GetParError(2);
 
      Gauss_Mean.push_back(gaus_mean);
      Gauss_Sigma.push_back(gaus_sigma);
      Gauss_Mean_Error.push_back(mean_error);      
      Gauss_Sigma_Error.push_back(sigma_error);      
    }

    std::vector< std::vector<float> > v_Results;
    v_Results.push_back(Gauss_Mean);
    v_Results.push_back(Numerical_Mean);
    v_Results.push_back(Gauss_Sigma);
    v_Results.push_back(Numerical_Sigma);
    v_Results.push_back(Gauss_Mean_Error);
    v_Results.push_back(Gauss_Sigma_Error);

    return v_Results;
  }

std::vector<TGraphErrors*> Make_TGraphs(int N_pT_Bins, float *pT_Centers, TString cut_string,TString r_cut_string,std::vector< std::vector<float> > v_Results, TString root_name, TString title)
{

  //not necesary, but redifined for clarity
    std::vector<float> Gaus_Mean = v_Results[0];
    std::vector<float> Numerical_Mean = v_Results[1];
    std::vector<float> Gaus_Sigma = v_Results[2];
    std::vector<float> Numerical_Sigma = v_Results[3];
    std::vector<float> Gaus_Mean_Error = v_Results[4];
    std::vector<float> Gaus_Sigma_Error = v_Results[5];

    float pT_Width = pT_Centers[1] - pT_Centers[0];

    std::vector <TGraphErrors*> v_Graphs;
    v_Graphs.push_back(new TGraphErrors(N_pT_Bins,pT_Centers,&Gaus_Mean[0],0, &Gaus_Mean_Error[0]));
    v_Graphs[0]->SetNameTitle(Form("Gaus_%s_Mean_%s",root_name.Data(),r_cut_string.Data()),Form("Mean %s %s",title.Data(),cut_string.Data()));
    v_Graphs[0]->GetXaxis()->SetTitle("p_{T}^{Truth}");
    v_Graphs[0]->SetMarkerStyle(8);
    v_Graphs[0]->SetMarkerColor(2);

    v_Graphs.push_back(new TGraphErrors(N_pT_Bins,pT_Centers,&Numerical_Mean[0],0,&Numerical_Sigma[0]));
    v_Graphs[1]->SetNameTitle(Form("Numerical_%s_Mean_%s",root_name.Data(),r_cut_string.Data()),Form("Mean %s",title.Data()));
    v_Graphs[1]->GetXaxis()->SetTitle("p_{T}^{Truth}");
    v_Graphs[1]->SetMarkerStyle(8);

    v_Graphs.push_back(new TGraphErrors(N_pT_Bins,pT_Centers,&Gaus_Sigma[0],0,&Gaus_Sigma_Error[0]));
    v_Graphs[2]->SetNameTitle(Form("Gaus_%s_Sigma_%s",root_name.Data(),r_cut_string.Data()),Form("#sigma(%s)",title.Data()));
    v_Graphs[2]->GetXaxis()->SetTitle("p_{T}^{Truth}");
    v_Graphs[2]->SetMarkerStyle(8);
    //    v_Graphs[2]->GetYaxis()->SetRangeUser(0,0.4);

    v_Graphs.push_back(new TGraphErrors(N_pT_Bins,pT_Centers,&Numerical_Sigma[0],0,0));
    v_Graphs[3]->SetNameTitle(Form("Numerical_%s_Sigma_%s",root_name.Data(),r_cut_string.Data()),Form("#sigma(%s)",title.Data()));
    v_Graphs[3]->GetXaxis()->SetTitle("p_{T}^{Truth}");
    v_Graphs[3]->SetMarkerStyle(8);
    //v_Graphs[3]->GetYaxis()->SetRangeUser(0,0.4);

    for (auto g_it = v_Graphs.begin(); g_it != v_Graphs.end(); ++g_it)
      (*g_it)->Write();

    return v_Graphs;
}


void Plot_Gaus(int N_pT_Bins, TString r_cut_string, std::vector<TH1F*> Histo_Array,TString save_descr)
{
  TCanvas canv("Histo_Array", "Histo_Array", 1600, 1200);                                                                                                  
  int n_div = int(TMath::Sqrt(N_pT_Bins)+0.999); //defaults to extra division
  canv.Divide(n_div,n_div);

  for (int ipt = 0; ipt < N_pT_Bins;ipt++){
    canv.cd(ipt+1);
    Histo_Array[ipt]->Draw();
    Histo_Array[ipt]->SetTitleOffset(1.6);
  }

  canv.Write();
  canv.SaveAs(Form("%s_%s.pdf",save_descr.Data(),r_cut_string.Data()));
}

void Num_Gaus_Overlay(TH2F *H2, std::vector<TGraphErrors*> v_Graphs, TString description, TString cut_string, TString r_cut_string)
  {

    //Vector: 
    // 0) Gaus Mean, 
    // 1) Numerical Mean, 
    // 2) Gaus Sigma, 
    // 3) Numerical Sigma (sqrt(RMS))

    TCanvas mean_canv(Form("mean_%s_%s",description.Data(),r_cut_string.Data()),
		      Form("mean_%s_%s",description.Data(),r_cut_string.Data()));
    H2->Draw("COLZ");
    v_Graphs[0]->SetLineColor(2); //Gauss
    v_Graphs[0]->SetLineWidth(3);
    v_Graphs[0]->Draw("SAME"); 
    v_Graphs[1]->SetLineColor(1); //Numerical
    v_Graphs[1]->SetMarkerStyle(8);
    v_Graphs[1]->Draw("SAME");

    mean_canv.Write();
    mean_canv.SaveAs(Form("%s_Means_%s.pdf",description.Data(),r_cut_string.Data()));

    TCanvas sigma_canv(Form("sigma_%s_%s",description.Data(),r_cut_string.Data()),
		       Form("sigma_%s_%s",description.Data(),r_cut_string.Data()));    
    v_Graphs[2]->SetLineColor(2);
    v_Graphs[2]->Draw();
    v_Graphs[3]->SetLineColor(4);
    v_Graphs[3]->Draw("SAME");
    auto legend = new TLegend(0.3,0.7,0.58,0.8);
    legend->AddEntry(v_Graphs[2], "Gaus Fit","lp");
    legend->AddEntry(v_Graphs[3], "Numerical","lp");
    legend->SetFillStyle(0);
    legend->Draw("SAME");

    sigma_canv.Write();
    sigma_canv.SaveAs(Form("%s_Sigma_%s.pdf",description.Data()));
  }


void plot_resolution(int N_pT_Bins, float *pT_Centers, TString cut_string, TString r_cut_string,std::vector<std::vector<float > > v_Results, TString description,bool corrected=false) {
  //resolution defined as sigma<R>/mu<R>
  //Alternative: sigma<dpT>/<pT_Truth>

  TString Corr_String = "Reco";

  if (corrected)
    Corr_String = "Corrected";

  std::vector<float> Gaus_Mean = v_Results[0];
  std::vector<float> Gaus_Sigma = v_Results[2];
  std::vector<float> Gaus_Mean_Error = v_Results[4];
  std::vector<float> Gaus_Sigma_Error = v_Results[5];
 
  std::vector<float> resolution;
  std::vector<float> resolution_error;

  std::vector<float> pT_Error;

  for (int w = 0; w < N_pT_Bins; w++)
    pT_Error.push_back((pT_Centers[1] - pT_Centers[0])/2); //FIXME: Genaralize for non-constant binning

  for (int i = 0; i < Gaus_Sigma.size(); i++){

    resolution.push_back(Gaus_Sigma[i]/TMath::Abs(Gaus_Mean[i]) );

    resolution_error.push_back( sqrt(pow( (Gaus_Sigma_Error[i]/Gaus_Sigma[i]),2 ) +
    				     pow( (Gaus_Mean_Error[i]/Gaus_Mean[i]),2 )) * resolution[i] );


    // resolution.push_back(Gaus_Sigma[i]/pT_Centers[i] );
    // resolution_error.push_back( sqrt(pow( (Gaus_Sigma_Error[i]/Gaus_Sigma[i]),2 ) +
    //                                  pow( (pT_Error/pT_Centers[i]),2 )) * resolution[i] );

    // fprintf(stderr, "\n %s \n",description.Data());
    // fprintf(stderr,"%d: Sigma = %f \n",__LINE__,Gaus_Sigma[i]);
    // fprintf(stderr,"%d: Mean = %f \n",__LINE__,Gaus_Mean[i]);
    // fprintf(stderr,"%d: Resolution = %f \n",__LINE__,resolution[i]);
  }

  TGraphErrors *Resolution_Graph =  new TGraphErrors(Gaus_Sigma.size(),pT_Centers,&resolution[0],&pT_Error[0],&resolution_error[0]);
  Resolution_Graph->SetNameTitle(Form("%s_Jet_%s_Resolution_%s",Corr_String.Data(),description.Data(),r_cut_string.Data()),
				 Form("%s Jet %s Resolution %s",Corr_String.Data(),description.Data(),cut_string.Data()));
  Resolution_Graph->GetXaxis()->SetTitle("p_{T}^{Truth}");
  Resolution_Graph->GetYaxis()->SetTitle("#sigma/#mu");

  TCanvas resolution_canv(Form("%s_%s_canv",Corr_String.Data(),description.Data()),
			  Form("%s_%s_canv",Corr_String.Data(),description.Data()),1600,1000);
  Resolution_Graph->SetLineWidth(3);
  Resolution_Graph->SetMarkerStyle(10);
  Resolution_Graph->Draw();
  resolution_canv.Write();
  resolution_canv.SaveAs(Form("%s_Jet_%s_Resolution_%s.pdf",Corr_String.Data(),description.Data(),r_cut_string.Data()));

}

TF1* do_gaus_fitting(int N_pT_Bins, float* pT_Centers, TString cut_string,TString r_cut_string,std::vector<std::vector<TH1F*> >vv_TH1F, std::vector<TH2F*> v_TH2F, bool is_corrected = false){

  TString r_c_string = "Reco";

  if (is_corrected)
    r_c_string = "Corrected";

  std::vector<TH1F*> pT_Differences = vv_TH1F[0];
  std::vector<TH1F*> pT_Ratios = vv_TH1F[1];
  std::vector<TH1F*> Reco_pT_Slices = vv_TH1F[2];
  std::vector<TH1F*> Phi_Deltas = vv_TH1F[3];
  std::vector<TH1F*> Eta_Deltas = vv_TH1F[4];


  TH2F* H_Diff_Tpt = v_TH2F[0];
  TH2F* H_Ratio_Tpt =  v_TH2F[1];
  TH2F* H_Reco_Tpt =  v_TH2F[2];

  std::vector< std::vector<float> > v_Diffs = Fit_Gauss(N_pT_Bins,pT_Differences);
  std::vector< std::vector<float> > v_Ratios = Fit_Gauss(N_pT_Bins,pT_Ratios);
  std::vector< std::vector<float> > v_Reco_Slices = Fit_Gauss(N_pT_Bins,Reco_pT_Slices);
  std::vector< std::vector<float> > v_dPhi = Fit_Gauss(N_pT_Bins,Phi_Deltas);
  std::vector< std::vector<float> > v_dEta = Fit_Gauss(N_pT_Bins,Eta_Deltas);
  
  std::vector<TGraphErrors*> v_Diff_Graph =  Make_TGraphs(N_pT_Bins, pT_Centers, cut_string, r_cut_string, v_Diffs, TString("Diff_Graph"), TString(Form("p_{T}^{%s} - p_{T}^{Truth}",r_c_string.Data() )) );
  std::vector<TGraphErrors*> v_Ratio_Graph = Make_TGraphs(N_pT_Bins, pT_Centers, cut_string, r_cut_string,v_Ratios, TString("Ratio_Graph"), TString(Form("p_{T}^{%s}/p_{T}^{Truth}",r_c_string.Data() )) );
  std::vector<TGraphErrors*> v_Reco_Slices_Graph = Make_TGraphs(N_pT_Bins, pT_Centers, cut_string, r_cut_string, v_Reco_Slices, TString(Form("pT_%s_Slice_Graph",r_c_string.Data())), TString(Form("p_{T}^{%s}",r_c_string.Data() )) );
  std::vector<TGraphErrors*> v_dPhi_Graph =  Make_TGraphs(N_pT_Bins, pT_Centers, cut_string, r_cut_string,v_dPhi, TString("dPhi_Graph"), TString(Form("#Delta#varphi %s Jets",r_c_string.Data() )) );
  std::vector<TGraphErrors*> v_dEta_Graph =  Make_TGraphs(N_pT_Bins, pT_Centers, cut_string, r_cut_string, v_dEta ,TString("dEta_Graph"), TString(Form("#Delta#eta %s Jets",r_c_string.Data() )) );


  Plot_Gaus(N_pT_Bins,r_cut_string,pT_Differences,TString(Form("%s_pT_Differences",r_c_string.Data()) ) );
  Plot_Gaus(N_pT_Bins,r_cut_string,pT_Ratios,TString(Form("%s_pT_Ratios",r_c_string.Data()) ) );
  Plot_Gaus(N_pT_Bins,r_cut_string,Reco_pT_Slices,TString(Form("%s_pT_Slices",r_c_string.Data()) ) );
  Plot_Gaus(N_pT_Bins,r_cut_string,Phi_Deltas,TString(Form("%s_dPhi",r_c_string.Data()) ) );
  Plot_Gaus(N_pT_Bins,r_cut_string,Eta_Deltas,TString(Form("%s_dEta",r_c_string.Data()) ) );

  Num_Gaus_Overlay(H_Diff_Tpt, v_Diff_Graph, TString(Form("%s_TwoD_pT_Differences",r_c_string.Data())),cut_string,r_cut_string);
  Num_Gaus_Overlay(H_Ratio_Tpt, v_Ratio_Graph, TString(Form("%s_TwoD_pT_Ratios",r_c_string.Data())),cut_string,r_cut_string);
  Num_Gaus_Overlay(H_Reco_Tpt, v_Reco_Slices_Graph, TString(Form("%s_TwoD_pT_RecovTruth",r_c_string.Data())),cut_string,r_cut_string);

  plot_resolution(N_pT_Bins,pT_Centers,cut_string,r_cut_string,v_Ratios,TString("Energy"),is_corrected);
  plot_resolution(N_pT_Bins,pT_Centers,cut_string,r_cut_string,v_dPhi,TString("dPhi"),is_corrected);
  plot_resolution(N_pT_Bins,pT_Centers,cut_string,r_cut_string,v_dEta,TString("dEta"),is_corrected);

  TF1* Fi;

  if (not(is_corrected)){
    TCanvas canv("Histo_Array", "Histo_Array", 1600, 1000);
    TGraphErrors *ratio_v_reco = new TGraphErrors(N_pT_Bins,&v_Reco_Slices[0][0],&v_Ratios[0][0],&v_Reco_Slices[4][0],&v_Ratios[4][0]);
    ratio_v_reco->GetXaxis()->SetTitle(" p_{T}^{Reco} [GeV/c]");
    ratio_v_reco->GetXaxis()->CenterTitle();
    ratio_v_reco->GetYaxis()->SetTitle(" p_{T}^{Reco}/p_{T}^{Truth} ");
    ratio_v_reco->SetTitle(Form("R vs p_{T}^{Reco} %s",cut_string.Data()));
    ratio_v_reco->GetYaxis()->CenterTitle();

    Fi = new TF1("Fi","[0]+[1]*(TMath::Log(x))+[2]*(TMath::Log(x))^2",15,55); //Fit Range
    gStyle->SetStatY(0.5);

    gStyle->SetOptFit(10111);
    ratio_v_reco->Fit(Fi,"R");

    //Save Graph
    ratio_v_reco->Draw();

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.DrawLatex(-12.0,0.9," #it{f} = a+b#upoint ln(p_{T}^{reco}) + c#upoint ln(p_{T}^{reco})^{2}");

    TString save_descr = TString("ratio_v_reco");
    canv.SaveAs(Form("%s_%s.pdf",save_descr.Data(),r_cut_string.Data()));

    TCanvas canv2(Form("Correction_%s",r_cut_string.Data()), "Correction", 1600,1000);
    std::vector<float> v_Cor_Factors;
    std::vector<float> v_Cor_pT;
    //for (auto ipt = pT_Centers.begin(); ipt != pT_Centers.end(); ++ipt){
    for (int ipt = 0; ipt < N_pT_Bins; ipt++){
      float F = Fi->Eval(pT_Centers[ipt]);
      v_Cor_Factors.push_back(1.0/F);
      v_Cor_pT.push_back(v_Reco_Slices[0][ipt]);
    }
    TGraph *Cor = new TGraph(N_pT_Bins,&v_Cor_pT[0],&v_Cor_Factors[0]);
    Cor->GetXaxis()->SetTitle("p_{T}^{Corrected} [GeV/c]");
    Cor->GetXaxis()->CenterTitle();
    Cor->GetXaxis()->SetRangeUser(10,60);
    Cor->GetYaxis()->SetTitle(" JES Correction ");
    Cor->SetTitle(Form("Jet p_{T} Correction  %s",cut_string.Data()));
    Cor->GetYaxis()->CenterTitle();

    Cor->Draw();
    save_descr = TString("Correction_Plot");
    canv2.SaveAs(Form("%s_%s.pdf",save_descr.Data(),r_cut_string.Data()));
    canv2.Write();

  }//Determine Fi

  return Fi;
}

std::vector<float> Get_pT_Widths(std::vector<float> pT_Bins){

  std::vector<float> pT_Widths;

  for (int ipt = 0; ipt < pT_Bins.size()-1; ipt++) {
    float iwidth = (pT_Bins[ipt+1] - pT_Bins[ipt])/2;
    pT_Widths.push_back(iwidth);
  }

  return pT_Widths;
}

void Closure_Test(int N_pT_Bins,float* pT_Centers,TString cut_string,TString r_cut_string,std::vector< std::vector<float> > v, TString Descr, TString Y_Title){
 
  float width = (pT_Centers[1] - pT_Centers[0])/2;

  float pT_Widths[N_pT_Bins];

  for (int iw = 0; iw < N_pT_Bins-1; iw++)
    pT_Widths[iw] = (pT_Centers[iw+1] - pT_Centers[iw])/2;
  pT_Widths[N_pT_Bins-1] = pT_Widths[N_pT_Bins-2];

  TCanvas check_canv(Form("%s_Check_%s",Descr.Data(),r_cut_string.Data()), Form("%s_Check",Descr.Data()), 1600, 1000);
    TGraphErrors *Correction_Check = new TGraphErrors(N_pT_Bins,pT_Centers,&v[0][0],pT_Widths,&v[4][0]);
    Correction_Check->SetTitle(Form("%s Closure Test   %s",Descr.Data(),cut_string.Data()));
    Correction_Check->GetXaxis()->SetTitle("p_{T}^{Truth}");
    Correction_Check->GetYaxis()->CenterTitle();
    Correction_Check->GetXaxis()->CenterTitle();
    Correction_Check->Fit("pol0","","",20,70);
    Correction_Check->GetXaxis()->SetRangeUser(10,70);
    Correction_Check->GetYaxis()->SetRangeUser(-1,2);
    Correction_Check->Draw();

    check_canv.Write();
    check_canv.SaveAs(Form("%s_Closure_Test_%s.pdf",Descr.Data(),r_cut_string.Data()));
    return;
  }

void resolution_overlay(int N_pT_Bins, float *pT_Centers, TString cut_string, TString r_cut_string, std::vector<TH1F*> v_H_1, TString name1, std::vector<TH1F*> v_H_2, TString name2, TString Title){
  std::vector< std::vector<float> > v_dPhis = Fit_Gauss(N_pT_Bins,v_H_1);
  std::vector< std::vector<float> > v_dEtas = Fit_Gauss(N_pT_Bins,v_H_2);

  std::vector<TGraphErrors*>  v_TG_dPhis = Make_TGraphs(N_pT_Bins, pT_Centers, cut_string, r_cut_string,v_dPhis, name1, Title);
  std::vector<TGraphErrors*>  v_TG_dEtas = Make_TGraphs(N_pT_Bins, pT_Centers,cut_string, r_cut_string, v_dEtas, name2, Title);
  TCanvas C_dPhidEta(Form("%s_%s_Overlay",name1.Data(),name2.Data()),Form("%s_%s_Overlay",name1.Data(),name2.Data()),1600,1000);

  v_TG_dPhis[2]->SetLineColor(2);
  v_TG_dPhis[2]->SetMarkerColor(2);
  v_TG_dPhis[2]->SetTitle(Title.Data());
  v_TG_dPhis[2]->Draw("");

  v_TG_dEtas[2]->SetLineColor(4);
  v_TG_dEtas[2]->SetMarkerColor(4);
  v_TG_dEtas[2]->Draw("SAME");

  auto legend = new TLegend(0.65,0.6,0.8,0.8);
  legend->AddEntry(v_TG_dPhis[2],name1.Data(),"lp");
  legend->AddEntry(v_TG_dEtas[2],name2.Data(),"lp");
  legend->Draw("SAME");
  C_dPhidEta.Write();
  C_dPhidEta.SaveAs(Form("resolution_%s_%s_Overlay_%s.pdf",name1.Data(),name2.Data(),r_cut_string.Data()));
}


void JER_Parametrization(int N_pT_Bins, float *pT_Centers, float Fit_Min,TString cut_string, TString r_cut_string, std::vector< std::vector<float> > v_dRatios, TString name, TString Title){

  gStyle->SetOptFit(111);

  std::vector<TGraphErrors*>  v_TG_dPhis = Make_TGraphs(N_pT_Bins, pT_Centers, cut_string, r_cut_string,v_dRatios, name, Title);

  //resolution of ratios
  TGraphErrors* Sigma_Ratio = v_TG_dPhis[2];
  TCanvas C(Form("%s_Parametrization",name.Data()),Form("%s_Parametrization",name.Data()),1600,1000);
TF1 func = TF1("func","TMath::Sqrt(TMath::Power(0.743423,2)/x + TMath::Power([0]/x,2) + TMath::Power(0.0763649,2))",Fit_Min,72);  
//TF1 func = TF1("func","TMath::Sqrt(TMath::Power([0],2)/x + TMath::Power([1]/x,2) + TMath::Power([2],2))",Fit_Min,72);
  //TF1 func = TF1("func","TMath::Sqrt(TMath::Power([0],2)/x + TMath::Power([1]/x,2))",16,72);
  Sigma_Ratio->Draw("");
  Sigma_Ratio->Fit("func","R");

  //   func.SetParameter(0,0.55)
  //   func.SetParameter(1,8.0)
  //   func.SetParameter(2,12.0)

  //   Purity_Graph.Fit("func","S")

  //   Fit_Band = ROOT.TGraphErrors(len(Purity));
  // pGeV_Edges_3 = np.linspace(12,40,40)
  //   for i in range(len(pGeV_Edges_3)):
  // Fit_Band.SetPoint(i, pGeV_Edges_3[i], 0)
  // 	    (ROOT.TVirtualFitter.GetFitter()).GetConfidenceIntervals(Fit_Band,0.68)
  C.Write();
  C.SaveAs(Form("resolution_%s_Parametrization.pdf",name.Data()));
}
int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cout<<"Syntax: [Command] [File]"<<std::endl;
    exit(EXIT_FAILURE);
  }
  //for (int iarg = 1; iarg < argc; iarg++) {
  int iarg = 1;
  TString root_file = (TString)argv[iarg];
  
  std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
  
  TFile *file = TFile::Open(root_file);
  
  if (file == NULL) {
    std::cout << " File Fail" << std::endl;
    exit(EXIT_FAILURE); 
  } 
  
  file->Print();
  
  TTree *_tree_event = dynamic_cast<TTree *>(file->Get("ttree"));

  if (_tree_event == NULL) {
    std::cout << " Tree Fail " << std::endl;
    exit(EXIT_FAILURE);
  }

  int embed_id;

  float Psi2;

  float HijingPt_Psi2;
  float HijingN_Psi2;
  float EPDPt_Psi2;
  float EPDN_Psi2;
  float fEMCPt_Psi2;
  float fEMCN_Psi2;

  float Hijing_v2;
  float b;

  int particle_n;
  float particle_pt[1000];
  float particle_eta[1000];
  float particle_phi[1000];
  int particle_pid[1000];
  int particle_embed[1000];

  int jet4_n;
  Float_t jet4_pt[200];
  Float_t jet4_E[200];
  Float_t jet4_m[200];
  int jet4_ncomp[250];
  Float_t jet4_eta[200];
  Float_t jet4_phi[200];
  float jet4_E_layer0[200];
  float jet4_E_layer1[200];
  float jet4_E_layer2[200];

  int jet4mod_n;
  float jet4mod_pt[200];
  float jet4mod_m[200];
  float jet4mod_eta[200];
  float jet4mod_phi[200];
  float jet4mod_E_layer0[200];
  float jet4mod_E_layer1[200];
  float jet4mod_E_layer2[200];

  int truthjet4_n;
  Float_t truthjet4_pt[200];
  Float_t truthjet4_m[200];
  Float_t truthjet4_eta[200];
  Float_t truthjet4_phi[200];

  int jet2_n;
  float jet2_pt[250];
  float jet2_E[250];
  float jet2_m[250];
  int jet2_ncomp[250];
  float jet2_eta[250];
  float jet2_phi[250];
  float jet2_E_layer0[250];
  float jet2_E_layer1[250];
  float jet2_E_layer2[250];

  int jet3_n;
  float jet3_pt[200];
  float jet3_E[200];
  float jet3_m[200];
  float jet3_eta[200];
  float jet3_phi[200];

  int jet5_n;
  float jet5_pt[200];
  float jet5_E[200];
  float jet5_m[200];
  float jet5_eta[200];
  float jet5_phi[200];


  int jet2mod_n;
  float jet2mod_pt[200];
  float jet2mod_m[200];
  float jet2mod_eta[200];
  float jet2mod_phi[200];
  float jet2mod_E_layer0[200];
  float jet2mod_E_layer1[200];
  float jet2mod_E_layer2[200];

  int truthjet2_n;
  float truthjet2_pt[200];
  float truthjet2_m[200];
  float truthjet2_eta[200];
  float truthjet2_phi[200];

  int cluster_n;
  float cluster_E[100];
  float cluster_ecore[100];
  float cluster_eta[100];
  float cluster_phi[100];

  int clustersub_n;
  float clustersub_E[100];
  float clustersub_ecore[100];
  float clustersub_eta[100];
  float clustersub_phi[100];

  int calo_n;
  int calo_layer[5000];
  float calo_eta[5000];
  float calo_phi[5000];
  float calo_E[5000];


  //  _tree_event->SetBranchAddress("embed_id",&embed_id);
  //  _tree_event->SetBranchAddress("particle_embed", particle_embed);

  // _tree_event->SetBranchAddress("jet4mod_n", &jet4mod_n);
  // _tree_event->SetBranchAddress("jet4mod_pt", jet4mod_pt);
  // _tree_event->SetBranchAddress("jet4mod_m", jet4mod_m);
  // _tree_event->SetBranchAddress("jet4mod_eta", jet4mod_eta);
  // _tree_event->SetBranchAddress("jet4mod_phi", jet4mod_phi);
  // _tree_event->SetBranchAddress("jet4mod_E_layer0", jet4mod_E_layer0);
  // _tree_event->SetBranchAddress("jet4mod_E_layer1", jet4mod_E_layer1);
  // _tree_event->SetBranchAddress("jet4mod_E_layer2", jet4mod_E_layer2);


  // _tree_event->SetBranchAddress("cluster_n", &cluster_n);
  // _tree_event->SetBranchAddress("cluster_E", cluster_E);
  // _tree_event->SetBranchAddress("cluster_ecore", cluster_ecore);
  // _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
  // _tree_event->SetBranchAddress("cluster_phi", cluster_phi);

  // _tree_event->SetBranchAddress("clustersub_n", &clustersub_n);
  // _tree_event->SetBranchAddress("clustersub_E", clustersub_E);
  // _tree_event->SetBranchAddress("clustersub_ecore", clustersub_ecore);
  // _tree_event->SetBranchAddress("clustersub_eta", clustersub_eta);
  // _tree_event->SetBranchAddress("clustersub_phi", clustersub_phi);

  // _tree_event->SetBranchAddress("calo_n", &calo_n);
  // _tree_event->SetBranchAddress("calo_layer", calo_layer);
  // _tree_event->SetBranchAddress("calo_eta", calo_eta);
  // _tree_event->SetBranchAddress("calo_phi", calo_phi);
  // _tree_event->SetBranchAddress("calo_E", calo_E);

  _tree_event->SetBranchAddress("Psi2", &Psi2);

  _tree_event->SetBranchAddress("HijingPt_Psi2", &HijingPt_Psi2);
  _tree_event->SetBranchAddress("HijingN_Psi2", &HijingN_Psi2);
  _tree_event->SetBranchAddress("EPDPt_Psi2", &EPDPt_Psi2);
  _tree_event->SetBranchAddress("EPDN_Psi2", &EPDN_Psi2);
  _tree_event->SetBranchAddress("fEMCPt_Psi2",&fEMCPt_Psi2);
  _tree_event->SetBranchAddress("fEMCN_Psi2", &fEMCN_Psi2);

  _tree_event->SetBranchAddress("Hijing_v2", &Hijing_v2);
  _tree_event->SetBranchAddress("b",&b);

  _tree_event->SetBranchAddress("particle_n", &particle_n);
  _tree_event->SetBranchAddress("particle_pt", particle_pt);
  _tree_event->SetBranchAddress("particle_eta", particle_eta);
  _tree_event->SetBranchAddress("particle_phi", particle_phi);
  _tree_event->SetBranchAddress("particle_pid", particle_pid);


  _tree_event->SetBranchAddress("jet4_n", &jet4_n);
  _tree_event->SetBranchAddress("jet4_pt", jet4_pt);
  _tree_event->SetBranchAddress("jet4_E", jet4_E);
  _tree_event->SetBranchAddress("jet4_m", jet4_m);
  _tree_event->SetBranchAddress("jet4_ncomp", jet4_ncomp);
  _tree_event->SetBranchAddress("jet4_eta", jet4_eta);
  _tree_event->SetBranchAddress("jet4_phi", jet4_phi);
  _tree_event->SetBranchAddress("jet4_E_layer0", jet4_E_layer0);
  _tree_event->SetBranchAddress("jet4_E_layer1", jet4_E_layer1);
  _tree_event->SetBranchAddress("jet4_E_layer2", jet4_E_layer2);

  _tree_event->SetBranchAddress("truthjet4_n", &truthjet4_n);
  _tree_event->SetBranchAddress("truthjet4_pt", truthjet4_pt);
  _tree_event->SetBranchAddress("truthjet4_m", truthjet4_m);
  _tree_event->SetBranchAddress("truthjet4_eta", truthjet4_eta);
  _tree_event->SetBranchAddress("truthjet4_phi", truthjet4_phi);

  _tree_event->SetBranchAddress("jet2_n", &jet2_n);
  _tree_event->SetBranchAddress("jet2_pt", jet2_pt);
  _tree_event->SetBranchAddress("jet2_E", jet2_E);
  _tree_event->SetBranchAddress("jet2_m", jet2_m);
  _tree_event->SetBranchAddress("jet2_ncomp", jet2_ncomp);
  _tree_event->SetBranchAddress("jet2_eta", jet2_eta);
  _tree_event->SetBranchAddress("jet2_phi", jet2_phi);
  _tree_event->SetBranchAddress("jet2_E_layer0", jet2_E_layer0);
  _tree_event->SetBranchAddress("jet2_E_layer1", jet2_E_layer1);
  _tree_event->SetBranchAddress("jet2_E_layer2", jet2_E_layer2);

  _tree_event->SetBranchAddress("truthjet2_n", &truthjet2_n);
  _tree_event->SetBranchAddress("truthjet2_pt", truthjet2_pt);
  _tree_event->SetBranchAddress("truthjet2_m", truthjet2_m);
  _tree_event->SetBranchAddress("truthjet2_eta", truthjet2_eta);
  _tree_event->SetBranchAddress("truthjet2_phi", truthjet2_phi);


  // gROOT->LoadMacro("/gpfs/mnt/gpfs02/sphenix/user/fernando/jet_analysis/macros-rec/macros/sPHENIXStyle/sPhenixStyle.C");
  // SetsPhenixStyle();
  gStyle->SetOptStat("emr");
  gStyle->SetStatY(0.85);
  gStyle->SetStatX(0.87);
  gStyle->SetStatW(0.15);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes  
  TH1::SetDefaultSumw2(true);
  
  //True pT Binning
  int N_pT_Bins = 26;
  float bin_width = 2.5;
  float pT_min = 10.0;
  float Truth_pT_Array[N_pT_Bins+1]; //edges
  float pT_Centers[N_pT_Bins];
  //std::vector<float> pT_Centers;

  float pT_Widths[N_pT_Bins];

  Truth_pT_Array[0] = pT_min;

  for (int ipt = 0; ipt < N_pT_Bins; ipt++){
    float pT_Truth = ipt*bin_width+pT_min;
    float pT_Truth_plus = (ipt+1)*bin_width+pT_min;
    pT_Centers[ipt] = pT_Truth + (pT_Truth_plus - pT_Truth)/2;
    Truth_pT_Array[ipt+1] = pT_Truth_plus;
    pT_Widths[ipt] = bin_width;
  }

  //|Eta| Binning
  int N_Eta_Bins = 4;
  std::vector<float> Eta_Bins;
  std::vector<float> Eta_Centers;
  std::vector<float> Eta_Widths;
  float EtaMax = 0.45; //|eta|
  float EtaMin = 0;

  //pT and Eta Binning. Creat N-dim vector for more dependancies: [eta][obsv.][pt] 
  std::vector<std::vector<std::vector<TH1F*> > > v_Eta_pT_TH1F;
  std::vector<std::vector<TH2F*> > v_Eta_pT_TH2F;

  std::vector<TString> root_cut_string; //naming/writing unique root objects
  std::vector<TString> title_cut_string; //Additional title string with cuts (eta)

  for (int ieta = 0; ieta < N_Eta_Bins; ieta++){
    Eta_Bins.push_back(ieta*EtaMax/N_Eta_Bins);
    Eta_Widths.push_back((ieta+1)*EtaMax/N_Eta_Bins - ieta*EtaMax/N_Eta_Bins);
    Eta_Centers.push_back(Eta_Bins[ieta]+Eta_Widths[ieta]/2);
    auto Eta_Range = std::make_pair(Eta_Bins[ieta],(ieta+1)*EtaMax/N_Eta_Bins);

    root_cut_string.push_back(TString(Form("Eta_%1.0f_%1.0f",Eta_Range.first*10,Eta_Range.second*10))); //10X avoid period in savenames
    title_cut_string.push_back(TString(Form("%1.2f < |#eta| < %1.2f",Eta_Range.first,Eta_Range.second))); //add dimensionality here and above (or larger flat vector)

    std::vector<std::vector<TH1F*> > v_TH1F_Calib = make_TH1F_Calib(N_pT_Bins, Truth_pT_Array,title_cut_string[ieta],root_cut_string[ieta]);
    std::vector<TH2F*> v_TH2F = make_TH2F(title_cut_string[ieta],root_cut_string[ieta]); 
    v_Eta_pT_TH1F.push_back(v_TH1F_Calib);    
    v_Eta_pT_TH2F.push_back(v_TH2F);  
  }
  Eta_Bins.push_back(EtaMax);

  //2D Histos

  //Supplementary Plots
  TH1F* H_NExtra_Matches  = new TH1F("H_NMatches","Number of True Jets Near Reco",4,0,4);
  TH1F* H_dR = new TH1F("H_dR","#Delta R Distribution R < 0.2",200,0.0,.2);

  Long64_t nentries = _tree_event->GetEntries();

  for(Long64_t ievent = 0; ievent < nentries ; ievent++){
    _tree_event->GetEntry(ievent);
    fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);

    for (Long64_t i_reco4 = 0; i_reco4 < jet4_n; i_reco4++){

      TVector3 v3_reco;
      v3_reco.SetPtEtaPhi(jet4_pt[i_reco4],jet4_eta[i_reco4],jet4_phi[i_reco4]);

      if (not(TMath::Abs(v3_reco.Eta()) <= EtaMax) and (TMath::Abs(v3_reco.Eta()) >= EtaMin)) continue;
      // if (jet4_pt[i_reco4] < 10) continue;

      Float_t reco4_pT = v3_reco.Pt();
      bool already_matched = false;
      int Extra_Match_Count = 0;
      
      for (Long64_t i_truth4 = 0; i_truth4 < truthjet4_n; i_truth4++){

	TVector3 v3_truth;
	v3_truth.SetPtEtaPhi(truthjet4_pt[i_truth4],truthjet4_eta[i_truth4],truthjet4_phi[i_truth4]);

	Float_t truth4_pT = v3_truth.Pt();
	Float_t dphi = v3_truth.DeltaPhi(v3_reco);
	Float_t deta = (v3_truth.Eta() - v3_reco.Eta());
	Float_t dR = v3_truth.DeltaR(v3_reco);
	if (dR <= 0.2){

	  H_dR->Fill(dR);

	  float pT_Difference = (reco4_pT - truth4_pT);
	  float pT_Ratio = (reco4_pT/truth4_pT);


	  for (int ieta = 0; ieta < N_Eta_Bins; ++ieta){
	    if (v3_reco.Eta() >= Eta_Bins[ieta] && v3_reco.Eta() < Eta_Bins[ieta+1]){	    
	      
	      v_Eta_pT_TH2F[ieta][0]->Fill(truth4_pT,pT_Difference);
	      v_Eta_pT_TH2F[ieta][1]->Fill(truth4_pT,pT_Ratio);
	      v_Eta_pT_TH2F[ieta][2]->Fill(truth4_pT,reco4_pT);
	      
	      for (int ipt = 0; ipt < N_pT_Bins; ipt++){
		if (truth4_pT >= Truth_pT_Array[ipt] && truth4_pT < Truth_pT_Array[ipt+1]){
		  v_Eta_pT_TH1F[ieta][0][ipt]->Fill(pT_Difference);
		  v_Eta_pT_TH1F[ieta][1][ipt]->Fill(pT_Ratio);
		  v_Eta_pT_TH1F[ieta][2][ipt]->Fill(reco4_pT);
		  v_Eta_pT_TH1F[ieta][3][ipt]->Fill(dphi);
		  v_Eta_pT_TH1F[ieta][4][ipt]->Fill(deta);
		  
		  // pT_Differences[ipt]->Fill(pT_Difference);
		  // pT_Ratios[ipt]->Fill(pT_Ratio);
		  // Reco_pT_Slices[ipt]->Fill(reco4_pT);
		  // Phi_Deltas[ipt]->Fill(dphi);
		  // Eta_Deltas[ipt]->Fill(deta);
		}
	      }//pt
	    }
	  }//eta

	  if (already_matched)
	    Extra_Match_Count ++;
	  already_matched = true;

	}//dR
      
	H_NExtra_Matches->Fill(Extra_Match_Count);
	
      }//truth jet loop
    }// reco jet loop
  } //event loop



//Write to new root file
  TFile* fout = new TFile("Histograms_Jet_Callibration.root","RECREATE");

  H_dR->Write();
  H_NExtra_Matches->Write();

  bool write_indv_histos = false;
  if (write_indv_histos){
    for (int ieta = 0; ieta < N_Eta_Bins; ++ieta){
      v_Eta_pT_TH2F[ieta][0]->Write();
      v_Eta_pT_TH2F[ieta][1]->Write();
      v_Eta_pT_TH2F[ieta][2]->Write();

      for (int ipt = 0; ipt < N_pT_Bins; ipt++){
	v_Eta_pT_TH1F[ieta][0][ipt]->Write();
	v_Eta_pT_TH1F[ieta][1][ipt]->Write();
	v_Eta_pT_TH1F[ieta][2][ipt]->Write();
      }
    }
  }

  //Obtain Correction Fuction

  std::vector<TF1*> v_Corrections;

  for (int ieta = 0; ieta < N_Eta_Bins; ++ieta){

    auto Eta_Range = std::make_pair(Eta_Bins[ieta],Eta_Bins[ieta+1]);

    TF1* Fi = do_gaus_fitting(N_pT_Bins,pT_Centers,title_cut_string[ieta],root_cut_string[ieta],v_Eta_pT_TH1F[ieta],v_Eta_pT_TH2F[ieta]);
    v_Corrections.push_back(Fi);

    resolution_overlay(N_pT_Bins,pT_Centers,title_cut_string[ieta],root_cut_string[ieta],v_Eta_pT_TH1F[ieta][3],TString("#Delta#varphi"),v_Eta_pT_TH1F[ieta][4],TString("#Delta#eta"),
		       TString("#sigma(#Delta#varphi) #sigma(#Delta#eta) Overlay"));
  }

  //APPLY CORRECTION
  //Using pointers allows for readable variable names here, and easy fill during the loop  

  std::vector<std::vector<std::vector<TH1F*> > > v_Eta_pT_TH1F_Cor;
  std::vector<std::vector<TH2F*> > v_Eta_pT_TH2F_Cor;
  for (int ieta = 0; ieta < N_Eta_Bins; ieta++){
    auto Eta_Range = std::make_pair(Eta_Bins[ieta],Eta_Bins[ieta+1]);
    std::vector<std::vector<TH1F*> > v_TH1F_Calib = make_TH1F_Calib(N_pT_Bins, pT_Centers,title_cut_string[ieta],true);
    std::vector<TH2F*> v_TH2F = make_TH2F(title_cut_string[ieta],true);
    v_Eta_pT_TH1F_Cor.push_back(v_TH1F_Calib);
    v_Eta_pT_TH2F_Cor.push_back(v_TH2F); 
  }

  
  for(Long64_t ievent = 0; ievent < nentries ; ievent++){
    fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);

    _tree_event->GetEntry(ievent);

    for (Long64_t i_reco4 = 0; i_reco4 < jet4_n; i_reco4++){
      
      TVector3 v3_reco;
      v3_reco.SetPtEtaPhi(jet4_pt[i_reco4],jet4_eta[i_reco4],jet4_phi[i_reco4]);

      // if (jet4_pt[i_reco4] < 10) continue;
      if (not(TMath::Abs(v3_reco.Eta()) <= EtaMax) and (TMath::Abs(v3_reco.Eta()) >= EtaMin)) continue;

      for (Long64_t i_truth4 = 0; i_truth4 < truthjet4_n; i_truth4++){

	TVector3 v3_truth;
	v3_truth.SetPtEtaPhi(truthjet4_pt[i_truth4],truthjet4_eta[i_truth4],truthjet4_phi[i_truth4]);

	Float_t reco4_pT = v3_reco.Pt();
	Float_t truth4_pT = v3_truth.Pt();
	Float_t dphi = v3_truth.DeltaPhi(v3_reco);
	Float_t deta = v3_truth.Eta() - v3_reco.Eta();
	Float_t dR = v3_truth.DeltaR(v3_reco);

	if (dR <= 0.2){


          for (int ieta = 0; ieta < N_Eta_Bins; ++ieta){
            if (v3_reco.Eta() >= Eta_Bins[ieta] && v3_reco.Eta() < Eta_Bins[ieta+1]){

	      float F = v_Corrections[ieta]->Eval(reco4_pT);	      
	      float Corr4_pT = reco4_pT/F;
	      float pT_Difference = (Corr4_pT - truth4_pT);
	      float pT_Ratio = (Corr4_pT/truth4_pT);


              v_Eta_pT_TH2F_Cor[ieta][0]->Fill(truth4_pT,pT_Difference);
              v_Eta_pT_TH2F_Cor[ieta][1]->Fill(truth4_pT,pT_Ratio);
              v_Eta_pT_TH2F_Cor[ieta][2]->Fill(truth4_pT,Corr4_pT);

              for (int ipt = 0; ipt < N_pT_Bins; ipt++){
                if (truth4_pT >= Truth_pT_Array[ipt] && truth4_pT < Truth_pT_Array[ipt+1]){
                  v_Eta_pT_TH1F_Cor[ieta][0][ipt]->Fill(pT_Difference);
                  v_Eta_pT_TH1F_Cor[ieta][1][ipt]->Fill(pT_Ratio);
                  v_Eta_pT_TH1F_Cor[ieta][2][ipt]->Fill(Corr4_pT);
                  v_Eta_pT_TH1F_Cor[ieta][3][ipt]->Fill(dphi);
                  v_Eta_pT_TH1F_Cor[ieta][4][ipt]->Fill(deta); 
                }
              }
            }
          }
	}//dR
      }//truth jet
    }//reco-corrected jet    
  }//event

  for (int ieta =0; ieta < N_Eta_Bins; ++ieta){
    auto Eta_Range = std::make_pair(Eta_Bins[ieta],(ieta+1)*EtaMax/N_Eta_Bins);

    TF1* F_null = do_gaus_fitting(N_pT_Bins,pT_Centers,title_cut_string[ieta],root_cut_string[ieta],v_Eta_pT_TH1F_Cor[ieta],v_Eta_pT_TH2F_Cor[ieta],true);

    std::vector< std::vector<float> > v_Diffs_Corrected = Fit_Gauss(N_pT_Bins,v_Eta_pT_TH1F_Cor[ieta][0]);
    std::vector< std::vector<float> > v_Ratios_Corrected = Fit_Gauss(N_pT_Bins,v_Eta_pT_TH1F_Cor[ieta][1]);
    std::vector< std::vector<float> > v_Reco_Slices_Corrected = Fit_Gauss(N_pT_Bins,v_Eta_pT_TH1F_Cor[ieta][2]);
  
    Closure_Test(N_pT_Bins,pT_Centers,title_cut_string[ieta],root_cut_string[ieta],v_Diffs_Corrected,TString("Differences"), TString("p_{T}^{Corrected} - p_{T}^{Truth}"));
    Closure_Test(N_pT_Bins,pT_Centers,title_cut_string[ieta],root_cut_string[ieta],v_Ratios_Corrected,TString("Ratio"), TString("p_{T}^{Corrected} /  p_{T}^{Truth}"));

    JER_Parametrization(N_pT_Bins,pT_Centers,15,title_cut_string[ieta],root_cut_string[ieta],v_Ratios_Corrected,TString("sigma_R"),TString("R"));
    JER_Parametrization(N_pT_Bins,pT_Centers,15,title_cut_string[ieta],root_cut_string[ieta],v_Ratios_Corrected,TString("sigma_R"),TString("R"));
  }
}
