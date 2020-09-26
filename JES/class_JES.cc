#include "class_JES.h"
#include <iostream>

void MyClass::set_binning(std::vector <float> &edges, std::vector <float> &centers,
			  float low, float high, float width)
{
  if (fmod((high-low),width) != 0)
    fprintf(stderr,"WARNING: Choice of binning yields remainder \n");

  int n_bins = int((high-low)/width);
  for (int i; i < n_bins+1; ++i)
    edges.push_back(low+i*width);
  for (int i; i < n_bins; ++i)
    centers.push_back(low + width*(0.5+i));

} //pass vectors by reference instead of editing member object to allow for locally defined 

void MyClass::set_E_binning()
{//// FIXME: just parameters in header. Use function to initialize vector of binning
  set_binning(E_Bin_Edges,E_Bin_Centers,E_Min,E_Max,E_Bin_Width);
}

void MyClass::set_eta_binning()
{
  set_binning(eta_Bin_Edges,eta_Bin_Centers,eta_Min,eta_Max,eta_Bin_Width);
}

void MyClass::set_E_FitRange(float low, float high)
{//range from which JES/Correction should be obtained
  //often low E slices don't have good gaussian fits.
  Fit_E_Min = low;
  Fit_E_Max = high;
}

TString MyClass::set_reco_or_corr(TString input_string)
{
  reco_or_corr = input_string;
}

void MyClass::reconstructed_loop(TTreeReader* Tree){

  TTreeReaderValue<int> njets(*Tree,"njets");
  TTreeReaderArray<Int_t> RecoNConst(*Tree,"nComponent");
  TTreeReaderArray<Float_t> RecoE(*Tree,"e");
  TTreeReaderArray<Float_t> RecoEta(*Tree,"eta");
  TTreeReaderArray<Float_t> RecoPhi(*Tree,"phi");
  TTreeReaderArray<Float_t> RecoPt(*Tree,"pt");

  TTreeReaderArray<Int_t> MatchedTruthNConst(*Tree,"matched_truthNComponent");
  TTreeReaderArray<Float_t> MatchedTruthE(*Tree,"matched_truthE");
  TTreeReaderArray<Float_t> MatchedTruthEta(*Tree,"matched_truthEta");
  TTreeReaderArray<Float_t> MatchedTruthPhi(*Tree,"matched_truthPhi");
  TTreeReaderArray<Float_t> MatchedTruthPt(*Tree,"matched_truthPt");

  while (Tree->Next()){
    
    for (int n = 0; n < *njets; ++n) {
      
      //CUTS
      if (RecoNConst[n] < min_N) continue;
      if (RecoE[n] < E_Min) continue;
      if (RecoEta[n] > eta_Max) continue;
      //if (isnan(MatchedTruthE[n])) continue;
      
      for (int i_eta = 0; i_eta < eta_Bin_Edges.size(); ++i_eta) {
  	bool eta_leftedge = RecoE[n] >= E_Bin_Edges[i_eta];
	bool eta_rightedge = RecoE[n] < E_Bin_Edges[i_eta+1];
	if (not((eta_leftedge) and (eta_rightedge))) continue;	

	for (int i_E= 0; i_E< E_Bin_Edges.size(); ++i_E) {
	  bool E_leftedge = RecoE[n] >= E_Bin_Edges[i_E];
	  bool E_rightedge = RecoE[n] < E_Bin_Edges[i_E+1];
	  if (not((E_leftedge) and (E_rightedge))) continue;
	  std::cout<<E_Bin_Edges[i_E]<<" < "<<RecoE[n]<<" < "<<E_Bin_Edges[i_E+1]<<std::endl;
	  
	  //FILLING
	  Float_t E_diff = abs(MatchedTruthE[n]-RecoE[n]);
	  Float_t E_ratio = (RecoE[n]/MatchedTruthEta[n]);
	  TLorentzVector RecoLorentz;
	  RecoLorentz.SetPtEtaPhiE(RecoPt[n],RecoEta[n],RecoPhi[n],RecoE[n]);
	  TLorentzVector TruthLorentz;
	  TruthLorentz.SetPtEtaPhiE(MatchedTruthPt[n],MatchedTruthEta[n],
				    MatchedTruthPhi[n],MatchedTruthE[n]);
	  
	  fill_histograms(RecoLorentz, TruthLorentz, i_eta, i_E);
	
	}
      }
    }
  }
}

void MyClass::fill_histograms(TLorentzVector reco, TLorentzVector truth, int etabin, int Ebin)
{//TLorentz used for DeltaPhi function
  
  int index = Ebin + etabin*Ebin;
  Float_t E_diff = abs( reco.E() - truth.E() );
  Float_t E_ratio = reco.E()/truth.E();

  std::cout<<"Filling Histos"<<std::endl;
  std::cout<<"size of diffs = "<<E_Differences.size()<<std::endl;
  E_Differences[index]->Fill(E_diff);
    std::cout<<"Filling Histos"<<std::endl;
  E_Ratios[index]->Fill(E_ratio);
  E_Slices[index]->Fill(reco.E());
  Phi_Deltas[index]->Fill(reco.DeltaPhi(truth));
  Eta_Deltas[index]->Fill(abs(reco.Eta() - truth.Eta()));
  // Ratio_TH2F_v[index]->Fill(truth.E(),E_ratio);
  // Diff_TH2F_v[index]->Fill(truth.E(),E_diff);
  // Slice_TH2F_v[index]->Fill(truth.E(),reco.E());
  //FIXME: Add Delta R
}

void MyClass::fit_histograms()
{  
  //for (auto it = all_TH1F_vecs.begin(); it != all_TH1F_vecs.end(); ++it) {
  //for (auto TH1F_vec = all_TH1F_vecs.begin(); TH1F_vec != all_TH1F_vecs.end(); ++TH1F_vec){

  for (int i = 0; i < all_TH1F_vecs.size(); ++i) {

    //for (auto it = TH1F_vec->begin(); it != TH1F_vec->end(); ++it){
    std::pair<float,float> mean_sigma = {NAN, NAN};
    for (auto it = all_TH1F_vecs[i]->begin(); it != all_TH1F_vecs[i]->end(); ++it){
      
      float num_mean = (*it)->GetMean();
      //std::cout<<"mean ="<<num_mean<<std::endl;
      float num_stdev = (*it)->GetStdDev();
      float fit_min = num_mean-(stdv_range*num_stdev);
      float fit_max = num_mean+(stdv_range*num_stdev);
      TF1 *gaus_fit = new TF1("gaus_fit","gaus",fit_min,fit_max);
      //(*it)->Fit(gaus_fit,"QR");
      mean_sigma = {gaus_fit->GetParameters()[1],gaus_fit->GetParameters()[2]};
      all_mean_vecs[i].push_back(mean_sigma);
    }
  }
  // std::cout<<all_mean_vecs[0].size()<<std::endl;
  // for (int i = 0; i < all_mean_vecs[0].size(); ++i){
  //   std::cout<<"Diff_means = "<<all_mean_vecs[0][i].first<<std::endl;
  //   std::cout<<"Diff_means = "<<E_Differences[i]->GetMean()<<std::endl;
  //}
  /*for (int ipt = 0; ipt < N_pT_Bins; ipt++){
      float c = 1.5;
      float mean = Histo[ipt]->GetMean();
      float stdev = Histo[ipt]->GetStdDev();
      Numerical_Mean.push_back(mean);
      Numerical_Sigma.push_back(stdev);

      float num_min = mean-(stdv_range*stdev);
      float num_max = mean+(stdv_range*stdev);

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
      Gauss_Sigma_Error.push_back(sigma_error);*/
      
      //stdv_range


  // for (auto it = E_Differences.begin(); it != E_Differences.end(); ++it) (*it)->Write();
  // for (auto it = E_Ratios.begin(); it != E_Ratios.end(); ++it) (*it)->Write();
  // for (auto it = E_Slices.begin(); it != E_Slices.end(); ++it) (*it)->Write();
  // for (auto it = Phi_Deltas.begin(); it != Phi_Deltas.end(); ++it) (*it)->Write();
  // for (auto it = Eta_Deltas.begin(); it != Eta_Deltas.end(); ++it) (*it)->Write();
  
}

std::vector<TH1F*> MyClass::create_TH1F(TString root_name,
					TString title)
{
  std::vector<TH1F*> TH1F_vector;

  root_name = reco_or_corr + root_name;  
  title = title.Insert(title.First("^")+2,reco_or_corr); //for latex

  //Eta loop here
  for (auto i_eta = eta_Bin_Edges.begin(); i_eta != eta_Bin_Edges.end(); ++i_eta) {    
    for (auto i_E = E_Bin_Edges.begin(); std::next(i_E,1) != E_Bin_Edges.end(); ++i_E){

      TString root_range = Form("eta_%1.0f_%1.0X10__E_%1.0f_%1.0fX10",
				*i_eta*19,*std::next(i_eta,1)*10,*i_E*10,*std::next(i_E,1)*10);

      TString title_range = TString(Form("(%1.2f < E^{Truth} < %1.2f GeV, %1.2f < #eta %1.2fB)",
					 *i_eta,*std::next(i_eta,1), *i_E,*std::next(i_E,1)));
      float binning[3];
      if(strstr(root_name, TString("Slices")) == NULL)
	memcpy(&binning, &th1_binning, sizeof(th1_binning));
      else
	binning[0] = n_slice_bins; binning[1] = *i_E, binning[2] = *std::next(i_E,1);
      //energy slice TH1F gets binning from E bin loop.
    
      TH1F* Single_Histo = new TH1F();
      Single_Histo->SetName(root_name + root_range);
      Single_Histo->SetTitle(title + title_range);
      Single_Histo->SetBins(binning[0],binning[1],binning[2]);
    
      TH1F_vector.push_back(Single_Histo);
    }
  }
  return TH1F_vector;
}

std::vector<TH2F*> MyClass::create_TH2F(TString root_name,
					TString title)
{
  std::vector<TH2F*> TH2F_vector;
  TString vs_E = " vs. E^{Truth}"; //common x-axis
  
  if (E_Bin_Edges.empty())
    fprintf(stderr,"%s:%d: Need to call set_E_binning first \n",__func__,__LINE__);
  
  root_name = reco_or_corr + root_name;
  title = title.Insert(title.First("^")+2,reco_or_corr);
  title = title + vs_E;

  TString root_range = Form("E_%1.0f_%1.0fX10",*E_Bin_Edges.begin()*10,*E_Bin_Edges.end());

  TString title_range = TString(Form("(%1.2f < E^{Truth} < %1.2f GeV)",*E_Bin_Edges.begin(), *E_Bin_Edges.end()));

  float y_binning[3];
  if(strstr(root_name, TString("Slices")) == NULL)
    memcpy(&y_binning, &th1_binning, sizeof(th1_binning));
  else{
    y_binning[0] = E_Bin_Centers.size();
    y_binning[1] = *E_Bin_Edges.begin();
    y_binning[2] = *E_Bin_Edges.end();
  }
  
  float x_binning[3];
  if (strstr(root_name, TString("E_")) != NULL){
    x_binning[0] = E_Bin_Centers.size();
    x_binning[1] = *E_Bin_Edges.begin();
    x_binning[2] = *E_Bin_Edges.end();
  }
  
  else if (strstr(root_name, TString("eta_")) != NULL){
    x_binning[0] = eta_Bin_Centers.size();
    x_binning[1] = *eta_Bin_Edges.begin();
    x_binning[2] = *eta_Bin_Edges.end();
  }
  
  TH2F* Single_Histo = new TH2F();
  Single_Histo->SetName(root_name + root_range);
  Single_Histo->SetTitle(title + title_range);
  Single_Histo->SetBins(x_binning[0],x_binning[1],x_binning[3],
			y_binning[0],y_binning[1],y_binning[3]);
  
  TH2F_vector.push_back(Single_Histo);
  
  return TH2F_vector;
}

void MyClass::initialize_histograms()
{
  E_Differences = create_TH1F("_E_Difference_","E^{} - E^{Truth} ");
  E_Ratios = create_TH1F("_E_Ratio_","E^{} / E^{Truth} ");
  E_Slices = create_TH1F("_E_Slices","E^{} ");
  
  Phi_Deltas = create_TH1F("_delta_phi_","#varphi^{} - #varphi^{Truth} ");
  Eta_Deltas = create_TH1F("_delta_eta_","#eta^{} - #eta^{Truth} ");
  
  Ratio_TH2F_v = create_TH2F("_E_Ratio_","E^{} / E^{Truth} ");
  Diff_TH2F_v = create_TH2F("_E_Difference_","E^{} - E^{Truth} ");
  Slice_TH2F_v = create_TH2F("_E_Slices","E^{} ");
  //written s.t. strings specified at  function call here
} 

void MyClass::write_histograms(TFile *out_file)
{
  for (auto TH1F_vec = all_TH1F_vecs.begin(); TH1F_vec != all_TH1F_vecs.end(); ++TH1F_vec)
    for (auto TH1_histo = (*TH1F_vec)->begin(); TH1_histo != (*TH1F_vec)->end(); ++TH1_histo)
      (*TH1_histo)->Write();
  
  for (auto TH2F_pair = TH2F_map.begin(); TH2F_pair != TH2F_map.end(); ++TH2F_pair) //{
    for (auto iTH2F = TH2F_pair->second.begin(); iTH2F != TH2F_pair->second.begin(); iTH2F++)
      (*iTH2F)->Write();

  // for (auto it = E_Differences.begin(); it != E_Differences.end(); ++it) (*it)->Write();
  // for (auto it = E_Ratios.begin(); it != E_Ratios.end(); ++it) (*it)->Write();
  // for (auto it = E_Slices.begin(); it != E_Slices.end(); ++it) (*it)->Write();
  // for (auto it = Phi_Deltas.begin(); it != Phi_Deltas.end(); ++it) (*it)->Write();
  // for (auto it = Eta_Deltas.begin(); it != Eta_Deltas.end(); ++it) (*it)->Write();

  // for (auto it = Ratio_TH2F_v.begin(); it != Ratio_TH2F_v.end(); ++it) (*it)->Write();
  // for (auto it = Diff_TH2F_v.begin(); it != Diff_TH2F_v.end(); ++it) (*it)->Write();
  // for (auto it = E_Slices.begin(); it != E_Slices.end(); ++it) (*it)->Write();
}


int main(int argc, char *argv[])
{
  if (argc < 2)
    {
      std::cout<<"Please Specify the ROOT File: [command] [file]"<<std::endl;
      exit(EXIT_FAILURE);
    }
  
  TFile *file = TFile::Open(TString(argv[1]));
  TTreeReader Tree("T",file); //T is the Tree name within the file
  
  MyClass myobj;
  //myobj.set_E_binning(0,80,1);
  //myobj.set_eta_binning(0,10,1);
  myobj.set_E_FitRange(2,80);
  myobj.reconstructed_loop(&Tree);
  //FIXME: Add TEnv capabilit
  myobj.set_reco_or_corr(TString("Reco"));
  myobj.initialize_histograms();  

  TFile outfile("EIC_JetEnergy_Scale+Resolution.root","recreate");

  myobj.fit_histograms();
  
  myobj.write_histograms(&outfile);
  return 0;
}
