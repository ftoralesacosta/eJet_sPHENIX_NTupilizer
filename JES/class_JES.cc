#include "class_JES.h"
#include <iostream>

void MyClass::set_binning(std::vector <float> &edges, std::vector <float> &centers,
			  float low, float high, float width)
{
  if (fmod((high-low),width) != 0)
    fprintf(stderr,"WARNING: Choice of binning leaves remainder \n");

  int n_bins = int((high-low)/width);
  for (int i; i < n_bins+1; ++i)
    edges.push_back(low+i*width);
  for (int i; i < n_bins; ++i)
    centers.push_back(low + width*(0.5+i));

} //pass vectors by reference instead of editing member object to allow for locally defined 

void MyClass::set_E_binning(float low, float high, float width)
{
  set_binning(E_Bin_Edges,E_Bin_Centers,low,high,width);
}

void MyClass::set_eta_binning(float low, float high, float width)
{
  set_binning(eta_bins,eta_Centers,low,high,width);
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
      if (RecoE[n] < min_E) continue;
      if (RecoEta[n] > max_Eta) continue;
      
      for (int i = 0; i < E_Bin_Edges.size(); ++i) {
	bool leftedge = RecoE[n] >= E_Bin_Edges[i];
	bool rightedge = RecoE[n] < E_Bin_Edges[i+1];
	if (not((leftedge) and (rightedge))) continue;
	
	//FILLING
	Float_t E_diff = abs(MatchedTruthE[n]-RecoE[n]);
	Float_t E_ratio = (RecoE[n]/MatchedTruthEta[n]);
	TLorentzVector *RecoLorentz;

	std::cout<<"Line: "<<__LINE__<<std::endl;
	std::cout<<RecoPt[n]<<std::endl;
	std::cout<<RecoEta[n]<<std::endl;
	std::cout<<RecoPhi[n]<<std::endl;
	std::cout<<RecoE[n]<<std::endl;
	//FIXME: Lorentz broken
	//RecoLorentz->SetPtEtaPhiE(RecoPt[n],RecoEta[n],RecoPhi[n],RecoE[n]);
	std::cout<<"Line: "<<__LINE__<<std::endl;
	if (isnan(MatchedTruthE[n])) continue;
	TLorentzVector *TruthLorentz;
	//TruthLorentz->SetPtEtaPhiE(MatchedTruthPt[n],MatchedTruthEta[n],
	//MatchedTruthPhi[n],MatchedTruthE[n]);

	//fill_histograms(RecoLorentz, TruthLorentz, i);
	
      }
    }
  }
}

void MyClass::fill_histograms(TLorentzVector *reco, TLorentzVector *truth, int Ebin)
{
  Float_t E_diff = abs( reco->E() - truth->E() );
  Float_t E_ratio = reco->E()/truth->E();
  
  E_Differences[Ebin]->Fill(E_diff);
  E_Ratios[Ebin]->Fill(E_ratio);
  E_Slices[Ebin]->Fill(reco->E());
  Phi_Deltas[Ebin]->Fill(reco->DeltaPhi(*truth));
  Eta_Deltas[Ebin]->Fill(abs(reco->Eta() - truth->Eta()));
  Ratio_TH2F_v[Ebin]->Fill(truth->E(),E_ratio);
  Diff_TH2F_v[Ebin]->Fill(truth->E(),E_diff);
  Slice_TH2F_v[Ebin]->Fill(truth->E(),reco->E());
  //FIXME: Add Delta R
}



void MyClass::fit_histograms()
{//Get Numerical Mean+Sigma. Fit gaus. Get Gauss mean+sigma
  //Can pass output to map of vectors using map keys
  

  
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
  // for (auto TH1F_vec = TH1F_2DVec.begin(); TH1F_vec != TH1F_2DVec.end(); ++TH1F_vec)
  //   for (auto TH1_histo = TH1F_vec->begin(); TH1_histo != TH1F_vec->end(); ++TH1_histo)
  //     (*TH1_histo)->Write();


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
  
  for (auto it = E_Bin_Edges.begin(); std::next(it,1) != E_Bin_Edges.end(); ++it){

    TString root_range = Form("E%1.0f_%1.0f",*it*10,*std::next(it,1)*10);
    TString title_range = TString(Form("(%1.2f < E^{Truth} < %1.2f)",*it,*std::next(it,1)));
    float binning[3];
    
    if(strstr(root_name, TString("Slices")) == NULL)
      memcpy(&binning, &th1_binning, sizeof(th1_binning));
    else
      binning[0] = n_slice_bins; binning[1] = *it, binning[2] = *std::next(it,1);
    //energy slcie TH1F gets binning from E bin loop.
    
    TH1F* Single_Histo = new TH1F();
    Single_Histo->SetName(root_name + root_range);
    Single_Histo->SetTitle(title + title_range);
    Single_Histo->SetBins(binning[0],binning[1],binning[2]);
    
    TH1F_vector.push_back(Single_Histo);
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

  //outer eta loop soon
  for (auto it = E_Bin_Edges.begin(); std::next(it,1) != E_Bin_Edges.end(); ++it){

    TString root_range = Form("E%1.0f_%1.0f",*it*10,*std::next(it,1)*10);
    TString title_range = TString(Form("(%1.2f < E^{Truth} < %1.2f)",*it,*std::next(it,1)));
    float binning[3];

    if(strstr(root_name, TString("Slices")) == NULL)
      memcpy(&binning, &th1_binning, sizeof(th1_binning));
    else
      binning[0] = n_slice_bins; binning[1] = *it, binning[2] =	*std::next(it,1);

    TH2F* Single_Histo = new TH2F();
    Single_Histo->SetName(root_name + root_range);
    Single_Histo->SetTitle(title + title_range);
    Single_Histo->SetBins(E_Bin_Edges.size(),E_Bin_Edges.front(),E_Bin_Edges.back(),
			  binning[0], binning[1], binning[2]);

    TH2F_vector.push_back(Single_Histo);
  }
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
  //written s.t. strings specified at function call here
} 

void MyClass::write_histograms(TFile *out_file)
{

  for (auto TH1F_vec = TH1F_2DVec.begin(); TH1F_vec != TH1F_2DVec.end(); ++TH1F_vec)
    for (auto TH1_histo = TH1F_vec->begin(); TH1_histo != TH1F_vec->end(); ++TH1_histo)
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
  
  MyClass myobj;\
  myobj.set_E_binning(0,10,2);
  myobj.set_eta_binning(0,10,1);
  myobj.set_E_FitRange(2,20);
  std::cout<<"before"<<std::endl;
  myobj.reconstructed_loop(&Tree);
  //FIXME: Add TEnv capabilit
  std::cout<<"after"<<std::endl;
  myobj.set_reco_or_corr(TString("Reco"));
  myobj.initialize_histograms();  

  TFile outfile("EIC_JetEnergy_Scale+Resolution.root","recreate");

  myobj.write_histograms(&outfile);
  return 0;
}
