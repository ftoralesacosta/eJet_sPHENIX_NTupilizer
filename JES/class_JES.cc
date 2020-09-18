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

void MyClass::apply_cuts(TTreeReader* Tree){

  TTreeReaderValue<int> njets(*Tree,"njets");
  TTreeReaderArray<Int_t> JetRecoNConst(*Tree,"nComponent");
  TTreeReaderArray<Float_t> JetRecoE(*Tree,"e");
  TTreeReaderArray<Float_t> JetRecoEta(*Tree,"eta");

  TTreeReaderArray<Int_t> Track_JetMatchedTruthNConst(*Tree,"matched_truthNComponent");
  TTreeReaderArray<Float_t> Track_JetMatchedTruthE(*Tree,"matched_truthE");
  TTreeReaderArray<Float_t> Track_JetMatchedTruthEta(*Tree,"matched_truthEta");
  
  while (Tree->Next()){
    for (int n = 0; n < *njets; ++n) {
      if (JetRecoNConst[n] < min_N) continue;
      if (JetRecoE[n] < min_E) continue;
      if (JetRecoEta[n] < min_Eta) continue;

      
    }
  }
}

// void MyClass::reco_loop(TTreeReader Tree)
// {
//   TTreeReaderArray<Int_t> Track_JetRecoNConst(Tree,"nComponent");
//   TTreeReaderArray<Float_t> Track_JetRecoE(Tree,"e");
//   TTreeReaderArray<Float_t> Track_JetRecoE(Tree,"pt");

//   TTreeReaderArray<Int_t> Track_JetMatchedTruthNConst(Tree,"matched_truthNComponent");
//   TTreeReaderArray<Float_t> Track_JetMatchedTruthE(Tree,"matched_truthE");
//   TTreeReaderArray<Float_t> Track_JetMatchedTruthE(Tree,"matched_truthPt");

//   for (auto event = jet_list.begin(); event != jet_list.end(); ++event) {
//     for (auto jet = event->begin(); jet != event->end(); ++jet) {
      
//     }
//   }
// }

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
    //energy th1f get binning from loop
    
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
  //written s.t. strings specified at function call
} 

void MyClass::write_histograms(TFile *out_file)
{
  for (auto it = E_Differences.begin(); it != E_Differences.end(); ++it)
    (*it)->Write();
  //FIXME: write for other TH1Vecs and TH2Vecs
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
  myobj.apply_cuts(&Tree);
  myobj.set_E_binning(0,10,2);
  myobj.set_eta_binning(0,10,1);
  myobj.set_E_FitRange(2,20);
  //FIXME: Add TEnv capabilit

  myobj.set_reco_or_corr(TString("Reco"));
  myobj.initialize_histograms();  

  TFile outfile("EIC_JetEnergy_Scale_Resolution.root","recreate");

  myobj.write_histograms(&outfile);
  //calll myobj.create_TH1F_Calib the same way
  return 0;
}
