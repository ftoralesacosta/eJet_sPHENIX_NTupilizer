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

  //myobj.get_binning(myobj.pT_bins,myobj.pT_Centers,0,10,2);
} //pass vectors by reference instead of editing member object to allow for locally defined 

void MyClass::set_pT_binning(float low, float high, float width)
{
  set_binning(pT_bins,pT_Centers,low,high,width);
}

void MyClass::set_eta_binning(float low, float high, float width)
{
  set_binning(eta_bins,eta_Centers,low,high,width);
}

void MyClass::set_pT_FitRange(float low, float high)
{
  Fit_pT_Min = low;
  Fit_pT_Max = high;
}

TString MyClass::set_reco_or_corr(TString input_string)
{
  reco_or_corr = input_string;
}

// TTreeReader MyClass::get_TTree(TString file_name)
// {//FIXME: make a TTreeReader member variable
//   std::cout << "Opening: " << file_name << std::endl;
//   TFile *file = TFile::Open(file_name);

//   if (file == NULL) {
//     std::cout << " Failed to open "<<file_name<< std::endl;
//     exit(EXIT_FAILURE);
//   }

//   TString tree_name = "T";
//   TTree *_tree_event = dynamic_cast<TTree *>(file->Get("T"));

//   if (_tree_event == NULL) {
//     std::cout << " Failed to load TTree \""<< tree_name<<"\""<<std::endl;
//     exit(EXIT_FAILURE);
//   }

//   TTreeReader The_Tree(_tree_event);
//   return The_Tree;
//   //Would like to return pointer, but don't want to make a member variable of TTReader as I want to change it in loops, and one should never return a pointer to a local variable (local in the sense that the scope is within the function).

//   //This is because the variable the pointer points to is out of scope after the function, and this memeroy is freed. But as soon as memory is used, the pointer will point to whatever occupies that memory now.

  
//   // TTreeReaderValue<int> njets(The_Tree,"njets");
//   // TTreeReaderArray<Int_t> JetRecoNConst(The_Tree,"nComponent");
//   // TTreeReaderArray<Float_t> JetRecoE(The_Tree,"e");
//   // TTreeReaderArray<Float_t> JetRecopT(The_Tree,"pt");

//   // TTreeReaderValue<int> nAlltruthjets(The_Tree,"nAlltruthjets");
//   // TTreeReaderArray<Int_t> JetAllTruthNConst(The_Tree,"all_truthNComponent");
//   // TTreeReaderArray<Float_t> JetAllTruthE(The_Tree,"all_truthE");
//   // TTreeReaderArray<Float_t> JetAllTruthpT(The_Tree,"all_truthPt");
// }

std::vector<TH1F*> MyClass::create_TH1F(TString reco_or_corr,TString root_name,
					  TString title,float *binning)
{
  std::vector<TH1F*> TH1F_vector;

  //Internal Root Name.
  root_name = reco_or_corr + root_name;

  //Histogram title
  title = title.Insert(title.First("^")+2,reco_or_corr); //insert into latex
  
  for (auto it = pT_bins.begin(); std::next(it,1) != pT_bins.end(); ++it){

    TString root_range = Form("pT%1.0f_%1.0f",*it*10,*std::next(it,1)*10);
    TString title_range = TString(Form("(%1.2f < p_{T}^{Truth} < %1.2f)",*it,*std::next(it,1)));

    TH1F* Single_Histo = new TH1F();
    Single_Histo->SetName(root_name + root_range);
    Single_Histo->SetTitle(title + title_range);
    Single_Histo->SetBins(binning[0],binning[1],binning[2]);

    TH1F_vector.push_back(Single_Histo);
  }
  return TH1F_vector;
}

std::vector<TH2F*> MyClass::create_TH2F(TString reco_or_corr,TString root_name,
					TString title,float *binning)
{
  std::vector<TH2F*> TH2F_vector;

  TString vs_pT = " vs. p_{T}^{Truth}";
  
  //Internal Root Name.
  root_name = reco_or_corr + root_name;

  //Histogram title, insert reco/corr for latex
  title = title.Insert(title.First("^")+2,reco_or_corr);
  title = title + vs_pT;
  
  for (auto it = pT_bins.begin(); std::next(it,1) != pT_bins.end(); ++it){

    TString root_range = Form("pT%1.0f_%1.0f",*it*10,*std::next(it,1)*10);
    TString title_range = TString(Form("(%1.2f < p_{T}^{Truth} < %1.2f)",*it,*std::next(it,1)));

    TH2F* Single_Histo = new TH2F();
    Single_Histo->SetName(root_name + root_range);
    Single_Histo->SetTitle(title + title_range);
    Single_Histo->SetBins(binning[0],binning[1],binning[2],
			  binning[0], binning[1], binning[2]);

    TH2F_vector.push_back(Single_Histo);
  }
  return TH2F_vector;
}

//void get_calibration()

void MyClass::initialize_histograms(){

  float binning[3] = {16,-2,2}; //{nbins,min,max}
  //FIXME: Pass this as argument. For pT_Slices copy pT truth binning member variable
  pT_Differences = create_TH1F("Reco","_pT_Difference_","p_{T}^{} - p_{T}^{Truth} ",binning);
  pT_Ratios = create_TH1F("Reco","_pT_Ratio_","p_{T}^{} / p_{T}^{Truth} ",binning);
  pT_Slices = create_TH1F("Reco","_pT_Slices","p_{T}^{} ",binning);
  
  Phi_Deltas = create_TH1F("Reco","_delta_phi_","#varphi^{} - #varphi^{Truth} ",binning);
  Eta_Deltas = create_TH1F("Reco","_delta_eta_","#eta^{} - #eta^{Truth} ",binning);
  
  Ratio_TH2F_v = create_TH2F("Reco","_pT_Ratio_","p_{T}^{} / p_{T}^{Truth} ",binning);
  Diff_TH2F_v = create_TH2F("Reco","_pT_Difference_","p_{T}^{} - p_{T}^{Truth} ",binning);
  Slice_TH2F_v = create_TH2F("Reco","_pT_Slices","p_{T}^{} ",binning);
}

int main(int argc, char *argv[])
{
  if (argc < 2)
    {
      std::cout<<"Please Specify the ROOT File: [command] [file]"<<std::endl;
      exit(EXIT_FAILURE);
    }
  
  MyClass myobj;
  //myobj.get_TTree(TString(argv[1]));
  //Can loop through argv[i] for multiple files

  TFile *file = TFile::Open(TString(argv[1]));
  TTreeReader The_Tree("T",file); //T is the Tree name within the file
  
  myobj.set_pT_binning(0,10,2);
  myobj.set_eta_binning(0,10,1);
  myobj.set_pT_FitRange(2,20);
  //FIXME: Add TEnv capabilit

  myobj.initialize_histograms();
  
  
  //float 2d binning = {idk idk idk idk}
  //calll myobj.create_TH1F_Calib the same way
  return 0;
}
