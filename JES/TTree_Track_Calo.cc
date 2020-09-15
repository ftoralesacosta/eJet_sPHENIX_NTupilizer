#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"
#include <TLegend.h>
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
#include <TMath.h>
#include <vector>
#include <math.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

float calc_Q_square(float inE, TLorentzVector v)
{
  return 2.*inE*v.E()*(1-TMath::Abs(v.CosTheta()));
}

int main(int argc, char *argv[])
{
  if (argc < 3) {
    std::cout<<"Syntax: [Command] [Track File] [Calo File]"<<std::endl;
    exit(EXIT_FAILURE);
  }

  //FILE 1
  int iarg = 1;
  TString root_file = (TString)argv[iarg];
  
  std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
  
  TFile *file = TFile::Open(root_file);
  
  if (file == NULL) {
    std::cout << " File Fail" << std::endl;
    exit(EXIT_FAILURE); 
  } 
  
  file->Print();
  
  TTree *_tree_event = dynamic_cast<TTree *>(file->Get("T"));
  
  if (_tree_event == NULL) {
    std::cout << " Tree Fail " << std::endl;
    exit(EXIT_FAILURE);
  }

  //FILE 2
  iarg = 2;
  TString root_file2 = (TString)argv[iarg];
  
  std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
  
  TFile *file2 = TFile::Open(root_file2);
  
  if (file2 == NULL) {
    std::cout << " File 2 Fail" << std::endl;
    exit(EXIT_FAILURE); 
  } 
  
  file->Print();
  
  TTree *_tree_event2 = dynamic_cast<TTree *>(file2->Get("T"));

  if (_tree_event2 == NULL) {
    std::cout << " Tree Fail " << std::endl;
    exit(EXIT_FAILURE);
  }


  
  //gStyle Plotting
  gStyle->SetOptStat("emr");
  gStyle->SetStatY(0.85);
  gStyle->SetStatX(0.87);
  gStyle->SetStatW(0.15);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes  
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  //Track Tracking
  TTreeReader Tree_Track(_tree_event);
  TTreeReaderValue<int> Track_njets(Tree_Track,"njets");
  TTreeReaderArray<Int_t> Track_JetRecoNConst(Tree_Track,"nComponent");
  TTreeReaderArray<Float_t> Track_JetRecoE(Tree_Track,"e");
  TTreeReaderArray<Float_t> Track_JetRecopT(Tree_Track,"pt");
  
  TTreeReaderArray<Int_t> Track_JetMatchedTruthNConst(Tree_Track,"matched_truthNComponent");
  TTreeReaderArray<Float_t> Track_JetMatchedTruthE(Tree_Track,"matched_truthE");
  TTreeReaderArray<Float_t> Track_JetMatchedTruthpT(Tree_Track,"matched_truthPt");

  TTreeReaderValue<int> Track_nAlltruthjets(Tree_Track,"nAlltruthjets");
  TTreeReaderArray<Int_t> Track_JetAllTruthNConst(Tree_Track,"all_truthNComponent");
  TTreeReaderArray<Float_t> Track_JetAllTruthE(Tree_Track,"all_truthE");
  TTreeReaderArray<Float_t> Track_JetAllTruthpT(Tree_Track,"all_truthPt");


  //Calo
  TTreeReader Tree_Calo(_tree_event2);
  TTreeReaderValue<int> Calo_njets(Tree_Calo,"njets");
  TTreeReaderArray<Int_t> Calo_JetRecoNConst(Tree_Calo,"nComponent");
  TTreeReaderArray<Float_t> Calo_JetRecoE(Tree_Calo,"e");
  TTreeReaderArray<Float_t> Calo_JetRecopT(Tree_Calo,"pt"); 
  
  TTreeReaderArray<Int_t> Calo_JetMatchedTruthNConst(Tree_Calo,"matched_truthNComponent");
  TTreeReaderArray<Float_t> Calo_JetMatchedTruthE(Tree_Calo,"matched_truthE");
  TTreeReaderArray<Float_t> Calo_JetMatchedTruthpT(Tree_Calo,"matched_truthPt");

  TTreeReaderValue<int> Calo_nAlltruthjets(Tree_Calo,"nAlltruthjets");
  TTreeReaderArray<Int_t> Calo_JetAllTruthNConst(Tree_Calo,"all_truthNComponent");
  TTreeReaderArray<Float_t> Calo_JetAllTruthE(Tree_Calo,"all_truthE");
  TTreeReaderArray<Float_t> Calo_JetAllTruthpT(Tree_Calo,"all_truthPt");

 //2D Histos

  TH2F * Calo_RecovTruth = new TH2F("Calo_Reco_v_Truth", "E^{Reco}_{Calo} vs E^{True}_{Calo}",100,0,50,100,0,50);
  TH2F * Track_RecovTruth = new TH2F("Track_Reco_v_Truth", "E^{Reco}_{Track} vs E^{True}_{Track}",100,0,50,100,0,50);
  
  TH2F * Calo_v_Hybdrid_Truth_Energy = new TH2F("TRUTH_Energy_CalovHybdrid", "E^{True}_{Calo} vs E^{True}_{Track}",100,0,50,100,0,50);
  Calo_v_Hybdrid_Truth_Energy->GetXaxis()->SetTitle("E^{True}_{Calo}");
  Calo_v_Hybdrid_Truth_Energy->GetXaxis()->CenterTitle(true);
  Calo_v_Hybdrid_Truth_Energy->GetYaxis()->SetTitle("E^{True}_{Track}");
  Calo_v_Hybdrid_Truth_Energy->GetYaxis()->CenterTitle(true);
  TH2F * Calo_v_Hybdrid_Reco_Energy = new TH2F("RECO_Energy_CalovHybdrid", "E^{Reco}_{Calo} vs E^{Reco}_{Track}",100,0,50,100,0,50);
  
  TH2F * Calo_v_Hybdrid_Truth_pT = new TH2F("TRUTH_pT_CalovHybdrid", "pT^{True}_{Calo} vs pT^{True}_{Track}",100,0,50,100,0,50);
  TH2F * Calo_v_Hybdrid_Reco_pT = new TH2F("RECO_pT_CalovHybdrid", "pT^{Reco}_{Calo} vs pT^{Reco}_{Track}",100,0,50,100,0,50);
  
  TH2F * Similar_Calo_v_Hybdrid_Reco_Energy = new TH2F("Similar_RECO_Energy_CalovHybdrid", "E^{Reco}_{Calo} vs E^{Reco}_{Track}",100,0,50,100,0,50);
  
  int counter = 0;
  Long64_t nentries = _tree_event->GetEntries();
  float MinE = 3.0; //Min Jet GeV
  int N_Const_Min = 3;
  while (Tree_Track.Next()){
    Tree_Calo.Next(); //TTrees should have same Nentries

   //fprintf(stderr,"\r%d: Event %i / %llud",__LINE__,counter,nentries);
   int Track_HardestN = 0;
   float Track_HardestE = 0;

   int Calo_HardestN = 0;
   float Calo_HardestE = 0;

   
   for (Int_t n = 0; n < *Track_njets; n++){
     if (Track_JetRecoNConst[n] < N_Const_Min) continue;
     if (isnan(Track_JetMatchedTruthE[n]))
       continue;

     Track_RecovTruth->Fill(Track_JetMatchedTruthE[n],Track_JetRecoE[n]);
     if (Track_JetRecoE[n] < MinE) continue;
     if (Track_JetRecoE[n] < Track_HardestE) continue;
     else {
       Track_HardestE = Track_JetRecoE[n];
       Track_HardestN = n;
     }
   }

   for (Int_t n = 0; n < *Calo_njets; n++){
     //fprintf(stderr,"%d: Calo Jet N = %i\n",__LINE__,n);
     if (Calo_JetRecoNConst[n] < N_Const_Min) continue;
     //if (not(isnan(Calo_JetMatchedTruthE[n])))
     if (isnan(Calo_JetMatchedTruthE[n]))
       continue;
     Calo_RecovTruth->Fill(Calo_JetMatchedTruthE[n],Calo_JetRecoE[n]);
     if (Calo_JetRecoE[n] < MinE) continue;
     if (Calo_JetRecoE[n] < Calo_HardestE) continue;
     else {       
       Calo_HardestE = Calo_JetRecoE[n];
       Calo_HardestN = n;      
     }
   }

   //Fill Histo with Hardest Jet
   //if ()
   fprintf(stderr,"%d: Track_HardestE = %1.3f, Calo_HardestE = %1.3f\n",__LINE__,Track_HardestE,Calo_HardestE);
   Calo_v_Hybdrid_Reco_Energy->Fill(Track_HardestE,Calo_HardestE);
   Calo_v_Hybdrid_Reco_pT->Fill(Track_JetRecopT[Track_HardestN],Calo_JetRecopT[Calo_HardestN]);


   //Truth Jet Loop
   for (Int_t n = 0; n < *Track_nAlltruthjets; n++){
     if (Track_JetAllTruthNConst[n] < N_Const_Min) continue;
     if (Track_JetAllTruthE[n] < MinE) continue;

     Calo_v_Hybdrid_Truth_Energy->Fill(Track_JetAllTruthE[n],Calo_JetAllTruthE[n]);
     Calo_v_Hybdrid_Truth_pT->Fill(Track_JetAllTruthpT[n],Calo_JetAllTruthpT[n]);
   }
   counter++;
 }

 Tree_Track.Restart();
 Tree_Calo.Restart();

 // while(Tree_Track.Next()){
 //   Tree_Calo.Next();
   
 //   fprintf(stderr,"\r%d: Event %i / %llud",__LINE__,counter,nentries);
 //   std::pair<int, int> Similar_Track_Calo = std::make_pair(0,0);// Track first, Calo Second
 //   float Delta_E = 1000.0;//reset each event, with the aim of finding the closest pair of jets per event
 //   for (int n_hyb = 0; n_hyb < *Track_njets; n_hyb++){
 //     if (Track_JetRecoNConst[n_hyb] < N_Const_Min) continue;
 //     if (Track_JetRecoE[n_hyb] < MinE) continue;
     
 //     //float Delta_hyb = 1000.0; //Delta for middle loop;
 //     float Track_E = Track_JetRecoE[n_hyb];

 //     float Track_pT = Track_JetRecopT[n_hyb];

 //     for (int n_emc = 0; n_emc < *Calo_njets; n_emc++){
 //       if (Calo_JetRecoNConst[n_emc] < N_Const_Min) continue;
 //       if (Calo_JetRecoE[n_emc] < MinE) continue;

 //       float Calo_E = Calo_JetRecoE[n_emc];

 //       if (abs(Calo_E-Track_E)/Track_E < Delta_E) continue;
 //   	 Delta_E = abs(Calo_E-Track_E)/Track_E;
 // 	 float Calo_pT = Calo_JetRecopT[n_emc];
 // 	 Similar_Track_Calo = std::make_pair(n_hyb,n_emc);
 //     }
 //     Similar_Calo_v_Hybdrid_Reco_Energy->Fill(Track_JetRecoE[Similar_Track_Calo.first],Calo_JetRecoE[Similar_Track_Calo.second]); 
 //   }
 // }
//Write to new root file
  TFile* fout = new TFile("Histograms_Jet_Callibration.root","RECREATE");

  Track_RecovTruth->GetXaxis()->SetTitle("E^{True}_{Jet} [GeV]");
  Track_RecovTruth->GetYaxis()->SetTitle("E^{Reco}_{Jet} [GeV]");
  Track_RecovTruth->Write();
  Calo_RecovTruth->GetXaxis()->SetTitle("E^{True}_{Jet} [GeV]");
  Calo_RecovTruth->GetYaxis()->SetTitle("E^{Reco}_{Jet} [GeV]");  
  Calo_RecovTruth->Write();
  
  Calo_v_Hybdrid_Truth_Energy->Write();
  Calo_v_Hybdrid_Reco_Energy->Write();
  Calo_v_Hybdrid_Truth_pT->Write();
  Calo_v_Hybdrid_Reco_pT->Write();
  Similar_Calo_v_Hybdrid_Reco_Energy->Write();
  
  fout->Close();
 
}

