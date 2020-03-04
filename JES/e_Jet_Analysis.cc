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

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

const int MAX_INPUT_LENGTH = 200;


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
  
  TTree *_tree_event = dynamic_cast<TTree *>(file->Get("T"));

  if (_tree_event == NULL) {
    std::cout << " Tree Fail " << std::endl;
    exit(EXIT_FAILURE);
  }

  //Declare Leaf Types
  Int_t           m_event;
  Int_t           id;
  Int_t           nComponent;
  Float_t         eta;
  Float_t         phi;
  Float_t         e;
  Float_t         pt;
  Int_t           truthID;
  Int_t           truthNComponent;
  Float_t         truthEta;
  Float_t         truthPhi;
  Float_t         truthE;
  Float_t         truthPt;
  Int_t           nMatchedTrack;
  //Electron Truth Variables
  Float_t         etruthEta;
  Float_t         etruthPhi;
  Float_t         etruthE;
  Float_t         etruthPt;
  Float_t         etruthpX;
  Float_t         etruthpY;
  Float_t         etruthpZ;
  Int_t           etruthPID;
  Int_t           etruthParentID;

  //Declare Branches
  TBranch        *b_event;
  TBranch        *b_id;
  TBranch        *b_nComponent;
  TBranch        *b_eta;
  TBranch        *b_phi;
  TBranch        *b_e;
  TBranch        *b_pt;
  TBranch        *b_truthID;
  TBranch        *b_truthNComponent;
  TBranch        *b_truthEta;
  TBranch        *b_truthPhi;
  TBranch        *b_truthE;
  TBranch        *b_truthPt;
  TBranch        *b_nMatchedTrack;

  TBranch        *b_etruthEta;
  TBranch        *b_etruthPhi;
  TBranch        *b_etruthE;
  TBranch        *b_etruthPt;
  TBranch        *b_etruthpX;
  TBranch        *b_etruthpY;
  TBranch        *b_etruthpZ;
  TBranch        *b_etruthPID;
  TBranch        *b_etruthParentID;

  
  _tree_event->SetBranchAddress("m_event", &m_event, &b_event);
  _tree_event->SetBranchAddress("id", &id, &b_id);
  _tree_event->SetBranchAddress("nComponent", &nComponent, &b_nComponent);
  _tree_event->SetBranchAddress("eta", &eta, &b_eta);
  _tree_event->SetBranchAddress("phi", &phi, &b_phi);
  _tree_event->SetBranchAddress("e", &e, &b_e);
  _tree_event->SetBranchAddress("pt", &pt, &b_pt);
  _tree_event->SetBranchAddress("truthID", &truthID, &b_truthID);
  _tree_event->SetBranchAddress("truthNComponent", &truthNComponent, &b_truthNComponent);
  _tree_event->SetBranchAddress("truthEta", &truthEta, &b_truthEta);
  _tree_event->SetBranchAddress("truthPhi", &truthPhi, &b_truthPhi);
  _tree_event->SetBranchAddress("truthE", &truthE, &b_truthE);
  _tree_event->SetBranchAddress("truthPt", &truthPt, &b_truthPt);
  _tree_event->SetBranchAddress("nMatchedTrack", &nMatchedTrack, &b_nMatchedTrack);
  
  _tree_event->SetBranchAddress("etruthEta", &etruthEta, &b_etruthEta);
  _tree_event->SetBranchAddress("etruthPhi", &etruthPhi, &b_etruthPhi);
  _tree_event->SetBranchAddress("etruthE", &etruthE, &b_etruthE);
  _tree_event->SetBranchAddress("etruthPt", &etruthPt, &b_etruthPt);
  _tree_event->SetBranchAddress("etruthpX", &etruthpX, &b_etruthpX);
  _tree_event->SetBranchAddress("etruthpY", &etruthpY, &b_etruthpY);
  _tree_event->SetBranchAddress("etruthpZ", &etruthpZ, &b_etruthpZ);
  _tree_event->SetBranchAddress("etruthPt", &etruthPt, &b_etruthPt);
  _tree_event->SetBranchAddress("etruthPID", &etruthPID, &b_etruthPID);
  _tree_event->SetBranchAddress("etruthParentID", &etruthParentID, &b_etruthParentID);

  int jet4_n;
  Float_t reco_jet_pt[200];
  Float_t reco_jet_E[200];
  Float_t reco_jet_m[200];
  int reco_jet_ncomp[250];
  Float_t reco_jet_eta[200];
  Float_t reco_jet_phi[200];
  float reco_jet_E_layer0[200];
  float reco_jet_E_layer1[200];
  float reco_jet_E_layer2[200];

  int truth_jet_n;
  Float_t truth_jet_pt[200];
  Float_t truth_jet_m[200];
  Float_t truth_jet_eta[200];
  Float_t truth_jet_phi[200];


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
  
  //2D Histos
  TH2F * Tjve = new TH2F("ETrueJet_vs_Eelectron", "E^{True}_{Jet} (|#eta^{Jet}|<0.7) vs. E_{e}^{True}",100,0,25,100,0,25);
  TH2F * Rjve = new TH2F("ERecoJet_vs_Eelectron", "E^{Reco}_{Jet} (|#eta^{Jet}|<0.7) vs. E_{e}^{True}",100,0,25,100,0,25);
  //Ratio Histos
  TH1F * eoTj = new TH1F("Eelectron_over_ETrueJet", "E_{e}^{True}/E^{True}_{Jet} (|#eta^{Jet}|<0.7)",30,0,3);
  TH1F * eoRj = new TH1F("Eelectron_over_ERecoJet", "E_{e}^{True}/E^{Reco}_{Jet} (|#eta^{Jet}|<0.7)",30,0,3);
  //Difference Histos
  TH1F * emTj = new TH1F("Eelectron_minus_ETrueJet", "E_{e}^{True} - E^{True}_{Jet} (|#eta^{Jet}|<0.7)",100,-20,30);
  TH1F * emRj = new TH1F("Eelectron_minus_ERecoJet", "E_{e}^{True} - E^{Reco}_{Jet} (|#eta^{Jet}|<0.7)",100,-20,30);
  //Detector Coordinate Histos
  TH1F * dPhiTj = new TH1F("dPhi_e_TrueJet", "|#Delta#varphi| (#varphi_{e} - #varphi^{True}_{Jet}) ", 32,0,M_PI);
  TH1F * dPhiRj = new TH1F("dPhi_e_RecoJet", "|#Delta#varphi| #varphi_{e} - #varphi(Jet^{Reco}_{Jet}) ", 32,0,M_PI);
  TH1F * dEtaTj = new TH1F("dEta_e_TrueJet", "|#Delta#eta| (#eta_{e} - #eta^{True}_{Jet})", 80,-10,10);
  TH1F * dEtaRj = new TH1F("dEta_e_RecoJet", "|#Delta#eta| (#eta_{e} - #eta^{Reco}_{Jet})", 80,-10,10);
  
  
  Long64_t nentries = _tree_event->GetEntries();
  for(Long64_t ie = 0; ie < nentries ; ie++){
    _tree_event->GetEntry(ie); //each entry is a 5GeV Electron
    fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ie, nentries);

    //cuts
    if(isnan(truthEta)) continue;
    Float_t True_DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(etruthPhi - truthPhi));
    if(nMatchedTrack < 3) continue;
    if (True_DeltaPhi < M_PI/2) continue;
    
    dPhiTj->Fill(True_DeltaPhi);
    Float_t Reco_DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(etruthPhi - phi));
    dPhiRj->Fill(Reco_DeltaPhi);

    dEtaTj->Fill(etruthEta-truthEta);
    dEtaRj->Fill(etruthEta-eta);

    
    Rjve->Fill(e,etruthE);
    Tjve->Fill(truthE,etruthE);

    eoTj->Fill(etruthE/truthE);
    eoRj->Fill(etruthE/e);

    emTj->Fill(etruthE-truthE);
    emRj->Fill(etruthE-e);

    // if (already_matched)
    //   Extra_Match_Count ++;
    // already_matched = true;
    
    // H_NExtra_Matches->Fill(Extra_Match_Count);
	
  } //entry loop



//Write to new root file
  TFile* fout = new TFile("Histograms_Jet_Callibration.root","RECREATE");

  Tjve->Write();
  Rjve->Write();

  eoTj->Write();
  eoRj->Write();

  emTj->Write();
  emRj->Write();
  
  dPhiTj->Write();
  dPhiRj->Write();
  dEtaTj->Write();
  dEtaRj->Write();

  // H_dR->Write();
  // H_NExtra_Matches->Write();

}
