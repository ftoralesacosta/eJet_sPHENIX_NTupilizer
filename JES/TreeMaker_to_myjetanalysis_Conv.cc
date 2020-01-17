#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace std;


int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout<<"Syntax: [Command] [Old Input File] [Output File Name]"<<std::endl;
    exit(EXIT_FAILURE);
  }


  //Open Old Input File
  TString root_file = (TString)argv[1];
  std::cout << "Opening: " << (TString)argv[1] << std::endl;
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
  
  //OLD MonteCarlo
  int jet4_n;
  Float_t jet4_pt[200];
  Float_t jet4_E[200];
  Float_t jet4_m[200];
  Float_t jet4_eta[200];
  Float_t jet4_phi[200];

  int truthjet4_n;
  Float_t truthjet4_pt[200];
  Float_t truthjet4_m[200];
  Float_t truthjet4_eta[200];
  Float_t truthjet4_phi[200];

  _tree_event->SetBranchAddress("jet4_n", &jet4_n);
  _tree_event->SetBranchAddress("jet4_pt", jet4_pt);
  _tree_event->SetBranchAddress("jet4_E", jet4_E);
  _tree_event->SetBranchAddress("jet4_m", jet4_m);
  _tree_event->SetBranchAddress("jet4_eta", jet4_eta);
  _tree_event->SetBranchAddress("jet4_phi", jet4_phi);
  
  _tree_event->SetBranchAddress("truthjet4_n", &truthjet4_n);
  _tree_event->SetBranchAddress("truthjet4_pt", truthjet4_pt);
  _tree_event->SetBranchAddress("truthjet4_m", truthjet4_m);
  _tree_event->SetBranchAddress("truthjet4_eta", truthjet4_eta);
  _tree_event->SetBranchAddress("truthjet4_phi", truthjet4_phi);
  
  // NEW TFILE
  TString out_root_file = (TString)argv[2];
  std::cout << "Opening: " << (TString)argv[2] << std::endl;
  //TFile *out_file = new TFile::Open(out_root_file);
  TFile *out_file = new TFile(out_root_file,"RECREATE");
  if (out_file == NULL) {
    std::cout << " File Fail" << std::endl;
    exit(EXIT_FAILURE);
  }  

  //! Output Tree variables
  TTree *m_T = new TTree("T", "MyJetAnalysis Tree");
  int m_event;
  int m_id;
  int m_nComponent;
  float m_eta;
  float m_phi;
  float m_e;
  float m_pt;
  int m_truthID;
  int m_truthNComponent;
  float m_truthEta;
  float m_truthPhi;
  float m_truthE;
  float m_truthPt;
  //! number of matched tracks
  int m_nMatchedTrack;
  enum
  {
    //! max number of tracks
    kMaxMatchedTrack = 1000
  };
  std::array<float, kMaxMatchedTrack> m_trackdR;
  std::array<float, kMaxMatchedTrack> m_trackpT;

  // Histograms
  const char* m_recoJetName = "AntiKt_Tower_r04";
  TH1*  m_hInclusiveE = new TH1F(
				 "hInclusive_E",  //
      TString(m_recoJetName) + " inclusive jet E;Total jet energy (GeV)", 100, 0, 100);
 
  TH1*  m_hInclusiveEta =
      new TH1F(
          "hInclusive_eta",  //
          TString(m_recoJetName) + " inclusive jet #eta;#eta;Jet energy density", 50, -1, 1);
  TH1* m_hInclusivePhi =
      new TH1F(
          "hInclusive_phi",  //
          TString(m_recoJetName) + " inclusive jet #phi;#phi;Jet energy density", 50, -M_PI, M_PI);

  //New Trees
  m_T->Branch("m_event", &m_event, "event/I");
  m_T->Branch("id", &m_id, "id/I");
  m_T->Branch("nComponent", &m_nComponent, "nComponent/I");
  m_T->Branch("eta", &m_eta, "eta/F");
  m_T->Branch("phi", &m_phi, "phi/F");
  m_T->Branch("e", &m_e, "e/F");
  m_T->Branch("pt", &m_pt, "pt/F");
  m_T->Branch("truthID", &m_truthID, "truthID/I");
  m_T->Branch("truthNComponent", &m_truthNComponent, "truthNComponent/I");
  m_T->Branch("truthEta", &m_truthEta, "truthEta/F");
  m_T->Branch("truthPhi", &m_truthPhi, "truthPhi/F");
  m_T->Branch("truthE", &m_truthE, "truthE/F");
  m_T->Branch("truthPt", &m_truthPt, "truthPt/F");
  m_T->Branch("nMatchedTrack", &m_nMatchedTrack, "nMatchedTrack/I");
  //  m_T->Branch("id", m_trackdR.data(), "trackdR[nMatchedTrack]/F");
  //  m_T->Branch("id", m_trackpT.data(), "trackpT[nMatchedTrack]/F");


//loop through Dennis Events
  Long64_t nentries = _tree_event->GetEntries();
  for(Long64_t ievent = 0; ievent < nentries ; ievent++){
    _tree_event->GetEntry(ievent);
    fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);
    for (Long64_t i_reco4 = 0; i_reco4 < jet4_n; i_reco4++){
      ++m_event;
      TVector3 v3_reco;
      v3_reco.SetPtEtaPhi(jet4_pt[i_reco4],jet4_eta[i_reco4],jet4_phi[i_reco4]);
      
	// fill histograms
	assert(m_hInclusiveE);
	m_hInclusiveE->Fill(jet4_E[i_reco4]);
	assert(m_hInclusiveEta);
	m_hInclusiveEta->Fill(v3_reco.Eta());
	assert(m_hInclusivePhi);
	m_hInclusivePhi->Fill(v3_reco.Phi());

	// fill trees - jet spectrum
	m_id = 0;
	m_nComponent = 0;
	m_eta = v3_reco.Eta();
	m_phi = v3_reco.Phi();
	m_e = jet4_E[i_reco4];
	m_pt = v3_reco.Pt();	

	bool do_printing = false;
	if (do_printing){
	  cout<<"m_event = "<<m_event<<endl;
	  cout<<"m_id = "<<m_id<<endl;
	  cout<<"m_nComponent = "<<m_nComponent<<endl;
	  cout<<"m_eta = "<<m_eta<<endl;
	  cout<<"m_phi = "<<m_phi<<endl;
	  cout<<"m_e = "<<m_e<<endl;
	  cout<<"m_pt = "<<m_pt<<endl;
	}

	for (Long64_t i_truth4 = 0; i_truth4 < truthjet4_n; i_truth4++){
	    TVector3 v3_truth;
	    v3_truth.SetPtEtaPhi(truthjet4_pt[i_truth4],truthjet4_eta[i_truth4],truthjet4_phi[i_truth4]);

	    bool already_matched = false;
	    if (already_matched) continue;
	    
	    m_truthID = NAN;
	    m_truthNComponent = NAN;
	    m_truthEta = NAN;
	    m_truthPhi = NAN;
	    m_truthE = NAN;
	    m_truthPt = NAN;

	    Float_t dR = v3_truth.DeltaR(v3_reco);
	    bool truthjet =  (dR <= 0.2);
	    
	    m_trackdR[i_truth4] = dR;
	    m_trackpT[i_truth4] = truthjet4_pt[i_truth4];
	    ++m_nMatchedTrack;
	    if (truthjet)
	      {
		already_matched = true;
		m_truthID = 1;
		m_truthNComponent = 2;
		m_truthEta = v3_truth.Eta();
		m_truthPhi = v3_truth.Phi();
		m_truthE = jet4_E[i_reco4];//Old files don't have truth E
		m_truthPt = v3_truth.Pt();

		if(do_printing){
		  cout<<"m_TruthID = "<<m_truthID<<endl;
		  cout<<"m_truthNComponent = "<<m_truthNComponent<<endl;
		  cout<<"m_truthEta = "<<m_truthEta<<endl;
		  cout<<"m_truthPhi = "<<m_truthPhi<<endl;
		  cout<<"m_truthE = "<<m_truthE<<endl;
		  cout<<"m_truthPt = "<<m_truthPt<<endl;
		  cout<<"m_nMatchedTrack"<<m_nMatchedTrack<<endl;
		}

	      }//if truth dR Matched
	} //truthjet loop
	m_T->Fill();
    } //reco loop
  } //event loop
  m_hInclusiveE->Write();
  m_hInclusiveEta->Write();
  m_hInclusivePhi->Write();
  m_T->Write();

  file->Close();
  out_file->Close();
  
  return 0;
  }


