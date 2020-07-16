#include "MyJetAnalysis.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <phool/getClass.h>

#include <g4eval/JetEvalStack.h>

#include <trackbase_historic/SvtxTrackMap.h>

#include <g4jets/JetMap.h>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace std;

MyJetAnalysis::MyJetAnalysis(const std::string& recojetname, const std::string& truthjetname, const std::string& outputfilename)
  : SubsysReco("MyJetAnalysis_" + recojetname + "_" + truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange(-8, 8)//-0.7 0.7 for sPHENIX track_reco jets
  , m_ptRange(0.5, 500)
  , m_jet_R(get_jet_radius_from_string(recojetname))
  , m_electronJetMatchingRadius(m_jet_R/2.0)
  , m_eEmin(1.0)
  , m_trackJetMatchingRadius(.7)
  , m_hInclusiveE(nullptr)
  , m_hInclusiveEta(nullptr)
  , m_hInclusivePhi(nullptr)
  , m_hInclusiveNJets(nullptr)
  , m_T(nullptr)
  , m_event(-1)
  , m_electron_truthEta(numeric_limits<float>::signaling_NaN())
  , m_electron_truthPhi(numeric_limits<float>::signaling_NaN())
  , m_electron_truthE(numeric_limits<float>::signaling_NaN())
  , m_electron_truthPt(numeric_limits<float>::signaling_NaN())
  , m_electron_truthpX(numeric_limits<float>::signaling_NaN())
  , m_electron_truthpY(numeric_limits<float>::signaling_NaN())
  , m_electron_truthpZ(numeric_limits<float>::signaling_NaN())
  , m_electron_truthPID(-1)
  , m_njets(-1)
  , m_nAlltruthjets(-1)
{
  m_id.fill(-1);
  m_nComponent.fill(-1);
  m_eta.fill(numeric_limits<float>::signaling_NaN());
  m_phi.fill(numeric_limits<float>::signaling_NaN());
  m_e.fill(numeric_limits<float>::signaling_NaN());
  m_pt.fill(numeric_limits<float>::signaling_NaN());

  m_matched_truthID.fill(-1);
  m_matched_truthNComponent.fill(-1);
  m_matched_truthEta.fill(numeric_limits<float>::signaling_NaN());
  m_matched_truthPhi.fill(numeric_limits<float>::signaling_NaN());
  m_matched_truthE.fill(numeric_limits<float>::signaling_NaN());
  m_matched_truthPt.fill(numeric_limits<float>::signaling_NaN());

  m_all_truthID.fill(-1);
  m_all_truthNComponent.fill(-1);
  m_all_truthEta.fill(numeric_limits<float>::signaling_NaN());
  m_all_truthPhi.fill(numeric_limits<float>::signaling_NaN());
  m_all_truthE.fill(numeric_limits<float>::signaling_NaN());
  m_all_truthPt.fill(numeric_limits<float>::signaling_NaN());

  std::fill_n(&m_matched_Constituent_truthPID[0][0],MaxNumJets*kMaxConstituents,NaN);
  std::fill_n(&m_matched_Constituent_truthCharge[0][0],MaxNumJets*kMaxConstituents,NaN);
  std::fill_n(&m_matched_Constituent_truthEta[0][0],MaxNumJets*kMaxConstituents,NaN);
  std::fill_n(&m_matched_Constituent_truthPhi[0][0],MaxNumJets*kMaxConstituents,NaN);
  std::fill_n(&m_matched_Constituent_truthPt[0][0],MaxNumJets*kMaxConstituents,NaN);
  std::fill_n(&m_matched_Constituent_truthE[0][0],MaxNumJets*kMaxConstituents,NaN);

  std::fill_n(&m_all_Constituent_truthPID[0][0],MaxNumJets*kMaxConstituents,NaN);
  std::fill_n(&m_all_Constituent_truthCharge[0][0],MaxNumJets*kMaxConstituents,NaN);
  std::fill_n(&m_all_Constituent_truthEta[0][0],MaxNumJets*kMaxConstituents,NaN);
  std::fill_n(&m_all_Constituent_truthPhi[0][0],MaxNumJets*kMaxConstituents,NaN);
  std::fill_n(&m_all_Constituent_truthPt[0][0],MaxNumJets*kMaxConstituents,NaN);
  std::fill_n(&m_all_Constituent_truthE[0][0],MaxNumJets*kMaxConstituents,NaN);
}

MyJetAnalysis::~MyJetAnalysis()
{
}

int MyJetAnalysis::Init(PHCompositeNode* topNode)
{
  if (Verbosity() >= MyJetAnalysis::VERBOSITY_SOME)
    cout << "MyJetAnalysis::Init - Output to " << m_outputFileName << endl;
  PHTFileServer::get().open(m_outputFileName, "RECREATE");

  // Histograms
  m_hInclusiveE = new TH1F("hInclusive_E",
    TString(m_recoJetName) + " inclusive jet E;Total jet energy (GeV)", 100, 0, 100);

  m_hInclusiveEta = new TH1F("hInclusive_eta",
    TString(m_recoJetName) + " inclusive jet #eta;#eta;Jet energy density", 50, -1, 1);

  m_hInclusivePhi = new TH1F("hInclusive_phi",
    TString(m_recoJetName) + " inclusive jet #phi;#phi;Jet energy density", 50, -M_PI, M_PI);

  m_hInclusiveNJets = new TH1F("hInclusive_njets",
    TString(m_recoJetName) + " inclusive number of jets",10,0,10);

  //Event Branches
  m_T = new TTree("T", "MyJetAnalysis Tree");
  m_T->Branch("event", &m_event, "event/I");
  m_T->Branch("njets", &m_njets, "njets/I");
  m_T->Branch("nAlltruthjets", &m_nAlltruthjets, "nAlltruthjets/I");
  //Root name (sans "m_") is in dynamic array size
  
  //Reconstructed Jet Branches
  m_T->Branch("id", m_id.data(), "id[njets]/I");
  m_T->Branch("nComponent", m_nComponent.data(), "nComponent[njets]/I");
  m_T->Branch("eta", m_eta.data(), "eta[njets]/F");
  m_T->Branch("phi", m_phi.data(), "phi[njets]/F");
  m_T->Branch("e", m_e.data(), "e[njets]/F");
  m_T->Branch("pt", m_pt.data(), "pt[njets]/F");

  //Matched Truth Jet Branches
  m_T->Branch("matched_truthID", m_matched_truthID.data(), "matched_truthID[njets]/I");
  m_T->Branch("matched_truthNComponent", m_matched_truthNComponent.data(), "matched_truthNComponent[njets]/I");
  m_T->Branch("matched_truthEta", m_matched_truthEta.data(), "matched_truthEta[njets]/F");
  m_T->Branch("matched_truthPhi", m_matched_truthPhi.data(), "matched_truthPhi[njets]/F");
  m_T->Branch("matched_truthE", m_matched_truthE.data(), "matched_truthE[njets]/F");
  m_T->Branch("matched_truthPt", m_matched_truthPt.data(), "matched_truthPt[njets]/F");
  m_T->Branch("matched_Constituent_truthPID", &m_matched_Constituent_truthPID,"matched_Constituent_truthPID[njets][100]/F");
  m_T->Branch("matched_Constituent_truthCharge", &m_matched_Constituent_truthCharge,"matched_Constituent_truthCharge[njets][100]/F");
  m_T->Branch("matched_Constituent_truthEta", &m_matched_Constituent_truthEta,"matched_Constituent_truthEta[njets][100]/F");
  m_T->Branch("matched_Constituent_truthPhi", &m_matched_Constituent_truthPhi,"matched_Constituent_truthPhi[njets][100]/F");
  m_T->Branch("matched_Constituent_truthPt", &m_matched_Constituent_truthPt,"matched_Constituent_truthPt[njets][100]/F");
  m_T->Branch("matched_Constituent_truthE", &m_matched_Constituent_truthE,"matched_Constituent_truthE[njets][100]/F");
  //100 = max number of constituents. ROOT does not support branches of 2D dynamic arrays nor with 2+ alpha-numeric dimensions
  
  //Track Constituents
  // m_T->Branch("nMatchedTrack", m_nMatchedTrack.data(), "nMatchedTrack/I");
  // m_T->Branch("TrackdR", m_trackdR.data(), "trackdR[nMatchedTrack]/F");
  // m_T->Branch("trackpT", m_trackpT.data(), "trackpT[nMatchedTrack]/F");

  // ALL Truth Jet Branches
  m_T->Branch("all_truthID", m_all_truthID.data(), "all_truthID[nAlltruthjets]/I");
  m_T->Branch("all_truthNComponent", m_all_truthNComponent.data(), "all_truthNComponent[nAlltruthjets]/I");
  m_T->Branch("all_truthEta", m_all_truthEta.data(), "all_truthEta[nAlltruthjets]/F");
  m_T->Branch("all_truthPhi", m_all_truthPhi.data(), "all_truthPhi[nAlltruthjets]/F");
  m_T->Branch("all_truthE", m_all_truthE.data(), "all_truthE[nAlltruthjets]/F");
  m_T->Branch("all_truthPt", m_all_truthPt.data(), "all_truthPt[nAlltruthjets]/F");

  m_T->Branch("all_Constituent_truthPID", &m_all_Constituent_truthPID,"all_Constituent_truthPID[nAlltruthjets][100]/F");
  m_T->Branch("all_Constituent_truthCharge", &m_all_Constituent_truthCharge,"all_Constituent_truthCharge[nAlltruthjets][100]/F");
  m_T->Branch("all_Constituent_truthEta", &m_all_Constituent_truthEta,"all_Constituent_truthEta[nAlltruthjets][100]/F");
  m_T->Branch("all_Constituent_truthPhi", &m_all_Constituent_truthPhi,"all_Constituent_truthPhi[nAlltruthjets][100]/F");
  m_T->Branch("all_Constituent_truthPt", &m_all_Constituent_truthPt,"all_Constituent_truthPt[nAlltruthjets][100]/F");
  m_T->Branch("all_Constituent_truthE", &m_all_Constituent_truthE,"all_Constituent_truthE[nAlltruthjets][100]/F");
  //100 = max number of constituents. ROOT does not support branches of 2D dynamic arrays nor with 2+ alpha-numeric dimensions
  
  //Electron Branches
  m_T->Branch("electron_truthEta", &m_electron_truthEta, "electron_truthEta/F");
  m_T->Branch("electron_truthPhi", &m_electron_truthPhi, "electron_truthPhi/F");
  m_T->Branch("electron_truthE", &m_electron_truthE, "electron_truthE/F");
  m_T->Branch("electron_truthPt", &m_electron_truthPt, "electron_truthPt/F");
  m_T->Branch("electron_truthpX", &m_electron_truthpX, "electron_truthpX/F");
  m_T->Branch("electron_truthpY", &m_electron_truthpY, "electron_truthpY/F");
  m_T->Branch("electron_truthpZ", &m_electron_truthpZ, "electron_truthpZ/F");
  m_T->Branch("electron_truthPID", &m_electron_truthPID, "electron_truthPID/I");
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::End(PHCompositeNode* topNode)
{
  cout << "MyJetAnalysis::End - Output to " << m_outputFileName << endl;
  PHTFileServer::get().cd(m_outputFileName);

  m_hInclusiveE->Write();
  m_hInclusiveEta->Write();
  m_hInclusivePhi->Write();
  m_hInclusiveNJets->Write();
  m_T->Write();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::InitRun(PHCompositeNode* topNode)
{
  m_jetEvalStack = shared_ptr<JetEvalStack>(new JetEvalStack(topNode, m_recoJetName, m_truthJetName));
  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() >= MyJetAnalysis::VERBOSITY_SOME)
    cout << "MyJetAnalysis::process_event() entered" << endl;

  m_jetEvalStack->next_event(topNode);
  JetRecoEval* recoeval = m_jetEvalStack->get_reco_eval();
  ++m_event;

  // interface to jets
  JetMap* jets = findNode::getClass<JetMap>(topNode, m_recoJetName);
  if (!jets)
  {
    cout
        << "MyJetAnalysis::process_event - Error can not find DST JetMap node "
        << m_recoJetName << endl;
    exit(-1);
  }

  JetMap* all_truth_jets = findNode::getClass<JetMap>(topNode, m_truthJetName);
  if (!all_truth_jets)
  {
    cout
        << "MyJetAnalysis::process_event - Error can not find DST JetMap node "
        << m_truthJetName << endl;
    exit(-1);
  }
  
  // // interface to tracks
  // SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  // if (!trackmap)
  // {
  //   cout
  //       << "MyJetAnalysis::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
  //   exit(-1);
  // }

  //Interface to True Electrons
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if ( !truthinfo )
    {
      cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
      exit(-1);
    }
  
  float hardest_electron_E = 0;
  int hardest_electron_int = -1;
  
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  for ( PHG4TruthInfoContainer::ConstIterator eter = range.first; eter != range.second; ++eter )
  {
    PHG4Particle* g4particle = eter->second;
    int particleID = g4particle->get_pid();

    bool iselectron = fabs(particleID) == 11;
    if (!iselectron) continue;

    float temp_electron_truthE = g4particle->get_e();
    if (temp_electron_truthE < m_eEmin) continue;
    if (temp_electron_truthE < hardest_electron_E) continue;
    
    hardest_electron_E = temp_electron_truthE;
    hardest_electron_int = eter->first;
  }

  for ( PHG4TruthInfoContainer::ConstIterator eter = range.first; eter != range.second; ++eter )
  {
    int e = eter->first;
    if (not(e == hardest_electron_int)) continue;
    PHG4Particle* g4particle = eter->second;
      
    m_electron_truthE = g4particle->get_e();
    m_electron_truthpX = g4particle->get_px();
    m_electron_truthpY = g4particle->get_py();
    m_electron_truthpZ = g4particle->get_pz();
    //Scalars need not be manually reset
    
    TLorentzVector e_vec;
    e_vec.SetPxPyPzE(m_electron_truthpX,m_electron_truthpY,m_electron_truthpZ,m_electron_truthE);
    m_electron_truthEta = e_vec.Eta();
    m_electron_truthPhi = e_vec.Phi();
    m_electron_truthPt = e_vec.Pt();

    //Remove Previous Event Data
    {  
      m_njets=0;
      m_id.fill(-1);
      m_nComponent.fill(-1);
      m_eta.fill(numeric_limits<float>::signaling_NaN());
      m_phi.fill(numeric_limits<float>::signaling_NaN());
      m_e.fill(numeric_limits<float>::signaling_NaN());
      m_pt.fill(numeric_limits<float>::signaling_NaN());

      m_matched_truthID.fill(-1);
      m_matched_truthNComponent.fill(-1);
      m_matched_truthEta.fill(numeric_limits<float>::signaling_NaN());
      m_matched_truthPhi.fill(numeric_limits<float>::signaling_NaN());
      m_matched_truthE.fill(numeric_limits<float>::signaling_NaN());
      m_matched_truthPt.fill(numeric_limits<float>::signaling_NaN());

      std::fill_n(&m_matched_Constituent_truthPID[0][0],MaxNumJets*kMaxConstituents,NaN);
      std::fill_n(&m_matched_Constituent_truthCharge[0][0],MaxNumJets*kMaxConstituents,NaN);
      std::fill_n(&m_matched_Constituent_truthEta[0][0],MaxNumJets*kMaxConstituents,NaN);
      std::fill_n(&m_matched_Constituent_truthPhi[0][0],MaxNumJets*kMaxConstituents,NaN);
      std::fill_n(&m_matched_Constituent_truthPt[0][0],MaxNumJets*kMaxConstituents,NaN);
      std::fill_n(&m_matched_Constituent_truthE[0][0],MaxNumJets*kMaxConstituents,NaN);
      
      m_nAlltruthjets = 0;
      m_all_truthID.fill(-1);
      m_all_truthNComponent.fill(-1);
      m_all_truthPhi.fill(numeric_limits<float>::signaling_NaN());
      m_all_truthE.fill(numeric_limits<float>::signaling_NaN());
      m_all_truthPt.fill(numeric_limits<float>::signaling_NaN());

      std::fill_n(&m_all_Constituent_truthPID[0][0],MaxNumJets*kMaxConstituents,NaN);
      std::fill_n(&m_all_Constituent_truthCharge[0][0],MaxNumJets*kMaxConstituents,NaN);
      std::fill_n(&m_all_Constituent_truthEta[0][0],MaxNumJets*kMaxConstituents,NaN);
      std::fill_n(&m_all_Constituent_truthPhi[0][0],MaxNumJets*kMaxConstituents,NaN);
      std::fill_n(&m_all_Constituent_truthPt[0][0],MaxNumJets*kMaxConstituents,NaN);
      std::fill_n(&m_all_Constituent_truthE[0][0],MaxNumJets*kMaxConstituents,NaN);
    }

    int inc_jet_counter = 0;
    int j = 0; //Jet element index. Same index for reco and matched truth, but in separate arrays
    for (JetMap::Iter jter = jets->begin(); jter != jets->end(); ++jter)
      {
	Jet* jet = jter->second;
	assert(jet); //checks if not null. Aborts if null pointer.

	// fill inclusive histograms
	assert(m_hInclusiveE);
	m_hInclusiveE->Fill(jet->get_e());
	assert(m_hInclusiveEta);
	m_hInclusiveEta->Fill(jet->get_eta());
	assert(m_hInclusivePhi);
	m_hInclusivePhi->Fill(jet->get_phi());
	++inc_jet_counter;

	TLorentzVector JReco_vec(jet->get_px(), jet->get_py(), jet->get_pz(),jet->get_e());
	
 	//Apply cuts
 	bool eta_cut = (jet->get_eta() >= m_etaRange.first) and (jet->get_eta() <= m_etaRange.second); 
	bool pt_cut = (jet->get_pt() >= m_ptRange.first) and (jet->get_pt() <= m_ptRange.second);
	bool electron_cut = (JReco_vec.DeltaR(e_vec) > m_electronJetMatchingRadius);
	
 	if ((not eta_cut) or (not pt_cut) or (not electron_cut))
	  {
	    if (Verbosity() >= MyJetAnalysis::VERBOSITY_MORE)
	      {
		cout << "MyJetAnalysis::process_event() - jet failed acceptance cut: ";
		cout << "eta cut: " << eta_cut << ", ptcut: " << pt_cut << endl;
		cout << "jet eta: " << jet->get_eta() << ", jet pt: " << jet->get_pt() << endl;
		cout << "electron dR cut: " << m_electronJetMatchingRadius << ", electron-jet dR: "
		     << JReco_vec.DeltaR(e_vec) << endl;
		jet->identify();
	      }
	    continue;
	  }
  
	//Jet Arrays
	m_id[j] = jet->get_id();
	m_nComponent[j] = jet->size_comp();
	m_e[j] = jet->get_e();
	m_eta[j] = jet->get_eta();
	m_phi[j] = jet->get_phi();
	m_pt[j] = jet->get_pt();
	
	//Which unique truth jet contributed the most enery to this reco jet?
	Jet* truthjet = recoeval->unique_truth_jet_from_reco(jet);
	if (truthjet)
	      {
		TLorentzVector JTruth_vec(truthjet->get_px(), truthjet->get_py(), truthjet->get_pz(),truthjet->get_e());
		
		//Apply cuts
		bool truth_eta_cut = (truthjet->get_eta() >= m_etaRange.first) and (truthjet->get_eta() <= m_etaRange.second); 
		bool truth_pt_cut = (truthjet->get_pt() >= m_ptRange.first) and (truthjet->get_pt() <= m_ptRange.second);
		bool truth_electron_cut = (JTruth_vec.DeltaR(e_vec) > m_electronJetMatchingRadius);
		if ((not truth_eta_cut) or (not truth_pt_cut) or (not truth_electron_cut))
		  continue;

		m_matched_truthID[j] = truthjet->get_id();
		m_matched_truthNComponent[j] = truthjet->size_comp();
		m_matched_truthEta[j] = truthjet->get_eta();
		m_matched_truthPhi[j] = truthjet->get_phi();
		m_matched_truthE[j] = truthjet->get_e();
		m_matched_truthPt[j] = truthjet->get_pt();

		
		cout<<"EVENT NUMBER "<<m_event<<", Jet Number = "<<j
		    <<", # Constituents = "<<m_matched_truthNComponent[j]<<std::endl;
		  		
		std::set<PHG4Particle*> truthj_particle_set = m_jetEvalStack->get_truth_eval()->all_truth_particles(truthjet);
		int c_index = 0;
		for (auto i_t:truthj_particle_set)
		  {
		    TLorentzVector constituent_vec;
		    constituent_vec.SetPxPyPzE(i_t->get_px(),i_t->get_py(),i_t->get_pz(),i_t->get_e());
		    TParticlePDG * pdg_p = TDatabasePDG::Instance()->GetParticle(i_t->get_pid());
		    m_matched_Constituent_truthPID[j][c_index] = i_t->get_pid();
		    m_matched_Constituent_truthCharge[j][c_index] = (pdg_p->Charge()/3.);
		    m_matched_Constituent_truthEta[j][c_index] = constituent_vec.Eta();		  
		    m_matched_Constituent_truthPhi[j][c_index] = constituent_vec.Phi();
		    m_matched_Constituent_truthPt[j][c_index] = constituent_vec.Pt();
		    m_matched_Constituent_truthE[j][c_index] = constituent_vec.E();
		    c_index++;
		  }
	      }
	++j;
	m_njets=j;
	//j is incremented outside of the matching criteria,
	//s.t. matched_truth is filled NaN if no match.
	
	if (j >= MaxNumJets) break;

	////// TRACK-JET MATCHING //////

	// //fill trees - jet track matching
	// m_nMatchedTrack[j] = 0;
	
	// Float_t jet_eta = jet->get_eta();
	// Float_t jet_phi = jet->get_phi();
	// for (SvtxTrackMap::Iter iter = trackmap->begin();
	//      iter != trackmap->end(); ++iter)

	//   {
	//     SvtxTrack* track = iter->second;

	//     TVector3 v(track->get_px(), track->get_py(), track->get_pz());
	//     const double dEta = v.Eta() - jet_eta;
	//     const double dPhi = v.Phi() - jet_phi;
	//     const double dR = sqrt(dEta * dEta + dPhi * dPhi);

	//     if (dR < m_trackJetMatchingRadius)
	//       {
	// 	//matched track to jet
	// 	assert(m_nMatchedTrack[j] < kMaxMatchedTrack);
	// 	m_trackdR[j][m_nMatchedTrack[j]] = dR;
	// 	m_trackpT[j][m_nMatchedTrack[j]] = v.Perp();
	// 	++m_nMatchedTrack[j];
		
	//       }
	//     if (m_nMatchedTrack[j] >= kMaxMatchedTrack)
	//       {
	// 	cout << "MyJetAnalysis::process_event() - reached max track that matching a jet. Quit iterating tracks" << endl;
	// 	break;
	//       }
	    
	//   }  //    for (SvtxTrackMap::Iter iter = trackmap->begin();
	
      }  //   for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)      

    m_hInclusiveNJets->Fill(inc_jet_counter);
    assert(m_hInclusiveNJets);

    //All Truth Jet Loop    
    int i_alltruth = 0;
    for (JetMap::Iter tter = all_truth_jets->begin(); tter != all_truth_jets->end(); ++tter)
      {
	Jet* all_truthjet = tter->second;
	assert(all_truthjet); //Check if null pointer.

 	TLorentzVector JAllTruth_vec(all_truthjet->get_px(), all_truthjet->get_py(),
				  all_truthjet->get_pz(),all_truthjet->get_e());

     	//Apply cuts
 	bool eta_cut = (all_truthjet->get_eta() >= m_etaRange.first) and (all_truthjet->get_eta() <= m_etaRange.second); 
	bool pt_cut = (all_truthjet->get_pt() >= m_ptRange.first) and (all_truthjet->get_pt() <= m_ptRange.second);
	bool electron_cut = (JAllTruth_vec.DeltaR(e_vec) > m_electronJetMatchingRadius);
	
 	if ((not eta_cut) or (not pt_cut) or (not electron_cut))
	  {
	    if (Verbosity() >= MyJetAnalysis::VERBOSITY_MORE)
	      {
		cout << "MyJetAnalysis::process_event() - jet failed acceptance cut: ";
		cout << "eta cut: " << eta_cut << ", ptcut: " << pt_cut << endl;
		cout << "jet eta: " << all_truthjet->get_eta() << ", jet pt: " << all_truthjet->get_pt() << endl;
		cout << "electron dR cut: " << m_electronJetMatchingRadius << ", electron-jet dR: "
		     << JAllTruth_vec.DeltaR(e_vec) << endl;
		all_truthjet->identify();
	      }
	    continue;
	  }
	
	m_all_truthID[i_alltruth] = all_truthjet->get_id();
	m_all_truthNComponent[i_alltruth] = all_truthjet->size_comp();
	m_all_truthEta[i_alltruth] = all_truthjet->get_eta();
	m_all_truthPhi[i_alltruth] = all_truthjet->get_phi();
	m_all_truthE[i_alltruth] = all_truthjet->get_e();
	m_all_truthPt[i_alltruth] = all_truthjet->get_pt();

	std::set<PHG4Particle*> all_truthj_particle_set = m_jetEvalStack->get_truth_eval()->all_truth_particles(all_truthjet);
	int c_Allindex = 0;
	for (auto i_At:all_truthj_particle_set)
	  {
	    TLorentzVector constituent_vec;
	    constituent_vec.SetPxPyPzE(i_At->get_px(),i_At->get_py(),i_At->get_pz(),i_At->get_e());
	    TParticlePDG * pdg_p = TDatabasePDG::Instance()->GetParticle(i_At->get_pid());
	    m_all_Constituent_truthPID[i_alltruth][c_Allindex] = i_At->get_pid();
	    m_all_Constituent_truthCharge[i_alltruth][c_Allindex] = (pdg_p->Charge()/3.);
	    m_all_Constituent_truthEta[i_alltruth][c_Allindex] = constituent_vec.Eta();		  
	    m_all_Constituent_truthPhi[i_alltruth][c_Allindex] = constituent_vec.Phi();
	    m_all_Constituent_truthPt[i_alltruth][c_Allindex] = constituent_vec.Pt();
	    m_all_Constituent_truthE[i_alltruth][c_Allindex] = constituent_vec.E();
	    cout<<"All Constituent Number = "<<c_Allindex;
	    cout<<", PID = "<<m_all_Constituent_truthPID[i_alltruth][c_Allindex];
	    cout<<", Parent ID = "<<i_At->get_parent_id()<<endl;
	    cout<<", C = "<<m_all_Constituent_truthCharge[i_alltruth][c_Allindex];
	    cout<<", Eta = "<<m_all_Constituent_truthEta[i_alltruth][c_Allindex];
	    cout<<", Phi = "<<m_all_Constituent_truthPhi[i_alltruth][c_Allindex];
	    cout<<", Pt = "<<m_all_Constituent_truthPt[i_alltruth][c_Allindex];
	    cout<<", E = "<<m_all_Constituent_truthE[i_alltruth][c_Allindex]<<std::endl;
	    c_Allindex++;
	  }	
	++i_alltruth;
	m_nAlltruthjets = i_alltruth;
      }
    m_T->Fill(); //Fill Tree inside electron Loop, after reco&truth loops
  }//electron Loop  
  return Fun4AllReturnCodes::EVENT_OK;
}
