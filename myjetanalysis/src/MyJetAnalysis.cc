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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace std;

MyJetAnalysis::MyJetAnalysis(const std::string& recojetname, const std::string& truthjetname, const std::string& outputfilename)
  : SubsysReco("MyJetAnalysis_" + recojetname + "_" + truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange(-5, 5)//-0.7 0.7 for track_reco jets
  , m_ptRange(0.5, 500)
  , m_eEmin(1.0)
  , m_trackJetMatchingRadius(.7)
  , m_hInclusiveE(nullptr)
  , m_hInclusiveEta(nullptr)
  , m_hInclusivePhi(nullptr)
  , m_hInclusiveNJets(nullptr)
  , m_T(nullptr)
  , m_event(-1)
  , m_id(-1)
  , m_nComponent(-1)
  , m_eta(numeric_limits<float>::signaling_NaN())
  , m_phi(numeric_limits<float>::signaling_NaN())
  , m_e(numeric_limits<float>::signaling_NaN())
  , m_pt(numeric_limits<float>::signaling_NaN())
  , m_subleading_id(-1)
  , m_subleading_nComponent(-1)
  , m_subleading_eta(numeric_limits<float>::signaling_NaN())
  , m_subleading_phi(numeric_limits<float>::signaling_NaN())
  , m_subleading_e(numeric_limits<float>::signaling_NaN())
  , m_subleading_pt(numeric_limits<float>::signaling_NaN())
  , m_truthID(-1)
  , m_truthNComponent(-1)
  , m_truthEta(numeric_limits<float>::signaling_NaN())
  , m_truthPhi(numeric_limits<float>::signaling_NaN())
  , m_truthE(numeric_limits<float>::signaling_NaN())
  , m_truthPt(numeric_limits<float>::signaling_NaN())
  , m_subleading_truthID(-1)
  , m_subleading_truthNComponent(-1)
  , m_subleading_truthEta(numeric_limits<float>::signaling_NaN())
  , m_subleading_truthPhi(numeric_limits<float>::signaling_NaN())
  , m_subleading_truthE(numeric_limits<float>::signaling_NaN())
  , m_subleading_truthPt(numeric_limits<float>::signaling_NaN())
  , m_nMatchedTrack(-1)
  , m_subleading_nMatchedTrack(-1)
  , m_etruthEta(numeric_limits<float>::signaling_NaN())
  , m_etruthPhi(numeric_limits<float>::signaling_NaN())
  , m_etruthE(numeric_limits<float>::signaling_NaN())
  , m_etruthPt(numeric_limits<float>::signaling_NaN())
  , m_etruthpX(numeric_limits<float>::signaling_NaN())
  , m_etruthpY(numeric_limits<float>::signaling_NaN())
  , m_etruthpZ(numeric_limits<float>::signaling_NaN())
  , m_etruthPID(-1)
  , m_etruthParentID(-1)
{
  m_trackdR.fill(numeric_limits<float>::signaling_NaN());
  m_trackpT.fill(numeric_limits<float>::signaling_NaN());
  m_subleading_trackdR.fill(numeric_limits<float>::signaling_NaN());
  m_subleading_trackpT.fill(numeric_limits<float>::signaling_NaN());

}

MyJetAnalysis::~MyJetAnalysis()
{
}

int MyJetAnalysis::Init(PHCompositeNode* topNode)
{
  if (Verbosity() >= MyJetAnalysis::VERBOSITY_SOME)
    cout << "MyJetAnalysis::Init - Outoput to " << m_outputFileName << endl;

  PHTFileServer::get().open(m_outputFileName, "RECREATE");

  // Histograms
  m_hInclusiveE = new TH1F(
      "hInclusive_E",  //
      TString(m_recoJetName) + " inclusive jet E;Total jet energy (GeV)", 100, 0, 100);

  m_hInclusiveEta =
      new TH1F(
          "hInclusive_eta",  //
          TString(m_recoJetName) + " inclusive jet #eta;#eta;Jet energy density", 50, -1, 1);
  m_hInclusivePhi =
      new TH1F(
          "hInclusive_phi",  //
          TString(m_recoJetName) + " inclusive jet #phi;#phi;Jet energy density", 50, -M_PI, M_PI);

  m_hInclusiveNJets =
      new TH1F(
	   "hInclusive_njets",  //
	   TString(m_recoJetName) + " inclusive number of jets",10,0,10);
	       
  
  //Trees
  m_T = new TTree("T", "MyJetAnalysis Tree");

  //Reconstructed Branches
  m_T->Branch("m_event", &m_event, "event/I");
  m_T->Branch("id", &m_id, "id/I");
  m_T->Branch("nComponent", &m_nComponent, "nComponent/I");
  m_T->Branch("eta", &m_eta, "eta/F");
  m_T->Branch("phi", &m_phi, "phi/F");
  m_T->Branch("e", &m_e, "e/F");
  m_T->Branch("pt", &m_pt, "pt/F");

  m_T->Branch("subleading_id", &m_subleading_id, "subleading_id/I");
  m_T->Branch("subleading_nComponent", &m_subleading_nComponent, "subleading_nComponent/I");
  m_T->Branch("subleading_eta", &m_subleading_eta, "subleading_eta/F");
  m_T->Branch("subleading_phi", &m_subleading_phi, "subleading_phi/F");
  m_T->Branch("subleading_e", &m_subleading_e, "subleading_e/F");
  m_T->Branch("subleading_pt", &m_subleading_pt, "subleading_pt/F");

  //Truth Branches
  m_T->Branch("truthID", &m_truthID, "truthID/I");
  m_T->Branch("truthNComponent", &m_truthNComponent, "truthNComponent/I");
  m_T->Branch("truthEta", &m_truthEta, "truthEta/F");
  m_T->Branch("truthPhi", &m_truthPhi, "truthPhi/F");
  m_T->Branch("truthE", &m_truthE, "truthE/F");
  m_T->Branch("truthPt", &m_truthPt, "truthPt/F");
  m_T->Branch("nMatchedTrack", &m_nMatchedTrack, "nMatchedTrack/I");
  m_T->Branch("TrackdR", m_trackdR.data(), "trackdR[nMatchedTrack]/F");
  m_T->Branch("trackpT", m_trackpT.data(), "trackpT[nMatchedTrack]/F");

  m_T->Branch("subleading_truthID", &m_subleading_truthID, "subleading_truthID/I");
  m_T->Branch("subleading_truthNComponent", &m_subleading_truthNComponent, "subleading_truthNComponent/I");
  m_T->Branch("subleading_truthEta", &m_subleading_truthEta, "subleading_truthEta/F");
  m_T->Branch("subleading_truthPhi", &m_subleading_truthPhi, "subleading_truthPhi/F");
  m_T->Branch("subleading_truthE", &m_subleading_truthE, "subleading_truthE/F");
  m_T->Branch("subleading_truthPt", &m_subleading_truthPt, "subleading_truthPt/F");
  m_T->Branch("subleading_nMatchedTrack", &m_subleading_nMatchedTrack, "subleading_nMatchedTrack/I");
  m_T->Branch("subleading_TrackdR", m_subleading_trackdR.data(), "subleading_trackdR[nMatchedTrack]/F");
  m_T->Branch("subleading_trackpT", m_subleading_trackpT.data(), "subleading_trackpT[nMatchedTrack]/F");

  //Electron Branches
  m_T->Branch("etruthEta", &m_etruthEta, "etruthEta/F");
  m_T->Branch("etruthPhi", &m_etruthPhi, "etruthPhi/F");
  m_T->Branch("etruthE", &m_etruthE, "etruthE/F");
  m_T->Branch("etruthPt", &m_etruthPt, "etruthPt/F");
  m_T->Branch("etruthpX", &m_etruthpX, "etruthpX/F");
  m_T->Branch("etruthpY", &m_etruthpY, "etruthpY/F");
  m_T->Branch("etruthpZ", &m_etruthpZ, "etruthpZ/F");
  m_T->Branch("etruthPt", &m_etruthPt, "etruthPt/F");
  m_T->Branch("etruthPID", &m_etruthPID, "etruthPID/I");
  m_T->Branch("etruthParentID", &m_etruthParentID, "etruthParentID/I");
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::End(PHCompositeNode* topNode)
{
  cout << "MyJetAnalysis::End - Outoput to " << m_outputFileName << endl;
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

  // interface to tracks
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
  {
    cout
        << "MyJetAnalysis::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
    exit(-1);
  }


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

    //int parent_ID = g4particle->get_parent_id(); 
    //bool iselectron = (fabs(particleID) == 11 && fabs(parent_ID) == 11 );
    //Cut on parent electrons. Nope. PID=NaN if hardscattered e
    
    bool iselectron = fabs(particleID) == 11;
    if (!iselectron) continue;

    m_etruthE = g4particle->get_e();
    if (m_etruthE < m_eEmin) continue;

    if (m_etruthE < hardest_electron_E) continue;
    hardest_electron_E = m_etruthE;
    hardest_electron_int = eter->first;
  }


  for ( PHG4TruthInfoContainer::ConstIterator eter = range.first; eter != range.second; ++eter )
  {
    if (not(eter->first == hardest_electron_int)) continue;
    PHG4Particle* g4particle = eter->second;

    m_etruthE = g4particle->get_e();
    m_etruthpX = g4particle->get_px();
    m_etruthpY = g4particle->get_py();
    m_etruthpZ = g4particle->get_pz();
	
    TLorentzVector e;
    e.SetPxPyPzE(m_etruthpX,m_etruthpY,m_etruthpZ,m_etruthE);
    m_etruthEta = e.Eta();
    m_etruthPhi = e.Phi();
    m_etruthPt = e.Pt();
    
    //Outer Jet Variables
    int jetcounter = 0;
    unsigned int hardest_jet_int = 0;
    float hardest_jet_energy = 0;
    unsigned int subleading_jet_int = 0;
    float subleading_jet_energy = 0;
    
    for (JetMap::Iter jter = jets->begin(); jter != jets->end(); ++jter)
      {
	Jet* jet = jter->second;
	assert(jet);
	++jetcounter;
	bool eta_cut = (jet->get_eta() >= m_etaRange.first) and (jet->get_eta() <= m_etaRange.second); 
	bool pt_cut = (jet->get_pt() >= m_ptRange.first) and (jet->get_pt() <= m_ptRange.second);

	if ((not eta_cut) or (not pt_cut))
	  {
	    if (Verbosity() >= MyJetAnalysis::VERBOSITY_MORE)
	      {
		cout << "MyJetAnalysis::process_event() - jet failed acceptance cut: ";
		cout << "eta cut: " << eta_cut << ", ptcut: " << pt_cut << endl;
		cout << "jet eta: " << jet->get_eta() << ", jet pt: " << jet->get_pt() << endl;
		jet->identify();
	      }
	    continue;
	  }

	// fill inclusive histograms
	assert(m_hInclusiveE);
	m_hInclusiveE->Fill(jet->get_e());
	assert(m_hInclusiveEta);
	m_hInclusiveEta->Fill(jet->get_eta());
	assert(m_hInclusivePhi);
	m_hInclusivePhi->Fill(jet->get_phi());

	//determine hardest and subleading jets
	float jet_e = jet->get_e();
	if (jet_e < subleading_jet_energy) continue;

	if ((jet_e >= subleading_jet_energy) && (jet_e < hardest_jet_energy))
	  {
	    subleading_jet_energy = jet_e;
	    subleading_jet_int = jter->first;
	  }

	if (jet_e >= hardest_jet_energy)
	  {
	    hardest_jet_energy = jet_e;
	    hardest_jet_int = jter->first;
	  }	
      }
    assert(m_hInclusiveNJets);
    m_hInclusiveNJets->Fill(jetcounter);
    
    for (JetMap::Iter jter = jets->begin(); jter != jets->end(); ++jter)
      {

	bool hardest_jet = (jter->first == hardest_jet_int);
	bool subleading_jet = (jter->first == subleading_jet_int);
	if (not(hardest_jet) || (subleading_jet)) continue;


	Jet* jet = jter->second;
	assert(jet);

	//LEADING JET
	if (hardest_jet)
	  {
	    m_id = jet->get_id();
	    m_nComponent = jet->size_comp();
	    m_e = jet->get_e();
	    m_eta = jet->get_eta();
	    m_phi = jet->get_phi();
	    m_pt = jet->get_pt();
	
	    m_truthID = NAN;
	    m_truthNComponent = NAN;
	    m_truthEta = NAN;
	    m_truthPhi = NAN;
	    m_truthE = NAN;
	    m_truthPt = NAN;

	    Jet* truthjet = recoeval->max_truth_jet_by_energy(jet);
	    if (truthjet)
	      {
		m_truthID = truthjet->get_id();
		m_truthNComponent = truthjet->size_comp();
		m_truthEta = truthjet->get_eta();
		m_truthPhi = truthjet->get_phi();
		m_truthE = truthjet->get_e();
		m_truthPt = truthjet->get_pt();
	      }
	  }
	
	//SUBLEADING JET
	if (subleading_jet)
	  {
	    m_subleading_id = jet->get_id();
	    m_subleading_nComponent = jet->size_comp();
	    m_subleading_e = jet->get_e();
	    m_subleading_eta = jet->get_eta();
	    m_subleading_phi = jet->get_phi();
	    m_subleading_pt = jet->get_pt();	    
    
	    Jet* truthjet = recoeval->max_truth_jet_by_energy(jet);
	    if (truthjet)
	      {
		m_subleading_truthID = truthjet->get_id();
		m_subleading_truthNComponent = truthjet->size_comp();
		m_subleading_truthEta = truthjet->get_eta();
		m_subleading_truthPhi = truthjet->get_phi();
		m_subleading_truthE = truthjet->get_e();
		m_subleading_truthPt = truthjet->get_pt();
	      }
	  }
	
	// fill trees - jet track matching
	m_nMatchedTrack = 0;
	m_subleading_nMatchedTrack = 0;
	
	Float_t jet_eta = jet->get_eta();
	Float_t jet_phi = jet->get_phi();
	for (SvtxTrackMap::Iter iter = trackmap->begin();
	     iter != trackmap->end(); ++iter)

	  {
	    SvtxTrack* track = iter->second;

	    TVector3 v(track->get_px(), track->get_py(), track->get_pz());
	    const double dEta = v.Eta() - jet_eta;
	    const double dPhi = v.Phi() - jet_phi;
	    const double dR = sqrt(dEta * dEta + dPhi * dPhi);

	    if (dR < m_trackJetMatchingRadius)
	      {
		//matched track to jet

		if (hardest_jet)
		  {
		    assert(m_nMatchedTrack < kMaxMatchedTrack);
		    m_trackdR[m_nMatchedTrack] = dR;
		    m_trackpT[m_nMatchedTrack] = v.Perp();
		    ++m_nMatchedTrack;
		  }
		if (subleading_jet)
		  {
		    assert(m_subleading_nMatchedTrack < kMaxMatchedTrack);
		    m_subleading_trackdR[m_subleading_nMatchedTrack] = dR;
		    m_subleading_trackpT[m_subleading_nMatchedTrack] = v.Perp();
		    ++m_subleading_nMatchedTrack;
		  }
	      }

	    if ((m_nMatchedTrack >= kMaxMatchedTrack) || (m_nMatchedTrack >= kMaxMatchedTrack))
	      {
		cout << "MyJetAnalysis::process_event() - reached max track that matching a jet. Quit iterating tracks" << endl;
		break;
	      }

	  }  //    for (SvtxTrackMap::Iter iter = trackmap->begin();
	
      }  //   for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)      
    m_T->Fill(); //Fill Tree inside electron Loop
  }//electron Loop  
  return Fun4AllReturnCodes::EVENT_OK;
}
