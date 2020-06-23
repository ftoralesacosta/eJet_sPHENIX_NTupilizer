#include "MyJetAnalysis.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/getClass.h>

#include <g4eval/JetEvalStack.h>
#include <g4main/PHG4Particle.h>
#include <g4eval/SvtxTrackEval.h>

#include <trackbase_historic/SvtxTrackMap.h>

#include <g4jets/JetMap.h>

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
#include <stdio.h>
#include <math.h>
using namespace std;

MyJetAnalysis::MyJetAnalysis(const std::string& recojetname, const std::string& truthjetname, const std::string& outputfilename)
  : SubsysReco("MyJetAnalysis_" + recojetname + "_" + truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange(-1, 1)
  , m_ptRange(5, 100)
  , m_trackJetMatchingRadius(.7)
  , m_hInclusiveE(nullptr)
  , m_hInclusiveEta(nullptr)
  , m_hInclusivePhi(nullptr)
  , m_T(nullptr)
  , m_event(-1)
  , m_id(-1)
  , m_nComponent(-1)
  , m_eta(numeric_limits<float>::signaling_NaN())
  , m_phi(numeric_limits<float>::signaling_NaN())
  , m_e(numeric_limits<float>::signaling_NaN())
  , m_pt(numeric_limits<float>::signaling_NaN())
  , m_truthID(-1)
  , m_truthNComponent(-1)
  , m_truthEta(numeric_limits<float>::signaling_NaN())
  , m_truthPhi(numeric_limits<float>::signaling_NaN())
  , m_truthE(numeric_limits<float>::signaling_NaN())
  , m_truthPt(numeric_limits<float>::signaling_NaN())
  , m_nMatchedTrack(-1)
{
  m_trackdR.fill(numeric_limits<float>::signaling_NaN());
  m_trackpT.fill(numeric_limits<float>::signaling_NaN());
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

  //Trees
  m_T = new TTree("T", "MyJetAnalysis Tree");

  //      int m_event;
  m_T->Branch("m_event", &m_event, "event/I");
  //      int m_id;
  m_T->Branch("id", &m_id, "id/I");
  //      int m_nComponent;
  m_T->Branch("nComponent", &m_nComponent, "nComponent/I");
  //      float m_eta;
  m_T->Branch("eta", &m_eta, "eta/F");
  //      float m_phi;
  m_T->Branch("phi", &m_phi, "phi/F");
  //      float m_e;
  m_T->Branch("e", &m_e, "e/F");
  //      float m_pt;
  m_T->Branch("pt", &m_pt, "pt/F");
  //
  //      int m_truthID;
  m_T->Branch("truthID", &m_truthID, "truthID/I");
  //      int m_truthNComponent;
  m_T->Branch("truthNComponent", &m_truthNComponent, "truthNComponent/I");
  //      float m_truthEta;
  m_T->Branch("truthEta", &m_truthEta, "truthEta/F");
  //      float m_truthPhi;
  m_T->Branch("truthPhi", &m_truthPhi, "truthPhi/F");
  //      float m_truthE;
  m_T->Branch("truthE", &m_truthE, "truthE/F");
  //      float m_truthPt;
  m_T->Branch("truthPt", &m_truthPt, "truthPt/F");
  //
  //      //! number of matched tracks
  //      int m_nMatchedTrack;
  m_T->Branch("nMatchedTrack", &m_nMatchedTrack, "nMatchedTrack/I");
  //      std::array<float, kMaxMatchedTrack> m_trackdR;
  m_T->Branch("id", m_trackdR.data(), "trackdR[nMatchedTrack]/F");
  //      std::array<float, kMaxMatchedTrack> m_trackpT;
  m_T->Branch("id", m_trackpT.data(), "trackpT[nMatchedTrack]/F");

  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::End(PHCompositeNode* topNode)
{
  cout << "MyJetAnalysis::End - Outoput to " << m_outputFileName << endl;
  PHTFileServer::get().cd(m_outputFileName);

  m_hInclusiveE->Write();
  m_hInclusiveEta->Write();
  m_hInclusivePhi->Write();
  m_T->Write();

  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::InitRun(PHCompositeNode* topNode)
{
  m_jetEvalStack = shared_ptr<JetEvalStack>(new JetEvalStack(topNode, m_recoJetName, m_truthJetName));
  
  //m_jetTruthEval = shared_ptr<JetTruthEval>(new JetTruthEval(topNode,m_truthJetName));
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

  //interface to truth particles
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if ( !truthinfo )
    {
      cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
      exit(-1);
    }


  //Debugging File. Print all particles and Tracks.
  FILE * fp = fopen ("./test_particles.txt","a+");
  fprintf(fp,"EVENT NUMBER %i\n",m_event);
  fprintf(fp,"List of all truth Particles:\n");

  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  for ( PHG4TruthInfoContainer::ConstIterator eter = range.first; eter != range.second; ++eter )
  //for ( auto eter = range.first; eter != range.second; ++eter ))
    {
      PHG4Particle* g4particle = eter->second;
      int particleID = g4particle->get_pid();
      float particleE = g4particle->get_e();
      fprintf(fp,"Particle ID: %i, Energy: %f\n",particleID,particleE);
    }
  fprintf(fp,"\n List of all Reco Tracks:\n");
  for (SvtxTrackMap::Iter iter = trackmap->begin();iter != trackmap->end();++iter)
    {
      SvtxTrack* track = iter->second;
      fprintf(fp,"Reco Track pT: %f \n",track->get_pt());
    }
   fprintf(fp,"\n");
  

   //Main Jet Loop
   for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)
  {
    Jet* jet = iter->second;
    assert(jet);

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

    // fill histograms
    assert(m_hInclusiveE);
    m_hInclusiveE->Fill(jet->get_e());
    assert(m_hInclusiveEta);
    m_hInclusiveEta->Fill(jet->get_eta());
    assert(m_hInclusivePhi);
    m_hInclusivePhi->Fill(jet->get_phi());

    // fill trees - jet spectrum
    Jet* truthjet = recoeval->max_truth_jet_by_energy(jet);
 
    m_id = jet->get_id();
    m_nComponent = jet->size_comp();
    m_eta = jet->get_eta();
    m_phi = jet->get_phi();
    m_e = jet->get_e();
    m_pt = jet->get_pt();  
      
    m_truthID = -1;
    m_truthNComponent = -1;
    m_truthEta = NAN;
    m_truthPhi = NAN;
    m_truthE = NAN;
    m_truthPt = NAN;

    if (truthjet)
    {
      m_truthID = truthjet->get_id();
      m_truthNComponent = truthjet->size_comp();
      m_truthEta = truthjet->get_eta();
      m_truthPhi = truthjet->get_phi();
      m_truthE = truthjet->get_e();
      m_truthPt = truthjet->get_pt();

      //Print Truth Jet Energy & Constituent Info
      fprintf(fp,"Jet TruthE: %f, \n",m_truthE); 
      std::set<PHG4Particle*> truthj_particle_set = m_jetEvalStack->get_truth_eval()->all_truth_particles(truthjet);
      for (auto i_t:truthj_particle_set)
	{
	  // float i_px = pow(i_t->get_px(),2.0);
	  // float i_py = pow(i_t->get_py(),2.0);
	  // float i_pt = sqrt(pow(i_px,2.0),pow(i_py,2.0));
	  fprintf(fp,"Truth Jet Constituent PID: %i, True Energy: %f\n",i_t->get_pid(),i_t->get_e());
	}
       fprintf(fp,"\n\n");
    }

    //Print Reco Jet Energy & Constituent Info
    // fprintf(fp,"Jet RecoE:%f \n",m_e);
    // std::set<PHG4Particle*> particle_set = m_jetEvalStack->get_reco_eval()->all_truth_particles(jet);
    // for (auto i_r:particle_set)
    //   {
    // 	fprintf(fp,"Reco Jet Truth Constituents PID: %i, True Energy: %f\n",i_r->get_pid(),i_r->get_e());
    //   }
    // fprintf(fp,"\n\n");

    // fill trees - jet track matching


    //Print Reco Jet Energy & Constituent Info
    fprintf(fp,"Jet RecoE:%f \n",m_e);

    m_nMatchedTrack = 0;
  
     for (SvtxTrackMap::Iter iter = trackmap->begin();
         iter != trackmap->end();
         ++iter)
    {
      SvtxTrack* track = iter->second;

      TVector3 v(track->get_px(), track->get_py(), track->get_pz());
      const double dEta = v.Eta() - m_eta;
      const double dPhi = v.Phi() - m_phi;
      const double dR = sqrt(dEta * dEta + dPhi * dPhi);

      if (dR < m_trackJetMatchingRadius)
      {
        //matched track to jet
	
        assert(m_nMatchedTrack < kMaxMatchedTrack);

        m_trackdR[m_nMatchedTrack] = dR;
        m_trackpT[m_nMatchedTrack] = v.Perp();

	PHG4Particle* particle_dr = m_jetEvalStack->get_stvx_eval_stack()->get_track_eval()->max_truth_particle_by_nclusters(track);
	
	fprintf(fp,"Reco Track Constituent; PID: %i  track pT: = %f, track eta = %f, track phi = %f\n",particle_dr->get_pid(),track->get_pt(),track->get_eta(),track->get_phi());
	
        ++m_nMatchedTrack;
      }

      if (m_nMatchedTrack >= kMaxMatchedTrack)
      {
        cout << "MyJetAnalysis::process_event() - reached max track that matching a jet. Quit iterating tracks" << endl;
        break;
      }

    }  //reco track loop
     fprintf(fp,"\n\n");
     m_T->Fill();
  }  //Jet Map loop
   fclose(fp);
  return Fun4AllReturnCodes::EVENT_OK;
}