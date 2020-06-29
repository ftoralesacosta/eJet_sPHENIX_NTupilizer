#ifndef MYJETANALYSIS_ALLSI_H
#define MYJETANALYSIS_ALLSI_H

#include <fun4all/SubsysReco.h>

#include <g4eval/JetTruthEval.h>

#include <set>

#include <memory>
#include <string>
#include <utility>  // std::pair, std::make_pair

#if ! defined(__CINT__) || defined(__CLING__)
#include <array>
#endif  // #ifndef __CINT__

#define NaN numeric_limits<float>::signaling_NaN()

class PHCompositeNode;
class JetEvalStack;
class TTree;
class TH1;

class Jet;
class JetMap;
class SvtxTrackMap;
class PHG4TruthInfoContainer;

/// \class MyJetAnalysis_AllSi
class MyJetAnalysis_AllSi : public SubsysReco
{
	public:
		MyJetAnalysis_AllSi(
				const std::string &recojetname = "AntiKt_Track_r10",
				const std::string &truthjetname = "AntiKt_Truth_r10",
				const std::string &outputfilename = "myjetanalysis.root");

		virtual ~MyJetAnalysis_AllSi();

		//! set eta range
		void setEtaRange(double low, double high)
		{
			m_etaRange.first = low;
			m_etaRange.second = high;
		}
		//! set eta range
		void setPtRange(double low, double high)
		{
			m_ptRange.first = low;
			m_ptRange.second = high;
		}

		void setelectronEmin(double Emin)
		{
			m_eEmin = Emin;  
		}

		int Init(PHCompositeNode *topNode);
		int InitRun(PHCompositeNode *topNode);
		int process_event(PHCompositeNode *topNode);
		int End(PHCompositeNode *topNode);

		float get_jet_radius_from_string( string jetname ){ // Assumming jetname is a string that ends, for instance, in "_r04" if the radius is 0.4 
			return .1*stof(jetname.substr(jetname.find("_r")+2));
		}

	private:

#if ! defined(__CINT__) || defined(__CLING__)

		//! cache the jet evaluation modules
		std::shared_ptr<JetEvalStack> m_jetEvalStack;

		Jet* max_truth_jet_by_track_fastsim(Jet*);
		int get_track_fastsim_contribution(Jet* recojet, Jet* truthjet);
		std::set<PHG4Particle*> all_truth_particles(Jet* truthjet);

		JetMap* all_truth_jets;
		SvtxTrackMap* _trackmap;
		PHG4TruthInfoContainer* _truthinfo;

		std::string m_recoJetName;
		std::string m_truthJetName;
		std::string m_outputFileName;

		//! eta range
		std::pair<double, double> m_etaRange;

		//! pT range
		std::pair<double, double> m_ptRange;

		//electron Energy min
		double m_eEmin;

		//electron-jet matching radius for veto
		double m_electronJetMatchingRadius;

		//! max track-jet matching radius
		double m_trackJetMatchingRadius;

		//! Output histograms
		TH1 *m_hInclusiveE;
		TH1 *m_hInclusiveEta;
		TH1 *m_hInclusivePhi;
		TH1 *m_hInclusiveNJets;

		//! Output Tree variables
		TTree *m_T;

		int m_event;

		//Electron Truth Variables
		float m_electron_truthEta;
		float m_electron_truthPhi;
		float m_electron_truthE;
		float m_electron_truthPt;
		float m_electron_truthpX;
		float m_electron_truthpY;
		float m_electron_truthpZ;
		int m_electron_truthPID;

		int m_njets;
		int m_ntruthjets;
		int m_nAlltruthjets;

		enum {MaxNumJets = 20, kMaxConstituents = 100};
		
		// Reconstructed jets
		std::array<int,MaxNumJets> m_id;
		std::array<int,MaxNumJets> m_nComponent;
		std::array<float,MaxNumJets> m_eta;
		std::array<float,MaxNumJets> m_phi;
		std::array<float,MaxNumJets> m_e;
		std::array<float,MaxNumJets> m_pt;

		// Truth jets matched to a reconstructed jet
		std::array<int,MaxNumJets> m_matched_truthID;
		std::array<int,MaxNumJets> m_matched_truthNComponent;
		std::array<float,MaxNumJets> m_matched_truthEta;
		std::array<float,MaxNumJets> m_matched_truthPhi;
		std::array<float,MaxNumJets> m_matched_truthE;
		std::array<float,MaxNumJets> m_matched_truthPt;  

		// All truth jets
		std::array<int,MaxNumJets> m_all_truthID;
		std::array<int,MaxNumJets> m_all_truthNComponent;
		std::array<float,MaxNumJets> m_all_truthEta;
		std::array<float,MaxNumJets> m_all_truthPhi;
		std::array<float,MaxNumJets> m_all_truthE;
		std::array<float,MaxNumJets> m_all_truthPt;

		// Constituent for all matched jets
		std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_matched_Constituent_truthPID; //PID float in order to fill branches with NaN correctly.
		std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_matched_Constituent_truthCharge;
		std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_matched_Constituent_truthEta;
		std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_matched_Constituent_truthPhi;
		std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_matched_Constituent_truthPt;
		std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_matched_Constituent_truthE;

		// Charged (neutral-subtracted) jets
		std::array<int,MaxNumJets> m_matched_charged_truthNComponent; 
		std::array<float,MaxNumJets> m_matched_charged_truthEta;
		std::array<float,MaxNumJets> m_matched_charged_truthPhi;
		std::array<float,MaxNumJets> m_matched_charged_truthE;
		std::array<float,MaxNumJets> m_matched_charged_truthPt;

		// Constituents for all truth jets
		std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_All_Constituent_truthPID;
		std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_All_Constituent_truthCharge;
		std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_All_Constituent_truthEta;
		std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_All_Constituent_truthPhi;
		std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_All_Constituent_truthPt;
		std::array<std::array<float, kMaxConstituents >, MaxNumJets > m_All_Constituent_truthE;

#endif  // #ifndef __CINT__
};

#endif  // MYJETANALYSIS_ALLSI_H
