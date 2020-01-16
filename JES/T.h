//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 19 13:56:59 2019 by ROOT version 6.16/00
// from TTree T/MyJetAnalysis Tree
// found on file: myjetanalysis.root
//////////////////////////////////////////////////////////

#ifndef T_h
#define T_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class T {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           m_event;
   Int_t           id;
   Int_t           nComponent;
   Int_t           m_eta;
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
   Float_t         id[4];   //[nMatchedTrack]
   Float_t         id[4];   //[nMatchedTrack]

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_id;   //!
   TBranch        *b_nComponent;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_e;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_truthID;   //!
   TBranch        *b_truthNComponent;   //!
   TBranch        *b_truthEta;   //!
   TBranch        *b_truthPhi;   //!
   TBranch        *b_truthE;   //!
   TBranch        *b_truthPt;   //!
   TBranch        *b_nMatchedTrack;   //!
   TBranch        *b_id;   //!
   TBranch        *b_id;   //!

   T(TTree *tree=0);
   virtual ~T();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef T_cxx
T::T(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("myjetanalysis.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("myjetanalysis.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

T::~T()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t T::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t T::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void T::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("m_event", &m_event, &b_event);
   fChain->SetBranchAddress("id", &id, &b_id);
   fChain->SetBranchAddress("nComponent", &nComponent, &b_nComponent);
   fChain->SetBranchAddress("m_eta", &m_eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("e", &e, &b_e);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("truthID", &truthID, &b_truthID);
   fChain->SetBranchAddress("truthNComponent", &truthNComponent, &b_truthNComponent);
   fChain->SetBranchAddress("truthEta", &truthEta, &b_truthEta);
   fChain->SetBranchAddress("truthPhi", &truthPhi, &b_truthPhi);
   fChain->SetBranchAddress("truthE", &truthE, &b_truthE);
   fChain->SetBranchAddress("truthPt", &truthPt, &b_truthPt);
   fChain->SetBranchAddress("nMatchedTrack", &nMatchedTrack, &b_nMatchedTrack);
//    fChain->SetBranchAddress("id", id, &b_id);
//    fChain->SetBranchAddress("id", id, &b_id);
   Notify();
}

Bool_t T::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void T::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t T::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef T_cxx
