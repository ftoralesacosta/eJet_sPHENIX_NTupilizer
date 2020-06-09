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
#include <TH1D.h>
#include <TObjArray.h>
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

int main(int argc, char *argv[])
{
  
  if (argc < 2) {
    std::cout<<"Syntax: [Command] [File]"<<std::endl;
    exit(EXIT_FAILURE);
  }

  int iarg = 1;
  TString root_file = (TString)argv[iarg];
  std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
  TFile *file = TFile::Open(root_file);
  
  if (file == NULL) {
    std::cout << " File Fail" << std::endl;
    exit(EXIT_FAILURE); 
  } 
  
  file->Print();
  TTree *_tree_event = dynamic_cast<TTree *>(file->Get("ntp_track"));

  if (_tree_event == NULL) {
    std::cout << " Tree Fail " << std::endl;
    exit(EXIT_FAILURE);
  }

  //Declare Leaf Types
  Float_t         gpx;
  Float_t         gpy;
  Float_t         gpz;
  Float_t         gp;
  Float_t         px;
  Float_t         py;
  Float_t         pz;
  Float_t         p;

  //Declare Branches
  TBranch        *b_gpx;
  TBranch        *b_gpy;
  TBranch        *b_gpz;
  TBranch        *b_px;
  TBranch        *b_py;
  TBranch        *b_pz;
 
  _tree_event->SetBranchAddress("gpx", &gpx, &b_gpx);
  _tree_event->SetBranchAddress("gpy", &gpy, &b_gpy);
  _tree_event->SetBranchAddress("gpz", &gpz, &b_gpz);
  _tree_event->SetBranchAddress("px", &px, &b_px);
  _tree_event->SetBranchAddress("py", &py, &b_py);
  _tree_event->SetBranchAddress("pz", &pz, &b_pz);

  TObjArray proj;
  TObjArray fitArr;
  float max1 = 0.05;
  float max2 = 0.10;
  int nbin = 50;
  int bin[5][2] = {{0,1}, {1,2}, {2,5}, {5,10}, {10,20}};
  for (int i = 0; i < 5; i++){
    if (i < 3){
	proj.Add(new TH1D(Form("%d < p^tru < %d GeV", bin[i][0], bin[i][1]),"", nbin, -max1, max1));
        fitArr.Add(new TF1(Form("fit%lu", i), "[0]/sqrt(2*TMath::Pi()*[2]^2)*exp(-0.5*((x-[1])/[2])**2) + [3]/sqrt(2*TMath::Pi()*[4]^2)*exp(-0.5*((x-[1])/[4])**2)", -max1, max1));
    }
    else {
	proj.Add(new TH1D(Form("%d < p^tru < %d GeV", bin[i][0], bin[i][1]),"", nbin, -max2, max2));
        fitArr.Add(new TF1(Form("fit%lu", i), "[0]/sqrt(2*TMath::Pi()*[2]^2)*exp(-0.5*((x-[1])/[2])**2) + [3]/sqrt(2*TMath::Pi()*[4]^2)*exp(-0.5*((x-[1])/[4])**2)", -max2, max2));
    }
  }

  Long64_t nentries = _tree_event->GetEntries();
  for(Long64_t ie = 0; ie < nentries ; ie++){
    
    _tree_event->GetEntry(ie); //each entry is a 5 GeV Electron
    fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ie, nentries);

    gp = sqrt( gpx*gpx + gpy*gpy + gpz*gpz );
    p = sqrt( px*px + py*py + pz*pz );

    int binNum = -1;
    for (int i = 0; i < 5; i++){
      if (gp > bin[i][0] && gp < bin[i][1]){
        binNum = i;
      }
    }
    if (binNum > -1){
      double dp_p = (p - gp) / gp;
      ((TH1D *)proj.At(binNum))->Fill(dp_p);  
    }
  } //entry loop
  
  //Write to new file
  
  TFile f("output/default.root","recreate");
  gStyle->SetOptStat(111);
  gStyle->SetOptFit(11);
  gStyle->SetStatW(0.155);
  gStyle->SetStatH(0.255);
  TCanvas *canv = new TCanvas("canv", "canv", 1600, 1200);
  canv->Divide(3,2);
  
  for (int i = 0; i < 5; i++){
    canv->cd(i+1);
    //dist = (TH1D *)proj.At(i);
    //fit = ((TF1 *)fitArr.At(i));
    ((TF1 *)fitArr.At(i))->SetParameter(4, -0.1);
    ((TF1 *)fitArr.At(i))->SetParameter(2, -0.1);
    ((TH1D *)proj.At(i))->Fit((TF1 *)fitArr.At(i), "rll");
    ((TH1D *)proj.At(i))->GetXaxis()->SetTitle("dp/p");
    ((TH1D *)proj.At(i))->GetYaxis()->SetTitle("counts");
    ((TH1D *)proj.At(i))->SetMarkerStyle(20);
    ((TH1D *)proj.At(i))->SetLineColor(kBlack);
    ((TH1D *)proj.At(i))->Write();
    ((TF1 *)fitArr.At(i))->Write();
    ((TH1D *)proj.At(i))->Draw();
  }
  canv->Write();

}
