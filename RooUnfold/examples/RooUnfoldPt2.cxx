//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
//#include "RooUnfoldIds.h"
#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;
double pt2bins_truth[]={0.,0.05,0.1,0.15,0.3,0.5,1.2,3.0};
double pt2bins_measu[]={0.,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.3,0.4,0.5,0.85,1.2,2.1,3.0};
//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

Double_t smear (Double_t xt)
{
  Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Double_t x= gRandom->Rndm();
  if (x>xeff) return cutdummy;
  Double_t xsmear= gRandom->Gaus(-2.5,0.2);     // bias and smear
  return xt+xsmear;
}

//==============================================================================
// Example Unfolding
//==============================================================================

void RooUnfoldPt2()
{
  TFile* file_emb = new TFile("../../macros/upc-dst-tinyTree-emb.root");
  TTree* tree = (TTree*) file_emb->Get("tinyTree");

  Int_t mMCnOS_tiny;
  Int_t mMCnTrack_tiny;
  Double32_t mMC_px_tiny[100];
  Double32_t mMC_py_tiny[100];
  Double32_t mMC_pz_tiny[100];
  Double32_t mMC_E_tiny[100];

  Int_t eventPass_tiny;
  Int_t mRECnTracks_tiny;
  Int_t mRECnSS_tiny;
  Int_t mRECnOS_tiny;
  Double32_t mREC_SS_px_tiny[100];
  Double32_t mREC_SS_py_tiny[100];
  Double32_t mREC_SS_pz_tiny[100];
  Double32_t mREC_SS_E_tiny[100];
  Double32_t mREC_OS_px_tiny[100];
  Double32_t mREC_OS_py_tiny[100];
  Double32_t mREC_OS_pz_tiny[100];
  Double32_t mREC_OS_E_tiny[100];

  tree->SetBranchAddress("mMCnTrack_tiny",&mMCnTrack_tiny);
  tree->SetBranchAddress("mMCnOS_tiny",&mMCnOS_tiny);
  tree->SetBranchAddress("mMC_px_tiny",&mMC_px_tiny);
  tree->SetBranchAddress("mMC_py_tiny",&mMC_py_tiny);
  tree->SetBranchAddress("mMC_pz_tiny",&mMC_pz_tiny);
  tree->SetBranchAddress("mMC_E_tiny",&mMC_E_tiny);
  
  tree->SetBranchAddress("eventPass_tiny",&eventPass_tiny);
  tree->SetBranchAddress("mRECnTracks_tiny",&mRECnTracks_tiny);
  tree->SetBranchAddress("mRECnSS_tiny",&mRECnSS_tiny);
  tree->SetBranchAddress("mRECnOS_tiny",&mRECnOS_tiny);
  tree->SetBranchAddress("mREC_SS_px_tiny",&mREC_SS_px_tiny);
  tree->SetBranchAddress("mREC_SS_py_tiny",&mREC_SS_py_tiny);
  tree->SetBranchAddress("mREC_SS_pz_tiny",&mREC_SS_pz_tiny);
  tree->SetBranchAddress("mREC_SS_E_tiny",&mREC_SS_E_tiny);
  tree->SetBranchAddress("mREC_OS_px_tiny",&mREC_OS_px_tiny);
  tree->SetBranchAddress("mREC_OS_py_tiny",&mREC_OS_py_tiny);
  tree->SetBranchAddress("mREC_OS_pz_tiny",&mREC_OS_pz_tiny);
  tree->SetBranchAddress("mREC_OS_E_tiny",&mREC_OS_E_tiny);

  cout << "==================================== TRAIN ====================================" << endl;
  int pt2TruthNbins = sizeof(pt2bins_truth)/sizeof(pt2bins_truth[0]) - 1;
  int pt2MeasuNbins = sizeof(pt2bins_measu)/sizeof(pt2bins_measu[0]) - 1;

  TH1D* hTruth = new TH1D("hTruth","hTruth",pt2TruthNbins,pt2bins_truth);
  TH1D* hMeasu = new TH1D("hMeasu","hMeasu",pt2MeasuNbins,pt2bins_measu);

  RooUnfoldResponse response(hMeasu, hTruth);
  
  for(int i=0;i<tree->GetEntries();i++){
    tree->GetEntry(i);

    if(eventPass_tiny!=1) continue;

  }

  // Train with a Breit-Wigner, mean 0.3 and width 2.5.
  for (Int_t i= 0; i<100000; i++) {
    Double_t xt= gRandom->BreitWigner (0.3, 2.5);
    Double_t x= smear (xt);
    if (x!=cutdummy)
      response.Fill (x, xt);
    else
      response.Miss (xt);
  }

  cout << "==================================== TEST =====================================" << endl;
  TH1D* hTrue= new TH1D ("true", "Test Truth",    40, -10.0, 10.0);
  TH1D* hMeas= new TH1D ("meas", "Test Measured", 40, -10.0, 10.0);
  // Test with a Gaussian, mean 0 and width 2.
  // for (Int_t i=0; i<10000; i++) {
  //   Double_t xt= gRandom->Gaus (0.0, 2.0), x= smear (xt);
  //   hTrue->Fill(xt);
  //   if (x!=cutdummy) hMeas->Fill(x);
  // }

  cout << "==================================== UNFOLD ===================================" << endl;
  // RooUnfoldBayes   unfold (&response, hMeas, 4);    // OR
//RooUnfoldSvd     unfold (&response, hMeas, 20);   // OR
//RooUnfoldTUnfold unfold (&response, hMeas);       // OR
//RooUnfoldIds     unfold (&response, hMeas, 1);

  // TH1D* hReco= (TH1D*) unfold.Hreco();

  // TCanvas* c1= new TCanvas("canvas","canvas");

  // unfold.PrintTable (cout, hTrue);
  // hReco->Draw();
  // hMeas->Draw("SAME");
  // hTrue->SetLineColor(8);
  // hTrue->Draw("SAME");

  // c1->SaveAs("RooUnfoldExample.pdf");

}

#ifndef __CINT__
int main () { RooUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif
