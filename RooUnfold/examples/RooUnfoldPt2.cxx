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
using namespace std;

#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <vector>
#include "TLorentzVector.h"
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"

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
double pt2bins_measu[]={0.,0.05,0.1,0.15,0.3,0.5,1.2,3.0};
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
  //input
  TFile* file_emb = new TFile("../macros/upc-dst-tinyTree-emb.root");
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

  RooUnfoldResponse response (hMeasu, hTruth);
  
  for(int i=0;i<tree->GetEntries();i++){
    tree->GetEntry(i);
    
    if( eventPass_tiny != 1 ) continue;

    TLorentzVector pMC(0,0,0,0);
    TLorentzVector pREC(0,0,0,0);
    TLorentzVector pREC_max(0,0,0,0);
    TLorentzVector pMC_max(0,0,0,0);
    double ptmax = 0.;
    for(int imc=0;imc<mMCnOS_tiny;imc++){
      pMC.Clear();
      pMC.SetPxPyPzE(mMC_px_tiny[imc],mMC_py_tiny[imc],mMC_pz_tiny[imc],mMC_E_tiny[imc]);
      if( pMC.Pt() > ptmax ) {
        ptmax = pMC.Pt();
        pMC_max = pMC;
      }
    }
    ptmax = 0.;
    for(int irec=0;irec<mRECnOS_tiny;irec++){
      pREC.Clear();
      pREC.SetPxPyPzE(mREC_OS_px_tiny[irec],mREC_OS_py_tiny[irec],mREC_OS_pz_tiny[irec],mREC_OS_E_tiny[irec]);
      if( pREC.Pt() > ptmax ) {
        ptmax = pREC.Pt();
        pREC_max = pREC;
      }
    }
    
    if( pREC_max.Pt() != 0 ){
      double pt2REC = pREC_max.Pt()*pREC_max.Pt();
      double pt2MC  = pMC_max.Pt()*pMC_max.Pt();
      response.Fill(pt2REC,pt2MC);
    }
    else{
      if( pMC_max.Pt() != 0){
        double pt2MC  = pMC_max.Pt()*pMC_max.Pt();
        response.Miss(pt2MC);
      }
    }
  }

  cout << "==================================== TEST =====================================" << endl;
  TFile* file_test = new TFile("../macros/upc-dst-histo-jpsi-simple_emb_zerobias_cohNincoh.root");
  TH2D* hJpsiMass_Pt2 = (TH2D*) file_test->Get("hJpsiMass_Pt2");
  TH1D* hMeasured = (TH1D*) hJpsiMass_Pt2->ProjectionY("hMeasured",1,100);
  TH1D* hTrue = (TH1D*) file_test->Get("hMCDielectronPt2");

  cout << "==================================== UNFOLD ===================================" << endl;
  // RooUnfoldBayes   unfold (&response, hMeasured, 10);    // OR
//RooUnfoldSvd     unfold (&response, hMeas, 20);   // OR
  RooUnfoldTUnfold unfold (&response, hMeasured);       // OR
//RooUnfoldIds     unfold (&response, hMeas, 1);

  TH1D* hReco= (TH1D*) unfold.Hreco();
  TH2D* h2D = (TH2D*) response.HresponseNoOverflow();

  TCanvas* c1= new TCanvas("canvas","canvas");
  
  unfold.PrintTable (cout, hTrue);
  for(int j=0;j<hMeasured->GetNbinsX();j++){
    hReco->SetBinContent(j+1, hReco->GetBinContent(j+1)/(hReco->GetBinWidth(j+1)) );
    hReco->SetBinError(j+1, hReco->GetBinError(j+1)/(hReco->GetBinWidth(j+1)) );

    hMeasured->SetBinContent(j+1, hMeasured->GetBinContent(j+1)/(hMeasured->GetBinWidth(j+1)) );
    hMeasured->SetBinError(j+1, hMeasured->GetBinError(j+1)/(hMeasured->GetBinWidth(j+1)) );

    hTrue->SetBinContent(j+1, hTrue->GetBinContent(j+1)/(hTrue->GetBinWidth(j+1)) );
    hTrue->SetBinError(j+1, hTrue->GetBinError(j+1)/(hTrue->GetBinWidth(j+1)) );
  }
  hReco->SetMarkerStyle(20);
  hReco->Draw("P");
  hMeasured->Draw("SAME");
  hTrue->SetMarkerStyle(25);
  hTrue->Draw("PSAME");

  TCanvas* c2 = new TCanvas("c2","c2");
  h2D->Draw("colz");

  // c1->SaveAs("RooUnfoldExample.pdf");

}

#ifndef __CINT__
int main () { RooUnfoldPt2(); return 0; }  // Main program when run stand-alone
#endif
