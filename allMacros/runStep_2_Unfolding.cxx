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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TSystem.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
//#include "RooUnfoldIds.h"
#endif

//==============================================================================
// Global definitions
//==============================================================================
#include "RiceStyle.h"
#include "inputRootFile.h"
const Double_t cutdummy= -99999.0;

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

Double_t fitNonExp(Double_t *x, Double_t *par){

  return par[0]*TMath::Power((1+(par[1]/par[2])*x[0]), -par[2]);

}
Double_t fitExpPlusDisso(Double_t *x, Double_t *par){
  
  double coherent = par[0]*TMath::Exp(par[1]*x[0]);
  double incoherent = par[2]*TMath::Exp(par[3]*x[0]);
  double dissoc = par[4]*TMath::Power((1+(par[5]/par[6])*x[0]), -par[6]);

  return coherent+incoherent+dissoc;

}
//==============================================================================
// Example Unfolding
//==============================================================================

void runStep_2_Unfolding( const bool doSys_ = false )
{

  gStyle->SetErrorX(0);

  //input
  TFile* file_emb = new TFile("output-PreStep_4-embedding.root");
  TTree* tree = (TTree*) file_emb->Get("tinyTree");
  TFile* file_binbybin = new TFile("output-PreStep_2-embedding.root");
  TH1D* hMCDielectronPt2 = (TH1D*) file_binbybin->Get("hMCDielectronPt2");
  TH1D* hDielectronPt2 = (TH1D*) file_binbybin->Get("hDielectronPt2");
  TH1D* hist_binBybin = (TH1D*) make_systematicRatio(hDielectronPt2,hMCDielectronPt2);//cannot be used at the moment

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
  Double32_t weight_evt;
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
  tree->SetBranchAddress("weight_evt",&weight_evt);
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
    vector< TLorentzVector> pREC_collection, pMC_collection;
   
    for(int imc=0;imc<mMCnOS_tiny;imc++){
      pMC.Clear();
      pMC.SetPxPyPzE(mMC_px_tiny[imc],mMC_py_tiny[imc],mMC_pz_tiny[imc],mMC_E_tiny[imc]);
      pMC_collection.push_back( pMC );
    }

    for(int irec=0;irec<mRECnOS_tiny;irec++){
      pREC.Clear();
      pREC.SetPxPyPzE(mREC_OS_px_tiny[irec],mREC_OS_py_tiny[irec],mREC_OS_pz_tiny[irec],mREC_OS_E_tiny[irec]);
      pREC_collection.push_back( pREC );
    }
    
    //case 1 to find fake
    if( pMC_collection.size() < pREC_collection.size() ){

      double deltaR_min = 999.;
      int real_index = -1;

      for( int j=0;j<pREC_collection.size();j++){
        for(int k = 0; k<pMC_collection.size();k++){
          TVector3 rec3 = pREC_collection[j].Vect();
          TVector3 mc3 = pMC_collection[k].Vect();
          if(rec3.DeltaR( mc3 ) < deltaR_min){
            deltaR_min = rec3.DeltaR( mc3 );
            real_index = j;
          }
        }
      }

      for(int ifak=0;ifak<pREC_collection.size();ifak++){
        if( ifak == real_index ) continue;
        double pt2REC = pREC_collection[ ifak ].Pt()*pREC_collection[ ifak ].Pt();
        response.Fake(pt2REC, weight_evt);
      }
    }
    //case 2, should be mostly 1-1
    else if( pMC_collection.size() == pREC_collection.size() ){

      double deltaR_min = 999.;
      int real_index_rec = -1;
      int real_index_mc = -1;
      for( int j=0;j<pREC_collection.size();j++){
        for(int k = 0; k<pMC_collection.size();k++){
          TVector3 rec3 = pREC_collection[j].Vect();
          TVector3 mc3 = pMC_collection[k].Vect();
          if(rec3.DeltaR( mc3 ) < deltaR_min){
            deltaR_min = rec3.DeltaR( mc3 );
            real_index_rec = j;
            real_index_mc = k;
          }
            
        }
      }
      if( real_index_rec < 0  ) continue;
      double pt2REC = pREC_collection[real_index_rec].Pt()*pREC_collection[real_index_rec].Pt();
      double pt2MC  = pMC_collection[real_index_mc].Pt()*pMC_collection[real_index_mc].Pt();
      response.Fill(pt2REC,pt2MC, weight_evt);
    }
    //case to find miss
    else{
      double deltaR_min = 999.;
      int real_index_rec = -1;
      int real_index_mc = -1;
      for( int j=0;j<pREC_collection.size();j++){
        for(int k = 0; k<pMC_collection.size();k++){
          TVector3 rec3 = pREC_collection[j].Vect();
          TVector3 mc3 = pMC_collection[k].Vect();
          if(rec3.DeltaR( mc3 ) < deltaR_min){
            deltaR_min = rec3.DeltaR( mc3 );
            real_index_rec = j;
            real_index_mc = k;
          }
        }
      }

      for(int imiss=0;imiss<pMC_collection.size();imiss++){
        if( imiss == real_index_mc ) {
          double pt2REC = pREC_collection[real_index_rec].Pt()*pREC_collection[real_index_rec].Pt();
          double pt2MC  = pMC_collection[real_index_mc].Pt()*pMC_collection[real_index_mc].Pt();
          response.Fill(pt2REC, pt2MC, weight_evt);
          continue;
        }
        double pt2MC = pMC_collection[ imiss ].Pt()*pMC_collection[ imiss ].Pt();
        response.Fake(pt2MC, weight_evt);
      }
    }
  }

  cout << "==================================== TEST =====================================" << endl;
  
  TFile* file_data = 0;
  TString sys_filename = "systematics-input/electron/output-Step_1-electron_loose.root";
  //for systematic uncertainty, use
  if(doSys_){file_data = new TFile( sys_filename );}
  //Default production
  else{file_data = new TFile("output-Step_1.root");}
  
  TH1D* hMeasured = (TH1D*) file_data->Get("JpsiPt2");

  cout << "==================================== UNFOLD ===================================" << endl;
  RooUnfoldBayes   unfold (&response, hMeasured, 20);    // OR
// RooUnfoldSvd     unfold (&response, hMeasured, 20);   // OR
// RooUnfoldBinByBin     unfold (&response, hMeasured);   // OR
  // RooUnfoldTUnfold unfold (&response, hMeasured);       // OR
//RooUnfoldIds     unfold (&response, hMeas, 1);

  TH1D* hReco= (TH1D*) unfold.Hreco();

  TCanvas* c1= new TCanvas("canvas","canvas");
  
  // unfold.PrintTable (cout, hTrue);
  for(int j=0;j<hMeasured->GetNbinsX();j++){
    hReco->SetBinContent(j+1, hReco->GetBinContent(j+1)/(hReco->GetBinWidth(j+1)) );
    hReco->SetBinError(j+1, hReco->GetBinError(j+1)/(hReco->GetBinWidth(j+1)) );

    double eff = hist_binBybin->GetBinContent(j+1);
    hMeasured->SetBinContent(j+1, hMeasured->GetBinContent(j+1)/(hMeasured->GetBinWidth(j+1)) );
    hMeasured->SetBinError(j+1, hMeasured->GetBinError(j+1)/(hMeasured->GetBinWidth(j+1)) );
  }

  hReco->SetMarkerStyle(20);
  hReco->Draw("P");
  hMeasured->SetMarkerStyle(24);
  hMeasured->Draw("PSAME");

  TFile outfile("output-Step_2.root","RECREATE");
  hReco->Write();
  response.Hresponse()->Write();

}

#ifndef __CINT__
int main () { runStep_2_Unfolding(); return 0; }  // Main program when run stand-alone
#endif
