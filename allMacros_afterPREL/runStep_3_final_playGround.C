#include "RiceStyle.h"
#include "inputRootFile.h"
using namespace std;
#define PI 3.1415926

/*
- d\sigma/dt 
cross section relation:
(1/flux) * Nobs / [ Lint * Br * dt * (A x e) * trig * event * dy ]
*/

double Lint = 93.418*0.912;// nb**-1 .. 0.934 fraction TPCvsVPD in hadronic events
double flux = 11.78; //given by starlight, k(dn/dk), where k is photon energy ~ M/2 Exp(-y) 
double Br = 0.0593; // ee
double trigBBCTOF = 0.8*0.98*0.78; //only TOFmult < 2 and > 6 efficiency and BBC veto efficiency
double event = 1.0; //0.2156; //includes vertex finding efficiency and event selection efficiency
double dy = 2.0; //should be normalized by [dy]
double phase_space_factor = 0.4540; //|y|<1,|Eta|<1 / |y|<1 starlight

double ratioIncohToDissco = 3.4;
Double_t fitExpPlusDisso(Double_t *x, Double_t *par){
  
  //par[2] leaves it blank
  double coherent = par[0]*TMath::Exp(par[1]*x[0]);
  double incoherent = (ratioIncohToDissco*par[4])*TMath::Exp(par[3]*x[0]);
  double dissoc = par[4]*TMath::Power((1+(par[5]/par[6])*x[0]), -par[6]);

  return coherent+incoherent+dissoc;

}

Double_t getChargeFFofD(Double_t *x, Double_t *par){

	double gamma_c = par[0];
	double theta_c = par[1];
	double alpha_c = par[2];
	double beta_c  = par[3];
	double momega = 0.782;
	double mphi = 1.012;
	double value = 1. /( TMath::Power(1.+gamma_c*x[0], theta_c) ) * ( 1. -  alpha_c - beta_c + alpha_c*(momega*momega/(momega*momega + x[0])) + beta_c*( mphi*mphi /(mphi*mphi+x[0]) ) ); 

	return value*value;
}

void runStep_3_final_playGround( const bool doSys_ = false ){

	gStyle->SetErrorX(0);

	TFile* file_data1 = new TFile("output-Step_2.root");
	TH1D* hInclus = (TH1D*) file_data1->Get("jpsi_total");
	TH1D* hZDC = (TH1D*) file_data1->Get("jpsi_zdc");
	
	TFile* file_raw = new TFile("output-Step_1.root");
	TH1D* RawJpsiPt2 = (TH1D*) file_raw->Get("JpsiPt2;1");
	TH1D* RawJpsiPt2_zdc = (TH1D*) file_raw->Get("JpsiPt2_zdc;1");

	TH1D* hTotal = (TH1D*) hInclus->Clone("hTotal");
	TH1D* hIncoh = (TH1D*) hZDC->Clone("hZDC");
	for(int i=0;i<hInclus->GetNbinsX();i++){

		//not use relative error for now..
		double rawyield = RawJpsiPt2->GetBinContent( i+1 );
		double rawerror = RawJpsiPt2->GetBinError( i+1 );
		double relative_err = rawerror / rawyield;
		//end relative error
		double value = hInclus->GetBinContent( i+1 );
		double error = hInclus->GetBinError( i+1 );

		error = error / (Lint*flux*Br*trigBBCTOF*dy*phase_space_factor);
		value = value / (Lint*flux*Br*trigBBCTOF*dy*phase_space_factor);

		hTotal->SetBinContent( i+1, value );
		hTotal->SetBinError( i+1, error );

		rawyield = RawJpsiPt2_zdc->GetBinContent( i+1 );
		rawerror = RawJpsiPt2_zdc->GetBinError( i+1 );
		relative_err = rawerror / rawyield;
		//end relative error
		value = hZDC->GetBinContent( i+1 );
		error = hZDC->GetBinError( i+1 );

		error = (value*relative_err) / (Lint*flux*Br*trigBBCTOF*dy*phase_space_factor);
		value = value / (Lint*flux*Br*trigBBCTOF*dy*phase_space_factor);


		hIncoh->SetBinContent( i+1, value );
		hIncoh->SetBinError( i+1, error );
		if(i==6){
			hIncoh->SetBinContent( i+1, 0 );
			hIncoh->SetBinError( i+1, 0 );
		}
	}

	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gPad->SetLogy(1);
	
	TH1D* base1 = makeHist("base1", "", "-t #approx p^{2}_{T, J/#psi} (GeV^{2})", "d#sigma/dtdy (nb/GeV^{2})", 100,0,2.5,kBlack);
	base1->GetYaxis()->SetRangeUser(0.1, 5e3);
	base1->GetXaxis()->SetRangeUser(0, 2);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.1,1.25);
	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,6,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);
	base1->Draw();

	drawAllTLatex();

	TFile* file_beagle = new TFile("eD_photo_main_Beagle.root");
	TH1D* h_truth = (TH1D*) file_beagle->Get("t_truth");
	TH1D* h_ZDC = (TH1D*) file_beagle->Get("t_ZDC");
	TH1D* h_ratio = (TH1D*) h_ZDC->Clone("h_ratio");
	h_ratio->Divide( h_truth );

	TFile* file_cgc = new TFile("/Users/kong/google_drive/BNL_folder/Work/STAR/DeuteronCGC/gamma_d_data/CGC_gammaAu.root","READ");
	TGraphErrors* coh = (TGraphErrors*) file_cgc->Get("coh_cgc");
	TGraphErrors* coh_hulthen = (TGraphErrors*) file_cgc->Get("coh_hulthen");
	TGraphErrors* incoh = (TGraphErrors*) file_cgc->Get("incoh_fluc_cgc");
	TGraphErrors* incoh_nofluc = (TGraphErrors*) file_cgc->Get("incoh_nofluc_cgc");
	TGraphErrors* total_VMCfluc = (TGraphErrors*) file_cgc->Get("total_VMCfluc");
	TGraphErrors* total_VMCnofluc = (TGraphErrors*) file_cgc->Get("total_VMCnofluc");
	TGraphErrors* total_HulthenNofluc = (TGraphErrors*) file_cgc->Get("total_HulthenNofluc");
	TGraphErrors* total_Hulthenfluc = (TGraphErrors*) file_cgc->Get("total_Hulthenfluc");
	coh->SetLineStyle(2);
	coh->SetLineWidth(4);
	coh->SetLineColor(kBlue);

	TGraphErrors* incoh_acceptanceZDC = (TGraphErrors*) incoh->Clone("incoh_acceptanceZDC");
	for(int i=0;i<incoh_acceptanceZDC->GetN();i++){
		double x,y;
		incoh_acceptanceZDC->GetPoint(i,x,y);
		int bin = h_ratio->FindBin(x);
		double accept = h_ratio->GetBinContent(bin);
		incoh_acceptanceZDC->SetPoint(i,x,y*accept);
	}

	coh_hulthen->SetLineStyle(1);
	coh_hulthen->SetLineWidth(2);
	coh_hulthen->SetLineColor(kBlue);

	incoh->SetLineStyle(1);
	incoh->SetLineColor(kYellow-2);
	incoh->SetLineWidth(2);

	incoh_nofluc->SetLineStyle(1);
	incoh_nofluc->SetLineColor(kGreen-2);
	incoh_nofluc->SetLineWidth(2);

	incoh_acceptanceZDC->SetLineStyle(2);
	incoh_acceptanceZDC->SetLineColor(kGreen-2);
	incoh_acceptanceZDC->SetLineWidth(2);

	total_VMCfluc->SetLineStyle(2);
	total_VMCfluc->SetLineColor(kBlack);
	total_VMCfluc->SetLineWidth(2);

	total_VMCnofluc->SetLineStyle(2);
	total_VMCnofluc->SetLineColor(kRed);
	total_VMCnofluc->SetLineWidth(2);

	total_HulthenNofluc->SetLineStyle(2);
	total_HulthenNofluc->SetLineColor(kBlack);
	total_HulthenNofluc->SetLineWidth(2);

	total_Hulthenfluc->SetLineStyle(2);
	total_Hulthenfluc->SetLineColor(kRed);
	total_Hulthenfluc->SetLineWidth(2);

	// coh->Draw("Lsame");
	// total_VMCfluc->Draw("Lsame");
	// total_VMCnofluc->Draw("Lsame");
	coh_hulthen->Draw("Lsame");
	incoh->Draw("Lsame");
	// incoh_nofluc->Draw("Lsame");
	// incoh_acceptanceZDC->Draw("Lsame");
	total_Hulthenfluc->Draw("Lsame");
	// total_HulthenNofluc->Draw("Lsame");

	hIncoh->SetMarkerStyle(25);
	hIncoh->SetMarkerSize(1.2);
	hIncoh->Draw("Psame");
	hTotal->SetMarkerSize(1.2);
	hTotal->Draw("Psame");

	TLegend *w5 = new TLegend(0.45,0.61,0.77,0.83);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(20);
	w5->SetTextFont(45);
	w5->AddEntry(hTotal, "Total data ", "P");
	w5->AddEntry(hIncoh, "n-tagged data ", "P");
	w5->AddEntry(coh_hulthen, "Coherent: Hulthen ", "L");
	w5->AddEntry(incoh, "Incoherent: Q_{s}/shape fluc. ", "L");
	w5->AddEntry(total_Hulthenfluc, "Total = Coh. + Incoh. ", "L");
	// w5->AddEntry(incoh_acceptanceZDC, "Incoherent * ZDC acceptance", "L");
	// w5->AddEntry(total_HulthenNofluc, "Total no shape/Q_{s} fluc. ", "L");

	// w5->AddEntry(coh, "Coherent CGC ", "L");
	// w5->AddEntry(total_VMCfluc, "Total with shape/Q_{s} fluc. ", "L");
	// w5->AddEntry(total_VMCnofluc, "Total no shape/Q_{s} fluc. ", "L");

	w5->Draw("same");

	// TCanvas* c2 = new TCanvas("c2","c2",1,1,600,600);
	// h_ratio->SetTitle("BeAGLE eD incoh. J/#psi 18x100 GeV, Q^{2} < 1");
	// h_ratio->GetYaxis()->SetTitleOffset(1.2);
	// h_ratio->GetYaxis()->SetTitle("ZDC neutron acceptance");
	// h_ratio->Draw("hist");
	// h_ratio->SetStats(kFALSE);

	// c1->Print("Preliminary/Figure_03_cgc_hulthen.pdf");



}