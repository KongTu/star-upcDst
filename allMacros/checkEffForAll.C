#include "RiceStyle.h"
#include "inputRootFile.h"
#include "TEfficiency.h"
using namespace std;
#define PI            3.1415926
#define MASS_ELECTRON 0.00051

void checkEffForAll(){

	TFile* file = new TFile("./output-PreStep_2-embedding.root");
	TH1D* hSingleTrackREC = (TH1D*) file->Get("hSingleTrackREC");
	TH1D* hSingleTrackGEN = (TH1D*) file->Get("hSingleTrackGEN");

	TH1D* h_Pm_BEMC = (TH1D*) file->Get("h_Pm_BEMC");
	TH1D* h_Pm = (TH1D*) file->Get("h_Pm");

	TH1D* hDielectronPt = (TH1D*) file->Get("hDielectronPt");
	TH1D* hMCDielectronPt = (TH1D*) file->Get("hMCDielectronPt");

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);

	TH1D* base1 = makeHist("base1", "", "p_{T}", "single track eff", 100,0,10.0,kBlack);
	base1->GetYaxis()->SetRangeUser(0., 1.0);
	base1->GetXaxis()->SetRangeUser(0.55, 2.0);
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

	TGraphAsymmErrors* gr = new TGraphAsymmErrors();
	gr->BayesDivide(hSingleTrackREC,hSingleTrackGEN);
	gr->SetMarkerStyle(20);
	gr->Draw("Psame");

	TGraphAsymmErrors* gr1 = new TGraphAsymmErrors();
	gr1->BayesDivide(h_Pm_BEMC,h_Pm);
	gr1->SetMarkerStyle(24);
	gr1->Draw("Psame");

	TLegend *w5 = new TLegend(0.28,0.21,0.53,0.3);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(20);
	w5->SetTextFont(45);
	w5->AddEntry(gr, "(Tracking + BEMC matching) eff ", "P");
	w5->AddEntry(gr1, "BEMC matching eff only ", "P");

	w5->Draw("same");

	TCanvas* c2 = new TCanvas("c2","c2",1,1,600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);

	TH1D* base2 = (TH1D*) base1->Clone("base1");
	base2->GetYaxis()->SetTitle("J/#psi eff");
	base2->GetXaxis()->SetRangeUser(0,1.5);
	base2->GetYaxis()->SetRangeUser(0,0.5);
	base2->Draw("");
	TGraphAsymmErrors* gr2 = new TGraphAsymmErrors();
	gr2->BayesDivide(hDielectronPt,hMCDielectronPt);
	gr2->SetMarkerStyle(24);
	gr2->Draw("Psame");

	TFile* file1 = new TFile("./output-Step_2.root");
	TH1D* h_pt_REC = (TH1D*) file1->Get("h_pt_REC");
	TH1D* h_pt_GEN = (TH1D*) file1->Get("h_pt_GEN");
	TH1D* h_eff = (TH1D*) h_pt_GEN->Clone("h_eff");
	for(int i=0;i<h_pt_REC->GetNbinsX();i++){
		double nn = h_pt_REC->GetBinContent(i+1);
		double dd = h_pt_GEN->GetBinContent(i+1);
		if( dd <= 0 ) continue;
		double nn_error = h_pt_REC->GetBinError(i+1);
		double dd_error = h_pt_GEN->GetBinError(i+1);

		h_eff->SetBinContent(i+1, nn/dd );
		h_eff->SetBinError(i+1, calColError(nn,dd,nn_error,dd_error) );

	}

	h_eff->SetMarkerStyle(20);
	h_eff->Draw("Psame");

	TLegend *w1 = new TLegend(0.28,0.21,0.53,0.3);
	w1->SetLineColor(kWhite);
	w1->SetFillColor(0);
	w1->SetTextSize(20);
	w1->SetTextFont(45);
	w1->AddEntry(gr2, "J/#psi eff only ", "P");
	w1->AddEntry(h_eff, "(J/#psi + vertex finding) eff ", "P");
	w1->Draw("same");


	c1->Print("~/Desktop/singleEff.pdf");
	c2->Print("~/Desktop/jpsiEff.pdf");
}