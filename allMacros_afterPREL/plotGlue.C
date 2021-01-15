#include "RiceStyle.h"

using namespace std;

void plotGlue(){

	TFile* file_glue[3];
	file_glue[0] = new TFile("glue_default.root");
	file_glue[1] = new TFile("glue_plus.root");
	file_glue[2] = new TFile("glue_minus.root");

	TH1D* hGlue[3];
	for(int i=0;i<3;i++){
		hGlue[i]=(TH1D*)file_glue[i]->Get("sourceHisto");
	}
	TH1D* hCharge = (TH1D*) file_glue[0]->Get("globalLookupTable");

	TGraphAsymmErrors* gr = new TGraphAsymmErrors();
	for(int j=0;j<hGlue[0]->GetNbinsX();j++){
		double value = hGlue[0]->GetBinContent(j+1);
		double error_plus = hGlue[1]->GetBinContent(j+1);
		double error_minus = hGlue[2]->GetBinContent(j+1);
		gr->SetPoint(j, hGlue[0]->GetBinCenter(j+1), value);
		gr->SetPointEYhigh(j, fabs(error_plus-value) );
		gr->SetPointEYlow(j, fabs(value-error_minus) );

	}

	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	
	TH1D* base1 = makeHist("base1", "", "b (fm)", "F(b)", 100,-3,3,kBlack);
	base1->GetYaxis()->SetRangeUser(0.0, 1);
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
	gr->SetMarkerStyle(20);
	gr->SetFillStyle(1001);
    gr->SetFillColorAlpha(2,0.4);
    gr->Draw("e3 same");
	hCharge->Draw("Psame");

    TLegend *w4 = new TLegend(0.16,0.75,0.45,0.86);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(20);
    w4->SetTextFont(45);
    w4->AddEntry(gr, "Gluon density  ", "F");
    w4->AddEntry(hCharge, "Charge density ", "P");
    w4->Draw("same");

    TLatex* r42 = new TLatex(0.18, 0.91, "Deuteron");
    r42->SetNDC();
    r42->SetTextSize(22);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);

    TLatex* r43 = new TLatex(0.6,0.91, "STAR");
    r43->SetNDC();
    r43->SetTextFont(62);
    r43->SetTextSize(0.04);

    // TLatex* r44 = new TLatex(0.78,0.91, "Internal");
    TLatex* r44 = new TLatex(0.72,0.91, "Preliminary");
    r44->SetNDC();
    r44->SetTextSize(21);
    r44->SetTextFont(53);



    r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");

    c1->Print("Preliminary/Figure_04.pdf");


}