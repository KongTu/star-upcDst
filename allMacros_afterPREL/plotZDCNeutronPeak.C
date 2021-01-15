#include "RiceStyle.h"
#include "inputRootFile.h"
using namespace std;
#define PI 3.1415926

void plotZDCNeutronPeak(){

	TFile* file_data = new TFile("output-PreStep_2-data.root");
	TH1D* h_ZDCEast_3 = (TH1D*) file_data->Get("h_ZDCEast_3");
	TH1D* h_ZDCWest_3 = (TH1D*) file_data->Get("h_ZDCWest_3");

	//neutron peak
	TCanvas* c4 = new TCanvas("c4","c4",1,1,600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);

	TH1D* base2 = makeHist("base2", "", "ZDC ADC", "counts", 1000,0,260,kBlack);
	base2->GetYaxis()->SetRangeUser(0.1, 2.5e4);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,1.1,1.25);

	base2->GetYaxis()->SetTitleOffset(1.5);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.4);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetNdivisions(4,6,0);
	base2->GetYaxis()->SetNdivisions(3,6,0);
	TGaxis::SetMaxDigits(2);

	base2->Draw("P");
	h_ZDCWest_3->SetLineColor(kBlack);
	h_ZDCWest_3->SetLineWidth(2);
	h_ZDCEast_3->SetLineColor(kRed);
	h_ZDCEast_3->SetLineWidth(2);
	h_ZDCWest_3->Draw("HIST same");
	h_ZDCEast_3->Draw("HIST same");

	TH1D* h_ZDCWest_3_shade = (TH1D*) h_ZDCWest_3->Clone("h_ZDCWest_3_shade");
	for(int i=0;i<h_ZDCWest_3_shade->GetNbinsX();i++){
		double value = h_ZDCWest_3_shade->GetBinContent(i+1);
		double error = h_ZDCWest_3_shade->GetBinError(i+1);
		double bincenter = h_ZDCWest_3_shade->GetBinCenter(i+1);

		if( bincenter < 40 ){
			h_ZDCWest_3_shade->SetBinContent(i+1, 0);
			h_ZDCWest_3_shade->SetBinError(i+1, 0);
		}
	}
	h_ZDCWest_3_shade->SetFillStyle(1001);
    h_ZDCWest_3_shade->SetFillColorAlpha(kBlue,0.3);
    h_ZDCWest_3_shade->Draw("HIST F same");

	TLatex* r42 = new TLatex(0.20, 0.85, "d+Au #sqrt{s_{_{NN}}} = 200 GeV");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);

	TLatex* r43 = new TLatex(0.8,0.91, "STAR");
	r43->SetNDC();
	r43->SetTextSize(0.04);

	TLatex* r44 = new TLatex(0.72,0.91, "Preliminary");
	r44->SetNDC();
	r44->SetTextSize(21);
	r44->SetTextFont(53);

	TLatex* r45 = new TLatex(0.31,0.31, "Neutron peak");
	r45->SetNDC();
	r45->SetTextSize(21);
	r45->SetTextFont(43);

	r42->Draw("same");
	r43->Draw("same");
	// r44->Draw("same");
	r45->Draw("same");

	TLegend *w5 = new TLegend(0.6,0.38,0.78,0.48);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(20);
	w5->SetTextFont(45);
	w5->AddEntry(h_ZDCWest_3, "d-going ZDC ", "L");
	// w5->AddEntry(h_ZDCWest_3_shade, "Neutron peak region ", "F");
	w5->AddEntry(h_ZDCEast_3, "Au-going ZDC ", "L");
	

	w5->Draw("same");

	c4->Print("../paper/figures/Figure_02_b.pdf");
	c4->Print("../paper/figures/Figure_02_b.png");

}