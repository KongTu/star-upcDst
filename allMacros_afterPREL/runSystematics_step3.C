#include "RiceStyle.h"
using namespace std;
#define PI 3.1415926

void runSystematics_step3( TString sys_name = "unfolding", bool draw_ = true){

	TFile* file[3];

	TString default_name = "systematics-input/"+sys_name+"/output-Step_2_default.root";
	TString sys_name_1   = "systematics-input/"+sys_name+"/output-Step_2_TUnfold.root";
	TString sys_name_2   = "systematics-input/"+sys_name+"/output-Step_2_weight.root";

	file[0] = new TFile( default_name );//default
	file[1] = new TFile( sys_name_1 );//lower
	file[2] = new TFile( sys_name_2 );//higher

	double y_min = 0.2;
	double y_max = 1.8;

	TH1D* histo[3];
	for(int ifile=0;ifile<1;ifile++){
		histo[ifile] = (TH1D*) file[ifile]->Get("jpsi_total");
	}
	for(int ifile=1;ifile<3;ifile++){
		histo[ifile] = (TH1D*) file[ifile]->Get("jpsi_total");
	}

	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	// gPad->SetLogy(1);
	
	TH1D* base1 = makeHist("base1", "", "-t #approx p^{2}_{T, J/#psi} (GeV^{2})", "ratio of dN/dt", 100,0,2,kBlack);
	base1->GetYaxis()->SetRangeUser(y_min, y_max);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.1,1.25);
	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,6,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	TH1D* h_ratio_sys_1 = make_systematicRatio(histo[1],histo[0]);
	TH1D* h_ratio_sys_2 = make_systematicRatio(histo[2],histo[0]);

	h_ratio_sys_1->SetMarkerStyle(24);
	h_ratio_sys_2->SetMarkerStyle(25);
	h_ratio_sys_1->Fit("pol0");
	h_ratio_sys_2->Fit("pol0");
	TF1* fit1 = (TF1*) h_ratio_sys_1->GetFunction("pol0");
	fit1->SetLineColor(kRed);
	TF1* fit2 = (TF1*) h_ratio_sys_2->GetFunction("pol0");
	fit2->SetLineColor(kBlue);

	base1->Draw("");
	fit1->Draw("same");
	fit2->Draw("same");
	h_ratio_sys_1->Draw("Psame");
	h_ratio_sys_2->Draw("Psame");

 	TLegend *w4 = new TLegend(0.6,0.67,0.8,0.82);
	w4->SetLineColor(kWhite);
	w4->SetFillColor(0);
	w4->SetTextSize(20);
	w4->SetTextFont(45);
	w4->AddEntry(h_ratio_sys_1, "sys_1/default");
	w4->AddEntry(fit1, "pol0 fit");
	w4->AddEntry(h_ratio_sys_2, "sys_2/default");
	w4->AddEntry(fit2, "pol0 fit");
	w4->Draw("same");

	TLatex* r42 = new TLatex(0.18, 0.91, "dAu #sqrt{s_{_{NN}}} = 200 GeV");
	r42->SetNDC();
	r42->SetTextSize(23);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);

	TLatex* r43 = new TLatex(0.67,0.91, "STAR");
	r43->SetNDC();
	r43->SetTextSize(0.04);

	TLatex* r44 = new TLatex(0.78,0.91, "Internal");
	r44->SetNDC();
	r44->SetTextSize(21);
	r44->SetTextFont(53);

	TLatex* r45 = new TLatex(0.21, 0.84, "#gamma*d #rightarrow J/#psi + X");
    r45->SetNDC();
    r45->SetTextSize(22);
    r45->SetTextFont(43);
    r45->SetTextColor(kBlack);

	r42->Draw("same");
	r43->Draw("same");
	r44->Draw("same");
	r45->Draw("same");

	if(draw_) c1->Print("systematics-figures/"+sys_name+"/ratio-figure-systematics.pdf");

}