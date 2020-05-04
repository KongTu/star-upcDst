#include "RiceStyle.h"
using namespace std;
#define PI 3.1415926

void runSystematics( bool draw_ = false){

	TFile* file[3];

	TString default_name = "systematics-output/electron/output-Step_3-electron_default.root";
	TString sys_name_1   = "systematics-output/electron/output-Step_3-electron_tight.root";
	TString sys_name_2   = "systematics-output/electron/output-Step_3-electron_loose.root";

	file[0] = new TFile( default_name );//default
	file[1] = new TFile( sys_name_1 );//lower
	file[2] = new TFile( sys_name_2 );//higher

	double y_min = 0.2;
	double y_max = 1.8;

	TH1D* histo[3];
	TF1* fTotal[3];
	TF1* fCoh[3];
	TF1* fIncoh[3];
	TF1* fDisso[3];
	TFitResult* rFit[3];
	TMatrixD covMatrix_coh[3];
	TMatrixD covMatrix_disso[3];
	TGraphErrors* gr[3];

	for(int ifile=0;ifile<3;ifile++){
		histo[ifile] = (TH1D*) file[ifile]->Get("hTotal");
		fTotal[ifile] = (TF1*) file[ifile]->Get("func2");
		fCoh[ifile] = (TF1*) file[ifile]->Get("fCoh");
		fIncoh[ifile] = (TF1*) file[ifile]->Get("fIncoh");
		fDisso[ifile] = (TF1*) file[ifile]->Get("fDisso");
		rFit[ifile] = (TFitResult*) file[ifile]->Get("TFitResult--func2;1");
		gr[ifile] = (TGraphErrors*) file[ifile]->Get("gr_totalIntegral");
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

	vector<double> coh_collection, incoh_collection, disso_collection;
	vector<double> coh_err_collection, incoh_err_collection, disso_err_collection;

	for(int ifile=0;ifile<3;ifile++){
		TMatrixD COV = rFit[ifile]->GetCovarianceMatrix();
	    TMatrixD covMatrix_coh(2,2);
	    covMatrix_coh[0][0] = COV[0][0];
	    covMatrix_coh[0][1] = COV[0][1];
	    covMatrix_coh[1][0] = COV[1][0];
	    covMatrix_coh[1][1] = COV[1][1];
		double params_coh[2] = {fCoh[ifile]->GetParameter(0),fCoh[ifile]->GetParameter(1)};
		
		TMatrixD covMatrix_disso(3,3);
	    covMatrix_disso[0][0] = COV[4][4];
	    covMatrix_disso[0][1] = COV[4][5];
	    covMatrix_disso[0][2] = COV[4][6];
	    covMatrix_disso[1][0] = COV[5][4];
	    covMatrix_disso[1][1] = COV[5][5];
	    covMatrix_disso[1][2] = COV[5][6];
	    covMatrix_disso[2][0] = COV[6][4];
	    covMatrix_disso[2][1] = COV[6][5];
	    covMatrix_disso[2][2] = COV[6][6];
		double params_disso[3] = {fDisso[ifile]->GetParameter(0),fDisso[ifile]->GetParameter(1),fDisso[ifile]->GetParameter(2)};
		double total_cross_section, total_cross_section_error;
		gr[ifile]->GetPoint(0,total_cross_section,total_cross_section_error);

		cout << "File = " << ifile << endl;
		cout << "Integral coherent ~ " << fCoh[ifile]->Integral(0,2.0) << " +/- " << fCoh[ifile]->IntegralError(0,2.0,params_coh,covMatrix_coh.GetMatrixArray(),1e-2) << endl;
	    cout << "Integral incoherent ~ " << fIncoh[ifile]->Integral(0,2.0) << endl;
	    cout << "Integral dissoc ~ " << fDisso[ifile]->Integral(0,2.0) << " +/- " << fDisso[ifile]->IntegralError(0,2.0,params_disso,covMatrix_disso.GetMatrixArray(),1e-2) << endl;
	    cout << "coherent fraction ~ " << fCoh[ifile]->Integral(0,2.0) / (fCoh[ifile]->Integral(0,2.0) + fIncoh[ifile]->Integral(0,2.0) + fDisso[ifile]->Integral(0,2.0)) << endl;
	    cout << "total ~ " << total_cross_section << " +/- " << total_cross_section_error << endl;
	    cout << "Slope and Intercept" << endl;
	    cout << "Coherent slope ~ " << fCoh[ifile]->GetParameter(1) << " +/- " << fCoh[ifile]->GetParError(1) << endl;
	    cout << "Coherent Intercept ~ " << fCoh[ifile]->GetParameter(0) << " +/- " << fCoh[ifile]->GetParError(0) << endl;
	    cout << "Dissoc Intercept ~ " << fDisso[ifile]->GetParameter(0) << " +/- " << fDisso[ifile]->GetParError(0) << endl;
	    cout << "Incoherent Intercept ~ " << fIncoh[ifile]->GetParameter(0) << " +/- " << fIncoh[ifile]->GetParError(0) << endl;
	    cout << "-------------------------" << endl;
	}

	double coh_integral_systematics = (fabs(fCoh[0]->Integral(0,2.0) - fCoh[1]->Integral(0,2.0))
	 + fabs(fCoh[0]->Integral(0,2.0) - fCoh[2]->Integral(0,2.0)))/2.;
	double coh_slope_systematics = (fabs(fCoh[0]->GetParameter(1) - fCoh[1]->GetParameter(1))
	 + fabs(fCoh[0]->GetParameter(1) - fCoh[2]->GetParameter(1)))/2.;
	double coh_intercept_systematics = (fabs(fCoh[0]->GetParameter(0) - fCoh[1]->GetParameter(0))
	 + fabs(fCoh[0]->GetParameter(0) - fCoh[2]->GetParameter(0)))/2.;

	cout << "Integral coherent systematics ~ " << coh_integral_systematics << " , percentage ~ " << coh_integral_systematics*0.5 / fabs(fCoh[0]->Integral(0,2.0)) << endl;
	cout << "Slope coherent systematics ~ " << coh_slope_systematics << " , percentagle ~ " << coh_slope_systematics*0.5 / -fCoh[0]->GetParameter(1) << endl;
	cout << "Intercept coherent systematics ~ " << coh_intercept_systematics << " , percentagle ~ " <<  coh_intercept_systematics*0.5 / fCoh[0]->GetParameter(0) << endl;

	double incoh_integral_systematics = (fabs(fIncoh[0]->Integral(0,2.0) - fIncoh[1]->Integral(0,2.0))
	 + fabs(fIncoh[0]->Integral(0,2.0) - fIncoh[2]->Integral(0,2.0)))/2.;
	double incoh_intercept_systematics = (fabs(fIncoh[0]->GetParameter(0) - fIncoh[1]->GetParameter(0))
	 + fabs(fIncoh[0]->GetParameter(0) - fIncoh[2]->GetParameter(0)))/2.;
	
	cout << "Integral incoherent systematics ~ " << incoh_integral_systematics << " , percentage ~ " << incoh_integral_systematics*0.5 / fabs(fIncoh[0]->Integral(0,2.0)) << endl;
	cout << "Intercept incoherent systematics ~ " << incoh_intercept_systematics << " , percentagle ~ " <<  incoh_intercept_systematics*0.5 / fIncoh[0]->GetParameter(0) << endl;

	double disso_integral_systematics = (fabs(fDisso[0]->Integral(0,2.0) - fDisso[1]->Integral(0,2.0))
	 + fabs(fDisso[0]->Integral(0,2.0) - fDisso[2]->Integral(0,2.0)))/2.;
	double disso_intercept_systematics = (fabs(fDisso[0]->GetParameter(0) - fDisso[1]->GetParameter(0))
	 + fabs(fDisso[0]->GetParameter(0) - fDisso[2]->GetParameter(0)))/2.;

	cout << "Integral dissociative systematics ~ " << disso_integral_systematics << " , percentage ~ " << disso_integral_systematics*0.5 / fabs(fDisso[0]->Integral(0,2.0)) << endl;
	cout << "Intercept dissociative systematics ~ " << disso_intercept_systematics << " , percentagle ~ " <<  disso_intercept_systematics*0.5 / fDisso[0]->GetParameter(0) << endl;

	TLatex* rsys1 = new TLatex(0.21, 0.24, Form("Coherent Integral systemtaics ~ %.3f",coh_integral_systematics*0.5 / fabs(fCoh[0]->Integral(0,2.0))) );
    rsys1->SetNDC();
    rsys1->SetTextSize(20);
    rsys1->SetTextFont(43);
    rsys1->SetTextColor(kBlack);
    rsys1->Draw("same");
    
    TLatex* rsys2 = new TLatex(0.21, 0.20, Form("Coherent Slope systemtaics ~ %.3f",coh_slope_systematics*0.5 / -fCoh[0]->GetParameter(1) ) );
    rsys2->SetNDC();
    rsys2->SetTextSize(20);
    rsys2->SetTextFont(43);
    rsys2->SetTextColor(kBlack);
    rsys2->Draw("same");

    TLatex* rsys3 = new TLatex(0.21, 0.16, Form("Coherent Intercept systemtaics ~ %.3f",coh_intercept_systematics*0.5 / fCoh[0]->GetParameter(0) ) );
    rsys3->SetNDC();
    rsys3->SetTextSize(20);
    rsys3->SetTextFont(43);
    rsys3->SetTextColor(kBlack);
    rsys3->Draw("same");

	c1->Print("systematics-figures/electron/ratio-figure-systematics.pdf");

}