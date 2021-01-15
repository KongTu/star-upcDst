#include "RiceStyle.h"
#include "inputRootFile.h"
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include <fstream>
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

TH1D* hist = new TH1D("hist","hist",44,pt2bins_cgc);
TH1D* hist_incoh = new TH1D("hist_incoh","hist_incoh",44,pt2bins_cgc);

Double_t fitTotalCrosSection(Double_t *x, Double_t *par){
  
  	int bin = hist->FindBin( x[0] );
	double histBinContent = hist->GetBinContent( bin );
	double histBinWidth = hist->GetBinWidth( bin );
	histBinContent = histBinContent / histBinWidth;
	double coherent = par[0] * histBinContent;

	double incoherent = hist_incoh->GetBinContent( bin );
	incoherent = par[1] * (incoherent / hist_incoh->GetBinWidth( bin ));

	return coherent+incoherent;

}

Double_t fitCoh(Double_t *x, Double_t *par){
	
	int bin = hist->FindBin( x[0] );
	double histBinContent = hist->GetBinContent( bin );
	double histBinWidth = hist->GetBinWidth( bin );
	histBinContent = histBinContent / histBinWidth;
	double coherent = par[0] * histBinContent;

	return coherent;
}
Double_t fitIncoh(Double_t *x, Double_t *par){
	
	int bin = hist_incoh->FindBin( x[0] );
	double histBinContent = hist_incoh->GetBinContent( bin );
	double histBinWidth = hist_incoh->GetBinWidth( bin );
	histBinContent = histBinContent / histBinWidth;
	double incoherent = par[0] * histBinContent;

	return incoherent;
}


void runStep_3_final_useCGCcohAndIncoh( const int which_coh = 0, const int which_incoh = 0 ){

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
	base1->GetXaxis()->SetRangeUser(0.,2);
	base1->GetYaxis()->SetRangeUser(0.1, 5e3);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.1,1.25);
	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,6,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	//template comes from here;
	TFile* file_cgc = new TFile("/Users/kong/google_drive/BNL_folder/Work/STAR/DeuteronCGC/gamma_d_data/CGC_gammaAu.root");
	TH1D* coh_template = 0;
	TH1D* incoh_template = 0;
	if( which_coh == 0 ) coh_template = (TH1D*) file_cgc->Get("coh_template_bin");
	if( which_coh == 1 ) coh_template = (TH1D*) file_cgc->Get("coh_hulthen_template_bin");
	if( which_incoh == 0 ) incoh_template = (TH1D*) file_cgc->Get("incoh_template_bin");
	if( which_incoh == 1 ) incoh_template = (TH1D*) file_cgc->Get("incoh_nofluc_template_bin");

	hist = (TH1D*) coh_template;
	hist_incoh = (TH1D*) incoh_template;
    TString coh_name = "VMC";
    TString incoh_name = "Q_{s} fluc.";
    if( which_coh == 1 ){
    	coh_name = "Hulthen";
    }
    if( which_incoh == 1){
    	incoh_name = "no fluc.";
    }

    //total fit
	TF1 *funcCGC = new TF1("funcCGC",fitTotalCrosSection,0,1.1,2);
	funcCGC->SetParameter(0, 0.1);
	funcCGC->SetParameter(1, 0.3);
	funcCGC->SetLineColor(kRed+1);
	funcCGC->SetLineWidth(2);

	double fit_low = 0.;
	double fit_high = 1.2;

	//bin center correction
	TGraphErrors* total_data = BinCenterCorrection(hTotal);
	//first fit
	TFitResultPtr r_fit = hTotal->Fit("funcCGC","RMES0+","",fit_low,fit_high);
	r_fit = hTotal->Fit("funcCGC","RMES0+","",fit_low,fit_high);
	r_fit = hTotal->Fit("funcCGC","RMES0+","",fit_low,fit_high);
	r_fit = hTotal->Fit("funcCGC","RMES0+","",fit_low,fit_high);

	base1->Draw();
	TF1* fCoh = new TF1("fCoh",fitCoh,0,1.1,1);
	fCoh->SetParameter(0, funcCGC->GetParameter(0));
	fCoh->SetLineColor(kBlue);
	fCoh->SetLineWidth(3);
	fCoh->Draw("Lsame");

	TF1* fIncoh = new TF1("fIncoh",fitIncoh,0,1.1,1);
	fIncoh->SetParameter(0, funcCGC->GetParameter(1));
	fIncoh->SetLineColor(kBlack);
	fIncoh->SetLineWidth(3);
	fIncoh->SetLineStyle(2);
	fIncoh->Draw("Lsame");

	TH1D *hint2 = new TH1D("hint2","Fitted gaussian with .68 conf.band", 100, 0.0,1.1);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint2, 0.68);
     //Now the "hint2" histogram has the fitted function values as the
     //bin contents and the confidence intervals as bin errors
    hint2->SetLineColor(kRed);
    hint2->SetStats(kFALSE);
    hint2->SetFillColor(kGreen-2);
    hint2->SetMarkerColor(2);
    hint2->SetFillStyle(1001);
    hint2->SetFillColorAlpha(kGreen-2,0.4);
    hint2->GetXaxis()->SetRangeUser(0,1.1);
    hint2->Draw("e3 same");

    // hTotal->Draw("Psame"); //before bin centering
    total_data->SetMarkerStyle(20);
    total_data->Draw("Psame");
	funcCGC->Draw("L][same");

	hIncoh->SetMarkerStyle(25);
	hIncoh->Draw("Psame");

	//all TLatex defined and drew.
	drawAllTLatex();

    TLegend *w5 = new TLegend(0.45,0.6,0.8,0.78);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(20);
	w5->SetTextFont(45);
	w5->AddEntry(hTotal, "Total data ", "P");
	w5->AddEntry(hIncoh, "n-tagged data ", "P");
	w5->AddEntry(hint2, "Total fit ", "FL");
	w5->AddEntry(fCoh, "Coherent "+coh_name, "L");
	w5->AddEntry(fIncoh, "Incoherent sum "+incoh_name, "L");
	w5->Draw("same");

	//Using the fit covariance matrix to properly get all the errors
	//on the fit parameters and integral.
	TMatrixD COV = r_fit->GetCovarianceMatrix();
    TMatrixD covMatrix_coh(1,1);
    covMatrix_coh[0][0] = COV[0][0];
	double params_coh[1] = {fCoh->GetParameter(0)};

	TMatrixD incovMatrix_coh(1,1);
    incovMatrix_coh[0][0] = COV[1][1];
	double params_incoh[1] = {fIncoh->GetParameter(0)};

	double integral_low = 0;
	double integral_high = 1.2;
	
    double x[1] = {0.};
	double err[1];
	r_fit->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, false);
    //printing out all cross sections;
    cout << "chi2 ** coherent ** incoherent ** dissoc ** ds/dt|t=0 coherent ** ds/dt|t=0 incoherent ** ds/dt|t=0 dissoc ** coherent fraction ** total cross section" << endl;
	cout << r_fit->Chi2()/r_fit->Ndf() << " ; ";
    cout << fCoh->Integral(integral_low,integral_high) << ";" << fCoh->IntegralError(integral_low,integral_high,params_coh,covMatrix_coh.GetMatrixArray(),1e-2) << " ; ";
    cout << fIncoh->Integral(integral_low,integral_high) <<  ";" << fIncoh->IntegralError(integral_low,integral_high,params_incoh,incovMatrix_coh.GetMatrixArray(),1e-2) << " ; ";
    cout << " " <<  ";" << " " << " ; ";

    cout << fCoh->Eval(0) << ";" << (err[0]/funcCGC->Eval(0)) * fCoh->Eval(0) <<  " ; ";
    cout << fIncoh->Eval(0) << ";" << (err[0]/funcCGC->Eval(0)) * fIncoh->Eval(0) <<  " ; ";
    cout << " " << ";" << " " <<  " ; ";
    cout << fCoh->Integral(integral_low,integral_high) / (fCoh->Integral(integral_low,integral_high) + fIncoh->Integral(integral_low,integral_high)) <<  " ; ";
    cout << funcCGC->Integral(integral_low,integral_high) << ";" << funcCGC->IntegralError(integral_low,integral_high) <<  " ; ";
    cout << endl;

    double total_xSection_value = funcCGC->Integral(integral_low,integral_high); double total_xSection_error = funcCGC->IntegralError(integral_low,integral_high);
    double total_xSection_coh = fCoh->Integral(integral_low,integral_high); double total_xSection_coh_error = fCoh->IntegralError(integral_low,integral_high,params_coh,covMatrix_coh.GetMatrixArray(),1e-2);
    double total_xSection_incoh = fIncoh->Integral(integral_low,integral_high); double total_xSection_incoh_error = fIncoh->IntegralError(integral_low,integral_high,params_incoh,incovMatrix_coh.GetMatrixArray(),1e-2);
    double coh_t0 = fCoh->Eval(0); double coh_t0_error = (err[0]/funcCGC->Eval(0)) * fCoh->Eval(0);
    double chi2 = r_fit->Chi2()/r_fit->Ndf();

    string string_total_xSection = std::to_string(total_xSection_value); string string_total_xSection_error = std::to_string(total_xSection_error);
    string string_total_xSection_coh = std::to_string(total_xSection_coh); string string_total_xSection_coh_error = std::to_string(total_xSection_coh_error);
    string string_total_xSection_incoh = std::to_string(total_xSection_incoh); string string_total_xSection_incoh_error = std::to_string(total_xSection_incoh_error);
    string string_coh_t0 = std::to_string(coh_t0); string string_coh_t0_error = std::to_string(coh_t0_error);
    string string_chi2 = std::to_string(chi2); 

    cout << endl;
    cout << "& " << std::setprecision(3) << total_xSection_value << " $\\pm$ " << total_xSection_error << " $\\pm$ "<<total_xSection_value*0.186<<"; & " << total_xSection_coh << " $\\pm$ " << total_xSection_coh_error << " $\\pm$ "<<total_xSection_coh*0.186<<"; & " << total_xSection_incoh << " $\\pm$ " << total_xSection_incoh_error << " $\\pm$ "<<total_xSection_incoh*0.186<<"; & " << coh_t0 << " $\\pm$ " << coh_t0_error << " $\\pm$ "<<coh_t0*0.186<<"; & " << chi2 << " \\\\" << endl;

    // c1->Print("Preliminary/Figure_03_main_new_cgc"+coh_name+"_"+incoh_name+".pdf");

}