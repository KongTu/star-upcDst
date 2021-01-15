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

TH1D* hist = new TH1D("hist","hist",44,pt2bins_cgc);

Double_t fitTotalCrosSection(Double_t *x, Double_t *par){
  
  	int bin = hist->FindBin( x[0] );
	double histBinContent = hist->GetBinContent( bin );
	double histBinWidth = hist->GetBinWidth( bin );
	histBinContent = histBinContent / histBinWidth;
	double incoherent = par[2] * histBinContent;
	double coherent = par[0] * TMath::Exp(par[1]*x[0]);

	return coherent+incoherent;

}

Double_t fitDFF(Double_t *x, Double_t *par){
	
	int bin = hist->FindBin( x[0] );
	double histBinContent = hist->GetBinContent( bin );
	double histBinWidth = hist->GetBinWidth( bin );
	histBinContent = histBinContent / histBinWidth;
	double incoherent = par[0] * histBinContent;

	return incoherent;
}

void runStep_3_final_useCGCincoh( const int which_incoh = 0 ){

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

		error = (value*relative_err) / (Lint*flux*Br*trigBBCTOF*dy*phase_space_factor);
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
	TH1D* incoh_template = 0;
	if( which_incoh == 0 ) incoh_template = (TH1D*) file_cgc->Get("incoh_template_bin");
	if( which_incoh == 1 ) incoh_template = (TH1D*) file_cgc->Get("incoh_nofluc_template_bin");
	hist = (TH1D*) incoh_template;
    TString incoh_name = "Q_{s} fluc.";
    if( which_incoh == 1 ){
    	incoh_name = "No fluc.";
    }

    //total fit
	TF1 *func2 = new TF1("func2",fitTotalCrosSection,0,2.5,3);
	func2->SetParameter(0, 6000);
	func2->SetParameter(1, -10);
	func2->SetParameter(2, 300);

	func2->SetLineColor(kRed+1);
	func2->SetLineWidth(2);

	double fit_low = 0.;
	double fit_high = 2.0;

	//first fit
	TFitResultPtr r = hTotal->Fit("func2","RMES+","",fit_low,fit_high);
	TGraphErrors* gr[5];  
	TFitResultPtr r2[5];

	//begin bin center correction. 2 iterations are enough and stable.
	for(int j=0;j<5;j++){
		gr[j] = new TGraphErrors();
		for(int i=0;i<hTotal->GetNbinsX();i++){
			double value = 0.;
			double error = 0.;
			double binwidth = 0.;
			if(j==0){
				value = hTotal->GetBinContent( i+1 );
				error = hTotal->GetBinError( i+1 );
				binwidth = hTotal->GetBinWidth( i+1 );
			}
			else{

				double x1,y1;
				gr[j-1]->GetPoint(i, x1, y1);
				value = y1;
				error = hTotal->GetBinError( i+1 );
				binwidth = hTotal->GetBinWidth( i+1 );
			}
			double integral = func2->Integral(pt2bins[i],pt2bins[i+1]);
			double y = integral/binwidth;
			double x = func2->GetX(y);

			gr[j]->SetPoint(i, x, value);
			gr[j]->SetPointError(i, 0, error);

		}
		r2[j] = gr[j]->Fit("func2","RMES0+","",fit_low,fit_high);
	}
	base1->Draw();
	TF1* fCoh = new TF1("fCoh","[0]*TMath::Exp([1]*x[0])",0,2.0);
	fCoh->SetParameter(0, func2->GetParameter(0));
	fCoh->SetParError(0, func2->GetParError(0));
	fCoh->SetParameter(1, func2->GetParameter(1));
	fCoh->SetParError(1, func2->GetParError(1));
	fCoh->SetLineColor(kBlue);
	fCoh->SetLineWidth(3);
	fCoh->Draw("Lsame");

	TF1* fIncoh = new TF1("fIncoh",fitDFF,0,2.0,1);
	fIncoh->SetParameter(0, func2->GetParameter(2) );
	fIncoh->SetParError(0, func2->GetParError(2) );
	fIncoh->SetLineStyle(3);
	fIncoh->SetLineWidth(3);
	fIncoh->SetLineColor(kBlack);
	fIncoh->Draw("Lsame");

	TH1D *hint2 = new TH1D("hint2","Fitted gaussian with .68 conf.band", 100, 0.0,3.14);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint2, 0.68);
     //Now the "hint2" histogram has the fitted function values as the
     //bin contents and the confidence intervals as bin errors
    hint2->SetLineColor(kRed);
    hint2->SetStats(kFALSE);
    hint2->SetFillColor(kGreen-2);
    hint2->SetMarkerColor(2);
    hint2->SetFillStyle(1001);
    hint2->SetFillColorAlpha(kGreen-2,0.4);
    hint2->Draw("e3 same");

    // hTotal->Draw("Psame"); //before bin centering
    gr[4]->SetMarkerStyle(20);
    gr[4]->Draw("Psame");
	func2->Draw("Lsame");

	hIncoh->SetMarkerStyle(25);
	hIncoh->Draw("Psame");

	//all TLatex defined and drew.
	drawAllTLatex();

    TLegend *w5 = new TLegend(0.45,0.55,0.75,0.78);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(20);
	w5->SetTextFont(45);
	w5->AddEntry(hTotal, "Total data ", "P");
	w5->AddEntry(hIncoh, "n-tagged data ", "P");
	w5->AddEntry(hint2, "Total fit ", "FL");
	w5->AddEntry(fCoh, "Coherent Exponential ", "L");
	w5->AddEntry(fIncoh, "Incoherent sum "+incoh_name, "L");
	w5->Draw("same");


	//Using the fit covariance matrix to properly get all the errors
	//on the fit parameters and integral.
	TMatrixD COV = r2[4]->GetCovarianceMatrix();
    TMatrixD covMatrix_coh(2,2);
    covMatrix_coh[0][0] = COV[0][0];
    covMatrix_coh[0][1] = COV[0][1];
    covMatrix_coh[1][0] = COV[1][0];
    covMatrix_coh[1][1] = COV[1][1];
	double params_coh[2] = {fCoh->GetParameter(0),fCoh->GetParameter(1)};

	TMatrixD covMatrix_incoh(1,1);
    covMatrix_incoh[0][0] = COV[2][2];
	double params_incoh[1] = {fIncoh->GetParameter(0)};
	
	double integral_low = 0.;
	double integral_high = 1.2;

	//printing out all cross sections;
    double x[1] = {0.};
	double err[1];
	r2[4]->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, false);
    cout << "chi2 ** coherent ** incoherent ** dissoc ** ds/dt|t=0 coherent ** ds/dt|t=0 incoherent ** ds/dt|t=0 dissoc ** coherent fraction ** total cross section" << endl;
	cout << r2[4]->Chi2()/r2[4]->Ndf() << " ; ";
    cout << fCoh->Integral(integral_low,integral_high) << ";" << fCoh->IntegralError(integral_low,integral_high,params_coh,covMatrix_coh.GetMatrixArray(),1e-2) << " ; ";
    cout << fIncoh->Integral(integral_low,integral_high) <<  ";" <<  fIncoh->IntegralError(integral_low,integral_high,params_incoh,covMatrix_incoh.GetMatrixArray(),1e-2)<< " ; ";
    cout << " " <<  ";" << " " <<  " ; ";

    cout << fCoh->Eval(0) << ";" << func2->GetParError(0) <<  " ; ";
    cout << fIncoh->Eval(0) << ";" << (err[0]/func2->Eval(0)) * fIncoh->Eval(0) <<  " ; ";
    cout << " " << ";" << " " <<  " ; ";
    cout << fCoh->Integral(integral_low,integral_high) / (fCoh->Integral(integral_low,integral_high) + fIncoh->Integral(integral_low,integral_high)) <<  " ; ";
    cout << func2->Integral(integral_low,integral_high) << ";" << func2->IntegralError(integral_low,integral_high) <<  " ; ";
    cout << endl;

    // c1->Print("Preliminary/Figure_03_main_cgcIncoh.pdf");

}