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
	double coherent = 0.;
	coherent = par[0]*TMath::Exp(par[1]*x[0]);
	double incoherent = (ratioIncohToDissco*par[4])*TMath::Exp(par[3]*x[0]);
	double dissoc = par[4]*TMath::Power((1+(par[5]/par[6])*x[0]), -par[6]);

	return coherent+incoherent+dissoc;

}

TH1D* hist = new TH1D("hist","hist",44,pt2bins_cgc);

Double_t fitExpTemplate(Double_t *x, Double_t *par){
  
  	int bin = hist->FindBin( x[0] );
	double histBinContent = hist->GetBinContent( bin );
	double histBinWidth = hist->GetBinWidth( bin );
	histBinContent = histBinContent;
	histBinContent = par[2] * histBinContent;

	double coherent = 0.;
	coherent = par[0]*TMath::Exp(par[1]*x[0]);

	return coherent+histBinContent;

}

Double_t incohTemp(Double_t *x, Double_t *par){

	double normalization = par[0];
	int bin = hist->FindBin( x[0] );
	double histBinContent = hist->GetBinContent( bin );
	double histBinWidth = hist->GetBinWidth( bin );
	histBinContent = histBinContent;

	return normalization * histBinContent;

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

void runStep_3_final( const bool doSys_ = false ){

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

	TH1D* hTotal_copy = (TH1D*) hTotal->Clone("hTotal_copy");

	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gPad->SetLogy(1);
	
	TH1D* base1 = makeHist("base1", "", "-t #approx p^{2}_{T, J/#psi} (GeV^{2})", "d#sigma/dtdy (nb/GeV^{2})", 100,0,2.5,kBlack);
	base1->GetYaxis()->SetRangeUser(0.1, 5e3);
	base1->GetXaxis()->SetRangeUser(0, 2);
	// base1->GetXaxis()->SetRangeUser(0.,1.15);
	// base1->GetYaxis()->SetRangeUser(5,1e3);
	// TH1D* base1 = makeHist("base1", "", "-t #approx p^{2}_{T, J/#psi} (GeV^{2})", "dN/dt (GeV^{-2})", 100,0,2.5,kBlack);
	// base1->GetYaxis()->SetRangeUser(5e0, 8e4);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,6,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	TF1 *func2 = new TF1("func2",fitExpPlusDisso,0.0,2.5,7);
	func2->SetParameter(0, 9000);
	func2->SetParameter(1, -15);
	func2->SetParameter(2, 27.6251);
	func2->FixParameter(3, -4.1);
	func2->SetParameter(4, 38.6);
	func2->FixParameter(5, 1.6);
	func2->FixParameter(6, 3.85);

	func2->SetLineColor(kRed+1);
	func2->SetLineWidth(2);

	// func2->SetParLimits(1, -30,-15);

	TFitResultPtr r = hTotal->Fit("func2","RMES0+","",0.0,2.5);
	r = hTotal->Fit("func2","RMES0+","",0.0,1.2);
	r = hTotal->Fit("func2","RMES0+","",0.0,1.2);
	r = hTotal->Fit("func2","RMES0+","",0.0,1.2);
	r = hTotal->Fit("func2","RMES+","",0.0,1.2);

	TGraphErrors* gr[5];  
	TFitResultPtr r2[5];

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
		r2[j] = gr[j]->Fit("func2","RMES0+","",0.0,1.2);
		r2[j] = gr[j]->Fit("func2","RMES0+","",0.0,1.2);
		r2[j] = gr[j]->Fit("func2","RMES0+","",0.0,1.2);
		r2[j] = gr[j]->Fit("func2","RMES0+","",0.0,1.2);
		r2[j] = gr[j]->Fit("func2","RMES+","",0.0,1.2);
	}

	base1->Draw();
	hIncoh->SetMarkerStyle(25);
	hIncoh->Draw("Psame");
	// hVetoZDC->Draw("Psame");
	// cout << "#Chi2 ~ " << r2[4]->Chi2()/r->Ndf() << endl;

	TF1* fCoh = new TF1("fCoh","[0]*TMath::Exp([1]*x[0])",0.0,2.0);
	fCoh->SetParameter(0, func2->GetParameter(0) );
	fCoh->SetParameter(1, func2->GetParameter(1) );
	// fCoh->SetParameter(2, func2->GetParameter(2) );
	// fCoh->SetParameter(3, func2->GetParameter(7) );
	fCoh->SetParError(0, func2->GetParError(0) );
	fCoh->SetParError(1, func2->GetParError(1) );
	// fCoh->SetParError(2, func2->GetParError(2) );
	// fCoh->SetParError(3, func2->GetParError(7) );
	fCoh->SetLineStyle(4);
	fCoh->SetLineWidth(4);
	fCoh->SetLineColor(kBlue);
	fCoh->Draw("Lsame");

	TF1* fIncoh = new TF1("fIncoh","[0]*TMath::Exp([1]*x[0])",0.0,2.0);
	
	fIncoh->SetParameter(0, func2->GetParameter(4)*ratioIncohToDissco );
	fIncoh->SetParameter(1, func2->GetParameter(3) );
	fIncoh->SetParError(0, func2->GetParError(4)*ratioIncohToDissco );
	fIncoh->SetParError(1, func2->GetParError(3) );
	fIncoh->SetLineStyle(3);
	fIncoh->SetLineWidth(4);
	fIncoh->SetLineColor(kRed);
	fIncoh->Draw("Lsame");

	TF1* fDisso = new TF1("fDisso","[0]*TMath::Power((1+([1]/[2])*x[0]), -[2])",0.0,2.0);
	fDisso->SetParameter(0, func2->GetParameter(4) );
	fDisso->SetParameter(1, func2->GetParameter(5) );
	fDisso->SetParameter(2, func2->GetParameter(6) );
	fDisso->SetParError(0, func2->GetParError(4) );
	fDisso->SetParError(1, func2->GetParError(5) );
	fDisso->SetParError(2, func2->GetParError(6) );
	fDisso->SetLineStyle(2);
	fDisso->SetLineWidth(4);
	fDisso->SetLineColor(kBlack);
	fDisso->Draw("Lsame");

	TH1D *hint2 = new TH1D("hint2","Fitted gaussian with .68 conf.band", 100, 0.0,3.14);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint2, 0.68);
     //Now the "hint2" histogram has the fitted function values as the
     //bin contents and the confidence intervals as bin errors
    hint2->SetStats(kFALSE);
    hint2->SetFillColor(kGreen-2);
    hint2->SetMarkerColor(2);
    hint2->SetFillStyle(1001);
    hint2->SetFillColorAlpha(kGreen-2,0.4);
    hint2->Draw("e3 same");

	gr[1]->SetMarkerStyle(20);
	gr[1]->Draw("Psame");

	TLatex* r42 = new TLatex(0.18, 0.91, "d+Au #sqrt{s_{_{NN}}} = 200 GeV");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);

	TLatex* r43 = new TLatex(0.6,0.91, "STAR");
	r43->SetNDC();
	r43->SetTextSize(0.04);

	// TLatex* r44 = new TLatex(0.78,0.91, "Internal");
	// r44->SetNDC();
	// r44->SetTextSize(21);
	// r44->SetTextFont(53);

	TLatex* r44 = new TLatex(0.72,0.91, "Preliminary");
	r44->SetNDC();
	r44->SetTextSize(21);
	r44->SetTextFont(53);

	TLatex* r45 = new TLatex(0.18, 0.85, "|#eta_{e}| < 1.0 ");
	r45->SetNDC();
	r45->SetTextSize(22);
	r45->SetTextFont(43);
	r45->SetTextColor(kBlack);

   	TLatex* r47 = new TLatex(0.18, 0.73, "1.0 < M_{e^{+}e^{-}} < 2.8 GeV/c^{2}");
    r47->SetNDC();
    r47->SetTextSize(22);
    r47->SetTextFont(43);
    r47->SetTextColor(kBlack);

    TLatex* r47_1 = new TLatex(0.18, 0.73, "3.3 < M_{e^{+}e^{-}} < 4.5 GeV/c^{2}");
    r47_1->SetNDC();
    r47_1->SetTextSize(22);
    r47_1->SetTextFont(43);
    r47_1->SetTextColor(kBlack);

    TLatex* r48 = new TLatex(0.6, 0.85, "#gamma*d #rightarrow J/#psi + X");
    r48->SetNDC();
    r48->SetTextSize(23);
    r48->SetTextFont(43);
    r48->SetTextColor(kBlack);

    TLatex* r49 = new TLatex(0.18, 0.85, "#LT W_{#gamma*p} #GT #approx 25 GeV");
    r49->SetNDC();
    r49->SetTextSize(23);
    r49->SetTextFont(43);
    r49->SetTextColor(kBlack);

    TLatex* r49_1 = new TLatex(0.18, 0.79, "|y_{J/#psi}| < 1");
    r49_1->SetNDC();
    r49_1->SetTextSize(23);
    r49_1->SetTextFont(43);
    r49_1->SetTextColor(kBlack);

    r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
   	// r45->Draw("same");
    // r47->Draw("same");
    r48->Draw("same");
    r49->Draw("same");
    r49_1->Draw("same");



	TLegend *w5 = new TLegend(0.6,0.55,0.85,0.78);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(20);
	w5->SetTextFont(45);
	w5->AddEntry(hTotal, "Total data ", "P");
	w5->AddEntry(func2, "Total fit ", "L");
	w5->AddEntry(fCoh, "Coherent ", "L");
	w5->AddEntry(fIncoh, "Elastic nucleon ", "L");
	w5->AddEntry(fDisso, "Nucleon disso. ", "L");
	w5->AddEntry(hIncoh, "n-tagged data ", "P");
	// w5->AddEntry(coh, "Coherent CGC ", "L");
	// w5->AddEntry(incoh, "Incoherent CGC ", "L");
	// w5->AddEntry(total_VMCfluc, "Total CGC ", "L");

	w5->Draw("same");

	TH1D* hCohData = (TH1D*) hTotal->Clone("hCohData");
	for(int i=0;i<hCohData->GetNbinsX();i++){
		double value = hCohData->GetBinContent(i+1);
		double bincenter = hCohData->GetBinCenter(i+1);
		double incoh = fIncoh->Eval( bincenter );
		double disso = fDisso->Eval( bincenter );
		double newvalue = value-incoh-disso;
		double ratio = newvalue/value;

		hCohData->SetBinContent(i+1, newvalue);
		hCohData->SetBinError(i+1, hCohData->GetBinError(i+1)*ratio);
	}

	// c1->Print("pt2-final.pdf");

    TMatrixD COV = r2[4]->GetCovarianceMatrix();
    TMatrixD covMatrix_coh(2,2);
    covMatrix_coh[0][0] = COV[0][0];
    covMatrix_coh[0][1] = COV[0][1];
    covMatrix_coh[1][0] = COV[1][0];
    covMatrix_coh[1][1] = COV[1][1];
	double params_coh[2] = {fCoh->GetParameter(0),fCoh->GetParameter(1)};
	
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
	double params_disso[3] = {fDisso->GetParameter(0),fDisso->GetParameter(1),fDisso->GetParameter(2)};

    cout << "Integral coherent ~ " << fCoh->Integral(0,1.2) << " +/- " << fCoh->IntegralError(0,1.2,params_coh,covMatrix_coh.GetMatrixArray(),1e-2) << endl;
    cout << "Integral incoherent ~ " << fIncoh->Integral(0,1.2) <<  endl;
    cout << "Integral dissoc ~ " << fDisso->Integral(0,1.2) <<  " +/- " << fDisso->IntegralError(0,1.2,params_disso,covMatrix_disso.GetMatrixArray(),1e-2) << endl;
    cout << "coherent fraction ~ " << fCoh->Integral(0,1.2) / (fCoh->Integral(0,1.2) + fIncoh->Integral(0,1.2) + fDisso->Integral(0,1.2)) << endl;
    cout << "total ~ " << func2->Integral(0,1.2) << " +/- " << func2->IntegralError(0.,2.0) << endl;

    double total_integral_error = func2->IntegralError(0,1.2);
    TFile* outputfile = 0;
    outputfile = new TFile("output-Step_3-final.root","RECREATE");

    hTotal->Write();
    func2->Write();
    fCoh->Write();
    fIncoh->Write();
    fDisso->Write();
    r2[4]->Write();
    TGraphErrors* gr1 = new TGraphErrors(1);
    gr1->SetName("gr_totalIntegral");
    gr1->SetPoint(0, func2->Integral(0,2.0), func2->IntegralError(0.,2.0));
    gr1->Write();

    // c1->Print("Preliminary/Figure_03_main_diffFit.pdf");
    c1->Print("Preliminary/Figure_03_main.pdf");
    //use another template
    TCanvas* c2 = new TCanvas("c2","c2",600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gPad->SetLogy(1);

	TFile* file_cgc = new TFile("/Users/kong/google_drive/BNL_folder/Work/STAR/DeuteronCGC/gamma_d_data/CGC_gammaAu.root");
	TH1D* incoh_template = (TH1D*) file_cgc->Get("incoh_template_bin");
	hist = (TH1D*) incoh_template;

	TF1 *func3 = new TF1("func3",fitExpTemplate,0.0,1.2,3);
	func3->SetParameter(0, 9000);
	func3->SetParameter(1, -15);
	func3->SetParameter(2, 300);
	func3->SetLineColor(kRed+1);
	func3->SetLineWidth(2);
	func3->SetParLimits(0,0.,20000);

	TFitResultPtr r1 = hTotal_copy->Fit("func3","RMES0+","",0.0,2.5);
	r1 = hTotal_copy->Fit("func3","RMES0+","",0.0,1.2);
	r1 = hTotal_copy->Fit("func3","RMES0+","",0.0,1.2);
	r1 = hTotal_copy->Fit("func3","RMES0+","",0.0,1.2);
	r1 = hTotal_copy->Fit("func3","RMES0+","",0.0,1.2);

	TH1D* base2 = (TH1D*) base1->Clone("base2");
	base2->GetXaxis()->SetRangeUser(0.,1.15);
	base2->GetYaxis()->SetRangeUser(5,1e3);
	base2->Draw();
	TF1* func3Draw = (TF1*) hTotal_copy->GetFunction("func3");
	func3Draw->SetLineColor(kRed);
	func3Draw->Draw("Lsame");
	hIncoh->Draw("Psame");
	hTotal_copy->Draw("Psame");
	TF1* fCoh_alt = new TF1("fCoh_alt","[0]*TMath::Exp([1]*x[0])",0.0,1.2);
	cout << "test " << func3->GetParameter(0) << endl;
	fCoh_alt->SetParameter(0, func3->GetParameter(0) );
	fCoh_alt->SetParameter(1, func3->GetParameter(1) );
	fCoh_alt->SetParError(0, func3->GetParError(0) );
	fCoh_alt->SetParError(1, func3->GetParError(1) );
	fCoh_alt->SetLineStyle(4);
	fCoh_alt->SetLineWidth(4);
	fCoh_alt->SetLineColor(kBlue);
	fCoh_alt->Draw("Lsame");

	TF1*fIncoh_alt = new TF1("fIncoh_alt",incohTemp,0,1.2,1);
	fIncoh_alt->SetParameter(0, func3->GetParameter(2));
	fIncoh_alt->SetParError(0, func3->GetParError(2) );
	fIncoh_alt->SetLineColor(kBlack);
	fIncoh_alt->Draw("hist ][ same");

	TH1D *hint3 = new TH1D("hint3","Fitted gaussian with .68 conf.band", 100, 0.0,2.0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint3, 0.68);
     //Now the "hint3" histogram has the fitted function values as the
     //bin contents and the confidence intervals as bin errors
    hint3->SetStats(kFALSE);
    hint3->SetFillColor(kGreen-2);
    hint3->SetMarkerColor(kGreen-2);
    hint3->SetFillStyle(1001);
    hint3->SetFillColorAlpha(kGreen-2,0.4);
    // int lastbin = hint3->FindLastBinAbove(1);
    // for(int ibin=0;ibin<hint3->GetNbinsX();ibin++){
    // 	if( ibin+1 > lastbin ){
    // 		hint3->SetBinContent(ibin+1, hint3->GetBinContent(ibin));
    // 		hint3->SetBinError(ibin+1, hint3->GetBinError(ibin));
    // 	}
    // }

    hint3->Draw("e3 same");

     r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
   	// r45->Draw("same");
    // r47->Draw("same");
    r48->Draw("same");
    r49->Draw("same");
    r49_1->Draw("same");

    TLegend *w6 = new TLegend(0.6,0.55,0.85,0.78);
	w6->SetLineColor(kWhite);
	w6->SetFillColor(0);
	w6->SetTextSize(20);
	w6->SetTextFont(45);
	w6->AddEntry(hTotal_copy, "Total data ", "P");
	w6->AddEntry(func3, "Total fit ", "L");
	w6->AddEntry(fCoh, "Coherent ", "L");
	w6->AddEntry(fIncoh_alt, "Incoherent sum ", "L");
	w6->AddEntry(hIncoh, "n-tagged data ", "P");
	// w6->AddEntry(coh, "Coherent CGC ", "L");
	// w6->AddEntry(incoh, "Incoherent CGC ", "L");
	// w6->AddEntry(total_VMCfluc, "Total CGC ", "L");

	w6->Draw("same");

	TMatrixD COV_alt = r1->GetCovarianceMatrix();
    TMatrixD covMatrix_coh_alt(2,2);
    covMatrix_coh_alt[0][0] = COV_alt[0][0];
    covMatrix_coh_alt[0][1] = COV_alt[0][1];
    covMatrix_coh_alt[1][0] = COV_alt[1][0];
    covMatrix_coh_alt[1][1] = COV_alt[1][1];
	double params_coh_alt[2] = {fCoh_alt->GetParameter(0),fCoh_alt->GetParameter(1)};
	
	cout << "chi2 ~ " << r1->Chi2()/r1->Ndf() << endl;
    cout << "Integral coherent ~ " << fCoh_alt->Integral(0,1.2) << " +/- " << fCoh_alt->IntegralError(0,1.2,params_coh_alt,covMatrix_coh_alt.GetMatrixArray(),1e-2) << endl;
    cout << "Integral incoherent ~ " << fIncoh_alt->Integral(0,1.2) <<  endl;
    cout << "coherent fraction ~ " << fCoh_alt->Integral(0,1.2) / (fCoh_alt->Integral(0,1.2) + fIncoh_alt->Integral(0,1.2)) << endl;
    cout << "total ~ " << func3->Integral(0.,1.2) << " +/- " << func3->IntegralError(0.,1.2) << endl;

    c2->Print("Preliminary/Figure_03_cgc_fit.pdf");



















    // c1->Print("Preliminary/Figure_03.pdf");

}