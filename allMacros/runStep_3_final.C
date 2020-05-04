#include "RiceStyle.h"
#include "inputRootFile.h"
using namespace std;
#define PI 3.1415926

/*
- d\sigma/dt 
cross section relation:
(1/flux) * Nobs / [ Lint * Br * dt * (A x e) * trig * event * dy ]
*/

double Lint = 93.418*0.934/2.5;// nb**-1
double flux = 11.78; //given by starlight, k(dn/dk), where k is photon energy ~ M/2 Exp(-y) 
double Br = 0.0593; // ee
double trigBBCTOF = 0.8*0.98*0.96; //only TOFmult < 2 and > 6 efficiency and BBC veto efficiency
double event = 0.2156; //0.2156; //includes vertex finding efficiency and event selection efficiency
double dy = 2.0; //should be normalized by [dy]

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

void runStep_3_final( const bool doSys_ = false ){

	gStyle->SetErrorX(0);

	TFile* file_data1 = new TFile("output-Step_2.root");
	TH1D* hInclus = (TH1D*) file_data1->Get("hMeasu_hTruth;1");
	
	TFile* file_raw = new TFile("output-Step_1.root");
	TH1D* RawJpsiPt2 = (TH1D*) file_raw->Get("JpsiPt2;1");

	TFile* file_emb = new TFile("output-PreStep_2-embedding.root");
	TH1D* hMCDielectronPt2 = (TH1D*) file_emb->Get("hMCDielectronPt2TrigEff");
	TH1D* hMCDielectronPt2Trig = (TH1D*) file_emb->Get("hMCDielectronPt2Trig");
	TH1D* trigBEMCpt2 = (TH1D*) hMCDielectronPt2Trig->Clone("trigBEMCpt2");
	trigBEMCpt2->Divide( hMCDielectronPt2 );
	TFile* output = new TFile("./trigBEMC.root","RECREATE");
	trigBEMCpt2->Write();

	TH1D* hTotal = (TH1D*) hInclus->Clone("hTotal");
	for(int i=0;i<hInclus->GetNbinsX();i++){

		//not use relative error for now..
		double rawyield = RawJpsiPt2->GetBinContent( i+1 );
		double rawerror = RawJpsiPt2->GetBinError( i+1 );
		double relative_err = rawerror / rawyield;
		//end relative error

		double value = hInclus->GetBinContent( i+1 );
		double error = hInclus->GetBinError( i+1 );
		double trigEffBEMCbinBybin = trigBEMCpt2->GetBinContent( i+1 );

		error = (error) / (Lint*flux*Br*trigBBCTOF*event*dy*trigEffBEMCbinBybin);
		value = value / (Lint*flux*Br*trigBBCTOF*event*dy*trigEffBEMCbinBybin);

		hTotal->SetBinContent( i+1, value );
		hTotal->SetBinError( i+1, error );

	}

	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gPad->SetLogy(1);
	
	TH1D* base1 = makeHist("base1", "", "-t #approx p^{2}_{T, J/#psi} (GeV^{2})", "d#sigma/dt (nb/GeV^{2})", 100,0,2.5,kBlack);
	base1->GetYaxis()->SetRangeUser(0.5, 5e3);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,6,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	TF1 *func2 = new TF1("func2",fitExpPlusDisso,0,2.5,7);
	func2->SetParameter(0, 9000);
	func2->SetParameter(1, -12);
	//func2->SetParameter(2, 100);
	func2->FixParameter(3, -4.1);
	func2->SetParameter(4, 38.6);
	func2->FixParameter(5, 1.6);
	func2->FixParameter(6, 3.85);

	func2->SetLineColor(kRed+1);
	func2->SetLineWidth(2);

	TFitResultPtr r = hTotal->Fit("func2","S0+","",0,2.5);
	r = hTotal->Fit("func2","S0+","",0,2.0);
	r = hTotal->Fit("func2","S0+","",0,2.0);
	r = hTotal->Fit("func2","S0+","",0,2.0);
	r = hTotal->Fit("func2","S0+","",0,2.0);

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
		r2[j] = gr[j]->Fit("func2","S0+","",0,2.0);
		r2[j] = gr[j]->Fit("func2","S0+","",0,2.0);
		r2[j] = gr[j]->Fit("func2","S0+","",0,2.0);
		r2[j] = gr[j]->Fit("func2","S0+","",0,2.0);
		r2[j] = gr[j]->Fit("func2","S+","",0,2.0);
	}

	base1->Draw();
	// hVetoZDC->Draw("Psame");
	cout << "#Chi2 ~ " << r2[4]->Chi2()/r->Ndf() << endl;

	TF1* fCoh = new TF1("fCoh","[0]*TMath::Exp([1]*x[0])",0,2.0);
	fCoh->SetParameter(0, func2->GetParameter(0) );
	fCoh->SetParameter(1, func2->GetParameter(1) );
	fCoh->SetParError(0, func2->GetParError(0) );
	fCoh->SetParError(1, func2->GetParError(1) );
	fCoh->SetLineStyle(4);
	fCoh->SetLineWidth(4);
	fCoh->SetLineColor(kBlue);
	fCoh->Draw("Lsame");

	TF1* fIncoh = new TF1("fIncoh","[0]*TMath::Exp([1]*x[0])",0,2.0);
	
	fIncoh->SetParameter(0, func2->GetParameter(4)*ratioIncohToDissco );
	fIncoh->SetParameter(1, func2->GetParameter(3) );
	fIncoh->SetParError(0, func2->GetParError(4)*ratioIncohToDissco );
	fIncoh->SetParError(1, func2->GetParError(3) );
	fIncoh->SetLineStyle(3);
	fIncoh->SetLineWidth(4);
	fIncoh->SetLineColor(kRed);
	fIncoh->Draw("Lsame");

	TF1* fDisso = new TF1("fDisso","[0]*TMath::Power((1+([1]/[2])*x[0]), -[2])",0,2.0);
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

	TLatex* r42 = new TLatex(0.18, 0.91, "dAu #sqrt{s_{_{NN}}} = 200 GeV");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);

	TLatex* r43 = new TLatex(0.67,0.91, "STAR");
	r43->SetNDC();
	r43->SetTextSize(0.04);

	TLatex* r44 = new TLatex(0.78,0.91, "Internal");
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

    TLatex* r48 = new TLatex(0.6, 0.85, "#gamma*d #rightarrow J/#psi + Y");
    r48->SetNDC();
    r48->SetTextSize(23);
    r48->SetTextFont(43);
    r48->SetTextColor(kBlack);

    TLatex* r49 = new TLatex(0.18, 0.85, "#LT W_{#gamma*d} #GT #approx 35 GeV");
    r49->SetNDC();
    r49->SetTextSize(23);
    r49->SetTextFont(43);
    r49->SetTextColor(kBlack);

    TLatex* r49_1 = new TLatex(0.18, 0.79, "#LT W_{#gamma*p} #GT #approx 25 GeV");
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

    TLegend *w5 = new TLegend(0.58,0.51,0.83,0.7);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(20);
	w5->SetTextFont(45);
	w5->AddEntry(hTotal, "Unfolded data ", "P");
	w5->AddEntry(func2, "Total fit ", "L");
	w5->AddEntry(fCoh, "Coherent ", "L");
	w5->AddEntry(fIncoh, "Incoherent ", "L");
	w5->AddEntry(fDisso, "Dissociative ", "L");
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

    cout << "Integral coherent ~ " << fCoh->Integral(0,2.0) << " +/- " << fCoh->IntegralError(0,2.0,params_coh,covMatrix_coh.GetMatrixArray(),1e-2) << endl;
    cout << "Integral incoherent ~ " << fIncoh->Integral(0,2.0) <<  endl;
    cout << "Integral dissoc ~ " << fDisso->Integral(0,2.0) <<  " +/- " << fDisso->IntegralError(0,2.0,params_disso,covMatrix_disso.GetMatrixArray(),1e-2) << endl;
    cout << "coherent fraction ~ " << fCoh->Integral(0,2.0) / (fCoh->Integral(0,2.0) + fIncoh->Integral(0,2.0) + fDisso->Integral(0,2.0)) << endl;
    cout << "total ~ " << func2->Integral(0,2.0) << " +/- " << func2->IntegralError(0.,2.0) << endl;

    double total_integral_error = func2->IntegralError(0.,2.0);
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


}