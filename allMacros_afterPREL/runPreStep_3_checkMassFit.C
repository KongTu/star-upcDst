#include "RiceStyle.h"
#include "TFitResult.h"
using namespace std;

#define PI            3.1415926
#define MASS_ELECTRON 0.00051

TH1D* hist = new TH1D("hist","hist",60,0.4,4);
Double_t fitJpsiMass(Double_t *x, Double_t *par){

	int bin = hist->FindBin( x[0] );
	double histBinContent = hist->GetBinContent( bin );
	double histBinWidth = hist->GetBinWidth( bin );
	histBinContent = histBinContent/histBinWidth;
	histBinContent = par[0] * histBinContent;
	/*simple exponential background*/
		// double bkgBinContent = par[1]*TMath::Exp(par[2]*x[0]);
	/* more complicated exponential background */
		double bkgBinContent = (x[0]-par[1])*TMath::Exp( par[2]*(x[0]-par[1])*(x[0]-par[3]) + par[3]*x[0]*x[0]*x[0] );
	return (histBinContent + par[4]*bkgBinContent);
}

Double_t JpsiBkg(Double_t *x, Double_t *par){

	double bkgBinContent = (x[0]-par[0])*TMath::Exp( par[1]*(x[0]-par[0])*(x[0]-par[2]) + par[2]*x[0]*x[0]*x[0] );
	return par[3]*bkgBinContent;
}

Double_t JpsiSig(Double_t *x, Double_t *par){

	double normalization = par[0];
	int bin = hist->FindBin( x[0] );
	double histBinContent = hist->GetBinContent( bin );
	double histBinWidth = hist->GetBinWidth( bin );
	histBinContent = histBinContent / histBinWidth;

	return normalization * histBinContent;

}

void runPreStep_3_checkMassFit(const bool draw_ = false){

	TFile* file_data = new TFile("output-PreStep_2-data.root");
	TFile* file_emb = new TFile("output-PreStep_2-embedding.root");
	TFile* file_mc = new TFile("output-PreStep_1.root");

	TH1D* hJpsiMass_emb = (TH1D*) file_emb->Get("hJpsiMass"); 
	TH1D* hJpsiMass_data = (TH1D*) file_data->Get("hJpsiMass");
	TH1D* hLikeSignMass_data = (TH1D*) file_data->Get("hLikeSignMass");

	for(int j=0;j<hJpsiMass_data->GetNbinsX();j++){
		double binwidth = hJpsiMass_data->GetBinWidth(j+1);
		double value = hJpsiMass_data->GetBinContent(j+1);
		double error = hJpsiMass_data->GetBinError(j+1);

		hJpsiMass_data->SetBinContent(j+1, value/binwidth);
		hJpsiMass_data->SetBinError(j+1, error/binwidth);

		binwidth = hJpsiMass_emb->GetBinWidth(j+1);
		value = hJpsiMass_emb->GetBinContent(j+1);
		error = hJpsiMass_emb->GetBinError(j+1);

		hJpsiMass_emb->SetBinContent(j+1, value/binwidth);
		hJpsiMass_emb->SetBinError(j+1, error/binwidth);
	}
	
	TCanvas* c3 = new TCanvas("c3","c3",600,600);
	int maxbin = hJpsiMass_data->FindBin(3.08);
	double maxvalue = hJpsiMass_data->GetBinContent(maxbin);
	int maxbin_emb = hJpsiMass_emb->FindBin(3.08);
	double maxvalue_emb = hJpsiMass_emb->GetBinContent(maxbin_emb);
	hJpsiMass_emb->Scale( maxvalue/maxvalue_emb );

	hJpsiMass_emb->SetTitle("mass");
	hJpsiMass_emb->SetStats(kFALSE);
	hJpsiMass_emb->GetXaxis()->SetRangeUser(2.0,3.5);
	hJpsiMass_emb->GetYaxis()->SetRangeUser(2.0,1.3*hJpsiMass_emb->GetMaximum());
	hJpsiMass_emb->GetXaxis()->SetTitle("mass (GeV/c^{2})");
	hJpsiMass_emb->SetMarkerStyle(24);
	hJpsiMass_emb->Draw("P");
	hJpsiMass_data->SetMarkerStyle(20);
	hJpsiMass_data->Draw("Psame");

	TLegend *w4 = new TLegend(0.2,0.75,0.55,0.85);
	w4->SetLineColor(kWhite);
	w4->SetFillColor(0);
	w4->SetTextSize(20);
	w4->SetTextFont(45);
	w4->AddEntry(hJpsiMass_data, "Data ", "P");
	w4->AddEntry(hJpsiMass_emb, "Embedding  ", "P");
	w4->Draw("same");

	// c3->Print("fitting/raw_mass.pdf");
	// return;

	TCanvas* c1[10][200];
	TH1D* hJpsiMass_mc[10][200];
	TF1 *fitJ[10][200];
	TF1* fitJSignal[10][200];
	TF1* fitJBkg[10][200];

	TGraphErrors* gr[10]; 
	for(int k=0;k<10;k++){
		gr[k] = new TGraphErrors();
		gr[k]->SetName(Form("gr_%d",k));
	}

	for(int j=0;j<1;j++){
	for(int i=0;i<40;i++){
		if(draw_) c1[j][i] = new TCanvas(Form("c1_%d_%d",i,j),"c1",1,1,600,600);

		hJpsiMass_mc[j][i] = (TH1D*) file_mc->Get(Form("hJpsiMass_%d_%d",i,j) );
		int mc_max = hJpsiMass_mc[j][i]->FindBin(3.08);
		int data_max = hJpsiMass_data->FindBin(3.08);
		// hJpsiMass_mc[j][i]->Scale( (hJpsiMass_data->GetBinContent(data_max)) / (1*(hJpsiMass_mc[j][i]->GetBinContent(mc_max))) );
		hJpsiMass_mc[j][i]->SetMarkerStyle(24);
		hJpsiMass_mc[j][i]->GetXaxis()->SetRangeUser(2.,3.5);
		hJpsiMass_mc[j][i]->GetYaxis()->SetRangeUser(0.1,2000);
		hJpsiMass_data->SetMarkerStyle(20);
		hJpsiMass_data->GetXaxis()->SetRangeUser(2.,3.5);
		hJpsiMass_data->GetYaxis()->SetRangeUser(0.1,2500);
		// hJpsiMass_mc[j][i]->Draw("PE");
		if( draw_) hJpsiMass_data->Draw("PE");

		hist = (TH1D*) hJpsiMass_mc[j][i]->Clone(Form("hJpsiMass_mc_template_%d_%d",i,j) );
		hist->Scale(1./(hist->Integral()));
		
		fitJ[j][i] = new TF1(Form("fitJ_%d_%d",i,j),fitJpsiMass,2.0,3.5,5);
		fitJ[j][i]->SetParameter(0,300.);
		fitJ[j][i]->SetParameter(1,1.0 );
		fitJ[j][i]->SetParameter(2,-1);
		fitJ[j][i]->SetParameter(3,0.12);
		fitJ[j][i]->SetParameter(4,1000);

		TFitResultPtr r = hJpsiMass_data->Fit(Form("fitJ_%d_%d",i,j),"S0+");
		r = hJpsiMass_data->Fit(Form("fitJ_%d_%d",i,j),"S0+");
		r = hJpsiMass_data->Fit(Form("fitJ_%d_%d",i,j),"S0+");
		r = hJpsiMass_data->Fit(Form("fitJ_%d_%d",i,j),"S0+");
		int count = 0;
		while( (!r->IsValid() || r->CovMatrixStatus()!= 3) && (count < 5) ) {
			fitJ[j][i]->SetParameter(0,r->Parameter(0));
			fitJ[j][i]->SetParameter(1,r->Parameter(1) );
			fitJ[j][i]->SetParameter(2,r->Parameter(2));
			fitJ[j][i]->SetParameter(3,r->Parameter(3));
			fitJ[j][i]->SetParameter(4,r->Parameter(4));
			r = hJpsiMass_data->Fit(Form("fitJ_%d_%d",i,j),"S0+");
			r = hJpsiMass_data->Fit(Form("fitJ_%d_%d",i,j),"S0+");
			r = hJpsiMass_data->Fit(Form("fitJ_%d_%d",i,j),"S0+");
			r = hJpsiMass_data->Fit(Form("fitJ_%d_%d",i,j),"S0+");
			count++;
		}

		TF1* drawFit = (TF1*) hJpsiMass_data->GetFunction(Form("fitJ_%d_%d",i,j));	
		drawFit->SetLineColor(kRed);
		if( draw_) drawFit->Draw("Lsame");
		double chi2 = r->Chi2();
		cout << "Fit IsValid ~ " << r->IsValid() << endl;
		cout << "Fit Error Status ~ " << r->CovMatrixStatus() << endl;
		int ndf = r->Ndf();

		fitJSignal[j][i] = new TF1(Form("fitJSignal_%d_%d",i,j),JpsiSig,2.0,3.5,1);
		fitJSignal[j][i]->SetParameter(0, fitJ[j][i]->GetParameter(0));
		fitJSignal[j][i]->SetParError(0, fitJ[j][i]->GetParError(0) );
		fitJSignal[j][i]->SetLineColor(kOrange-3);
	    fitJSignal[j][i]->SetLineWidth(1);
	    fitJSignal[j][i]->SetLineStyle(2);
	    fitJSignal[j][i]->SetFillColorAlpha(kOrange-3,0.3);
	    fitJSignal[j][i]->SetFillStyle(1001);
	    fitJSignal[j][i]->SetNpx(5000);
		if( draw_) fitJSignal[j][i]->Draw("Lsame");

		fitJBkg[j][i] = new TF1(Form("fitJBkg_%d_%d",i,j),JpsiBkg,2.0,3.5,4);
		fitJBkg[j][i]->SetParameter(0, fitJ[j][i]->GetParameter(1) );
		fitJBkg[j][i]->SetParameter(1, fitJ[j][i]->GetParameter(2) );
		fitJBkg[j][i]->SetParameter(2, fitJ[j][i]->GetParameter(3) );
		fitJBkg[j][i]->SetParameter(3, fitJ[j][i]->GetParameter(4) );
		fitJBkg[j][i]->SetLineStyle(2);
		fitJBkg[j][i]->SetLineColor(kBlue);
		if( draw_) fitJBkg[j][i]->Draw("Lsame");

		TLatex* r50 = new TLatex(0.23,0.87, Form("chi2/ndf ~ %f", chi2/ndf ));
		r50->SetNDC();
		r50->SetTextSize(18);
		r50->SetTextFont(44);
		if( draw_) r50->Draw("same");

		gr[j]->SetPoint(i,i,chi2/ndf);
		cout << "chi2/ndf ~ " << chi2/ndf << endl;
		cout << "signal yield ~ " << fitJSignal[j][i]->GetParameter(0) << " +- " << fitJSignal[j][i]->GetParError(0) << endl;
		cout << "total yield ~ " << hJpsiMass_data->Integral(hJpsiMass_data->FindBin(2.7), hJpsiMass_data->FindBin(3.3)) << endl;

		TLatex* r501 = new TLatex(0.23,0.82, Form("signal yield ~ %.1f #pm %.1f", fitJSignal[j][i]->GetParameter(0), fitJSignal[j][i]->GetParError(0)));
		r501->SetNDC();
		r501->SetTextSize(18);
		r501->SetTextFont(44);
		if( draw_) r501->Draw("same");

		if(draw_) c1[j][i]->Print(Form("massfit-scan/fit-02-12-2020_%d_%d.pdf",i,j));
	}
	}

	TCanvas* c2[10];
	TFile outputfit("massfit-scan-output.root","RECREATE");
	hJpsiMass_data->Write();
	hJpsiMass_emb->Write();
	for (int j = 0; j < 1; ++j)
	 {
	 	/* code */
	 	c2[j] = new TCanvas();
	 	gr[j]->SetMarkerStyle(20);
		gr[j]->Draw("");
		gr[j]->Write();
	 
	 } 

	

}