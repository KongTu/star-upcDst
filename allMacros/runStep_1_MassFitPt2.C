#include "RiceStyle.h"
#include "inputRootFile.h"
using namespace std;
#define PI 3.1415926

TH1D* hist = new TH1D("hist","hist",60,0.4,4);

Double_t fitJpsiMass(Double_t *x, Double_t *par){

	int bin = hist->FindBin( x[0] );
	double histBinContent = hist->GetBinContent( bin );
	double histBinWidth = hist->GetBinWidth( bin );
	histBinContent = histBinContent/histBinWidth;
	histBinContent = par[0] * histBinContent;
	/*simple exponential background*/
		// double bkgBinContent = (par[1]+par[3])*TMath::Exp(par[2]*x[0]);
	/* more complicated exponential background */
		double bkgBinContent = (x[0]-par[1])*TMath::Exp( par[2]*(x[0]-par[1])*(x[0]-par[3]) + par[3]*x[0]*x[0]*x[0] );
	return (histBinContent + par[4]*bkgBinContent);
}

Double_t JpsiBkg(Double_t *x, Double_t *par){

	double bkgBinContent = (x[0]-par[0])*TMath::Exp( par[1]*(x[0]-par[0])*(x[0]-par[2]) + par[2]*x[0]*x[0]*x[0] );
		//simple exponential
	// double bkgBinContent = (par[0]+par[2])*TMath::Exp(par[1]*x[0]);
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

void runStep_1_MassFitPt2(const bool doSys_ = false, const bool draw_ = true){

	TFile* file_data = new TFile("output-PreStep_2-data.root");
	TH2D* hJpsiMass_Pt2 = (TH2D*) file_data->Get("hJpsiMass_Pt2");
	TH2D* hLikeSignMass_Pt2 = (TH2D*) file_data->Get("hLikeSignMass_Pt2");
	TH1D* hJpsiMass = (TH1D*) file_data->Get("hJpsiMass");
	TH1D* hLikeSignMass = (TH1D*) file_data->Get("hLikeSignMass");

	//LOAD Jpsi embedding best template
	TFile* file_mc = new TFile("output-PreStep_1.root");
	//Chi2/ndf ~ 1.1, best fit
	//a = 0.002+86*0.0001 = 0.0106
	//b = 0.008
	//systematic varations 70, 100 as opposed to 86
	TString name_of_template = "79";
	TH1D* hJpsiMass_mc = (TH1D*) file_mc->Get("hJpsiMass_"+name_of_template+"_0");

	TCanvas* c2 = new TCanvas("c2","c2",800,800);
	c2->Divide(3,3,0.01,0.01);
	int Nbins = hJpsiMass_Pt2->GetNbinsY();
	TH1D* hJpsiMass_Pt2_1D[Nbins+1];
	TH1D* hLikeSignMass_Pt2_1D[Nbins+1];
	for(int i = 0; i < Nbins; i++){
		hJpsiMass_Pt2_1D[i] = (TH1D*) hJpsiMass_Pt2->ProjectionX(Form("mass_%d",i),i+1,i+1);
		hLikeSignMass_Pt2_1D[i] = (TH1D*) hLikeSignMass_Pt2->ProjectionX(Form("likeSignMass_%d",i),i+1,i+1);
	}
	hJpsiMass_Pt2_1D[Nbins] = (TH1D*) hJpsiMass_Pt2->ProjectionX(Form("mass_%d",Nbins),1,Nbins);
	hLikeSignMass_Pt2_1D[Nbins] = (TH1D*) hLikeSignMass_Pt2->ProjectionX(Form("likeSignMass_%d",Nbins),1,Nbins);

	TF1* func[Nbins+1];
	TF1* fitJSignal[Nbins+1];
	TF1* fitJBkg[Nbins+1];

	int ptNbins = sizeof(ptbins)/sizeof(ptbins[0]) - 1;
	int pt2Nbins = sizeof(pt2bins)/sizeof(pt2bins[0]) - 1;
	TH1D* JpsiPt = new TH1D("JpsiPt","JpsiPt",ptNbins,ptbins);
	TH1D* JpsiPt2 = new TH1D("JpsiPt2","JpsiPt2",pt2Nbins,pt2bins);
	TH1D* EEPt2 = new TH1D("EEPt2","EEPt2",pt2Nbins,pt2bins);

 //Jpsi efficiency from embedding. 
	TFile * file_embed = new TFile("output-PreStep_2-embedding.root");
	TH2D* hJpsiMass_Pt2_embed = (TH2D*) file_embed->Get("hJpsiMass_Pt2");
	TH2D* hLikeSignMass_Pt2_embed = (TH2D*) file_embed->Get("hLikeSignMass_Pt2");
	TH1D* hMCDielectronPt2 = (TH1D*) file_embed->Get("hMCDielectronPt2");
	
	TH1D* hMCDielectronPt2_match = (TH1D*) file_embed->Get("hMCDielectronPt2_match");
	TH1D* hJpsiPt2_match = (TH1D*) file_embed->Get("hJpsiPt2_match");
	for(int i = 0; i < Nbins; i++){
		double binwidth = hJpsiPt2_match->GetBinWidth(i+1);
		double value = hJpsiPt2_match->GetBinContent(i+1);
		double error = hJpsiPt2_match->GetBinError(i+1);

		hJpsiPt2_match->SetBinContent(i+1, value/binwidth);
		hJpsiPt2_match->SetBinError(i+1, error/binwidth);

		binwidth = hMCDielectronPt2_match->GetBinWidth(i+1);
		value = hMCDielectronPt2_match->GetBinContent(i+1);
		error = hMCDielectronPt2_match->GetBinError(i+1);

		hMCDielectronPt2_match->SetBinContent(i+1, value/binwidth);
		hMCDielectronPt2_match->SetBinError(i+1, error/binwidth);
	}

//jpsi resolution
	TCanvas *c4 = new TCanvas("c4","c4",1,1,600,600);
	gPad->SetLogy(1);
	hJpsiPt2_match->SetStats(kFALSE);
	hJpsiPt2_match->SetTitle("");
	hJpsiPt2_match->GetXaxis()->SetTitle("p^{2}_{T} (GeV^{2})");
	hJpsiPt2_match->SetMarkerStyle(20);
	hJpsiPt2_match->Draw();
	hMCDielectronPt2_match->SetMarkerStyle(24);
	hMCDielectronPt2_match->Draw("Psame");

	TLegend *w44 = new TLegend(0.4,0.65,0.75,0.85);
	w44->SetLineColor(kWhite);
	w44->SetFillColor(0);
	w44->SetTextSize(20);
	w44->SetTextFont(45);
	w44->AddEntry(hJpsiPt2_match, "RECO matched ", "P");
	w44->AddEntry(hMCDielectronPt2_match, "GEN  ", "P");
	w44->Draw("same");
//end resolution

//jpsi bin-by-bin efficiency
	TH1D* hJpsiMass_Pt2_1D_embed[Nbins+1];
	TH1D* hLikeSignMass_Pt2_1D_embed[Nbins+1];
	TH1D* jpsiEff = new TH1D("jpsiEff","jpsiEff",pt2Nbins,pt2bins);
	for(int i = 0; i < Nbins; i++){
		hJpsiMass_Pt2_1D_embed[i] = (TH1D*) hJpsiMass_Pt2_embed->ProjectionX(Form("MCmass_%d",i),i+1,i+1);
		hLikeSignMass_Pt2_1D_embed[i] = (TH1D*) hLikeSignMass_Pt2_embed->ProjectionX(Form("MCLikeSignmass_%d",i),i+1,i+1);
		double reco = hJpsiMass_Pt2_1D_embed[i]->Integral(); // no mass window cut
		double likesign = hLikeSignMass_Pt2_1D_embed[i]->Integral();
		reco = reco - likesign;
		double gen = hMCDielectronPt2->GetBinContent(i+1);
		if( gen <= 0 ) continue;
		jpsiEff->SetBinContent(i+1, reco/gen);
		jpsiEff->SetBinError(i+1, sqrt(reco)/gen);
	}	
	jpsiEff->Fit("pol0","ER0","",0.1,3);
	TF1* fitFunc = (TF1*) jpsiEff->GetFunction("pol0");
	double eff_high = fitFunc->GetParameter(0);

//total eff:
	hJpsiMass_Pt2_1D_embed[Nbins] = (TH1D*) hJpsiMass_Pt2_embed->ProjectionX(Form("MCmass_%d",Nbins),1,Nbins);
	double reco = hJpsiMass_Pt2_1D_embed[Nbins]->Integral(); // no mass window cut
	double gen = hMCDielectronPt2->Integral();
	double eff_total = reco/gen;
//end efficiency

	for(int i = 0; i < Nbins+1; i++){
		
		c2->cd(i+1);

		TH1D* JpsiMassForFit = (TH1D*) hJpsiMass_Pt2_1D[i]->Clone(Form("JpsiMassForFit_%d",i));
		JpsiMassForFit->SetMarkerStyle(20);
		JpsiMassForFit->GetXaxis()->SetTitle("mass (GeV/c^{2})");
		TH1D* LikeSignMassForFit = (TH1D*) hLikeSignMass_Pt2_1D[i]->Clone(Form("LikeSignMassForFit_%d",i));
		LikeSignMassForFit->SetMarkerStyle(24);
		TH1D* LikeSignMassForCount = (TH1D*) LikeSignMassForFit->Clone("LikeSignMassForCount");
		if(i==Nbins) {
			JpsiMassForFit = (TH1D*) hJpsiMass->Clone(Form("JpsiMassForFit_%d",Nbins));
			LikeSignMassForFit = (TH1D*) hLikeSignMass->Clone(Form("LikeSignMassForFit_%d",Nbins));
			LikeSignMassForFit->SetMarkerStyle(24);
		}

		for(int j=0;j<JpsiMassForFit->GetNbinsX();j++){
			double binwidth = JpsiMassForFit->GetBinWidth(j+1);
			double value = JpsiMassForFit->GetBinContent(j+1);
			double error = JpsiMassForFit->GetBinError(j+1);

			JpsiMassForFit->SetBinContent(j+1, value/binwidth);
			JpsiMassForFit->SetBinError(j+1, error/binwidth);

			binwidth = LikeSignMassForFit->GetBinWidth(j+1);
			value = LikeSignMassForFit->GetBinContent(j+1);
			error = LikeSignMassForFit->GetBinError(j+1);

			LikeSignMassForFit->SetBinContent(j+1, value/binwidth);
			LikeSignMassForFit->SetBinError(j+1, error/binwidth);
		}
		
		JpsiMassForFit->SetMarkerStyle(20);
		JpsiMassForFit->GetXaxis()->SetRangeUser(2.,3.5);
		if( i<Nbins) JpsiMassForFit->GetYaxis()->SetRangeUser(0.1,40/(0.06));
		else JpsiMassForFit->GetYaxis()->SetRangeUser(0.1,180/(0.06));
		JpsiMassForFit->SetStats(kFALSE);
		if( i<Nbins) JpsiMassForFit->SetTitle(Form("p^{2}_{T} bin (%.3f-%.3f) GeV^{2}",pt2bins[i],pt2bins[i+1]));
		else {
			JpsiMassForFit->SetTitle(Form("Total yield, ave eff=%.2f",eff_total) );
			JpsiMassForFit->GetXaxis()->SetTitle("mass (GeV/c^{2})");
		}
		if( draw_) {
			JpsiMassForFit->Draw("PE");
			LikeSignMassForFit->Draw("PEsame");
		}

		hist = (TH1D*) hJpsiMass_mc->Clone(Form("hJpsiMass_mc_template_%d",i) );
		hist->Scale(1./(hist->Integral()));

		func[i] = new TF1(Form("func_%d",i),fitJpsiMass,2.0,3.5,5);
		func[i]->SetParameter(0,300.);
		func[i]->SetParameter(1,1.1 );
		func[i]->SetParameter(2,-1);
		func[i]->SetParameter(3,0.12);
		func[i]->SetParameter(4,1000);

		double fitrange_low = 2.0;
		double fitrange_high = 3.4;
		TFitResultPtr r = JpsiMassForFit->Fit(Form("func_%d",i),"S0+");
		int count = 0;
		while( (!r->IsValid() || r->CovMatrixStatus()!= 3) && count < 10 ) {
			func[i]->SetParameter(0,r->Parameter(0));
			func[i]->SetParameter(1,r->Parameter(1) );
			func[i]->SetParameter(2,r->Parameter(2));
			func[i]->SetParameter(3,r->Parameter(3));
			func[i]->SetParameter(4,r->Parameter(4));
			r = JpsiMassForFit->Fit(Form("func_%d",i),"S0+");
			count++;
		}

		TF1* drawFit = (TF1*) JpsiMassForFit->GetFunction(Form("func_%d",i));	
		drawFit->SetLineColor(kRed);
		if( draw_) drawFit->Draw("Lsame");
		double chi2 = r->Chi2();
		cout << "Fit IsValid ~ " << r->IsValid() << endl;
		cout << "Fit Error Status ~ " << r->CovMatrixStatus() << endl;
		int ndf = r->Ndf();

		fitJSignal[i] = new TF1(Form("fitJSignal_%d",i),JpsiSig,2.0,3.5,1);
		fitJSignal[i]->SetParameter(0, drawFit->GetParameter(0));
		fitJSignal[i]->SetParError(0, drawFit->GetParError(0) );
		fitJSignal[i]->SetLineColor(kOrange-3);
	    fitJSignal[i]->SetLineWidth(1);
	    fitJSignal[i]->SetLineStyle(2);
	    fitJSignal[i]->SetFillColorAlpha(kOrange-3,0.3);
	    fitJSignal[i]->SetFillStyle(1001);
	    fitJSignal[i]->SetNpx(5000);
		if( draw_) fitJSignal[i]->Draw("Lsame");

		fitJBkg[i] = new TF1(Form("fitJBkg_%d",i),JpsiBkg,2.0,3.5,4);
		fitJBkg[i]->SetParameter(0, drawFit->GetParameter(1) );
		fitJBkg[i]->SetParameter(1, drawFit->GetParameter(2) );
		fitJBkg[i]->SetParameter(2, drawFit->GetParameter(3) );
		fitJBkg[i]->SetParameter(3, drawFit->GetParameter(4) );
		fitJBkg[i]->SetLineStyle(2);
		fitJBkg[i]->SetLineColor(kBlue);
		if( draw_) fitJBkg[i]->Draw("Lsame");

		TLatex* r50 = new TLatex(0.23,0.84, Form("#chi^{2}/ndf ~ %f", chi2/ndf ));
		r50->SetNDC();
		r50->SetTextSize(16);
		r50->SetTextFont(44);
		if( draw_) r50->Draw("same");

		cout << "chi2/ndf ~ " << chi2/ndf << endl;
		cout << "signal yield ~ " << fitJSignal[i]->GetParameter(0) << " +- " << fitJSignal[i]->GetParError(0) << endl;
		cout << "total yield ~ " << JpsiMassForFit->Integral(JpsiMassForFit->FindBin(2.7), JpsiMassForFit->FindBin(3.3)) << endl;
		cout << "background yield ~ " << fitJBkg[i]->Integral(2.0,3.5) << endl;
		cout << "likesign yield ~ " << LikeSignMassForCount->Integral(LikeSignMassForFit->FindBin(2.0),LikeSignMassForFit->FindBin(3.5)) << endl;
		double ee_yield = fitJBkg[i]->Integral(2.0,3.5) - LikeSignMassForCount->Integral(LikeSignMassForFit->FindBin(2.0),LikeSignMassForFit->FindBin(3.5));
		double ee_error = sqrt(fitJBkg[i]->IntegralError(2.0,3.5)*fitJBkg[i]->IntegralError(2.0,3.5));

		EEPt2->SetBinContent(i+1, ee_yield);
		EEPt2->SetBinError(i+1, ee_error);

		TLatex* r501 = new TLatex(0.23,0.78, Form("signal yield ~ %.2f #pm %.2f", fitJSignal[i]->GetParameter(0), fitJSignal[i]->GetParError(0)));
		r501->SetNDC();
		r501->SetTextSize(16);
		r501->SetTextFont(44);
		if( draw_) r501->Draw("same");


		//yield from fitting
		double yield = fitJSignal[i]->GetParameter(0);
		double error = fitJSignal[i]->GetParError(0);		
		
		//yield from bin-counting
		TH1D* JpsiMassForCounting = (TH1D*) hJpsiMass_Pt2_1D[i]->Clone(Form("JpsiMassForCounting_%d",i));
		double total_yield_counting = JpsiMassForCounting->Integral(JpsiMassForCounting->FindBin(2.0), JpsiMassForCounting->FindBin(3.5));
		double signal_yield = total_yield_counting - fitJBkg[i]->Integral(2.0,3.5);
		double signal_error = sqrt(total_yield_counting - fitJBkg[i]->GetParError(3)*fitJBkg[i]->GetParError(3));

		JpsiPt2->SetBinContent(i+1, yield );
		JpsiPt2->SetBinError(i+1, error );

		TLegend *w4 = new TLegend(0.2,0.55,0.55,0.85);
		w4->SetLineColor(kWhite);
		w4->SetFillColor(0);
		w4->SetTextSize(20);
		w4->SetTextFont(45);
		w4->AddEntry(JpsiPt2, "unlike sign ", "P");
		w4->AddEntry(LikeSignMassForFit, "like sign  ", "P");
		c2->cd(Nbins+2) ;
		w4->Draw("same");
	}

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

	TLatex* r48 = new TLatex(0.21, 0.84, "#gamma*d #rightarrow J/#psi + X");
    r48->SetNDC();
    r48->SetTextSize(22);
    r48->SetTextFont(43);
    r48->SetTextColor(kBlack);

	TCanvas* c1_2 = new TCanvas("c1_2","c1_2",600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);

	for(int i = Nbins; i < Nbins+1; i++){
		
		c1_2->cd();

		TH1D* JpsiMassForFit = (TH1D*) hJpsiMass_Pt2_1D[i]->Clone(Form("JpsiMassForFit_%d",i));
		JpsiMassForFit->SetMarkerStyle(20);
		JpsiMassForFit->GetXaxis()->SetTitle("mass (GeV/c^{2})");
		TH1D* LikeSignMassForFit = (TH1D*) hLikeSignMass_Pt2_1D[i]->Clone(Form("LikeSignMassForFit_%d",i));
		LikeSignMassForFit->SetMarkerStyle(24);
		TH1D* LikeSignMassForCount = (TH1D*) LikeSignMassForFit->Clone("LikeSignMassForCount");
		if(i==Nbins) {
			JpsiMassForFit = (TH1D*) hJpsiMass->Clone(Form("JpsiMassForFit_%d",Nbins));
			LikeSignMassForFit = (TH1D*) hLikeSignMass->Clone(Form("LikeSignMassForFit_%d",Nbins));
			LikeSignMassForFit->SetMarkerStyle(24);
		}

		for(int j=0;j<JpsiMassForFit->GetNbinsX();j++){
			double binwidth = JpsiMassForFit->GetBinWidth(j+1);
			double value = JpsiMassForFit->GetBinContent(j+1);
			double error = JpsiMassForFit->GetBinError(j+1);

			JpsiMassForFit->SetBinContent(j+1, value/binwidth);
			JpsiMassForFit->SetBinError(j+1, error/binwidth);

			binwidth = LikeSignMassForFit->GetBinWidth(j+1);
			value = LikeSignMassForFit->GetBinContent(j+1);
			error = LikeSignMassForFit->GetBinError(j+1);

			LikeSignMassForFit->SetBinContent(j+1, value/binwidth);
			LikeSignMassForFit->SetBinError(j+1, error/binwidth);
		}
		
		JpsiMassForFit->SetMarkerStyle(20);
		JpsiMassForFit->GetXaxis()->SetRangeUser(2.,3.5);
		if( i<Nbins) JpsiMassForFit->GetYaxis()->SetRangeUser(0.1,40/(0.06));
		else JpsiMassForFit->GetYaxis()->SetRangeUser(0.1,180/(0.06));
		JpsiMassForFit->SetStats(kFALSE);
		if( i<Nbins) JpsiMassForFit->SetTitle(Form("p^{2}_{T} bin (%.3f-%.3f) GeV^{2}",pt2bins[i],pt2bins[i+1]));
		else {
			JpsiMassForFit->SetTitle( " " );
			JpsiMassForFit->GetXaxis()->SetTitle("mass (GeV/c^{2})");
			JpsiMassForFit->GetXaxis()->CenterTitle();
			JpsiMassForFit->GetXaxis()->SetTitleOffset(1.2*(JpsiMassForFit->GetXaxis()->GetTitleOffset()));
			JpsiMassForFit->GetXaxis()->SetTitleSize(1.2*(JpsiMassForFit->GetXaxis()->GetTitleSize()));
			JpsiMassForFit->GetYaxis()->SetTitleOffset(1.2*(JpsiMassForFit->GetYaxis()->GetTitleOffset()));
			JpsiMassForFit->GetYaxis()->SetTitleSize(1.2*(JpsiMassForFit->GetYaxis()->GetTitleSize()));
			JpsiMassForFit->GetYaxis()->CenterTitle();
			JpsiMassForFit->GetYaxis()->SetTitle( Form("Candidates / %.2f GeV/c^{2}", JpsiMassForFit->GetBinWidth(1)));
		}
		if( draw_) {
			JpsiMassForFit->GetYaxis()->SetNdivisions(5,5,0);
			JpsiMassForFit->Draw("PE");
			LikeSignMassForFit->Draw("PEsame");
		}

		hist = (TH1D*) hJpsiMass_mc->Clone(Form("hJpsiMass_mc_template_%d",i) );
		hist->Scale(1./(hist->Integral()));

		func[i] = new TF1(Form("func_%d",i),fitJpsiMass,2.0,3.5,5);
		func[i]->SetParameter(0,300.);
		func[i]->SetParameter(1,1.1 );
		func[i]->SetParameter(2,-1);
		func[i]->SetParameter(3,0.12);
		func[i]->SetParameter(4,1000);

		double fitrange_low = 2.0;
		double fitrange_high = 3.4;
		TFitResultPtr r = JpsiMassForFit->Fit(Form("func_%d",i),"S0+");
		int count = 0;
		while( (!r->IsValid() || r->CovMatrixStatus()!= 3) && count < 10 ) {
			func[i]->SetParameter(0,r->Parameter(0));
			func[i]->SetParameter(1,r->Parameter(1) );
			func[i]->SetParameter(2,r->Parameter(2));
			func[i]->SetParameter(3,r->Parameter(3));
			func[i]->SetParameter(4,r->Parameter(4));
			r = JpsiMassForFit->Fit(Form("func_%d",i),"S0+");
			count++;
		}

		TF1* drawFit = (TF1*) JpsiMassForFit->GetFunction(Form("func_%d",i));	
		drawFit->SetLineColor(kRed);
		drawFit->SetLineWidth(2);
		if( draw_) drawFit->Draw("Lsame");
		double chi2 = r->Chi2();
		cout << "Fit IsValid ~ " << r->IsValid() << endl;
		cout << "Fit Error Status ~ " << r->CovMatrixStatus() << endl;
		int ndf = r->Ndf();

		fitJSignal[i] = new TF1(Form("fitJSignal_%d",i),JpsiSig,2.0,3.5,1);
		fitJSignal[i]->SetParameter(0, drawFit->GetParameter(0));
		fitJSignal[i]->SetParError(0, drawFit->GetParError(0) );
		fitJSignal[i]->SetLineColor(kOrange-3);
	    fitJSignal[i]->SetLineWidth(1);
	    fitJSignal[i]->SetLineStyle(2);
	    fitJSignal[i]->SetFillColorAlpha(kOrange-3,0.3);
	    fitJSignal[i]->SetFillStyle(1001);
	    fitJSignal[i]->SetNpx(5000);
		if( draw_) fitJSignal[i]->Draw("Lsame");

		fitJBkg[i] = new TF1(Form("fitJBkg_%d",i),JpsiBkg,2.0,3.5,4);
		fitJBkg[i]->SetParameter(0, drawFit->GetParameter(1) );
		fitJBkg[i]->SetParameter(1, drawFit->GetParameter(2) );
		fitJBkg[i]->SetParameter(2, drawFit->GetParameter(3) );
		fitJBkg[i]->SetParameter(3, drawFit->GetParameter(4) );
		fitJBkg[i]->SetLineStyle(2);
		fitJBkg[i]->SetLineWidth(2);
		fitJBkg[i]->SetLineColor(kBlue);
		if( draw_) fitJBkg[i]->Draw("Lsame");

		TLatex* r50 = new TLatex(0.57,0.84, Form("#chi^{2}/ndf = %.2f", chi2/ndf ));
		r50->SetNDC();
		r50->SetTextSize(20);
		r50->SetTextFont(44);
		if( draw_) r50->Draw("same");

		cout << "chi2/ndf ~ " << chi2/ndf << endl;
		cout << "signal yield ~ " << fitJSignal[i]->GetParameter(0) << " +- " << fitJSignal[i]->GetParError(0) << endl;
		cout << "total yield ~ " << JpsiMassForFit->Integral(JpsiMassForFit->FindBin(2.7), JpsiMassForFit->FindBin(3.3)) << endl;
		cout << "background yield ~ " << fitJBkg[i]->Integral(2.0,3.5) << endl;
		cout << "likesign yield ~ " << LikeSignMassForCount->Integral(LikeSignMassForFit->FindBin(2.0),LikeSignMassForFit->FindBin(3.5)) << endl;
		

		TLatex* r501 = new TLatex(0.57,0.78, Form("signal yield = %.0f #pm %.0f", fitJSignal[i]->GetParameter(0), fitJSignal[i]->GetParError(0)));
		r501->SetNDC();
		r501->SetTextSize(20);
		r501->SetTextFont(44);
		if( draw_) r501->Draw("same");

		TLegend *w4 = new TLegend(0.16,0.65,0.55,0.75);
		w4->SetLineColor(kWhite);
		w4->SetFillColor(0);
		w4->SetTextSize(20);
		w4->SetTextFont(45);
		w4->AddEntry(JpsiPt2, "unlike sign ", "P");
		w4->AddEntry(LikeSignMassForFit, "like sign  ", "P");
		c1_2->cd() ;
		w4->Draw("same");

		TLegend *w5 = new TLegend(0.18,0.45,0.55,0.59);
		w5->SetLineColor(kWhite);
		w5->SetFillColor(0);
		w5->SetTextSize(20);
		w5->SetTextFont(45);
		w5->AddEntry(drawFit, "Total fit ", "L");
		w5->AddEntry(fitJSignal[i], "Signal  ", "SL");
		w5->AddEntry(fitJBkg[i], "Background  ", "L");
		c1_2->cd() ;
		w5->Draw("same");

		r42->Draw("same");
		r43->Draw("same");
		r44->Draw("same");
		r48->Draw("same");
	}

	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gPad->SetLogy(1);
	
	TH1D* base1 = makeHist("base1", "", "-t #approx p^{2}_{T, J/#psi} (GeV^{2})", "dN/dt (GeV^{-2})", 100,0,3,kBlack);
	base1->GetYaxis()->SetRangeUser(0.1, 3e5);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,6,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	base1->Draw("");

	// make_dNdX(JpsiPt2);
	JpsiPt2->SetMarkerStyle(20);
	JpsiPt2->SetLineColor(kBlack);
	TF1* expoFit = new TF1("expoFit","[1]*TMath::Exp([0]*x[0])",0,0.2);
	JpsiPt2->Fit("expoFit","0","",0,0.2);
	JpsiPt2->Draw("Psame");

	EEPt2->SetMarkerStyle(24);
	EEPt2->Draw("Psame");

	TLatex* r45 = new TLatex(0.18, 0.80, "|#eta_{e}| < 1.0, p_{T,e} > 0.5 GeV/c");
	r45->SetNDC();
	r45->SetTextSize(22);
	r45->SetTextFont(43);
	r45->SetTextColor(kBlack);

   	TLatex* r47 = new TLatex(0.21, 0.78, "|y_{e^{+}e^{-}}| < 1.0");
    r47->SetNDC();
    r47->SetTextSize(22);
    r47->SetTextFont(43);
    r47->SetTextColor(kBlack);

    r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
   	// r45->Draw("same");
    r47->Draw("same");
    r48->Draw("same");
    // r49->Draw("same");
    // r50->Draw("same");

    TLegend *w4 = new TLegend(0.6,0.67,0.8,0.82);
	w4->SetLineColor(kWhite);
	w4->SetFillColor(0);
	w4->SetTextSize(20);
	w4->SetTextFont(45);

	TString outputname = "output-Step_1.root";
	if( doSys_ ) outputname = "systematics-input/electron/output-Step_1-electron_"+name_of_template+".root";
	TFile outputfile(outputname,"RECREATE");
	JpsiPt2->Write();

	c2->Print("systematics-figures/electron/fit-signal-template_"+name_of_template+".pdf");


}