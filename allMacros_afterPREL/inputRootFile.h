#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraphErrors.h"

TString input_embedding = "../examples/dstreader/output/output_defaultBEMC_Emb_zerobias_cohNincohVM_100k_v6_trigSimuNew_SL17f_bugfix.root";
TString input_data = "../examples/dstreader/output/output_defaultBEMC_new_trigSimu.root";
const double smearing_a_para = 0.005+20*0.0002;
const double smearing_b_para = 0.008;
const double pt_indep_const_sigma = 0.008;
double pt_wide_bins[]={1.0,1.2,1.4,1.6,1.8,2.0};
double pt_dep_sigma[]={0.01,0.01,0.01,0.01,0.01};
// double pt_dep_sigma[]={0.01131,0.0108,0.0108,0.010756,0.010756};
double BHrate = 0.15;
double cutNhitDedx = 15;//15
double cutNhitFit = 25;//25
double cutEtaDaug = 1.0;
double cutDCAxyz = 3.0;//3

TLorentzVector ePlus_cand;
TLorentzVector eMinus_cand;

TF1* getProb=new TF1("getProb","1",0,1);
double getBremPhotonEnergy(double elec_energy){
	//test with BH process, parameter [0] is electron energy;
	TF1* funcBH = new TF1("funcBH","100.*(1/x[0])*((4./3) - (4./3)*(x[0]/[0]) + TMath::Power(x[0]/[0],2))",0.001*elec_energy,elec_energy);
	funcBH->SetParameter(0,elec_energy);
	double energy = funcBH->GetRandom();
	return energy;
}

//default 10.
const double elecPIDcut_a = 10.;
const double elecPIDcut_b = 30.;

//define binning of histogram
double massbins[]={0.4,0.52,0.64,0.76,0.88,1,1.12,1.24,1.36,1.48,1.6,1.72,1.84,1.96,2.08,2.2,2.32,2.44,2.56,2.68,2.8,2.86,
2.92,2.98,3.04,3.10,3.16,3.22,3.28,3.4,3.52,3.64,3.76,3.88,4.0};

double ptbins[]={0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,4.0,6.0,10.0};//20 bins
double pt2binsREC[]={0.,0.05,0.1,0.15,0.3,0.5,1.2,2.0};
double pt2bins[]={0.,0.05,0.1,0.15,0.3,0.5,1.2,2.0};
double pt2bins_truth[]={0.,0.05,0.1,0.15,0.3,0.5,1.2,2.0};
double pt2bins_measu[]={0.,0.05,0.1,0.15,0.3,0.5,1.2,2.0};
double pt2bins_cgc[] = {-0.0125,0.0125,0.0375,0.0625,0.0875,0.1125,0.1375,0.1625,0.1875,0.2125,0.2375,0.2625,0.2875,0.3125,0.3375,0.3625,0.3875,0.4125,0.4375,0.4625,0.4875,0.5125,0.5375,0.5625,0.5875,0.6125,0.6375,0.6625,0.6875,0.7125,0.7375,0.7625,0.7875,0.8125,0.8375,0.8625,0.8875,0.9125,0.9375,0.9625,0.9875,1.0125,1.0375,1.0625,1.09875};

TLatex* r42[7];

Bool_t isMatch( TVector3 p1, TLorentzVector p2){

	TVector3 p2_v3 = p2.Vect();
	// if( fabs(p1.Pt() - p2.Pt())/p2.Pt() > 0.3 ) return false;
	if( p1.DeltaR( p2_v3 ) > 0.02 ) return false;
	return true;
}

void drawAllTLatex(){

	r42[0] = new TLatex(0.18, 0.91, "d+Au #sqrt{s_{_{NN}}} = 200 GeV");
	r42[0]->SetNDC();
	r42[0]->SetTextSize(22);
	r42[0]->SetTextFont(43);
	r42[0]->SetTextColor(kBlack);

	r42[1] = new TLatex(0.8,0.91, "STAR");
	r42[1]->SetNDC();
	r42[1]->SetTextSize(0.04);

	// r42[2] = new TLatex(0.72,0.91, "Preliminary");
	// r42[2]->SetNDC();
	// r42[2]->SetTextSize(21);
	// r42[2]->SetTextFont(53);

	r42[3] = new TLatex(0.63, 0.85, "#gamma* + d #rightarrow J/#psi + X");
	r42[3]->SetNDC();
	r42[3]->SetTextSize(23);
	r42[3]->SetTextFont(43);
	r42[3]->SetTextColor(kBlack);

	r42[4] = new TLatex(0.18, 0.85, "#LT W_{#gamma*p} #GT #approx 25 GeV");
	r42[4]->SetNDC();
	r42[4]->SetTextSize(23);
	r42[4]->SetTextFont(43);
	r42[4]->SetTextColor(kBlack);

	r42[5] = new TLatex(0.18, 0.79, "|y_{J/#psi}| < 1");
	r42[5]->SetNDC();
	r42[5]->SetTextSize(23);
	r42[5]->SetTextFont(43);
	r42[5]->SetTextColor(kBlack);

	for(int i=0; i<6; i++){
		if(i==2) continue;
		r42[i]->Draw("same");
	}

}

double ratioIncohToDissco = 3.4;
Double_t fitExpPlusDisso(Double_t *x, Double_t *par){
  
  //par[2] leaves it blank
	double coherent = 0.;
	coherent = par[0]*TMath::Exp(par[1]*x[0]);
	double incoherent = (ratioIncohToDissco*par[4])*TMath::Exp(par[3]*x[0]);
	double dissoc = par[4]*TMath::Power((1+(par[5]/par[6])*x[0]), -par[6]);

	return coherent+incoherent+dissoc;

}

TGraphErrors* BinCenterCorrection( TH1D* histForCorr ){

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

	TFitResultPtr r = histForCorr->Fit("func2","RMES0+","",0.0,2);

	TGraphErrors* gr = new TGraphErrors();
	for(int ibin=0;ibin<histForCorr->GetNbinsX();ibin++){

		double integral = func2->Integral(pt2binsREC[ibin],pt2binsREC[ibin+1]);
		double interval = pt2binsREC[ibin+1] - pt2binsREC[ibin];
		double average = integral / interval;
		double x = func2->GetX(average);

		gr->SetPoint(ibin, x, histForCorr->GetBinContent(ibin+1));
		gr->SetPointError(ibin, 0, histForCorr->GetBinError(ibin+1));
	}

	return gr;
}




















