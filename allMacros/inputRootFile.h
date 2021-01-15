#include "TLorentzVector.h"
#include "TVector3.h"
TString input_embedding = "../examples/dstreader/output/output_defaultBEMC_Emb_zerobias_cohNincohVM_100k_v6_trigSimuNew_SL17f_bugfix.root";
TString input_data = "../examples/dstreader/output/output_defaultBEMC_new_trigSimu.root";
const double smearing_a_para = 0.005+40*0.0001;
const double smearing_b_para = 0.008;
const double pt_indep_const_sigma = 0.0132;

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
double pt2bins_cgc[] = {0.0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.0,1.025,1.05,1.075,1.1};

Bool_t isMatch( TVector3 p1, TLorentzVector p2){

	TVector3 p2_v3 = p2.Vect();
	if( fabs(p1.Pt() - p2.Pt())/p2.Pt() > 0.3 ) return false;
	if( p1.DeltaR( p2_v3 ) > 0.1 ) return false;
	return true;
}