TString input_embedding = "../examples/dstreader/output/output_defaultBEMC_Emb_zerobias_cohNincohVM_100k_v6_trigSimuNew_SL17f_bugfix.root";
TString input_data = "../examples/dstreader/output/output_defaultBEMC_new_trigSimu.root";
const double smearing_a_para = 0.005+78*0.0001;
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
