#include "RiceStyle.h"
#include "inputRootFile.h"
using namespace std;
#define PI            3.1415926
#define MASS_ELECTRON 0.00051

void makeReweightVertex(){

	TFile* file_data = new TFile("output-PreStep_2-data.root");
	TFile* file_embd = new TFile("output-PreStep_2-embedding.root");

	TH1D* vzdata = (TH1D*) file_data->Get("hVtxZCut");
	TH1D* vzembd = (TH1D*) file_embd->Get("hVtxZCut");

	vzdata->Scale(1./vzdata->Integral());
	vzembd->Scale(1./vzembd->Integral());

	vzdata->Divide( vzembd );
	vzdata->SetMarkerStyle(24);
	vzdata->Draw("P");

	TFile* output = new TFile("reweight_vertex.root","RECREATE");
	vzdata->Write();



}