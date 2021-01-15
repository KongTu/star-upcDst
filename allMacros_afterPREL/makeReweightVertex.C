#include "RiceStyle.h"
#include "inputRootFile.h"
using namespace std;
#define PI            3.1415926
#define MASS_ELECTRON 0.00051

void makeReweightVertex(){

	TFile* file_data = new TFile("output-PreStep_2-data.root");
	TFile* file_embd = new TFile("output-PreStep_2-embedding.root");
	TFile* file_pt2 = new TFile("output-Step_3-final.root");

	TH1D* hTotal = (TH1D*) file_pt2->Get("hTotal");	
	TH1D* hMCpt2 = (TH1D*) file_embd->Get("hMCDielectronPt2");

	TH1D* vzdata = (TH1D*) file_data->Get("hVtxZCut");
	TH1D* vzembd = (TH1D*) file_embd->Get("hVtxZCut");

	vzdata->Scale(1./vzdata->Integral());
	vzembd->Scale(1./vzembd->Integral());

	vzdata->Divide( vzembd );
	vzdata->SetMarkerStyle(24);
	vzdata->Draw("P");

	//pt2
	for(int i=0;i<hMCpt2->GetNbinsX();i++){
		double value = hMCpt2->GetBinContent(i+1);
		double error = hMCpt2->GetBinError(i+1);
		double width = hMCpt2->GetBinWidth(i+1);

		hMCpt2->SetBinContent( i+1, value / width );
		hMCpt2->SetBinError( i+1, error / width );

	}

	hMCpt2->Scale( 1./hMCpt2->Integral() );
	hTotal->Scale( 1./hTotal->Integral() );

	hTotal->Divide( hMCpt2 );
	hTotal->Smooth(100);
	
	TFile* output = new TFile("reweight_vertex.root","RECREATE");
	vzdata->Write();
	hTotal->Write();



}