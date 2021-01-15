#include "RiceStyle.h"

using namespace std;

void plotEventAndTracks(){

	TFile* file_data = new TFile("output-PreStep_2-data.root");
	TFile* file_embd = new TFile("output-PreStep_2-embedding.root");

	TH1D* hist_data[10];
	TH1D* hist_embd[10];

	hist_data[0] = (TH1D*) file_data->Get("hVtxZ");
	hist_embd[0] = (TH1D*) file_embd->Get("hVtxZ");

	hist_data[1] = (TH1D*) file_data->Get("h_EoverP_BEMC");
	hist_embd[1] = (TH1D*) file_embd->Get("h_EoverP_BEMC");

	hist_data[2] = (TH1D*) file_data->Get("hElecPt");
	hist_embd[2] = (TH1D*) file_embd->Get("hElecPt");

	hist_data[3] = (TH1D*) file_data->Get("hElecEta");
	hist_embd[3] = (TH1D*) file_embd->Get("hElecEta");

	hist_data[4] = (TH1D*) file_data->Get("hElecPhi");
	hist_embd[4] = (TH1D*) file_embd->Get("hElecPhi");

	hist_data[5] = (TH1D*) file_data->Get("hElecDCAxy");
	hist_embd[5] = (TH1D*) file_embd->Get("hElecDCAxy");

	hist_data[6] = (TH1D*) file_data->Get("hElecDCAz");
	hist_embd[6] = (TH1D*) file_embd->Get("hElecDCAz");

	hist_data[7] = (TH1D*) file_data->Get("hElecNhitsFit");
	hist_embd[7] = (TH1D*) file_embd->Get("hElecNhitsFit");

	hist_data[8] = (TH1D*) file_data->Get("hElecNhitsDedx");
	hist_embd[8] = (TH1D*) file_embd->Get("hElecNhitsDedx");

	TH2D* hist_data_2D[5];
	TH2D* hist_embd_2D[5];

	hist_data_2D[0] = (TH2D*) file_data->Get("h_Nvertex");
	hist_embd_2D[0] = (TH2D*) file_embd->Get("h_Nvertex");

	hist_data_2D[1] = (TH2D*) file_data->Get("h_NsigmaPart");
	hist_embd_2D[1] = (TH2D*) file_embd->Get("h_NsigmaPart");

	TH1D* hist_data_nvertex = (TH1D*) hist_data_2D[0]->ProjectionX("hist_data_1",3,3);
	TH1D* hist_embd_nvertex = (TH1D*) hist_embd_2D[0]->ProjectionX("hist_embd_1",1,3);

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);

	hist_embd[0]->SetStats(kFALSE);
	hist_embd[0]->GetXaxis()->SetRangeUser(-200,200);
	hist_embd[0]->SetLineWidth(2);
	hist_embd[0]->DrawNormalized("");
	hist_data[0]->SetMarkerStyle(20);
	hist_data[0]->DrawNormalized("Psame");

	TLegend *w1 = new TLegend(0.6,0.75,0.78,0.86);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(20);
    w1->SetTextFont(45);
    w1->AddEntry(hist_data[0], "Data  ", "P");
    w1->AddEntry(hist_embd[0], "Embedding ", "L");
    w1->Draw("same");

    TLatex* r43 = new TLatex(0.65,0.91, "STAR");
    r43->SetNDC();
    r43->SetTextFont(62);
    r43->SetTextSize(0.04);

    TLatex* r44 = new TLatex(0.78,0.91, "Internal");
    r44->SetNDC();
    r44->SetTextSize(21);
    r44->SetTextFont(53);

    r43->Draw("same");
    r44->Draw("same");

	TCanvas* c2 = new TCanvas("c2","c2",1,1,600,600);
	gPad->SetLogy(1);
	hist_embd_nvertex->SetStats(kFALSE);
	hist_embd_nvertex->GetXaxis()->SetRangeUser(0,7);
	hist_embd_nvertex->SetLineWidth(2);
	hist_embd_nvertex->DrawNormalized("hist");
	hist_data_nvertex->SetMarkerStyle(20);
	hist_data_nvertex->DrawNormalized("Psame");

	w1->Draw("same");
	r43->Draw("same");
	r44->Draw("same");


	TCanvas* c3 = new TCanvas("c3","c3",1,1,600,600);
	gPad->SetLogy(0);
	hist_embd[1]->GetXaxis()->SetTitle("E/p ratio BEMC");
	hist_embd[1]->SetStats(kFALSE);
	hist_embd[1]->GetXaxis()->SetRangeUser(-200,200);
	hist_embd[1]->SetLineWidth(2);
	hist_embd[1]->DrawNormalized("");
	hist_data[1]->SetMarkerStyle(20);
	hist_data[1]->DrawNormalized("Psame");

	w1->Draw("same");
	r43->Draw("same");
	r44->Draw("same");

	TCanvas* c4 = new TCanvas("c4","c4",1,1,1200,600);
	c4->Divide(2,1,0.01,0.01);
	c4->cd(1);
	gPad->SetLogz(1);
	hist_data_2D[1]->SetTitle("data");
	hist_data_2D[1]->SetStats(kFALSE);
	hist_data_2D[1]->SetLineWidth(2);
	hist_data_2D[1]->Draw("colz");

	c4->cd(2);
	gPad->SetLogz(1);
	hist_embd_2D[1]->SetTitle("embedding");
	hist_embd_2D[1]->SetStats(kFALSE);
	hist_embd_2D[1]->SetLineWidth(2);
	hist_embd_2D[1]->Draw("colz");

	r43->Draw("same");
	r44->Draw("same");

	TCanvas* c5[7];
	for(int i=0;i<7;i++){
		c5[i] = new TCanvas(Form("c5_%d",i),"",1,1,600,600);
		double maxheight = hist_embd[2+i]->GetMaximum();
		hist_embd[2+i]->GetYaxis()->SetRangeUser(0,1.4* maxheight);
		hist_embd[2+i]->SetStats(kFALSE);
		hist_embd[2+i]->SetLineWidth(2);
		hist_embd[2+i]->DrawNormalized("");
		hist_data[2+i]->SetMarkerStyle(20);
		hist_data[2+i]->DrawNormalized("Psame");

		w1->Draw("same");
		r43->Draw("same");
		r44->Draw("same");
	}


	c1->Print("AN_figures/vertex_z.pdf");
	c2->Print("AN_figures/nvertex.pdf");
	c3->Print("AN_figures/EoverP.pdf");
	c4->Print("AN_figures/nsigma.pdf");
	for(int i=0;i<7;i++){
		c5[i]->Print(Form("AN_figures/tracks_%d.pdf",i));
	}























}