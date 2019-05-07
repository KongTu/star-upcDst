#include "RiceStyle.h"
#include "TLorentzVector.h"

using namespace std;

#define PI            3.1415926
#define MASS_ELECTRON 0.00051

void readUpcDstMiniTree(){

	TFile* file = new TFile("/Users/kong/google_drive/BNL_folder/Work/STAR/star-upcDst/examples/dstreader/output/output.root");
	if(!file) cout << "wrong name! Check input files" << endl;

	TH1D* hNtrk = new TH1D("hNtrk","hNtrk; Ntrk",10,0,10);
	TH1D* hJpsiMass = new TH1D("hJpsiMass","hJpsiMass; M_{e^{+}e^{-}} (GeV/c^{2})",60,0.4,4);
	TH1D* hDeltaPhi = new TH1D("hDeltaPhi","hDeltaPhi; #Delta#phi",300,-6.28,6.28);
	TH1D* hdielectronPt = new TH1D("hdielectronPt","hdielectronPt; #Dielectron p_{T}",50,0,5.0);

	TTree* tree = (TTree*) file->Get("myTree");

	Bool_t isTrigger_mini[3];
	Int_t mNumberOfTracks_mini;
	Float_t mPosZ_mini[10];
  	Int_t mNPrimaryTracks_mini[10];

	static const Int_t MAXTracks = 500;
	Bool_t mFlagBemc_mini[MAXTracks];
	Bool_t mFlagTof_mini[MAXTracks];
	Double32_t mPt_mini[MAXTracks]; // pT at point of dca to primary vertex
	Double32_t mEta_mini[MAXTracks]; // pseudorapidity at point of dca to primary vertex
	Double32_t mPhi_mini[MAXTracks]; // phi at point of dca to primary vertex

	Float_t mDcaXY_mini[MAXTracks]; // perpendicular dca to primary vertex of associated global track
	Float_t mDcaZ_mini[MAXTracks]; // longitudinal dca to primary vertex of associated global track

	Short_t mCharge_mini[MAXTracks]; // track electrical charge

	UShort_t mNhits_mini[MAXTracks]; // total number of hits on track
	UShort_t mNhitsFit_mini[MAXTracks]; // number of hits used in fit

	Double32_t mChi2_mini[MAXTracks]; // chi2 of fit

	UShort_t mNhitsDEdx_mini[MAXTracks]; // number of hits used for dE/dx measurement
	Double32_t mDEdxSignal_mini[MAXTracks]; // measured dE/dx value

	Float16_t mNSigmasTPCElectron_mini[MAXTracks]; // dE/dx n sigmas for particle species
	Float16_t mNSigmasTPCPion_mini[MAXTracks]; // dE/dx n sigmas for particle species

	tree->SetBranchAddress("isTrigger_mini",&isTrigger_mini);
	tree->SetBranchAddress("mNumberOfTracks_mini",&mNumberOfTracks_mini);
	tree->SetBranchAddress("mNPrimaryTracks_mini",&mNPrimaryTracks_mini);
	tree->SetBranchAddress("mPosZ_mini",&mPosZ_mini);
	tree->SetBranchAddress("mFlagTof_mini",&mFlagTof_mini);
	tree->SetBranchAddress("mPt_mini",&mPt_mini);
	tree->SetBranchAddress("mEta_mini",&mEta_mini);
	tree->SetBranchAddress("mPhi_mini",&mPhi_mini);
	tree->SetBranchAddress("mDcaXY_mini",&mDcaXY_mini);
	tree->SetBranchAddress("mDcaZ_mini",&mDcaZ_mini);
	tree->SetBranchAddress("mCharge_mini",&mCharge_mini);
	tree->SetBranchAddress("mNhits_mini",&mNhits_mini);
	tree->SetBranchAddress("mNhitsFit_mini",&mNhitsFit_mini);
	tree->SetBranchAddress("mChi2_mini",&mChi2_mini);
	tree->SetBranchAddress("mNhitsDEdx_mini",&mNhitsDEdx_mini);
	tree->SetBranchAddress("mDEdxSignal_mini",&mDEdxSignal_mini);
	tree->SetBranchAddress("mNSigmasTPCElectron_mini",&mNSigmasTPCElectron_mini);
	tree->SetBranchAddress("mNSigmasTPCElectron_mini",&mNSigmasTPCElectron_mini);

	cout << "========== Reading the upc-dst miniTree, total number of events = " << tree->GetEntries() << " =========== " << endl;
	
	vector< TLorentzVector> e_plus;
	vector< TLorentzVector> e_minus;
	vector<double> mNSigmaE1, mNSigmaE2;

	for(int ievt=0; ievt<tree->GetEntries(); ievt++){

		tree->GetEntry(ievt);

		// if( isTrigger_mini[0] != 1 && isTrigger_mini[1] != 1 && isTrigger_mini[2] != 1 ) continue;
		// if( mPosZ_mini[0] < 100 ) continue; //loop over vertex and cut on Z vertex

		TLorentzVector e1,e2;

		int Nparticles = 0;
		int Nparticles_plus = 0;
		int Nparticles_minus = 0;
		
		mNSigmaE1.clear();
		mNSigmaE2.clear();

		e_plus.clear();
		e_minus.clear();

		for(int itrk = 0; itrk < mNumberOfTracks_mini; itrk++){

			if( mPt_mini[itrk] < 1.1) continue;
			if( mNhitsFit_mini[itrk] < 14 ) continue;
			if( fabs(mEta_mini[itrk]) > 1.0 ) continue;
			if( fabs(mDcaXY_mini[itrk]) > 3.0 ) continue;
			if( fabs(mDcaZ_mini[itrk]) > 3.0 ) continue;
			if( mFlagBemc_mini[itrk] == 0 || mFlagTof_mini[itrk] == 0 ) continue;

			Nparticles++;

			if( mCharge_mini[itrk] == +1 ){
				Nparticles_plus++;
				e1.SetPtEtaPhiM(mPt_mini[itrk],mEta_mini[itrk],mPhi_mini[itrk],MASS_ELECTRON);
				
				e_plus.push_back( e1 );
				mNSigmaE1.push_back(mNSigmasTPCElectron_mini[itrk]);
			}

			if( mCharge_mini[itrk] == -1 ){
				Nparticles_minus++;
				e2.SetPtEtaPhiM(mPt_mini[itrk],mEta_mini[itrk],mPhi_mini[itrk],MASS_ELECTRON);
				
				e_minus.push_back( e2 );
				mNSigmaE2.push_back(mNSigmasTPCElectron_mini[itrk]);
			}

		}

		hNtrk->Fill( Nparticles );
		if( e_plus.size() < 1 || e_minus.size() < 1 ) continue;

		//unlike-sign pair
		for( unsigned itrk = 0; itrk < e_plus.size(); itrk++ ){
			for( unsigned jtrk = 0; jtrk < e_minus.size(); jtrk++ ){

				if( mNSigmaE1[itrk]*mNSigmaE1[itrk] + mNSigmaE2[jtrk]*mNSigmaE2[jtrk] > 9 ) continue;
				
				TLorentzVector e_inv;
				e_inv = e_plus[itrk]+e_minus[jtrk];
				hJpsiMass->Fill( e_inv.M() );
				
				if( e_inv.M() > 3.3 || e_inv.M() < 2.5 ) continue;
				hDeltaPhi->Fill( e_plus[itrk].Phi() - e_minus[jtrk].Phi() );
				hdielectronPt->Fill( e_inv.Pt() );

			}
		}
	}


	TFile output("upc-dst-histo.root","RECREATE");
	hJpsiMass->Write();
	hDeltaPhi->Write();
	hdielectronPt->Write();
	hNtrk->Write();


}