#include "RiceStyle.h"
#include "inputRootFile.h"
using namespace std;
#define PI            3.1415926
#define MASS_ELECTRON 0.00051

void runPreStep_1_smearMass(const bool doSmearing_ = true){

	int massNbins = sizeof(massbins)/sizeof(massbins[0]) - 1;
	
	TFile* file = 0;
	file = new TFile( input_embedding );
	if(!file) cout << "wrong name! Check input files" << endl;

	TH2D* hPtSmear2D = new TH2D("hPtSmear2D","hPtSmear2D",10,0,2,1000,-0.1,0.1);
	TH1D* hPtRes1D = new TH1D("hPtRes1D","hPtRes1D",1000,-0.1,0.1);

	TH1D* hJpsiMass[200][200];
	for(int j = 0; j < 200; j++){
		for(int i = 0; i < 200; i++){
			hJpsiMass[i][j] = new TH1D(Form("hJpsiMass_%d_%d",i,j),"hJpsiMass; M_{e^{+}e^{-}} (GeV/c^{2})",massNbins,massbins);
		} 
	}

	TTree* tree = (TTree*) file->Get("myTree");

	Bool_t isTrigger_mini[4];
	Int_t mNumberOfTracks_mini;
	Int_t mNvertex_mini;
	Float_t mPosZ_mini[30];
	Float_t mPosX_mini[30];
	Float_t mPosY_mini[30];
  	Int_t mNPrimaryTracks_mini[30];

	static const Int_t mNZdcPmt = 3;
	UShort_t mZdcEastADC_mini[mNZdcPmt]; // ZDC 3 PMT ADCs, east
	UShort_t mZdcWestADC_mini[mNZdcPmt]; // ZDC 3 PMT ADCs, west
	static const Int_t mNZdcSmd = 16; // number of ZDC SMD channels, east and west
	UShort_t mZdcSmdEast_mini[mNZdcSmd]; // ZDC SMD data east
	UShort_t mZdcSmdWest_mini[mNZdcSmd]; // ZDC SMD data west

	static const Int_t MAXTracks = 500;
	UInt_t mVtxId_mini[MAXTracks];
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
	Float_t mTofTime_mini[MAXTracks]; 
	Float16_t mNSigmasTPCElectron_mini[MAXTracks]; // dE/dx n sigmas for particle species
	Float16_t mNSigmasTPCPion_mini[MAXTracks]; // dE/dx n sigmas for particle species
	
	Double32_t mBemcPmag_mini[MAXTracks];
	Float_t mBemcEta_mini[MAXTracks];
	Float_t mBemcPhi_mini[MAXTracks];
	Float_t mBemcEnergy_mini[MAXTracks]; // energy of matched BEMC cluster
	Float_t mBemcHTEnergy_mini[MAXTracks]; // energy of matched hot tower
	Float_t mBemcClsEta_mini[MAXTracks]; // eta of matched BEMC cluster
	Float_t mBemcClsPhi_mini[MAXTracks]; // phi of matched BEMC cluster

	//MC particles
	
	Int_t mMCnTracks_mini = 0;
	Double32_t mMC_px_mini[MAXTracks];
	Double32_t mMC_py_mini[MAXTracks];
	Double32_t mMC_pz_mini[MAXTracks];
	Double32_t mMC_E_mini[MAXTracks];
  	Int_t mMC_pdg_mini[MAXTracks];

	tree->SetBranchAddress("isTrigger_mini",&isTrigger_mini);
	tree->SetBranchAddress("mNumberOfTracks_mini",&mNumberOfTracks_mini);
	tree->SetBranchAddress("mNvertex_mini",&mNvertex_mini);
	tree->SetBranchAddress("mNPrimaryTracks_mini",&mNPrimaryTracks_mini);
	tree->SetBranchAddress("mPosX_mini",&mPosX_mini);
	tree->SetBranchAddress("mPosY_mini",&mPosY_mini);
	tree->SetBranchAddress("mPosZ_mini",&mPosZ_mini);
	
	tree->SetBranchAddress("mZdcEastADC_mini",&mZdcEastADC_mini);
	tree->SetBranchAddress("mZdcWestADC_mini",&mZdcWestADC_mini);
	tree->SetBranchAddress("mZdcSmdEast_mini",&mZdcSmdEast_mini);
	tree->SetBranchAddress("mZdcSmdWest_mini",&mZdcSmdWest_mini);

	tree->SetBranchAddress("mFlagTof_mini",&mFlagTof_mini);
	tree->SetBranchAddress("mFlagBemc_mini",&mFlagBemc_mini);
	tree->SetBranchAddress("mVtxId_mini",&mVtxId_mini);
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
	tree->SetBranchAddress("mTofTime_mini",&mTofTime_mini);
	tree->SetBranchAddress("mNSigmasTPCElectron_mini",&mNSigmasTPCElectron_mini);
	tree->SetBranchAddress("mNSigmasTPCPion_mini",&mNSigmasTPCPion_mini);
	tree->SetBranchAddress("mBemcHTEnergy_mini",&mBemcHTEnergy_mini);
	tree->SetBranchAddress("mBemcEta_mini",&mBemcEta_mini);
	tree->SetBranchAddress("mBemcPhi_mini",&mBemcPhi_mini);
	tree->SetBranchAddress("mBemcPmag_mini",&mBemcPmag_mini);
	tree->SetBranchAddress("mBemcEnergy_mini",&mBemcEnergy_mini);
	tree->SetBranchAddress("mBemcClsEta_mini",&mBemcClsEta_mini);
	tree->SetBranchAddress("mBemcClsPhi_mini",&mBemcClsPhi_mini);
	
	tree->SetBranchAddress("mMCnTracks_mini",&mMCnTracks_mini);
	tree->SetBranchAddress("mMC_px_mini",&mMC_px_mini);
	tree->SetBranchAddress("mMC_py_mini",&mMC_py_mini);
	tree->SetBranchAddress("mMC_pz_mini",&mMC_pz_mini);
	tree->SetBranchAddress("mMC_E_mini",&mMC_E_mini);
	tree->SetBranchAddress("mMC_pdg_mini",&mMC_pdg_mini);

	cout << "========== Reading the upc-dst miniTree, total number of events = " << tree->GetEntries() << " =========== " << endl;
	
	vector< TLorentzVector> e_plus;
	vector< TLorentzVector> e_minus;
	vector<double> mNSigmaE1, mNSigmaE2, mNSigmaPi1, mNSigmaPi2;
	
	// double pt_wide_bins[]={1.0,1.2,1.4,1.6,1.8,2.0};
	// double pt_dep_sigma[]={0.01131,0.0108,0.0108,0.010756,0.010756};
	double pt_indep_sigma = pt_indep_const_sigma;//pt_independent smearing width for now. Change it with more Embedding events;

	for(int j = 0; j < 1; j++){
	for(int incre = 0; incre < 120; incre++){
		//smearing parameter:
		double a = 0.005+incre*0.0001;
		double b = 0.008;
		cout << "iteration ~ " << incre << endl;
		cout << "a ~ " << a << endl;
		cout << "b ~ " << b << endl;

		for(int ievt=0; ievt<tree->GetEntries(); ievt++){
			tree->GetEntry(ievt);
			for(int ivtx = 0; ivtx < mNvertex_mini; ivtx++){
				if(TMath::Abs(mPosX_mini[ivtx])<1.e-5 && TMath::Abs(mPosY_mini[ivtx])<1.e-5 && TMath::Abs(mPosZ_mini[ivtx])<1.e-5) continue;
				if(fabs(mPosZ_mini[ivtx]) > 100. ) continue; //loop over vertex and cut on Z vertex

				//MC particles:
				vector< TLorentzVector> e_MC_plus, e_MC_minus;
				e_MC_plus.clear();e_MC_minus.clear();
				for(int imc = 0; imc < mMCnTracks_mini; imc++){
					TLorentzVector eMCParticle(mMC_px_mini[imc], mMC_py_mini[imc],mMC_pz_mini[imc],mMC_E_mini[imc]);
					if( TMath::Abs(mMC_pdg_mini[imc]) != 11 ) continue;
					if( TMath::Abs(eMCParticle.Eta()) > 1.0 ) continue;
					if( mMC_pdg_mini[imc] == -11 ) e_MC_plus.push_back( eMCParticle );
					if( mMC_pdg_mini[imc] == +11 ) e_MC_minus.push_back( eMCParticle );
				}
				if( e_MC_plus.size() < 1 || e_MC_minus.size() < 1 ) continue;

				TLorentzVector e1,e2;
				mNSigmaE1.clear();
				mNSigmaE2.clear();
				mNSigmaPi1.clear();
				mNSigmaPi2.clear();
				e_plus.clear();
				e_minus.clear();

				int Nparticles = 0;
				for(int itrk = 0; itrk < mNumberOfTracks_mini; itrk++){

					if( (int) mVtxId_mini[itrk] != ivtx ) continue;
					if( mNhitsDEdx_mini[itrk] < 10 ) continue;
					if( mNhitsFit_mini[itrk] < 13 ) continue;
					if( fabs(mEta_mini[itrk]) > 1.0 ) continue;
					if( fabs(mDcaXY_mini[itrk]) > 3.0 ) continue;
					if( fabs(mDcaZ_mini[itrk]) > 3.0 ) continue;
					if( !mFlagBemc_mini[itrk] ) continue;				

					if( mCharge_mini[itrk] == +1 ){
						double pt_smearing = (mPt_mini[itrk] - e_MC_plus[0].Pt())/e_MC_plus[0].Pt();
						double Nsigma = pt_smearing/pt_indep_sigma;
						double sigma_new = sqrt( (a*e_MC_plus[0].Pt())*(a*e_MC_plus[0].Pt()) + b*b);
						if(doSmearing_) mPt_mini[itrk] = e_MC_plus[0].Pt()*(sigma_new*Nsigma+1);

					}
					if( mCharge_mini[itrk] == -1 ){
						double pt_smearing = (mPt_mini[itrk] - e_MC_minus[0].Pt())/e_MC_minus[0].Pt();
						double Nsigma = pt_smearing/pt_indep_sigma;
						double sigma_new = sqrt( (a*e_MC_minus[0].Pt())*(a*e_MC_minus[0].Pt()) + b*b);
						if(doSmearing_) mPt_mini[itrk] = e_MC_minus[0].Pt()*(sigma_new*Nsigma+1);

					}
					if( mPt_mini[itrk] < 0.5 ) continue;

					Nparticles++;

					if( mCharge_mini[itrk] == +1 ){
						e1.SetPtEtaPhiM(mPt_mini[itrk],mEta_mini[itrk],mPhi_mini[itrk],MASS_ELECTRON);
						e_plus.push_back( e1 );
						mNSigmaE1.push_back(mNSigmasTPCElectron_mini[itrk]);
						mNSigmaPi1.push_back(mNSigmasTPCPion_mini[itrk]);
					}

					if( mCharge_mini[itrk] == -1 ){
						e2.SetPtEtaPhiM(mPt_mini[itrk],mEta_mini[itrk],mPhi_mini[itrk],MASS_ELECTRON);
						e_minus.push_back( e2 );
						mNSigmaE2.push_back(mNSigmasTPCElectron_mini[itrk]);
						mNSigmaPi2.push_back(mNSigmasTPCPion_mini[itrk]);
					}
				}

				//at least 2 tracks events are required
				if( Nparticles < 2 ) continue;
				if( e_MC_plus.size() == 1 && e_plus.size() == 1 ){
					double reso = (e_plus[0].Pt() - e_MC_plus[0].Pt())/(e_MC_plus[0].Pt());
					hPtSmear2D->Fill(e_MC_plus[0].Pt(), reso );
					hPtRes1D->Fill( reso );
				}
				if( e_MC_minus.size() == 1 && e_minus.size() == 1 ){
					double reso = (e_minus[0].Pt() - e_MC_minus[0].Pt())/(e_MC_minus[0].Pt());
					hPtSmear2D->Fill(e_MC_minus[0].Pt(), reso );
					hPtRes1D->Fill( reso );
				}

				//unlike-sign pair loop
				for( unsigned itrk = 0; itrk < e_plus.size(); itrk++ ){
					for( unsigned jtrk = 0; jtrk < e_minus.size(); jtrk++ ){

						double chi2_e = mNSigmaE1[itrk]*mNSigmaE1[itrk] + mNSigmaE2[jtrk]*mNSigmaE2[jtrk];
						double chi2_p = mNSigmaPi1[itrk]*mNSigmaPi1[itrk] + mNSigmaPi2[jtrk]*mNSigmaPi2[jtrk];
				
						// if( chi2_p<30 && (chi2_e/chi2_p) > 1./3 ) continue;
						// if( chi2_p>30 && chi2_e > 10. ) continue;

						TLorentzVector e_inv;
						e_inv = e_plus[itrk]+e_minus[jtrk];
						if( fabs(e_inv.Rapidity()) > 1.0 ) continue;
						hJpsiMass[incre][j]->Fill( e_inv.M() );
						
					}
				}
			}
		}

	}
	}

	TString outname;
	outname = "output-PreStep_1.root";

	TFile output(outname,"RECREATE");
	hPtRes1D->Write();
	hPtSmear2D->Write();
	for(int j = 0; j < 1; j++){
		for(int i = 0; i < 120; i++){
			hJpsiMass[i][j]->Write();
		}
	}


}