#include "RiceStyle.h"
#include "inputRootFile.h"
using namespace std;
#define PI            3.1415926
#define MASS_ELECTRON 0.00051

struct MyEvent {

	enum {
	nMCtrack_MAX=100
	};
	// if there is no MC info, nMCtrack is set to zero
	Int_t mMCnOS_tiny;
	Int_t mMCnTrack_tiny;
	Double32_t mMC_px_tiny[nMCtrack_MAX];
	Double32_t mMC_py_tiny[nMCtrack_MAX];
	Double32_t mMC_pz_tiny[nMCtrack_MAX];
	Double32_t mMC_E_tiny[nMCtrack_MAX];

	enum {
	nRECtrack_MAX=100
	};
	Int_t eventPass_tiny;
	Int_t mRECnTracks_tiny;
	Int_t mRECnSS_tiny;
	Int_t mRECnOS_tiny;
	Double32_t mREC_SS_px_tiny[nRECtrack_MAX];
	Double32_t mREC_SS_py_tiny[nRECtrack_MAX];
	Double32_t mREC_SS_pz_tiny[nRECtrack_MAX];
	Double32_t mREC_SS_E_tiny[nRECtrack_MAX];
	Double32_t mREC_OS_px_tiny[nRECtrack_MAX];
	Double32_t mREC_OS_py_tiny[nRECtrack_MAX];
	Double32_t mREC_OS_pz_tiny[nRECtrack_MAX];
	Double32_t mREC_OS_E_tiny[nRECtrack_MAX];

};

void runPreStep_4_tinyTree( const bool doEmb_ = false ){

	TFile* file = 0;
	if(doEmb_){file = new TFile( input_embedding );}
	else{file = new TFile( input_data );}
	if(!file) cout << "wrong name! Check input files" << endl;

	TString outfile_name = "output-PreStep_4-data.root";
	if(doEmb_) outfile_name = "output-PreStep_4-embedding.root";

	TFile *outfile = new TFile( outfile_name, "RECREATE");
	TTree *tinyTree = new TTree("tinyTree","tinyTree");
	MyEvent myEvent;
	tinyTree->Branch("mMCnTrack_tiny",&myEvent.mMCnTrack_tiny,"mMCnTrack_tiny/I");
	tinyTree->Branch("mMCnOS_tiny",&myEvent.mMCnOS_tiny,"mMCnOS_tiny/I");
	tinyTree->Branch("mMC_px_tiny",myEvent.mMC_px_tiny,"mMC_px_tiny[mMCnOS_tiny]/D");
	tinyTree->Branch("mMC_py_tiny",myEvent.mMC_py_tiny,"mMC_py_tiny[mMCnOS_tiny]/D");
	tinyTree->Branch("mMC_pz_tiny",myEvent.mMC_pz_tiny,"mMC_pz_tiny[mMCnOS_tiny]/D");
	tinyTree->Branch("mMC_E_tiny",myEvent.mMC_E_tiny,"mMC_E_tiny[mMCnOS_tiny]/D");
	
	tinyTree->Branch("eventPass_tiny",&myEvent.eventPass_tiny,"eventPass_tiny/I");
	tinyTree->Branch("mRECnTracks_tiny",&myEvent.mRECnTracks_tiny,"mRECnTracks_tiny/I");
	tinyTree->Branch("mRECnSS_tiny",&myEvent.mRECnSS_tiny,"mRECnSS_tiny/I");
	tinyTree->Branch("mRECnOS_tiny",&myEvent.mRECnOS_tiny,"mRECnOS_tiny/I");
	tinyTree->Branch("mREC_SS_px_tiny",myEvent.mREC_SS_px_tiny,"mREC_SS_px_tiny[mRECnSS_tiny]/D");
	tinyTree->Branch("mREC_SS_py_tiny",myEvent.mREC_SS_py_tiny,"mREC_SS_py_tiny[mRECnSS_tiny]/D");
	tinyTree->Branch("mREC_SS_pz_tiny",myEvent.mREC_SS_pz_tiny,"mREC_SS_pz_tiny[mRECnSS_tiny]/D");
	tinyTree->Branch("mREC_SS_E_tiny",myEvent.mREC_SS_E_tiny,"mREC_SS_E_tiny[mRECnSS_tiny]/D");
	tinyTree->Branch("mREC_OS_px_tiny",myEvent.mREC_OS_px_tiny,"mREC_OS_px_tiny[mRECnOS_tiny]/D");
	tinyTree->Branch("mREC_OS_py_tiny",myEvent.mREC_OS_py_tiny,"mREC_OS_py_tiny[mRECnOS_tiny]/D");
	tinyTree->Branch("mREC_OS_pz_tiny",myEvent.mREC_OS_pz_tiny,"mREC_OS_pz_tiny[mRECnOS_tiny]/D");
	tinyTree->Branch("mREC_OS_E_tiny",myEvent.mREC_OS_E_tiny,"mREC_OS_E_tiny[mRECnOS_tiny]/D");

	TTree* tree = (TTree*) file->Get("myTree");
	Bool_t isTrigger_mini[4];
	Int_t mNumberOfTracks_mini;
	Int_t mNvertex_mini;
	Float_t mPosZ_mini[30];
	Float_t mPosX_mini[30];
	Float_t mPosY_mini[30];
  	Int_t mNPrimaryTracks_mini[30];

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

	if(doEmb_){
		tree->SetBranchAddress("mMCnTracks_mini",&mMCnTracks_mini);
		tree->SetBranchAddress("mMC_px_mini",&mMC_px_mini);
		tree->SetBranchAddress("mMC_py_mini",&mMC_py_mini);
		tree->SetBranchAddress("mMC_pz_mini",&mMC_pz_mini);
		tree->SetBranchAddress("mMC_E_mini",&mMC_E_mini);
		tree->SetBranchAddress("mMC_pdg_mini",&mMC_pdg_mini);
	}

	cout << "========== Reading the upc-dst miniTree, total number of events = " << tree->GetEntries() << " =========== " << endl;
	vector< TLorentzVector> e_plus;
	vector< TLorentzVector> e_minus;
	vector<double> mNSigmaE1, mNSigmaE2, mNSigmaPi1, mNSigmaPi2;

	for(int ievt=0; ievt<tree->GetEntries(); ievt++){

		tree->GetEntry(ievt);
		
		for(int ivtx = 0; ivtx < mNvertex_mini; ivtx++){

			myEvent.eventPass_tiny = 1;

			//vertex selections:
			if(!doEmb_) if( isTrigger_mini[2] != 1 ) myEvent.eventPass_tiny = 0;
			if(TMath::Abs(mPosX_mini[ivtx])<1.e-5 && TMath::Abs(mPosY_mini[ivtx])<1.e-5 && TMath::Abs(mPosZ_mini[ivtx])<1.e-5) myEvent.eventPass_tiny = 0;
			if(fabs(mPosZ_mini[ivtx]) > 100. ) myEvent.eventPass_tiny = 0; //loop over vertex and cut on Z vertex

			//event selections: at least two valid tracks that match BEMC and > 0.5 Gev in pT.
			int Nparticles = 0;
			for(int itrk = 0; itrk < mNumberOfTracks_mini; itrk++){

				if( (int) mVtxId_mini[itrk] != ivtx ) continue;

				if( mNhitsDEdx_mini[itrk] < 10 ) continue;
				if( mNhitsFit_mini[itrk] < 13 ) continue;
				if( fabs(mEta_mini[itrk]) > 1.0 ) continue;
				if( fabs(mDcaXY_mini[itrk]) > 3.0 ) continue;
				if( fabs(mDcaZ_mini[itrk]) > 3.0 ) continue;
				if( !mFlagBemc_mini[itrk] ) continue;
				if( mPt_mini[itrk] < 0.5 ) continue;
				
				Nparticles++;
			}
			if( Nparticles < 2 ) myEvent.eventPass_tiny = 0; //event selection end.
			
			//do embedding MC particle loop, fill the MC part
			vector< TLorentzVector> e_MC_plus, e_MC_minus;
			myEvent.mMCnTrack_tiny = -1;
			myEvent.mMCnOS_tiny = -1;
			int nOS_MC = 0;
			int Nparticles_MC = 0;
			myEvent.mMC_px_tiny[nOS_MC] = -999.;
			myEvent.mMC_py_tiny[nOS_MC] = -999.;
			myEvent.mMC_pz_tiny[nOS_MC] = -999.;
			myEvent.mMC_E_tiny[nOS_MC] = -999.;

			if(doEmb_){
				e_MC_plus.clear();e_MC_minus.clear();
				for(int imc = 0; imc < mMCnTracks_mini; imc++){
					TLorentzVector eMCParticle(mMC_px_mini[imc], mMC_py_mini[imc],mMC_pz_mini[imc],mMC_E_mini[imc]);
					if( fabs(eMCParticle.Eta()) > 1.0 ) continue;
					if( eMCParticle.Pt() < 0.5 ) continue;
					Nparticles_MC++;

					if( mMC_pdg_mini[imc] == -11 ) e_MC_plus.push_back( eMCParticle );
					if( mMC_pdg_mini[imc] == +11 ) e_MC_minus.push_back( eMCParticle );
				}
				if(e_MC_plus.size() < 1 || e_MC_minus.size() < 1 ) continue;
					for( unsigned i = 0; i < e_MC_plus.size(); i++){
						for( unsigned j = 0; j < e_MC_minus.size(); j++){
							TLorentzVector e_inv_MC = e_MC_plus[i] + e_MC_minus[j];
							if( TMath::Abs(e_inv_MC.Rapidity()) > 1.0 ) continue;
							if( e_inv_MC.Pt() < 0. ) continue;
							if( e_inv_MC.M() < 3.0 ) continue;

							myEvent.mMC_px_tiny[nOS_MC] = e_inv_MC.Px();
							myEvent.mMC_py_tiny[nOS_MC] = e_inv_MC.Py();
							myEvent.mMC_pz_tiny[nOS_MC] = e_inv_MC.Pz();
							myEvent.mMC_E_tiny[nOS_MC] = e_inv_MC.E();
							nOS_MC++;

						}
					}
			}
			myEvent.mMCnTrack_tiny = Nparticles_MC;
			myEvent.mMCnOS_tiny = nOS_MC;
			myEvent.mRECnTracks_tiny = Nparticles;

			TLorentzVector e1,e2;
			mNSigmaE1.clear();
			mNSigmaE2.clear();
			mNSigmaPi1.clear();
			mNSigmaPi2.clear();
			e_plus.clear();
			e_minus.clear();

			for(int itrk = 0; itrk < mNumberOfTracks_mini; itrk++){

				if( (int) mVtxId_mini[itrk] != ivtx ) continue;
				if( mNhitsDEdx_mini[itrk] < 10 ) continue;
				if( mNhitsFit_mini[itrk] < 13 ) continue;
				if( fabs(mEta_mini[itrk]) > 1.0 ) continue;
				if( fabs(mDcaXY_mini[itrk]) > 3.0 ) continue;
				if( fabs(mDcaZ_mini[itrk]) > 3.0 ) continue;
				if( !mFlagBemc_mini[itrk] ) continue;

				if( doEmb_ ){
					//already studied numbers:
					double a = smearing_a_para;
					double b = smearing_b_para;
					double pt_indep_sigma = pt_indep_const_sigma;
					if( mCharge_mini[itrk] == +1 ){
						double pt_smearing = (mPt_mini[itrk] - e_MC_plus[0].Pt())/e_MC_plus[0].Pt();
						double Nsigma = pt_smearing/pt_indep_sigma;
						double sigma_new = sqrt( (a*e_MC_plus[0].Pt())*(a*e_MC_plus[0].Pt()) + b*b);
						mPt_mini[itrk] = e_MC_plus[0].Pt()*(sigma_new*Nsigma+1);

					}
					if( mCharge_mini[itrk] == -1 ){
						double pt_smearing = (mPt_mini[itrk] - e_MC_minus[0].Pt())/e_MC_minus[0].Pt();
						double Nsigma = pt_smearing/pt_indep_sigma;
						double sigma_new = sqrt( (a*e_MC_minus[0].Pt())*(a*e_MC_minus[0].Pt()) + b*b);
						mPt_mini[itrk] = e_MC_minus[0].Pt()*(sigma_new*Nsigma+1);
					}
				}
				
				if( mPt_mini[itrk] < 0.5 ) continue;
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

			//like sign case
			int nSS = 0;
			if(e_plus.size() > 1){
				for(unsigned i = 0; i < e_plus.size(); i++){
					for(unsigned j = i+1; j < e_plus.size(); j++){

						double chi2_e = mNSigmaE1[i]*mNSigmaE1[i] + mNSigmaE1[j]*mNSigmaE1[j];
						double chi2_p = mNSigmaPi1[i]*mNSigmaPi1[i] + mNSigmaPi1[j]*mNSigmaPi1[j];
						
						// if( chi2_p<30 && (chi2_e/chi2_p) > elecPIDcut/30 ) continue;
						// if( chi2_p>30 && chi2_e > elecPIDcut ) continue;
						
						TLorentzVector e_inv;
						e_inv = e_plus[i]+e_plus[j];
						if( fabs(e_inv.Rapidity()) > 1.0 ) continue;
						myEvent.mREC_SS_px_tiny[nSS] = e_inv.Px();
						myEvent.mREC_SS_py_tiny[nSS] = e_inv.Py();
						myEvent.mREC_SS_pz_tiny[nSS] = e_inv.Pz();
						myEvent.mREC_SS_E_tiny[nSS] = e_inv.E();
						nSS++;
						
					}
				}
			}
			if(e_minus.size() > 1){
				for(unsigned i = 0; i < e_minus.size(); i++){
					for(unsigned j = i+1; j < e_minus.size(); j++){

						double chi2_e = mNSigmaE2[i]*mNSigmaE2[i] + mNSigmaE2[j]*mNSigmaE2[j];
						double chi2_p = mNSigmaPi2[i]*mNSigmaPi2[i] + mNSigmaPi2[j]*mNSigmaPi2[j];
						
						// if( chi2_p<30 && (chi2_e/chi2_p) > elecPIDcut/30 ) continue;
						// if( chi2_p>30 && chi2_e > elecPIDcut ) continue;
						
						TLorentzVector e_inv;
						e_inv = e_minus[i]+e_minus[j];
						if( fabs(e_inv.Rapidity()) > 1.0 ) continue;
						myEvent.mREC_SS_px_tiny[nSS] = e_inv.Px();
						myEvent.mREC_SS_py_tiny[nSS] = e_inv.Py();
						myEvent.mREC_SS_pz_tiny[nSS] = e_inv.Pz();
						myEvent.mREC_SS_E_tiny[nSS] = e_inv.E();
						nSS++;
						
					}
				}
			}

			myEvent.mRECnSS_tiny = nSS;

			//unlike sign case
			int nOS = 0;
			if( e_plus.size() >= 1 && e_minus.size() >= 1 ){
				//unlike-sign pair loop
				for( unsigned itrk = 0; itrk < e_plus.size(); itrk++ ){
					for( unsigned jtrk = 0; jtrk < e_minus.size(); jtrk++ ){

						double chi2_e = mNSigmaE1[itrk]*mNSigmaE1[itrk] + mNSigmaE2[jtrk]*mNSigmaE2[jtrk];
						double chi2_p = mNSigmaPi1[itrk]*mNSigmaPi1[itrk] + mNSigmaPi2[jtrk]*mNSigmaPi2[jtrk];
						
						// if( chi2_p<30 && (chi2_e/chi2_p) > elecPIDcut/30 ) continue;
						// if( chi2_p>30 && chi2_e > elecPIDcut ) continue;
						
						TLorentzVector e_inv;
						e_inv = e_plus[itrk]+e_minus[jtrk];
						if( fabs(e_inv.Rapidity()) > 1.0 ) continue;
						if( e_inv.Pt() < 0.0 ) continue;
						myEvent.mREC_OS_px_tiny[nOS] = e_inv.Px();
						myEvent.mREC_OS_py_tiny[nOS] = e_inv.Py();
						myEvent.mREC_OS_pz_tiny[nOS] = e_inv.Pz();
						myEvent.mREC_OS_E_tiny[nOS] = e_inv.E();
						nOS++;
					}
				}
			}
			myEvent.mRECnOS_tiny = nOS;
			
			tinyTree->Fill();

		}
	}

	
	
   	tinyTree->Write();




}