#include "RiceStyle.h"
#include "inputRootFile.h"
using namespace std;
#define PI            3.1415926
#define MASS_ELECTRON 0.00051

Double_t deltaPhi( double phi_1, double phi_2 ){

	double relAngle = 0.;
	if( phi_1 > phi_2 ){
		relAngle = phi_1 - phi_2;
		if( relAngle > PI ) relAngle = 2*PI - relAngle; 
	}
	else{
		relAngle = phi_1 - phi_2;
		if( relAngle > -PI ) relAngle = -relAngle;
		else relAngle = 2*PI + relAngle;
	}
	return relAngle;
}


void runPreStep_2_saveHisto( const bool doEmb_ = false, const bool doSmear_ = false ){

	/*
	Reweighting
	*/
	TFile * file_vtx = new TFile("reweight_vertex.root");
	TH1D* hVtxZ_data = (TH1D*) file_vtx->Get("hVtxZCut");

	/*
	Nsigma electron and pion in data
	*/
	TFile *  file_nsigma = new TFile("/Users/kong/google_drive/BNL_folder/Work/STAR/dileptonAnalyzer/macros/dAu200_hadronic.root");
	TH2D* hNSigmaEPion_pt[5];
	for(int j=0;j<5;j++){
		hNSigmaEPion_pt[j] = (TH2D*) file_nsigma->Get(Form("hNSigmaEPion_pt_%d",j));
	}
	double ptbins_nsigma[]={0.5,1.0,1.3,1.6,1.9,100.};

	//input
	TFile* file = 0;
	if(doEmb_) file = new TFile( input_embedding );
	else file = new TFile( input_data );
	if(!file) cout << "wrong name! Check input files" << endl;

	//output
	TString outname;
	if(!doEmb_){
		outname = "output-PreStep_2-data.root";
	}
	else{
		if( doSmear_ ){
			outname = "output-PreStep_2-embedding.root";
		}
		else{
			outname = "output-PreStep_2-embedding-unsmear.root";
		}
	}
	TFile *fout = new TFile(outname,"recreate");

	int ptNbins = sizeof(ptbins)/sizeof(ptbins[0]) - 1;
	int pt2Nbins = sizeof(pt2bins)/sizeof(pt2bins[0]) - 1;
	int pt2NbinsREC = sizeof(pt2binsREC)/sizeof(pt2binsREC[0]) - 1;
	int massNbins = sizeof(massbins)/sizeof(massbins[0]) - 1;

	TH2D* h_Nvertex = new TH2D("h_Nvertex",";n vertex;itrg",20,0,20,5,0,5);
	TH1D* hVtxZ = new TH1D("hVtxZ",";z vertex (cm)", 200,-500,500);
	TH1D* hVtxZCut = new TH1D("hVtxZCut",";z vertex (cm)", 200,-500,500);
	TH1D* hVtxZzdc = new TH1D("hVtxZzdc",";z vertex (cm)", 200,-1000,10000);
	TH2D* hVtxZzdcTPC = new TH2D("hVtxZzdcTPC",";zdc;tpc",200,-1000,1000,200,-1000,1000);
	TH1D* hVtxZzdcTPCcut = new TH1D("hVtxZzdcTPCcut",";zdc vertex",200,-1000,1000);
	TH1D* hTofMult = new TH1D("hTofMult","hTofMult",10,0,10);
	TH1D* hTofMultJpsi = new TH1D("hTofMultJpsi","hTofMultJpsi",10,0,10);
	TH1D* hNtrkValid = new TH1D("hNtrkValid","N",10,0,10);
	TH2D* hTrigSimu = new TH2D("hTrigSimu",";data;mc",2,0,2,2,0,2);
	TH2D* hTrigBemc = new TH2D("hTrigBemc",";adc0;#Delta#phi",300,0,300,100,0,PI);

	TH1D* hJpsiMass = new TH1D("hJpsiMass","hJpsiMass",massNbins,massbins);
	TH1D* hJpsiMassZDCveto = new TH1D("hJpsiMassZDCveto","hJpsiMassZDCveto",massNbins,massbins);
	TH2D* hJpsiMass_Pt2 = new TH2D("hJpsiMass_Pt2","hJpsiMass_Pt2",massNbins,massbins,pt2Nbins,pt2bins);
	TH1D* hJpsiPt2_match = new TH1D("hJpsiPt2_match","hJpsiPt2_match",pt2Nbins,pt2bins);
	TH2D* hJpsiMassZDCveto_Pt2 = new TH2D("hJpsiMassZDCveto_Pt2","hJpsiMassZDCveto_Pt2",massNbins,massbins,pt2Nbins,pt2bins);
	TH1D* hDielectronPt = new TH1D("hDielectronPt","hDielectronPt; #Dielectron p_{T} (GeV^{2})",100,0,3);
	TH1D* hDielectronPt2 = new TH1D("hDielectronPt2","hDielectronPt2; #Dielectron p^{2}_{T} (GeV^{2})",pt2NbinsREC,pt2binsREC);
	TH1D* hJpsiMassZDC = new TH1D("hJpsiMassZDC","hJpsiMassZDC",massNbins,massbins);
	TH2D* hJpsiMassZDC_Pt2 = new TH2D("hJpsiMassZDC_Pt2","hJpsiMassZDC_Pt2",massNbins,massbins,pt2Nbins,pt2bins);

	TH1D* hLikeSignMass = new TH1D("hLikeSignMass","hLikeSignMass",massNbins,massbins);
	TH1D* hLikeSignMassZDCveto = new TH1D("hLikeSignMassZDCveto","hLikeSignMassZDCveto",massNbins,massbins);
	TH2D* hLikeSignMass_Pt2 = new TH2D("hLikeSignMass_Pt2","hLikeSignMass_Pt2",massNbins,massbins,pt2Nbins,pt2bins);
	TH2D* hLikeSignMassZDCveto_Pt2 = new TH2D("hLikeSignMassZDCveto_Pt2","hLikeSignMassZDCveto_Pt2",massNbins,massbins,pt2Nbins,pt2bins);
	TH1D* hLikeSignMassZDC = new TH1D("hLikeSignMassZDC","hLikeSignMassZDC",massNbins,massbins);
	TH2D* hLikeSignMassZDC_Pt2 = new TH2D("hLikeSignMassZDC_Pt2","hLikeSignMassZDC_Pt2",massNbins,massbins,pt2Nbins,pt2bins);

	//MCParticle 
	TH1D* hSingleTrackGEN = new TH1D("hSingleTrackGEN",";pt",100,0,5);
	TH1D* hSingleTrackREC = new TH1D("hSingleTrackREC",";pt",100,0,5);
	TH1D* hMatchDeltaR = new TH1D("hMatchDeltaR","#DeltaR",1000,0,3);
	TH1D* hMatchVtxIndex = new TH1D("hMatchVtxIndex","vtx index",5,0,5);

	TH1D* hMCJpsiMass = new TH1D("hMCJpsiMass","hMCJpsiMass; M_{e^{+}e^{-}} (GeV/c^{2})",60,0.4,4);
	TH1D* hMCDielectronPt = new TH1D("hMCDielectronPt","hMCDielectronPt; #Dielectron p_{T} (GeV)",100,0,3);
	TH1D* hMCDielectronPt2 = new TH1D("hMCDielectronPt2","hMCDielectronPt2; #Dielectron p^{2}_{T} (GeV^{2})",pt2Nbins,pt2bins);
	
	//non physics distribution
	TH1D* h_EoverP_BEMC = new TH1D("h_EoverP_BEMC",";E/p",50,0,5.0);
	TH1D* h_PhiDist_BEMC = new TH1D("h_PhiDist_BEMC","h_PhiDist_BEMC; PhiDist",50,0,PI);
	TH1D* h_EtaDist_BEMC = new TH1D("h_EtaDist_BEMC","h_EtaDist_BEMC; EtaDist",50,-2,2);
	TH1D* h_E0_BEMC = new TH1D("h_E0_BEMC",";leading hot tower",100,0,10);
	TH2D* h_NsigmaElec = new TH2D("h_NsigmaElec","h_NsigmaElec",50,-10,10,50,-10,10);
	TH2D* h_NsigmaPion = new TH2D("h_NsigmaPion","h_NsigmaPion",50,-10,10,50,-10,10);
	TH2D* h_NsigmaPart = new TH2D("h_NsigmaPart",";#chi^{2}_{#pi#pi};#chi^{2}_{ee}",100,0,120,100,0,120);
	TH2D* h_NsigmaCorr = new TH2D("h_NsigmaCorr","h_NsigmaCorr",50,-10,10,50,-10,10);
	TH1D* h_Pm = new TH1D("h_Pm",";p (GeV/c)",50,0,3);
	TH1D* h_Pm_BEMC = new TH1D("h_Pm_BEMC",";p (GeV/c)",50,0,3);

	TH2D* h_deltaR = new TH2D("h_deltaR",";deltaR;pt",100,0,PI,100,-1,1);
	//single electron
	TH1D* hElecPt = new TH1D("hElecPt",";p_{T} (GeV/c)",50,0,3);
	TH1D* hElecEta = new TH1D("hElecEta",";#eta",10,-1.2,1.2);;
	TH1D* hElecPhi = new TH1D("hElecPhi",";#phi",10,-PI,PI);;
	TH1D* hElecDCAxy = new TH1D("hElecDCAxy",";DCA xy",50,0,10);;
	TH1D* hElecDCAz = new TH1D("hElecDCAz",";DCA z",100,-10,10);;
	TH1D* hElecNhitsFit = new TH1D("hElecNhitsFit",";NhitsFit",50,0,50);;
	TH1D* hElecNhitsDedx = new TH1D("hElecNhitsDedx",";NhitsDedx",50,0,50);;

	//ZDC and SMD 
	TH1D* h_ZDCEast[4];
	TH1D* h_ZDCWest[4];
	for(int i=0;i<4;i++){
 		h_ZDCEast[i]= new TH1D(Form("h_ZDCEast_%d",i),Form("h_ZDCEast_%d",i),600,0,1200);
 		h_ZDCWest[i]= new TH1D(Form("h_ZDCWest_%d",i),Form("h_ZDCWest_%d",i),600,0,1200);
	}

	TTree* tree = (TTree*) file->Get("myTree");

//tree branches
	Int_t mEvtNum_mini;
	Bool_t isTrigger_mini[4];
	Bool_t isSimuTrigger_mini[4];
	Int_t mNumberOfTracks_mini;
	Int_t mNvertex_mini;
	Float_t mPosZ_mini[30];
	Float_t mPosX_mini[30];
	Float_t mPosY_mini[30];
	Float_t mZdcVertexZ_mini;
  	Int_t mNPrimaryTracks_mini[30];
  	Short_t mTofMult_mini;

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
	Float16_t mNSigmasTPCElectron_mini[MAXTracks]; // dE/dx n sigmas for particle species
	Float16_t mNSigmasTPCPion_mini[MAXTracks]; // dE/dx n sigmas for particle species
	
	Double32_t mBemcPmag_mini[MAXTracks];
	Float_t mBemcEta_mini[MAXTracks];
	Float_t mBemcPhi_mini[MAXTracks];
	Float_t mBemcEnergy_mini[MAXTracks]; // energy of matched BEMC cluster
	Float_t mBemcHTEnergy_mini[MAXTracks]; // energy of matched hot tower
	Float_t mBemcClsEta_mini[MAXTracks]; // eta of matched BEMC cluster
	Float_t mBemcClsPhi_mini[MAXTracks]; // phi of matched BEMC cluster
	Float_t mBemcClsAdc0_mini[MAXTracks]; // adc0 of matched BEMC cluster
	Float_t mBemcClsDsmadc0_mini[MAXTracks]; // dsmadc0 of matched BEMC cluster

	//MC particles
	Int_t mMCnTracks_mini = 0;
	Double32_t mMC_px_mini[MAXTracks];
	Double32_t mMC_py_mini[MAXTracks];
	Double32_t mMC_pz_mini[MAXTracks];
	Double32_t mMC_E_mini[MAXTracks];
  	Int_t mMC_pdg_mini[MAXTracks];

	tree->SetBranchAddress("mEvtNum_mini",&mEvtNum_mini);
	tree->SetBranchAddress("isTrigger_mini",&isTrigger_mini);
	tree->SetBranchAddress("isSimuTrigger_mini",&isSimuTrigger_mini);
	tree->SetBranchAddress("mNumberOfTracks_mini",&mNumberOfTracks_mini);
	tree->SetBranchAddress("mNvertex_mini",&mNvertex_mini);
	tree->SetBranchAddress("mNPrimaryTracks_mini",&mNPrimaryTracks_mini);
	tree->SetBranchAddress("mPosX_mini",&mPosX_mini);
	tree->SetBranchAddress("mPosY_mini",&mPosY_mini);
	tree->SetBranchAddress("mPosZ_mini",&mPosZ_mini);
	tree->SetBranchAddress("mZdcVertexZ_mini",&mZdcVertexZ_mini);
	tree->SetBranchAddress("mTofMult_mini",&mTofMult_mini);
	
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
	tree->SetBranchAddress("mNSigmasTPCElectron_mini",&mNSigmasTPCElectron_mini);
	tree->SetBranchAddress("mNSigmasTPCPion_mini",&mNSigmasTPCPion_mini);
	tree->SetBranchAddress("mBemcHTEnergy_mini",&mBemcHTEnergy_mini);
	tree->SetBranchAddress("mBemcEta_mini",&mBemcEta_mini);
	tree->SetBranchAddress("mBemcPhi_mini",&mBemcPhi_mini);
	tree->SetBranchAddress("mBemcPmag_mini",&mBemcPmag_mini);
	tree->SetBranchAddress("mBemcEnergy_mini",&mBemcEnergy_mini);
	tree->SetBranchAddress("mBemcClsEta_mini",&mBemcClsEta_mini);
	tree->SetBranchAddress("mBemcClsPhi_mini",&mBemcClsPhi_mini);
	tree->SetBranchAddress("mBemcClsAdc0_mini",&mBemcClsAdc0_mini);
	tree->SetBranchAddress("mBemcClsDsmadc0_mini",&mBemcClsDsmadc0_mini);
	
	if(doEmb_){
		tree->SetBranchAddress("mMCnTracks_mini",&mMCnTracks_mini);
		tree->SetBranchAddress("mMC_px_mini",&mMC_px_mini);
		tree->SetBranchAddress("mMC_py_mini",&mMC_py_mini);
		tree->SetBranchAddress("mMC_pz_mini",&mMC_pz_mini);
		tree->SetBranchAddress("mMC_E_mini",&mMC_E_mini);
		tree->SetBranchAddress("mMC_pdg_mini",&mMC_pdg_mini);
	}
//end tree branch

	cout << "========== Reading the upc-dst miniTree, total number of events = " << tree->GetEntries() << " =========== " << endl;
	
	vector< TLorentzVector> e_plus;
	vector< TLorentzVector> e_minus;
	vector<double> mNSigmaE1, mNSigmaE2, mNSigmaPi1, mNSigmaPi2;
	vector<double> mTOFtime1, mTOFtime2;
	vector<double> mBEMCEnergy1,mBEMCEnergy2,mBEMCPmag1,mBEMCPmag2;

	for(int ievt=0; ievt<tree->GetEntries(); ievt++){

		tree->GetEntry(ievt);

		//data only, checking validity of trigger simulations.
		if(!doEmb_) {
			for(int itrg=0;itrg<4;itrg++){
				if( isTrigger_mini[itrg] == true ){
					h_Nvertex->Fill( mNvertex_mini, itrg);
				}
			}
			if(isTrigger_mini[2] == true) {
				hVtxZzdcTPC->Fill( mZdcVertexZ_mini, mPosZ_mini[0]);
				if( fabs(mPosZ_mini[0]) < 100 ){
					hVtxZzdcTPCcut->Fill( mZdcVertexZ_mini );
				}
			}
			hTrigSimu->Fill(isTrigger_mini[2], isSimuTrigger_mini[2] );
			vector< double> adc0,bemcPhi;
			for(int itrk = 0; itrk < mNumberOfTracks_mini; itrk++){
				adc0.push_back(mBemcClsAdc0_mini[itrk]);
				bemcPhi.push_back(mBemcClsPhi_mini[itrk]);
			}
			if( isTrigger_mini[2] && isSimuTrigger_mini[2] ){
				if(bemcPhi.size() < 1) continue;
				for(unsigned ipart=0;ipart<bemcPhi.size();ipart++){
					for(unsigned jpart=ipart+1;jpart<bemcPhi.size();jpart++){
						if( ipart == jpart ) continue;
						double deltaBemcPhi = deltaPhi(bemcPhi[ipart],bemcPhi[jpart]);
						hTrigBemc->Fill(adc0[ipart], deltaBemcPhi);
						hTrigBemc->Fill(adc0[jpart], deltaBemcPhi);
					}
				}
			}

			//trigger requirements
			if( isTrigger_mini[2] != 1 ) continue;
		}else{
			h_Nvertex->Fill( mNvertex_mini, isSimuTrigger_mini[2]);
		}
		
		//
		//No event without vertex will pass here
		//

		if( mNvertex_mini <= 0 ) continue;
		vector< TLorentzVector> e_MC_plus, e_MC_minus, j_MC;
		if(doEmb_){
			//MC single track
			for(int imc = 0; imc < mMCnTracks_mini; imc++){
				TLorentzVector eMCParticle(mMC_px_mini[imc], mMC_py_mini[imc],mMC_pz_mini[imc],mMC_E_mini[imc]);
				if( mMC_pdg_mini[imc] == -11 ) e_MC_plus.push_back( eMCParticle );
				if( mMC_pdg_mini[imc] == +11 ) e_MC_minus.push_back( eMCParticle );
				if( fabs(mMC_pdg_mini[imc]) != 11 ) continue;
				if( eMCParticle.Pt() < 0.5  ) continue;
				if( imc >= 2 ) continue;
				hSingleTrackGEN->Fill( eMCParticle.Pt() );
			}
			
			//MC J/psi
			for( unsigned i = 0; i < e_MC_plus.size(); i++){
				for( unsigned j = 0; j < e_MC_minus.size(); j++){
					TLorentzVector e_inv_MC = e_MC_plus[i] + e_MC_minus[j];
					if( TMath::Abs(e_inv_MC.Rapidity()) > 1.0 ) continue;
					if( e_inv_MC.Pt() < 0. ) continue;
					if( e_inv_MC.M() < 3.0 ) continue;
					hMCJpsiMass->Fill( e_inv_MC.M() );
					hMCDielectronPt->Fill( e_inv_MC.Pt() );
					hMCDielectronPt2->Fill( e_inv_MC.Pt()*e_inv_MC.Pt() );
				}
			}
		}


		//loop over different vertices
		for(int ivtx = 0; ivtx < mNvertex_mini; ivtx++){
			//fill vertex z before cuts
			hVtxZ->Fill( mPosZ_mini[ivtx] );
			//vertex selections:
			if(TMath::Abs(mPosX_mini[ivtx])<1.e-5 && TMath::Abs(mPosY_mini[ivtx])<1.e-5 && TMath::Abs(mPosZ_mini[ivtx])<1.e-5) continue;
			if(fabs(mPosZ_mini[ivtx]) > 100. ) continue; //loop over vertex and cut on Z vertex
			hVtxZzdc->Fill( mZdcVertexZ_mini );			
			//vertex weight
			double weight = 1.0;
			if( doEmb_ ){
				weight = hVtxZ_data->GetBinContent(hVtxZ_data->FindBin(mPosZ_mini[ivtx])); 
			}
			weight = 1.0;//use 1 for now

			if(doEmb_){			
				//single track efficiency
				for(int imc = 0; imc < mMCnTracks_mini; imc++){
					TLorentzVector eMCParticle(mMC_px_mini[imc], mMC_py_mini[imc],mMC_pz_mini[imc],mMC_E_mini[imc]);
					if( fabs(mMC_pdg_mini[imc]) != 11 ) continue;
					if( eMCParticle.Pt() < 0.5 ) continue;
					if( imc >= 2 ) continue;
			
					double deltaR_minPlus = 999.;
					double deltaR_minMinus = 999.;
					int bestIndex = -1;
					TVector3 bestMatchEplus(0,0,0);
					TVector3 bestMatchEminus(0,0,0);
					for(int itrk = 0; itrk < mNumberOfTracks_mini; itrk++){

						if( (int) mVtxId_mini[itrk] != ivtx ) continue;
						if( mNhitsDEdx_mini[itrk] < cutNhitDedx ) continue;
						if( mNhitsFit_mini[itrk] < cutNhitFit ) continue;
						if( fabs(mEta_mini[itrk]) > cutEtaDaug ) continue;
						if( fabs(mDcaXY_mini[itrk]) > cutDCAxyz ) continue;
						if( fabs(mDcaZ_mini[itrk]) > cutDCAxyz ) continue;
						if( !mFlagBemc_mini[itrk] ) continue;
						if( mPt_mini[itrk] < 0.5 ) continue;
						TVector3 part3;
						part3.SetPtEtaPhi(mPt_mini[itrk],mEta_mini[itrk],mPhi_mini[itrk]);
						
						if( mCharge_mini[itrk] == +1 && mMC_pdg_mini[imc] == -11 ){
							if( isMatch(part3,eMCParticle) ){
								TVector3 partmc3 = eMCParticle.Vect();
								hMatchDeltaR->Fill( part3.DeltaR(partmc3) );
								if( part3.DeltaR(partmc3) < deltaR_minPlus ) {
									deltaR_minPlus=part3.DeltaR(partmc3);
									bestMatchEplus = part3;
									bestIndex = mVtxId_mini[itrk];
									
								}
							}
						}
						if( mCharge_mini[itrk] == -1 && mMC_pdg_mini[imc] == +11 ){
							if( isMatch(part3,eMCParticle) ){
								TVector3 partmc3 = eMCParticle.Vect();
								hMatchDeltaR->Fill( part3.DeltaR(partmc3) );
								if( part3.DeltaR(partmc3) < deltaR_minMinus ) {
									deltaR_minMinus=part3.DeltaR(partmc3);
									bestMatchEminus = part3;
									bestIndex = mVtxId_mini[itrk];
								}
							}
						}
						
					}
					if( bestMatchEminus.Pt() !=0 ) {
						hSingleTrackREC->Fill( bestMatchEminus.Pt() );
					}
					if( bestMatchEplus.Pt() != 0 ){
						hSingleTrackREC->Fill( bestMatchEplus.Pt() );
					}
				}//single track efficiency

				//single track BEMC matching efficiency			
				for(int imc = 0; imc < mMCnTracks_mini; imc++){
					TLorentzVector eMCParticle(mMC_px_mini[imc], mMC_py_mini[imc],mMC_pz_mini[imc],mMC_E_mini[imc]);
					if( fabs(mMC_pdg_mini[imc]) != 11 ) continue;
					if( eMCParticle.Pt() < 0.5 ) continue;
					if( imc >= 2 ) continue;

					for(int itrk = 0; itrk < mNumberOfTracks_mini; itrk++){

						if( (int) mVtxId_mini[itrk] != ivtx ) continue;
						if( mNhitsDEdx_mini[itrk] < cutNhitDedx ) continue;
						if( mNhitsFit_mini[itrk] < cutNhitFit ) continue;
						if( fabs(mEta_mini[itrk]) > cutEtaDaug ) continue;
						if( fabs(mDcaXY_mini[itrk]) > cutDCAxyz ) continue;
						if( fabs(mDcaZ_mini[itrk]) > cutDCAxyz ) continue;
						if( mPt_mini[itrk] < 0.5 ) continue;

						TVector3 reco;
						reco.SetPtEtaPhi(mPt_mini[itrk],mEta_mini[itrk],mPhi_mini[itrk]);
						if( isMatch(reco, eMCParticle)  ) h_Pm->Fill( reco.Pt(), weight );
						if( isMatch(reco, eMCParticle) && mFlagBemc_mini[itrk] ) h_Pm_BEMC->Fill( reco.Pt(), weight );
					}
				}//bemc matching

			}
		
			//fill tofMult
			hTofMult->Fill( mTofMult_mini );
			double totalZDCEastADC = 0.;
			double totalZDCWestADC = 0.;
			for(int izdc=0;izdc<3;izdc++){	
				totalZDCEastADC += mZdcEastADC_mini[izdc];
				totalZDCWestADC += mZdcWestADC_mini[izdc];
			}
			h_ZDCEast[3]->Fill(totalZDCEastADC,weight);
			h_ZDCWest[3]->Fill(totalZDCWestADC,weight);
			
			TLorentzVector e1,e2;

			int Nparticles_plus = 0;
			int Nparticles_minus = 0;
			
			mNSigmaE1.clear();
			mNSigmaE2.clear();

			mNSigmaPi1.clear();
			mNSigmaPi2.clear();

			mBEMCEnergy1.clear();
			mBEMCEnergy2.clear();
			mBEMCPmag1.clear();
			mBEMCPmag2.clear();

			e_plus.clear();
			e_minus.clear();

			double deltaR_minPlus = 999.;
			double deltaR_minMinus = 999.;
			int bestIndex = -1;
			TVector3 bestMatchEplus(0,0,0);
			TVector3 bestMatchEminus(0,0,0);
			for(int itrk = 0; itrk < mNumberOfTracks_mini; itrk++){

				if( (int) mVtxId_mini[itrk] != ivtx ) continue;

				hElecNhitsDedx->Fill( mNhitsDEdx_mini[itrk], weight );
				hElecNhitsFit->Fill( mNhitsFit_mini[itrk], weight );
				hElecDCAxy->Fill( mDcaXY_mini[itrk], weight );
				hElecDCAz->Fill( mDcaZ_mini[itrk], weight );

				if( mNhitsDEdx_mini[itrk] < cutNhitDedx ) continue;
				if( mNhitsFit_mini[itrk] < cutNhitFit ) continue;
				if( fabs(mEta_mini[itrk]) > cutEtaDaug ) continue;
				if( fabs(mDcaXY_mini[itrk]) > cutDCAxyz ) continue;
				if( fabs(mDcaZ_mini[itrk]) > cutDCAxyz ) continue;
				if( !mFlagBemc_mini[itrk] ) continue;

				//Recalculate the nsigma electron and pion in embedding
				if( doEmb_ ){
					//already studied numbers:
					double a = smearing_a_para;
					double b = smearing_b_para;
					double pt_indep_sigma = pt_indep_const_sigma;
					//select matching tracks
					TVector3 part3;
					part3.SetPtEtaPhi(mPt_mini[itrk],mEta_mini[itrk],mPhi_mini[itrk]);
					TLorentzVector eMCParticle_plus = e_MC_plus[0];
					TLorentzVector eMCParticle_minus = e_MC_minus[0];
					if( mCharge_mini[itrk] == +1 ){
						if( isMatch(part3,eMCParticle_plus) ){
							TVector3 partmc3 = eMCParticle_plus.Vect();
							if( part3.DeltaR(partmc3) < deltaR_minPlus ) {
								deltaR_minPlus=part3.DeltaR(partmc3);
								bestMatchEplus = part3;

								int index_smear = 0;
								for(int ismear=0;ismear<5;ismear++){
									if( eMCParticle_plus.Pt()>pt_wide_bins[ismear] && eMCParticle_plus.Pt()<pt_wide_bins[ismear+1]) index_smear = ismear;
								}
								pt_indep_sigma = pt_dep_sigma[index_smear];

								double pt_smearing = (mPt_mini[itrk] - e_MC_plus[0].Pt())/e_MC_plus[0].Pt();
								double Nsigma = pt_smearing/pt_indep_sigma;
								double sigma_new = sqrt( (a*e_MC_plus[0].Pt())*(a*e_MC_plus[0].Pt()) + b*b);
								// if( doSmear_) mPt_mini[itrk] = e_MC_plus[0].Pt()*(sigma_new*Nsigma+1);
								if(doSmear_) {
									mPt_mini[itrk] = e_MC_plus[0].Pt()*(sigma_new*Nsigma+1);									
									ePlus_cand.SetPtEtaPhiM(mPt_mini[itrk],part3.Eta(),part3.Phi(),MASS_ELECTRON);
									double fraction = 0.;
									double gamma_energy = 0.;
									if( getProb->GetRandom() < BHrate ){
										gamma_energy = getBremPhotonEnergy(ePlus_cand.E());
										fraction = gamma_energy / ePlus_cand.E();
									}
									ePlus_cand = ePlus_cand*(1.-fraction);
									mPt_mini[itrk] = ePlus_cand.Pt();
								}	
								if( mPt_mini[itrk] < 0.5 ) bestMatchEplus.SetPtEtaPhi(0,0,0);
								else bestMatchEplus.SetPtEtaPhi(mPt_mini[itrk],part3.Eta(),part3.Phi() );
							}
						}
					}
					if( mCharge_mini[itrk] == -1 ){
						if( isMatch(part3,eMCParticle_minus) ){
							TVector3 partmc3 = eMCParticle_minus.Vect();
							if( part3.DeltaR(partmc3) < deltaR_minMinus ) {
								deltaR_minMinus=part3.DeltaR(partmc3);
								bestMatchEminus = part3;

								int index_smear = 0;
								for(int ismear=0;ismear<5;ismear++){
									if( eMCParticle_minus.Pt()>pt_wide_bins[ismear] && eMCParticle_minus.Pt()<pt_wide_bins[ismear+1]) index_smear = ismear;
								}
								pt_indep_sigma = pt_dep_sigma[index_smear];

								double pt_smearing = (mPt_mini[itrk] - e_MC_minus[0].Pt())/e_MC_minus[0].Pt();
								double Nsigma = pt_smearing/pt_indep_sigma;
								double sigma_new = sqrt( (a*e_MC_minus[0].Pt())*(a*e_MC_minus[0].Pt()) + b*b);
								// if( doSmear_ ) mPt_mini[itrk] = e_MC_minus[0].Pt()*(sigma_new*Nsigma+1);
								if(doSmear_) {
									mPt_mini[itrk] = e_MC_minus[0].Pt()*(sigma_new*Nsigma+1);
									eMinus_cand.SetPtEtaPhiM(mPt_mini[itrk],part3.Eta(),part3.Phi(),MASS_ELECTRON);
									double fraction = 0.;
									double gamma_energy = 0.;
									if( getProb->GetRandom() < BHrate ){
										gamma_energy = getBremPhotonEnergy(eMinus_cand.E());
										fraction = gamma_energy / eMinus_cand.E();
									}
									eMinus_cand = eMinus_cand*(1.-fraction);
									mPt_mini[itrk] = eMinus_cand.Pt();
								}
								if( mPt_mini[itrk] < 0.5 ) bestMatchEminus.SetPtEtaPhi(0,0,0);
								else bestMatchEminus.SetPtEtaPhi(mPt_mini[itrk],part3.Eta(),part3.Phi() );
							}
						}
					}

					for(int j=0;j<5;j++){
						if(mPt_mini[itrk] > ptbins_nsigma[j] && mPt_mini[itrk] < ptbins_nsigma[j+1]){
							double nsigma_e_data=0.;
							double nsigma_pi_data=0.;
							hNSigmaEPion_pt[j]->GetRandom2(nsigma_e_data,nsigma_pi_data);
							mNSigmasTPCElectron_mini[itrk] = nsigma_e_data;
							mNSigmasTPCPion_mini[itrk] = nsigma_pi_data;
						}
					}
				}
		
				if( mPt_mini[itrk] < 0.5 ) continue;
				
				h_E0_BEMC->Fill( mBemcHTEnergy_mini[itrk] , weight);
				h_EtaDist_BEMC->Fill(mBemcEta_mini[itrk]-mBemcClsEta_mini[itrk], weight);
				h_PhiDist_BEMC->Fill( deltaPhi(mBemcPhi_mini[itrk],mBemcClsPhi_mini[itrk]) , weight);				
				
				if( mCharge_mini[itrk] == +1 ){
					Nparticles_plus++;
					e1.SetPtEtaPhiM(mPt_mini[itrk],mEta_mini[itrk],mPhi_mini[itrk],MASS_ELECTRON);
					
					e_plus.push_back( e1 );
					mNSigmaE1.push_back(mNSigmasTPCElectron_mini[itrk]);
					mNSigmaPi1.push_back(mNSigmasTPCPion_mini[itrk]);
					mBEMCEnergy1.push_back(mBemcEnergy_mini[itrk]);
					mBEMCPmag1.push_back(mBemcPmag_mini[itrk]);
				}

				if( mCharge_mini[itrk] == -1 ){
					Nparticles_minus++;
					e2.SetPtEtaPhiM(mPt_mini[itrk],mEta_mini[itrk],mPhi_mini[itrk],MASS_ELECTRON);
					
					e_minus.push_back( e2 );
					mNSigmaE2.push_back(mNSigmasTPCElectron_mini[itrk]);
					mNSigmaPi2.push_back(mNSigmasTPCPion_mini[itrk]);
					mBEMCEnergy2.push_back(mBemcEnergy_mini[itrk]);
					mBEMCPmag2.push_back(mBemcPmag_mini[itrk]);
				}

			}//end of track loop

		//like sign case
			if(e_plus.size() > 1){
				for(unsigned i = 0; i < e_plus.size(); i++){
					for(unsigned j = i+1; j < e_plus.size(); j++){

						double chi2_e = mNSigmaE1[i]*mNSigmaE1[i] + mNSigmaE1[j]*mNSigmaE1[j];
						double chi2_p = mNSigmaPi1[i]*mNSigmaPi1[i] + mNSigmaPi1[j]*mNSigmaPi1[j];
						//embedding has correct nsigma e and pi
						if( chi2_p<elecPIDcut_b && (chi2_e/chi2_p) > elecPIDcut_a/elecPIDcut_b ) continue;
						if( chi2_p>elecPIDcut_b && chi2_e > elecPIDcut_a ) continue;
						
						TLorentzVector e_inv;
						e_inv = e_plus[i]+e_plus[j];
						if( fabs(e_inv.Rapidity()) > 1.0 ) continue;
						hLikeSignMass_Pt2->Fill( e_inv.M(), e_inv.Pt() * e_inv.Pt(), weight );
						hLikeSignMass->Fill( e_inv.M(), weight );
						//neutron coincidence or not
						if( totalZDCWestADC < 40 || (totalZDCWestADC > 9999.) ){
							hLikeSignMassZDCveto->Fill( e_inv.M() );
							hLikeSignMassZDCveto_Pt2->Fill( e_inv.M(), e_inv.Pt() * e_inv.Pt(), weight );
						}
						else{
							hLikeSignMassZDC->Fill( e_inv.M() );
							hLikeSignMassZDC_Pt2->Fill( e_inv.M(), e_inv.Pt() * e_inv.Pt(), weight );
						}
						
					}
				}
			}
			if(e_minus.size() > 1){
				for(unsigned i = 0; i < e_minus.size(); i++){
					for(unsigned j = i+1; j < e_minus.size(); j++){

						double chi2_e = mNSigmaE2[i]*mNSigmaE2[i] + mNSigmaE2[j]*mNSigmaE2[j];
						double chi2_p = mNSigmaPi2[i]*mNSigmaPi2[i] + mNSigmaPi2[j]*mNSigmaPi2[j];

						//embedding has correct nsigma e and pi
						if( chi2_p<elecPIDcut_b && (chi2_e/chi2_p) > elecPIDcut_a/elecPIDcut_b ) continue;
						if( chi2_p>elecPIDcut_b && chi2_e > elecPIDcut_a ) continue;
						
						TLorentzVector e_inv;
						e_inv = e_minus[i]+e_minus[j];
						if( fabs(e_inv.Rapidity()) > 1.0 ) continue;
						hLikeSignMass_Pt2->Fill( e_inv.M(), e_inv.Pt() * e_inv.Pt(), weight );
						hLikeSignMass->Fill( e_inv.M(), weight );
						//neutron coincidence or not
						if( totalZDCWestADC < 40 || (totalZDCWestADC > 9999.) ){
							hLikeSignMassZDCveto->Fill( e_inv.M() );
							hLikeSignMassZDCveto_Pt2->Fill( e_inv.M(), e_inv.Pt() * e_inv.Pt(), weight );
						}
						else{
							hLikeSignMassZDC->Fill( e_inv.M() );
							hLikeSignMassZDC_Pt2->Fill( e_inv.M(), e_inv.Pt() * e_inv.Pt(), weight );
						}
					}
				}
			}
		//end like sign

		//unlike sign case
			if( e_plus.size() < 1 || e_minus.size() < 1 ) continue;
			if( doEmb_ ){
				if( bestMatchEplus.Pt() == 0 || bestMatchEminus.Pt() == 0 ) continue;
			
				e_plus.clear();
				e_minus.clear();
				TLorentzVector bestMatchEplus4vect,bestMatchEminus4vect;
				bestMatchEplus4vect.SetPtEtaPhiM(bestMatchEplus.Pt(),bestMatchEplus.Eta(),bestMatchEplus.Phi(), MASS_ELECTRON);
				bestMatchEminus4vect.SetPtEtaPhiM(bestMatchEminus.Pt(),bestMatchEminus.Eta(),bestMatchEminus.Phi(), MASS_ELECTRON);
				e_plus.push_back( bestMatchEplus4vect );
				e_minus.push_back( bestMatchEminus4vect );
			}
			//after event selections 
			hVtxZCut->Fill( mPosZ_mini[ivtx] );
			h_ZDCEast[2]->Fill(totalZDCEastADC, weight);
			h_ZDCWest[2]->Fill(totalZDCWestADC, weight);
			
			//unlike-sign pair loop
			for( unsigned itrk = 0; itrk < e_plus.size(); itrk++ ){
				for( unsigned jtrk = 0; jtrk < e_minus.size(); jtrk++ ){

					double chi2_e = mNSigmaE1[itrk]*mNSigmaE1[itrk] + mNSigmaE2[jtrk]*mNSigmaE2[jtrk];
					double chi2_p = mNSigmaPi1[itrk]*mNSigmaPi1[itrk] + mNSigmaPi2[jtrk]*mNSigmaPi2[jtrk];
					h_NsigmaElec->Fill(mNSigmaE1[itrk], mNSigmaE2[jtrk], weight);
					h_NsigmaPion->Fill(mNSigmaPi1[itrk], mNSigmaPi2[jtrk], weight );
					h_NsigmaPart->Fill(chi2_p, chi2_e, weight);
					h_NsigmaCorr->Fill(mNSigmaE1[itrk], mNSigmaPi1[itrk], weight);
					h_NsigmaCorr->Fill(mNSigmaE2[jtrk], mNSigmaPi2[jtrk], weight);
					
					//embedding has correct nsigma e and pi
					if( chi2_p<elecPIDcut_b && (chi2_e/chi2_p) > elecPIDcut_a/elecPIDcut_b ) continue;
					if( chi2_p>elecPIDcut_b && chi2_e > elecPIDcut_a ) continue;
				
					h_EoverP_BEMC->Fill(mBEMCEnergy1[itrk]/mBEMCPmag1[itrk], weight);
					h_EoverP_BEMC->Fill(mBEMCEnergy2[jtrk]/mBEMCPmag2[jtrk], weight);

					TLorentzVector e_inv;
					e_inv = e_plus[itrk]+e_minus[jtrk];
					if( fabs(e_inv.Rapidity()) > 1.0 ) continue;
					if( e_inv.Pt() < 0.0 ) continue;
					
					hElecPt->Fill( e_plus[itrk].Pt(), weight );
					hElecPt->Fill( e_minus[jtrk].Pt(), weight );
					hElecEta->Fill( e_plus[itrk].Eta(), weight );
					hElecEta->Fill( e_minus[jtrk].Eta(), weight );
					hElecPhi->Fill( e_plus[itrk].Phi(), weight );
					hElecPhi->Fill( e_minus[jtrk].Phi(), weight );
					
					hJpsiMass->Fill( e_inv.M(), weight );
					hJpsiMass_Pt2->Fill(e_inv.M(), e_inv.Pt()*e_inv.Pt(), weight);
					hDielectronPt->Fill(e_inv.Pt(), weight);
					hDielectronPt2->Fill(e_inv.Pt()*e_inv.Pt(), weight);
					//within neutron peak;
					if( totalZDCWestADC < 40 || (totalZDCWestADC > 9999.) ){
						hJpsiMassZDCveto->Fill( e_inv.M(), weight );
						hJpsiMassZDCveto_Pt2->Fill(e_inv.M(), e_inv.Pt()*e_inv.Pt(), weight);
					}
					else{
						hJpsiMassZDC->Fill( e_inv.M(), weight );
						hJpsiMassZDC_Pt2->Fill(e_inv.M(), e_inv.Pt()*e_inv.Pt(), weight);
					}
				}
			}

		}//nvertex
	}

	
	fout->Write();
	fout->Close();


}