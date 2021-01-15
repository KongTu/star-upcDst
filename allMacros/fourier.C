// ==============================================================================
//  File: fourier.C
//
//  Macro to derive source distribution from dsig/dt distribution
//
//  Author: Thomas Ullrich
//  Last update: September 25, 2012
//
//  For details on method see M. Diehl's talk at INT workshop
//  Nov 1, 2010 "How well one needs to measure t for getting images in b space"
//  mainly page 3,4
// ==============================================================================
#include <iostream>
#include "TObject.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
#include "Math/GaussIntegrator.h"
#include "TPolyMarker.h"
#include "TMinuitMinimizer.h"
#include "TRandom3.h"

#include "RiceStyle.h"

#define PR(x) cout << #x << " = " << (x) << endl;
#define PRP(x) cout << #x << " = " << (void*)(x) << endl;

using namespace std;

//
//  Flags
//
bool PlotErrors = true;
bool MakeErrorBand = true;
bool SavePlots = false;
bool GaussianSmearing = false; // careful, when on peak finder goes nuts
bool OffsetCorrection = true; 

//
//  Constants
//
const double hbarc = 0.197327; // GeV*fm
const double hbarc2 = hbarc*hbarc;
const double twopi = 2*acos(-1.);
// Au Woods-Saxon parameter
const double mRho0 = 0.1693;   
const double mOmega = 0;  
const double mRadius = 6.38;   
const double mSurfaceThickness = 0.535;   
// Scale error bars here (e.g. to test different L)
const double scaleErrors = 1; // 3.16227766016838; // sqrt(10);

//
//  Prototypes
//
unsigned int readTable(const char*, double*, double*, double *, unsigned int);
double fun(double *x, double *par);
void giveMeTheSource(TH1D*&, TH1D*&, TH1D*&,const char*);
double rho(double b, double z);
double rhoForIntegration(double *x, double* par) ;
double T(double *x, double *);

//
//  Global variables
//
TH1D* globalHisto;
vector<double> globalPeaks;
TH1D* globalLookupTable;
double globalUpperTValue;
TRandom3 rndm;


//functional form


void fourier(double tmax = -1) 
{  
    //
    //  Arguments
    //
    globalUpperTValue = tmax;
    
    //
    //  General ROOT options
    //
    int ifont    = 42;     // text font
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    
    gStyle->SetTextFont(ifont);
    gStyle->SetLabelFont(ifont);
    gStyle->SetTitleTextColor(1);
    gStyle->SetFillColor(1);
    gStyle->SetStatColor(1);
    gStyle->SetTitleColor(1);
    gStyle->SetTitleBorderSize(1);
    

    unsigned int n_jpsi_upc = 50;
    double t_jpsi_upc[200];
    double y_jpsi_upc[200];
    double ey_jpsi_upc[200];

    TF1* coh = new TF1("coh","170*TMath::Exp(-6.8*x[0])",0,5);
    for(int i=0;i<50;i++){
        t_jpsi_upc[i] = 0+i*0.05;
        y_jpsi_upc[i] = coh->Eval(t_jpsi_upc[i]);
        ey_jpsi_upc[i] = 0.01*y_jpsi_upc[i];
    }

    //
    //  Fill histograms
    //
    double binsize = fabs(t_jpsi_upc[2]-t_jpsi_upc[1]); 
    double abstmin = 0.025;
    double abstmax = t_jpsi_upc[n_jpsi_upc-1]+binsize/2;
    int nbins = static_cast<int> (abstmax/binsize + 0.025);

    TH1D *h_jpsi_upc = new TH1D("h_jpsi_upc", "jpsi sat", nbins, abstmin, abstmax);
    h_jpsi_upc->GetXaxis()->SetTitle("|t| (GeV^{2})");

    for (unsigned int i=0; i<n_jpsi_upc; i++) {
        int k = h_jpsi_upc->FindBin(t_jpsi_upc[i]);
        h_jpsi_upc->SetBinContent(k, y_jpsi_upc[i]);
        h_jpsi_upc->SetBinError(k, ey_jpsi_upc[i]);
        double x = h_jpsi_upc->GetBinCenter(k);
        if (fabs(x-t_jpsi_upc[i]) > binsize/100) {
            cout << "fourier:: Warning, table entry and bin size are not matching for jpsi sat." << endl;
            exit(1);
        }
    }
    
    //
    //  Plot what we got
    //
    TCanvas *c1 = new TCanvas("c1","dsig/dt", 600, 600);
    c1->SetBorderSize(1);
    c1->SetFillColor(0);
    h_jpsi_upc->Draw("PE");
    gPad->SetLogy(1); 
    if (SavePlots) c1->SaveAs("plot1.eps");
    
    
    
    TH1D *source_jpsi_sat = 0;
    TH1D *source_jpsi_sat_band = 0;
    TH1D *source_jpsi_sat_upper = 0;
    TH1D *source_jpsi_sat_lower = 0;
    globalHisto = h_jpsi_upc;

    giveMeTheSource(source_jpsi_sat, source_jpsi_sat_upper, source_jpsi_sat_lower,  "J/psi dAu");
    source_jpsi_sat_band = new TH1D(*source_jpsi_sat);
    source_jpsi_sat_band = new TH1D("source_jpsi_sat_band", "",
                                    source_jpsi_sat->GetNbinsX(), 
                                    source_jpsi_sat->GetXaxis()->GetXmin(), 
                                    source_jpsi_sat->GetXaxis()->GetXmax());
    double xnorm = 0.;
    xnorm = source_jpsi_sat->Integral("width");
    cout << "Integral over jpsi sat source function: " << xnorm << endl;    
    source_jpsi_sat->Scale(1/xnorm);
    xnorm = source_jpsi_sat_upper->Integral("width");
    cout << "Integral over jpsi sat source function (upper): " << xnorm << endl;    
    source_jpsi_sat_upper->Scale(1/xnorm);
    xnorm = source_jpsi_sat_lower->Integral("width");
    cout << "Integral over jpsi sat source function (lower): " << xnorm << endl;    
    source_jpsi_sat_lower->Scale(1/xnorm);

    source_jpsi_sat_upper->SetLineStyle(1);
    source_jpsi_sat_upper->SetLineColor(1);
    source_jpsi_sat_lower->SetLineStyle(2);
    source_jpsi_sat_lower->SetLineColor(2);

    //  Generate lookup table for overlap integral T  
    //  (b and z in fm)   
    double range = 3.;  
    TF1* funcToIntegrate = new TF1("funcToIntegrate", rhoForIntegration, -10, 10, 1);  
    int numbins = 100; // size of lookup table  
    globalLookupTable = new TH1D("globalLookupTable","T lookup table", numbins, -3, 3); 
    for (int i=0; i<=numbins; i++) {  
        double b = globalLookupTable->GetBinCenter(i+1);  
        funcToIntegrate->SetNpx(1000);  
        funcToIntegrate->SetParameter(0, b);  
        double res = funcToIntegrate->Integral(-10, 10); 
        globalLookupTable->SetBinContent(i+1, res);      
    }  
    delete funcToIntegrate; 

    if (MakeErrorBand) {
        for (int i=1; i <= source_jpsi_sat->GetNbinsX(); i++) {
            double ymax = source_jpsi_sat->GetBinContent(i);
            double ymin = ymax;
            ymax = max(ymax, source_jpsi_sat_lower->GetBinContent(i));
            ymax = max(ymax, source_jpsi_sat_upper->GetBinContent(i));
            ymin = min(ymin, source_jpsi_sat_lower->GetBinContent(i));
            ymin = min(ymin, source_jpsi_sat_upper->GetBinContent(i));
            double ypos = (ymax+ymin)/2;
            source_jpsi_sat_band->SetBinContent(i, ypos);
            source_jpsi_sat_band->SetBinError(i, (ymax-ypos));
        }

        //kong's addition to surpress huge error:
        for(int i=0;i<source_jpsi_sat->GetNbinsX();i++){
            source_jpsi_sat->SetBinError(i+1, 0.);
        }
    }

    TCanvas* c2 = new TCanvas("c2","c2",600,600);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    
    TH1D* base1 = makeHist("base1", "", "b (fm)", "F(b)", 100,-3,3,kBlack);
    base1->GetYaxis()->SetRangeUser(0, 1);
    base1->GetXaxis()->SetTitleColor(kBlack);
    
    fixedFontHist1D(base1,1.2,1.25);

    base1->GetYaxis()->SetTitleOffset(1.3);
    base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.2);
    base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.2);
    base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
    base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);
    base1->GetXaxis()->SetNdivisions(4,6,0);
    base1->GetYaxis()->SetNdivisions(4,6,0);
    base1->Draw();

    double yMaximum = 1;
    // source_jpsi_sat->Scale(1./(source_jpsi_sat->Integral(source_jpsi_sat->FindBin(-2.4),source_jpsi_sat->FindBin(2.4))));
    source_jpsi_sat_band->Scale(1./(source_jpsi_sat_band->Integral(source_jpsi_sat_band->FindBin(-2.4),source_jpsi_sat_band->FindBin(2.4))));
    globalLookupTable->Scale((source_jpsi_sat->Integral(source_jpsi_sat->FindBin(-2.4),source_jpsi_sat->FindBin(2.4)))/(globalLookupTable->Integral(globalLookupTable->FindBin(-2.4),globalLookupTable->FindBin(2.4))));
    //kong's addition to supress huge error:
    for(int i=0;i<globalLookupTable->GetNbinsX();i++){
        globalLookupTable->SetBinError(i+1, 0.);
    }
    source_jpsi_sat->SetTitle("");
    source_jpsi_sat->SetLineColor(kBlue);
    source_jpsi_sat->GetYaxis()->SetNdivisions(5,5,0);
    source_jpsi_sat->GetYaxis()->SetTitle("F(b)");
    source_jpsi_sat->GetXaxis()->SetTitle("b (fm)");
    source_jpsi_sat->GetXaxis()->CenterTitle();
    source_jpsi_sat->GetYaxis()->CenterTitle();
    source_jpsi_sat->GetXaxis()->SetTitleOffset(1.5*source_jpsi_sat->GetXaxis()->GetTitleOffset());
    source_jpsi_sat->GetYaxis()->SetTitleOffset(1.5*source_jpsi_sat->GetYaxis()->GetTitleOffset());
    source_jpsi_sat->SetMarkerStyle(20);
    source_jpsi_sat->SetMarkerColor(kBlue);

    source_jpsi_sat->Draw("Psame");
    globalLookupTable->SetMarkerStyle(27);
    globalLookupTable->Draw("Psame");

    TLegend *w4 = new TLegend(0.16,0.75,0.45,0.86);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(20);
    w4->SetTextFont(45);
    w4->AddEntry(source_jpsi_sat, "Gluon density  ", "P");
    w4->AddEntry(globalLookupTable, "Charge density ", "P");
    w4->Draw("same");

    TLatex* r42 = new TLatex(0.18, 0.91, "Deuteron");
    r42->SetNDC();
    r42->SetTextSize(22);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);

    TLatex* r43 = new TLatex(0.67,0.91, "STAR");
    r43->SetNDC();
    r43->SetTextFont(62);
    r43->SetTextSize(0.04);

    TLatex* r44 = new TLatex(0.78,0.91, "Internal");
    r44->SetNDC();
    r44->SetTextSize(21);
    r44->SetTextFont(53);



    r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");

    TFile* output = new TFile("glue_minus.root","RECREATE");
    source_jpsi_sat->Write();
    globalLookupTable->Write();

    // sarteWS->Draw("");
    // if (PlotErrors) source_jpsi_sat_upper->Draw("same");
    // if (PlotErrors) source_jpsi_sat_lower->Draw("same");
    // if (MakeErrorBand) {
        // source_jpsi_sat_band->Draw("E2 same");
        // source_jpsi_sat_band->SetFillColor(kGray);
        // source_jpsi_sat_band->SetMarkerStyle(0);
    // }


}

void giveMeTheSource(TH1D*& sourceHisto, TH1D*& sourceHistoErrorUp, TH1D*& sourceHistoErrorLow, const char* txt)
{
    globalPeaks.clear();
    
    //
    //   Peak finding
    //  
    //   Tricky!
    //   First we invert the spectra so that the minima
    //   become maxima. Then we search the peaks and display them.
    //   Wouldn't trust them out of the box. Parameters to 
    //   ShowPeaks might need adjusment. 
    //
    TH1D *histoPeak = new TH1D("histoPeak","Peak finding", globalHisto->GetNbinsX(), globalHisto->GetXaxis()->GetXmin(), globalHisto->GetXaxis()->GetXmax());
    for (int i=1; i<=globalHisto->GetNbinsX(); i++) {
        double y = globalHisto->GetBinContent(i);
        if (! (y>0) ) continue;
        double newy = -log10(y);
        histoPeak->SetBinContent(i, newy);
    }
    histoPeak->Draw();

    //
    //  Getting the peak positions is ugly
    //
    //  ROOT says about ShowPeaks():
    //
    //  Int_t ShowPeaks(Double_t sigma = 2, Option_t* option = "", Double_t threshold = 0.05)
    //  Interface to TSpectrum::Search.
    //  The function finds peaks in this histogram where the width is > sigma
    //  and the peak maximum greater than threshold*maximum bin content of this.
    //  For more details see TSpectrum::Search.
    //  Note the difference in the default value for option compared to TSpectrum::Search
    //  option="" by default (instead of "goff").    //
    
    histoPeak->ShowPeaks(0.2, "", 3); // arguments might need to be adjusted
    
    TList *functions = histoPeak->GetListOfFunctions();
    TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
    double *tmppeaks = pm->GetX();
    for (int i=0; i<pm->GetN(); i++) 
        globalPeaks.push_back(tmppeaks[i]);
    globalPeaks.push_back(globalHisto->GetXaxis()->GetXmin());
    globalPeaks.push_back(globalHisto->GetXaxis()->GetXmax());
    sort( globalPeaks.begin(), globalPeaks.end() );  // absolutely necessary
    cout << "Found " << pm->GetN() << " peaks in dsig/dt of" << txt << ": ";
    for (unsigned int i=0; i<globalPeaks.size(); i++) cout << globalPeaks[i] << '\t';
    cout << endl;
    
    //
    //  Setup function for integration
    //
    TF1 ffun("ffun", fun, globalHisto->GetXaxis()->GetXmin(), globalHisto->GetXaxis()->GetXmax(), 2);
    
    //
    //  Histogram holding the result, F(b)
    //
    char htext[100];
    sprintf(htext, "F(b) %s", txt);
    sourceHisto = new TH1D("sourceHisto", htext, 100, -3, 3);
    sourceHistoErrorUp = new TH1D("sourceHistoErrorUp", htext, 100, -3, 3);
    sourceHistoErrorLow = new TH1D("sourceHistoErrorLow", htext, 100, -3, 3);
    
    //
    //  Generate F(b) - GSLIntegrator is probably the best one available (?)
    //
    for (int i=1; i<=sourceHisto->GetNbinsX(); i++) { // loop over all bins from h1
        double b = sourceHisto->GetBinCenter(i);
        ffun.SetParameter(0, fabs(b));
        ffun.SetParameter(1, 0);
        ROOT::Math::WrappedTF1 wf1(ffun);
        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
        ig.SetFunction(wf1);
        ig.SetRelTolerance(0.001);
        double F = ig.Integral(0, sqrt(globalHisto->GetXaxis()->GetXmax())); //?
        sourceHisto->SetBinContent(i, F);
    }
    
    //
    // Now the errors - take the two extremes
    //
    for (int i=1; i<=sourceHisto->GetNbinsX(); i++) {
        double b = sourceHisto->GetBinCenter(i);
        ffun.SetParameter(0, fabs(b));
        ffun.SetParameter(1, 1); // !
        ROOT::Math::WrappedTF1 wf1(ffun);
        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
        ig.SetFunction(wf1);
        ig.SetRelTolerance(0.001);
        double F = ig.Integral(0, sqrt(globalHisto->GetXaxis()->GetXmax()));
        sourceHistoErrorUp->SetBinContent(i, F);
    }
    for (int i=1; i<=sourceHisto->GetNbinsX(); i++) {
        double b = sourceHisto->GetBinCenter(i);
        ffun.SetParameter(0, fabs(b));
        ffun.SetParameter(1, 2);  // !
        ROOT::Math::WrappedTF1 wf1(ffun);
        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
        ig.SetFunction(wf1);
        ig.SetRelTolerance(0.001);
        double F = ig.Integral(0, sqrt(globalHisto->GetXaxis()->GetXmax()));
        sourceHistoErrorLow->SetBinContent(i, F);
    }
    
    //
    //  Cheat add offset to get rid of negative F(b) - effect is tiny
    //
    if (OffsetCorrection) {
        double ymin = sourceHisto->GetBinContent(sourceHisto->GetMinimumBin());
        if (ymin < 0) {
            cout << "giveMeTheSource: add offest of " << ymin << " to F(b)" << endl;
            for (int i=1; i<=sourceHisto->GetNbinsX(); i++) { // loop over all bins from h1
                double y = sourceHisto->GetBinContent(i);
                sourceHisto->SetBinContent(i, y+(fabs(ymin)));
            }
        }
        ymin = sourceHistoErrorUp->GetBinContent(sourceHistoErrorUp->GetMinimumBin());
        if (ymin < 0) {
            cout << "giveMeTheSource: add offest of " << ymin << " to F(b) upper" << endl;
            for (int i=1; i<=sourceHistoErrorUp->GetNbinsX(); i++) { // loop over all bins from h1
                double y = sourceHistoErrorUp->GetBinContent(i);
                sourceHistoErrorUp->SetBinContent(i, y+(fabs(ymin)));
            }
        }
        ymin = sourceHistoErrorLow->GetBinContent(sourceHistoErrorLow->GetMinimumBin());
        if (ymin < 0) {
            cout << "giveMeTheSource: add offest of " << ymin << " to F(b) lower" << endl;
            for (int i=1; i<=sourceHistoErrorLow->GetNbinsX(); i++) { // loop over all bins from h1
                double y = sourceHistoErrorLow->GetBinContent(i);
                sourceHistoErrorLow->SetBinContent(i, y+(fabs(ymin)));
            }
        }
    }
    
}

unsigned int readTable(const char* filename, double *t, double *y, double *ey, unsigned int n)
{
    ifstream ifs(filename);
    if (!ifs) {
        cout << "readTable: Cannot open file '" << filename << "'." << endl;
        return 0;
    }
    else {
        cout << "readTable: Reading file '" << filename << "'." << endl;
    }
    string line;
    unsigned int ncount = 0;
    
    while (ifs.good() && !ifs.eof()) {
        getline(ifs, line);
        size_t pos = line.find_first_of('#');
        if (pos <= line.size())
            line = line.erase(pos, line.length());
        if (line.empty()) continue;
        istringstream iss (line,istringstream::in);
        iss >> t[ncount] >> y[ncount] >> ey[ncount];
        // Manipulating points
        ey[ncount] *= scaleErrors; // scale errors by global scale (test different L)
        if (GaussianSmearing) {
            y[ncount] = rndm.Gaus(y[ncount], ey[ncount]);
        }
        ncount++;
        if (ncount >= n) {
            cout << "readTable: Warning: not enough room to store all tabel entries." << endl;
            break;
        }
    }
    cout << "readTable: All done. Entries read: " << ncount << endl;
    
    ifs.close();
    return ncount;
}

double fun(double *x, double *par)
{
    double Delta = *x;  // GeV
    double b = par[0]; // fm
    int iflag;  // 0 = data, 1=data+error, 2=data-error
    if (par[1] < 0.5) 
        iflag = 0;
    else if (par[1] > 0.5 && par[1] < 1.5)
        iflag = 1;
    else if (par[1] > 1.5 && par[1] < 2.5)
        iflag = 2;
    else {
        cout << "fun(): Warning, unknown iflag value (" << iflag << ")" << endl;
        exit(1);
    }
    double result = Delta/twopi;
    double t = Delta*Delta;
    
    //
    //  Cut everything above t > globalUpperTValue (used in systematic studies)
    //  globalUpperTValue = -1 means use full range
    //
    if (globalUpperTValue > 0 && t > globalUpperTValue) return 0;  // globalUpperTValue = -1 means no cuts
        
    double sigma = 0;
    if (t <= globalHisto->GetXaxis()->GetXmax()) {
        sigma = globalHisto->Interpolate(t);
        if (iflag != 0) {
            int k = globalHisto->FindBin(t);
            if (iflag == 1) 
                sigma += globalHisto->GetBinError(k);
            else if (iflag == 2)
                sigma -= globalHisto->GetBinError(k);
        }
        if (sigma < 0) sigma = 0;
    }
    
    int sign = 1;
    for (unsigned int i=1; i<globalPeaks.size(); i++) {
        if (t >= globalPeaks[i-1] && t < globalPeaks[i] ) break;
        sign *= -1;
    }
    
    result *= sign*sqrt(sigma)*TMath::BesselJ0(Delta*b/hbarc);
    return result;
}

double rho(double b, double z) // returns value in units of fm^-3  
{  
    // b and z in fm  
    double r = sqrt(b*b+z*z); 
    double A = 0.456;
    double B = 2.36; 
    return 0.39946312*TMath::Power( (TMath::Exp(-A*r)-TMath::Exp(-B*r))/r , 2 );
    // return mRho0*(1+mOmega*(r/mRadius)*(r/mRadius))/(1+exp((r-mRadius)/mSurfaceThickness));   
}  

double rhoForIntegration(double *x, double* par)  
{  
    // b and z in fm  
    double b = par[0];  
    double z = *x;  
    return rho(b, z);  
}  

double T(double *x, double *par)
{  
    //   
    //  Returns overlap integral   
    //  b in fm, overlap function T in GeV^2  
    //  HERE: scaled by par[0] !!! tu  
    //  
    double b = fabs(x[0]);
    int bin = globalLookupTable->FindBin(b);  
    double res = globalLookupTable->GetBinContent(bin);  
    return res*hbarc2*par[0];  
}  

