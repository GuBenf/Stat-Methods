bool USE_MC = true;
std::vector<double> x_var;

Double_t MASS_D_PLUS_PHI_PI_LOW = 1.84;
Double_t MASS_D_PLUS_PHI_PI_HIGH = 1.9;
int MASS_D_PLUS_PHI_PI_BINS = 120;

Double_t MASS_DS_PLUS_PHI_MU_NU_LOW = 1.6;
Double_t MASS_DS_PLUS_PHI_MU_NU_HIGH = 2.;
int MASS_DS_PLUS_PHI_MU_NU_BINS = 400;

Double_t MASS_DS_PLUS_PHI_PI_LOW = 1.94;
Double_t MASS_DS_PLUS_PHI_PI_HIGH = 2.;
int MASS_DS_PLUS_PHI_PI_BINS = 120;

Double_t MASS_DS_PLUS_TAU_NU_LOW = 1.75;
Double_t MASS_DS_PLUS_TAU_NU_HIGH = 1.81;
int MASS_DS_PLUS_TAU_NU_BINS = 120;

Double_t TOT_MASS_LOW = 1.6;
// Double_t TOT_MASS_HIGH = 2.6;
Double_t TOT_MASS_HIGH = 2.;

int TOT_MASS_BINS = 400;

//Double_t MASS_COMBINATORIA_LOW = 2.;
//Double_t MASS_COMBINATORIAL_HIGH = 3.;
//int MASS_COMBINATORIAL_BIN = 100.;

int ID_D_PLUS_PHI_PI = 34;
int ID_DS_PLUS_PHI_MU_NU = 42;
int ID_DS_PLUS_PHI_PI = 41;
int ID_DS_PLUS_TAU_NU = 44;

std::vector<double> vec_mass_Dplus_PhiPi_MC;
std::vector<double> vec_mass_DSplus_PhiMuNu_MC;
std::vector<double> vec_mass_DSplus_PhiPi_MC;
std::vector<double> vec_mass_DSplus_TauNu_MC;
std::vector<double> vec_tot_mass_MC;
//std::vector<double> vec_mass_combinatorial;
std::vector<double> vec_tot_mass_DATA;

//std::vector<double> x_var;

Double_t MU1_FIT_D_PLUS_PHI_PI = 1.86969e+00;
Double_t MU2_FIT_D_PLUS_PHI_PI = 1.87035e+00;
Double_t SIGMA1_FIT_D_PLUS_PHI_PI = 5.79484e-03;
Double_t SIGMA2_FIT_D_PLUS_PHI_PI = 8.15639e-03;
Double_t FRACTION_FIT_D_PLUS_PHI_PI = 8.83747e-01;

// Double_t MU1_FIT_D_PLUS_PHI_PI = 1.86970e+00;
// Double_t MU2_FIT_D_PLUS_PHI_PI = 1.87050e+00;
// Double_t SIGMA1_FIT_D_PLUS_PHI_PI = 5.84594e-03;
// Double_t SIGMA2_FIT_D_PLUS_PHI_PI = 8.60173e-03;
// Double_t FRACTION_FIT_D_PLUS_PHI_PI = 9.17802e-01;


Double_t M0_FIT_DS_PLUS_PHI_MU_NU = 1.97801e+00;
Double_t C_FIT_DS_PLUS_PHI_MU_NU = -3.82391e+00;
Double_t P_FIT_DS_PLUS_PHI_MU_NU = 1.41325e+00;

Double_t MU1_FIT_DS_PLUS_PHI_PI = 1.96843e+00;
Double_t MU2_FIT_DS_PLUS_PHI_PI = 1.96807e+00;
Double_t SIGMA1_FIT_DS_PLUS_PHI_PI = 6.08795e-03;
Double_t SIGMA2_FIT_DS_PLUS_PHI_PI = 7.87780e-03;
Double_t FRACTION_FIT_DS_PLUS_PHI_PI = 8.65940e-01;

Double_t MU1_FIT_MASS_DS_PLUS_TAU_NU = 1.77687e+00;
Double_t MU2_FIT_MASS_DS_PLUS_TAU_NU = 1.77748e+00;
Double_t SIGMA1_FIT_MASS_DS_PLUS_TAU_NU = 5.42707e-03;
Double_t SIGMA2_FIT_MASS_DS_PLUS_TAU_NU = 7.16697e-03;
Double_t FRACTION_FIT_MASS_DS_PLUS_TAU_NU = 7.88093e-01;

Double_t TAU_FIT_COMBINATORIAL = 1.0;


#define analysis_cxx
#include "analysis.h"
#include "fit_unbinned_InvariantMass.h"
#include "fit_unbinned_Argus.h"
#include "fit_unbinned_TotalSpectrum.h"
//#include "fit_unbinned_Combinatorial.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

Double_t modelFormula(double *xx, double *p){
  Double_t x = xx[0];

  // Signal
  Double_t sig = p[2] * TMath ::Gaus(x,p[3] ,p[4] ,1);

  // DS_PHI_PI
  Double_t ds_phipi = p[5]* (FRACTION_FIT_DS_PLUS_PHI_PI*TMath ::Gaus(x,MU1_FIT_DS_PLUS_PHI_PI,SIGMA1_FIT_DS_PLUS_PHI_PI,1) + (1-FRACTION_FIT_DS_PLUS_PHI_PI)* TMath ::Gaus(x,MU2_FIT_DS_PLUS_PHI_PI,SIGMA2_FIT_DS_PLUS_PHI_PI,1));

  // D_PHI_PI
  Double_t d_phipi = p[6]* (FRACTION_FIT_D_PLUS_PHI_PI*TMath ::Gaus(x,MU1_FIT_D_PLUS_PHI_PI,SIGMA1_FIT_D_PLUS_PHI_PI,1) + (1-FRACTION_FIT_D_PLUS_PHI_PI)* TMath ::Gaus(x,MU2_FIT_D_PLUS_PHI_PI,SIGMA2_FIT_D_PLUS_PHI_PI,1));

  // DS_PHI_MU_NU
  Double_t argus_par[2] = {C_FIT_DS_PLUS_PHI_MU_NU,P_FIT_DS_PLUS_PHI_MU_NU};
  Double_t ds_phimunu = p[7]* pdf_Argus(&x,argus_par);

  // Combinatorial
  Double_t comb = (1-p[2]-p[5]-p[6]-p[7])* _combinatorial_pdf(x,&TAU_FIT_COMBINATORIAL);

  // Total pdf
  //
  Double_t tot_pdf = p[0]*p[1]*(sig + ds_phipi + d_phipi + ds_phimunu + comb);

  return tot_pdf;
}


void analysis::Loop () 
{

  gRandom ->SetSeed (12345);
  const double mMin = 1.60;
  const double mMax = 2.10;
  const int nBins = 50;
  const double binWidth = (mMax - mMin) / nBins;
  const int nToys = 1000;

  const double nTotTrue = 500.0;
  const double fSigTrue = 0.05;
  const double meanTrue = 1.78;
  const double sigmaTrue = 0.02;
  const double fDsPhiPi = 0.04;
  const double fDPhiPi = 0.0215;
  const double fDsPhiMuNu = 0.653;

  // const char* modelFormula = "[0]*[5]*([1]* TMath ::Gaus(x,[2] ,[3] ,1) + (1.0 -[1]) *([4]* TMath::Exp([4]*(x -1.60))/(TMath::Exp ([4]*(2.10 -1.60)) -1.0)))";

  // Known model used to generate the pseudo - experiments .
  // TF1 :: GetRandom uses only the shape of the function.
  TF1 generatorModel("generatorModel", modelFormula , mMin , mMax, 8);
  generatorModel.SetParameters(nTotTrue, binWidth, fSigTrue, meanTrue, sigmaTrue, fDsPhiPi, fDPhiPi, fDsPhiMuNu);
  TH1D hfSigFit("hfSigFit", "Fitted signal yield;#hat{f}_{sig};Pseudo -experiments", 50, 0.0, 1.0);

  for (int iToy = 0; iToy < nToys; ++iToy) 
  {
    TH1D hToy("hToy", "Toy dataset;M;Candidates", nBins , mMin , mMax);
    // const int nObs = gRandom ->Poisson(nTotTrue); // Extended toy generation .
    const int nObs = nTotTrue;
    for (int i = 0; i < nObs; ++i) 
    {
      const double m = generatorModel.GetRandom(mMin , mMax);
      hToy.Fill(m);
    }
    // Fit model:
    // par [0] = Ntot
    // par [1] = fsig
    // par [2] = Gaussian mean
    // par [3] = Gaussian width
    // par [4] = exponential slope
    // par [5] = bin width , fixed
    
    TF1 fitModel("fitModel", modelFormula , mMin , mMax, 8);
    fitModel.SetParameters(nObs , binWidth, fSigTrue, meanTrue, sigmaTrue, fDsPhiPi, fDPhiPi, fDsPhiMuNu);
    fitModel.SetParNames("Ntot", "binWidth", "fsig", "mean", "sigma", "fDsPhiPi", "fDPhiPi", "fDsPhiMuNu");
    fitModel.FixParameter (0, nObs); // or comment to let it float
    fitModel.FixParameter (1, binWidth);
    hToy.Fit(&fitModel , "LI");
    const double nTotFit = fitModel.GetParameter (0);
    const double fSigFit = fitModel.GetParameter (2);
    const double nSigFit = fSigFit * nTotFit;
    hfSigFit.Fill(fSigFit);
  }
  TCanvas c("c", "Pseudo -experiment result", 800, 600);
  hfSigFit.Draw("HIST");
  c.SaveAs("toy_fitted_signal_yield.pdf");
  }
/*
void analysis::Loop()
{

  TFile * outfile = new TFile("outfile.root","RECREATE");

  TH1D * mass_Dplus_PhiPi_MC = new TH1D("mass_Dplus_PhiPi_MC","",MASS_D_PLUS_PHI_PI_BINS,MASS_D_PLUS_PHI_PI_LOW,MASS_D_PLUS_PHI_PI_HIGH);
  TH1D * mass_DSplus_PhiMuNu_MC = new TH1D("mass_DSplus_PhiMuNu_MC","",MASS_DS_PLUS_PHI_MU_NU_BINS,MASS_DS_PLUS_PHI_MU_NU_LOW,MASS_DS_PLUS_PHI_MU_NU_HIGH);
  TH1D * mass_DSplus_PhiPi_MC = new TH1D("mass_DSplus_PhiPi_MC","",MASS_DS_PLUS_PHI_PI_BINS,MASS_DS_PLUS_PHI_PI_LOW,MASS_DS_PLUS_PHI_PI_HIGH);
  TH1D * mass_DSplus_TauNu_MC = new TH1D("mass_DSplus_TauNu_MC","",MASS_DS_PLUS_TAU_NU_BINS,MASS_DS_PLUS_TAU_NU_LOW,MASS_DS_PLUS_TAU_NU_HIGH);
  TH1D * tot_mass_MC = new TH1D("tot_mass_MC","",TOT_MASS_BINS,TOT_MASS_LOW,TOT_MASS_HIGH);
  //TH1D * mass_combinatorial_DATA = new TH1D("mass_combinatorial_DATA","",MASS_COMBINATORIAL_BIN,MASS_COMBINATORIA_LOW,MASS_COMBINATORIAL_HIGH);
  TH1D * tot_mass_DATA = new TH1D("tot_mass_DATA","",TOT_MASS_BINS,TOT_MASS_LOW,TOT_MASS_HIGH);
  
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    Double_t m_tau_mc;
    Double_t m_sigma_tau_mc;
    m_tau_mc= 1.77699;
    m_sigma_tau_mc=  0.00581298;

    if(id == ID_D_PLUS_PHI_PI && D_M >= MASS_D_PLUS_PHI_PI_LOW && D_M <= MASS_D_PLUS_PHI_PI_HIGH)
    {
      mass_Dplus_PhiPi_MC->Fill(D_M);
      vec_mass_Dplus_PhiPi_MC.push_back(D_M);
    }
    if(id == ID_DS_PLUS_PHI_MU_NU && D_M >= MASS_DS_PLUS_PHI_MU_NU_LOW && D_M <= MASS_DS_PLUS_PHI_MU_NU_HIGH)
    {
      mass_DSplus_PhiMuNu_MC->Fill(D_M);
      vec_mass_DSplus_PhiMuNu_MC.push_back(D_M);
    }
    if(id == ID_DS_PLUS_PHI_PI && D_M >= MASS_DS_PLUS_PHI_PI_LOW && D_M <= MASS_DS_PLUS_PHI_PI_HIGH)
    {
      mass_DSplus_PhiPi_MC->Fill(D_M);
      vec_mass_DSplus_PhiPi_MC.push_back(D_M);
    }
    if(id == ID_DS_PLUS_TAU_NU && D_M >= MASS_DS_PLUS_TAU_NU_LOW && D_M <= MASS_DS_PLUS_TAU_NU_HIGH)
    {
      mass_DSplus_TauNu_MC->Fill(D_M);
      vec_mass_DSplus_TauNu_MC.push_back(D_M);
    }
    if(id != 0 && D_M >= TOT_MASS_LOW && D_M <= TOT_MASS_HIGH)
    {
      tot_mass_MC->Fill(D_M);
      vec_tot_mass_MC.push_back(D_M);
    }
     */

    // if(/*id == 0 &&*/ D_M >= TOT_MASS_LOW && D_M <= TOT_MASS_HIGH)
    // {
    //   vec_tot_mass_DATA.push_back(D_M);
    // }

    // if(/*id == 0 &&*/ TMath::Abs(D_M-m_tau_mc)>3*m_sigma_tau_mc && D_M >= TOT_MASS_LOW && D_M <= TOT_MASS_HIGH) //DON'T LOOK IN THE HISTO!!
    // {
    //   tot_mass_DATA->Fill(D_M);
    // }; //Mass blinding

    //if(id == 0 && TMath::Abs(D_M-m_tau_mc)>3*m_sigma_tau_mc && D_M >= MASS_COMBINATORIA_LOW && D_M <= MASS_COMBINATORIAL_HIGH) 
    //{
    //  vec_mass_combinatorial.push_back(D_M);
    //  mass_combinatorial_DATA->Fill(D_M);
    //}; //Mass blinding
     
  // }; //End Loop Events

  // cout << endl << "FIT INVARIANT MASS D_PLUS -> PHI PI" << endl << endl;
  // fit_unbinned_InvariantMass(vec_mass_Dplus_PhiPi_MC, 
  //                           {0.5, 1.86, 1.88, 0.005, 0.007}, 
  //                           {0.001, 0.001, 0.001, 0.001, 0.001}, 
  //                           {0., 0., 0., 0., 0.}, 
  //                           {1., 3., 3., 1., 1.}, 
  //                           {"fraction","mu1","mu2","sigma1","sigma2"}, 
  //                           mass_Dplus_PhiPi_MC, 
  //                           outfile, 
  //                           "D_Plus_Phi_Pi", 
  //                           MASS_D_PLUS_PHI_PI_LOW, 
  //                           MASS_D_PLUS_PHI_PI_HIGH
  //                         );

  // cout << endl << "FIT INVARIANT MASS DS_PLUS -> PHI PI" << endl << endl;
  // fit_unbinned_InvariantMass(vec_mass_DSplus_PhiPi_MC, 
  //                           {0.5, 1.966, 1.968, 0.005, 0.007}, 
  //                           {0.001, 0.001, 0.001, 0.001, 0.001}, 
  //                           {0., 0., 0., 0., 0.}, 
  //                           {1., 3., 3., 1., 1.}, 
  //                           {"fraction","mu1","mu2","sigma1","sigma2"}, 
  //                           mass_DSplus_PhiPi_MC, 
  //                           outfile, 
  //                           "D_Plus_Phi_Mu_Nu", 
  //                           MASS_DS_PLUS_PHI_PI_LOW, 
  //                           MASS_DS_PLUS_PHI_PI_HIGH
  //                         );

  // cout << endl << "FIT INVARIANT MASS DS_PLUS -> TAU NU" << endl << endl;
  // fit_unbinned_InvariantMass(vec_mass_DSplus_TauNu_MC, 
  //                           {0.5, 1.76, 1.78, 0.005, 0.007}, 
  //                           {0.001, 0.001, 0.001, 0.001, 0.001}, 
  //                           {0., 0., 0., 0., 0.}, 
  //                           {1., 3., 3., 1., 1.}, 
  //                           {"fraction","mu1","mu2","sigma1","sigma2"}, 
  //                           mass_DSplus_TauNu_MC, 
  //                           outfile, 
  //                           "D_Plus_Phi_Tau_Nu", 
  //                           MASS_DS_PLUS_TAU_NU_LOW, 
  //                           MASS_DS_PLUS_TAU_NU_HIGH
  //                         );

  // cout << endl << "FIT INVARIANT MASS DS_PLUS -> PHI MU NU" << endl << endl;
  // fit_unbinned_Argus(vec_mass_DSplus_PhiMuNu_MC, 
  //                   //{1.98, -4.0,  1.5}, 
  //                   //{0.001,  0.001,  0.001}, 
  //                   //{1.95, -10., 0.1}, 
  //                   //{2.05, -1., 2.5}, 
  //                   //{"m0", "c", "p"}, 
  //                   {-4.0,  1.5}, 
  //                   {0.001,  0.001}, 
  //                   {-10., 0.1}, 
  //                   {-1., 2.5}, 
  //                   {"c", "p"}, 
  //                   mass_DSplus_PhiMuNu_MC, 
  //                   outfile, 
  //                   "DS_Plus_Phi_Mu_Nu",
  //                   MASS_DS_PLUS_PHI_MU_NU_LOW,
  //                   MASS_DS_PLUS_PHI_MU_NU_HIGH
  //                 );

  //cout << endl << "FIT INVARIANT MASS COMBINATORIAL" << endl << endl;
  //fit_unbinned_Combinatorial(vec_mass_combinatorial, 
  //                  {1.0}, 
  //                  {0.001}, 
  //                  {0.}, 
  //                  {0.}, 
  //                  {"tau"}, 
  //                  mass_combinatorial_DATA, 
  //                  outfile, 
  //                  "Combinatorial",
  //                  MASS_COMBINATORIA_LOW,
  //                  MASS_COMBINATORIAL_HIGH
  //                );
/*
  cout << endl << "FIT INVARIANT MASS TOTAL SPECTRUM" << endl << endl;
  fit_unbinned_TotalSpectrum(vec_tot_mass_DATA, 
  // fit_unbinned_TotalSpectrum(vec_tot_mass_MC, 
                            {0.3, 0.15, 0.25, 0.01, 2.}, 
                            {0.01, 0.01, 0.01, 0.0001, 0.01}, 
                            {0., 0., 0., 0., 0.}, 
                            {1., 1., 1., 1., 10.}, 
                            {"f_D_PhiPi","f_DS_PhiMuNu","f_DS_PhiPi","f_DS_TauNu","tau"}, 
                            tot_mass_DATA, 
                            outfile, 
                            "Total_Invariant_Mass_Spectrum", 
                            TOT_MASS_LOW, 
                            TOT_MASS_HIGH
                          );

  // fit_unbinned_TotalSpectrum(vec_tot_mass_MC, 
  //                           {0.2, 0.1, 0.15, 0.3}, 
  //                           {0.01, 0.01, 0.01, 0.01}, 
  //                           {0., 0., 0., 0.}, 
  //                           {1., 1., 1., 1.}, 
  //                           {"f_D_PhiPi","f_DS_PhiMuNu","f_DS_PhiPi","f_DS_TauNu"}, 
  //                           tot_mass_DATA, 
  //                           outfile, 
  //                           "Total_Invariant_Mass_Spectrum", 
  //                           TOT_MASS_LOW, 
  //                           TOT_MASS_HIGH
  //                         );

  mass_Dplus_PhiPi_MC->Write();
  mass_DSplus_PhiMuNu_MC->Write();
  mass_DSplus_PhiPi_MC->Write();
  mass_DSplus_TauNu_MC->Write();
  tot_mass_MC->Write();
  outfile->Close();
}
*/

