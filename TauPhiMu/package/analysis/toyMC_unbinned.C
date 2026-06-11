Double_t MU1_FIT_D_PLUS_PHI_PI = 1.86969e+00;
Double_t MU2_FIT_D_PLUS_PHI_PI = 1.87035e+00;
Double_t SIGMA1_FIT_D_PLUS_PHI_PI = 5.79484e-03;
Double_t SIGMA2_FIT_D_PLUS_PHI_PI = 8.15639e-03;
Double_t FRACTION_FIT_D_PLUS_PHI_PI = 8.83747e-01;

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

Double_t TAU_FIT_COMBINATORIAL = 1.03105;

Double_t TOT_MASS_LOW = 1.6;
Double_t TOT_MASS_HIGH = 2.1;

//Double_t f_D_PhiPi = 0.0214727;
//Double_t f_Ds_PhiNuMu = 0.652866;
//Double_t f_Ds_PhiPi = 0.04305;

Double_t f_D_PhiPi = 0.0215021;
Double_t f_Ds_PhiMuNu = 0.653769;
Double_t f_Ds_PhiPi = 0.0431105;

Double_t fSigTrue = 0.01;

Double_t EFF_D_PhiPi = 53699./50e6;
Double_t EFF_Ds_PhiMuNu = 95217./5e6;
Double_t EFF_Ds_PhiPi = 32106./50e6;
Double_t EFF_Sig = 86069./1e6;

Double_t SIGMA_EFF_D_PhiPi = std::sqrt((EFF_D_PhiPi*(1-EFF_D_PhiPi))/50e6);
Double_t SIGMA_EFF_Ds_PhiMuNu = std::sqrt((EFF_Ds_PhiMuNu*(1-EFF_Ds_PhiMuNu))/5e6);
Double_t SIGMA_EFF_Ds_PhiPi = std::sqrt((EFF_Ds_PhiPi*(1-EFF_Ds_PhiPi))/50e6);
Double_t SIGMA_EFF_Sig = std::sqrt((EFF_Sig*(1-EFF_Sig))/1e6);

Double_t BR_D_PhiPi = 2.69e-3;
Double_t BR_Ds_PhiMuNu = 2.24e-2;
Double_t BR_Ds_PhiPi = 2.25e-2;
Double_t BR_Ds_TauNu = 5.39e-2;
Double_t BR_Phi_KK = 49.9e-2;

Double_t SIGMA_BR_D_PhiPi = 0.08e-3;
Double_t SIGMA_BR_Ds_PhiMuNu = 0.11e-2;
Double_t SIGMA_BR_Ds_PhiPi = 0.05e-2;
Double_t SIGMA_BR_Ds_TauNu = 0.09e-2;
Double_t SIGMA_BR_Phi_KK = 0.5e-2;

Double_t CS_pp_D = 834;
Double_t CS_pp_Ds = 353;

Double_t SIGMA_CS_pp_D = 2;
Double_t SIGMA_CS_pp_Ds = 9;

Double_t Gauss_pdf(Double_t *x, Double_t *par)
{
  return TMath::Gaus(x[0],par[0],par[1],1);
}

Double_t sample_val(Double_t mu, Double_t sigma)
{
  Double_t par[2] = {mu,sigma};

  TF1 *f = new TF1("f",Gauss_pdf,mu - 5 * sigma, mu + 5 * sigma, 2);
  f->SetParameters(mu,sigma);
  Double_t extracted_val = f->GetRandom(mu - 5 * sigma,mu + 5 * sigma);

  delete f;

  return extracted_val;
}

#define analysis_cxx
#include "fit_unbinned_TotalSpectrum_pseudo.h"
#include "fit_unbinned_Gauss.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;
std::vector<double> fsig;
std::vector<double> fdiffsig;
std::vector<double> fpullsig;
std::vector<double> fsig_over_fphiMuNu;

std::vector<double> BR_Sig;
std::vector<double> R_tau;


void toy() 
{

  TFile * outfile = TFile::Open("outfile.root","RECREATE");
  TFile * dummy = TFile::Open("dummy.root","RECREATE");

  gRandom ->SetSeed (12345);
  const double mMin = 1.60;
  const double mMax = 2.10;
  const int nBins = 100;
  const double binWidth = (mMax - mMin) / nBins;
  const int nToys = 10000;

  const double nTotTrue = 1000.0;

  // Known model used to generate the pseudo - experiments .
  // TF1 :: GetRandom uses only the shape of the function.
  TF1 generatorModel("generatorModel", total_mass_spectrum_pdf , mMin , mMax, 4);
  generatorModel.SetParameters(f_D_PhiPi, f_Ds_PhiMuNu, f_Ds_PhiPi, fSigTrue);
  TH1D hfSigFit("hfSigFit", "Fitted signal yield;#hat{f}_{sig};Pseudo experiments", 50, fSigTrue - 5*fSigTrue, fSigTrue + 5*fSigTrue);
  TH1D hdifffSigFit("hdifffSigFit", "Fitted signal yield - Simulated signal yield ; #hat{f_{sig}}-f_{sig} ; Pseudo experiments", 50, -0.04, 0.04);
  TH1D hpullSigFit("hpullSigFit", "(Fitted signal yield - Simulated signal yield)/Sigma fitted signal yield ; (#hat{f_{sig}}-f_{sig})/#sigma_{#hat{f_{sig}}} ; Pseudo experiments", 100, -7., 7.);
  TH1D hfSig_over_fPhiMuNu("hfSig_over_fPhiMuNu","Fitted signal yield / Fitted #var_phi#rightarrow#mu#nu yield ; hat{f_{sig}}/#hat{f_{#var_phi#rightarrow#mu#nu}} ; Pseudo experiments", 100, -0.04, 0.065);
  TH1D hfBRSig("hfBRSig","BR(#tau^+#rightarrow#var_phi#mu^+) ; BR(#tau^+#rightarrow#var_phi#mu^+) ; Pseudo experiments", 100, -0.025, 0.035);
  TH1D hfR_tau("hfR_tau","R_{#tau} ; R_{#tau} ; Pseudo experiments", 100, 0., 0.);

  for (int iToy = 0; iToy < nToys; ++iToy) 
  {
    if (iToy % 100 == 0){
      cout << "PROCESSING " << iToy << " TOY" << endl;

    }

    TH1D *hToy = new TH1D(Form("hToy_%d",iToy), "", nBins, mMin, mMax);
    // const int nObs = gRandom ->Poisson(nTotTrue); // Extended toy generation .
    const int nObs = nTotTrue;

    std::vector<double> vec_generated;
    for (int i = 0; i < nObs; ++i) 
    {
      const double m = generatorModel.GetRandom(mMin , mMax);
      hToy -> Fill(m);
      vec_generated.push_back(m);
    }

    std::tuple<Double_t,Double_t,Double_t,Double_t> tuple_fSigfromFIT = fit_unbinned_TotalSpectrum(
                            vec_generated, 
                            {0.02, 0.65, 0.04, 0.01}, 
                            {0.01, 0.01, 0.01, 0.0001}, 
                            {0., 0., 0., 0.}, 
                            {0., 0., 0., 0.}, 
                            {"f_D_PhiPi","f_DS_PhiMuNu","f_DS_PhiPi","f_DS_TauNu"}, 
                            hToy, 
                            dummy, 
                            Form("Generated_Invariant_Mass_Spectrum_%d",iToy), 
                            mMin, 
                            mMax
                          );

    Double_t fSigfromFIT     = std::get<0>(tuple_fSigfromFIT);
    Double_t err_fSigfromFIT = std::get<1>(tuple_fSigfromFIT);
    Double_t fPhiMuNufromFIT = std::get<2>(tuple_fSigfromFIT);
    Double_t err_fPhiMuNufromFIT = std::get<3>(tuple_fSigfromFIT);

    hfSigFit.Fill(fSigfromFIT);
    fsig.push_back(fSigfromFIT);
    hdifffSigFit.Fill(fSigfromFIT-fSigTrue);
    fdiffsig.push_back(fSigfromFIT-fSigTrue);
    hpullSigFit.Fill((fSigfromFIT-fSigTrue)/err_fSigfromFIT);
    fpullsig.push_back((fSigfromFIT-fSigTrue)/err_fSigfromFIT);
    hfSig_over_fPhiMuNu.Fill(fSigfromFIT/fPhiMuNufromFIT);
    fsig_over_fphiMuNu.push_back(fSigfromFIT/fPhiMuNufromFIT);
    hfBRSig.Fill( (fSigfromFIT) / ( sample_val(EFF_Sig,SIGMA_EFF_Sig) * sample_val(CS_pp_Ds,SIGMA_CS_pp_Ds) * sample_val(BR_Ds_TauNu,SIGMA_BR_Ds_TauNu) ) );
    BR_Sig.push_back( (fSigfromFIT) / ( sample_val(EFF_Sig,SIGMA_EFF_Sig) * sample_val(CS_pp_Ds,SIGMA_CS_pp_Ds) * sample_val(BR_Ds_TauNu,SIGMA_BR_Ds_TauNu) ) );
    hfR_tau.Fill( (fSigfromFIT * sample_val(EFF_Ds_PhiMuNu,SIGMA_EFF_Ds_PhiMuNu)) / (fPhiMuNufromFIT * sample_val(EFF_Sig,SIGMA_EFF_Sig)) );
    R_tau.push_back( (fSigfromFIT * sample_val(EFF_Ds_PhiMuNu,SIGMA_EFF_Ds_PhiMuNu)) / (fPhiMuNufromFIT * sample_val(EFF_Sig,SIGMA_EFF_Sig)) );

    delete hToy;

  }
    outfile -> cd();

  std::tuple <Double_t,Double_t> results_fit_fsig = fit_unbinned_Gauss( fsig, 
                                                                  {fSigTrue, 0.005}, 
                                                                  {0.001, 0.0001}, 
                                                                  {0., 0.}, 
                                                                  {0., 0.}, 
                                                                  {"mu", "sigma"}, 
                                                                  &hfSigFit, 
                                                                  outfile, 
                                                                  "fsig", 
                                                                  0.01 - 1, 
                                                                  0.01 + 1
                                                                  );
  
  std::tuple <Double_t,Double_t> results_fit_residuals = fit_unbinned_Gauss(  fdiffsig, 
                                                                        {0, 0.005}, 
                                                                        {0.001, 0.0001}, 
                                                                        {0., 0.}, 
                                                                        {0., 0.}, 
                                                                        {"mu", "sigma"}, 
                                                                        &hdifffSigFit, 
                                                                        outfile, 
                                                                        "residuals", 
                                                                        -0.04, 
                                                                        0.04
                                                                      );

  std::tuple <Double_t,Double_t> results_fit_pulls = fit_unbinned_Gauss(  fpullsig, 
                                                                    {0., 1.}, 
                                                                    {0.001, 0.001}, 
                                                                    {0., 0.}, 
                                                                    {0., 0.}, 
                                                                    {"mu", "sigma"},
                                                                    &hpullSigFit, 
                                                                    outfile, 
                                                                    "pulls", 
                                                                    -7, 
                                                                    7
                                                                  );  

  std::tuple <Double_t,Double_t> results_fit_fSig_over_fPhiMuNu = fit_unbinned_Gauss(  fsig_over_fphiMuNu, 
                                                                    {0.01, 0.01}, 
                                                                    {0.0001, 0.0001}, 
                                                                    {0., 0.}, 
                                                                    {0., 0.}, 
                                                                    {"mu", "sigma"},
                                                                    &hfSig_over_fPhiMuNu, 
                                                                    outfile, 
                                                                    "fSig_over_fPhiMuNu", 
                                                                    -0.04, 
                                                                    0.065
                                                                  ); 
  
  std::tuple <Double_t,Double_t> results_fit_BRSig = fit_unbinned_Gauss(  BR_Sig, 
                                                                    {0.006, 0.005}, 
                                                                    {0.0001, 0.0001}, 
                                                                    {0., 0.}, 
                                                                    {0., 0.}, 
                                                                    {"mu", "sigma"},
                                                                    &hfBRSig, 
                                                                    outfile, 
                                                                    "BRSig", 
                                                                    -0.025, 
                                                                    0.035
                                                                  );  
                                                                  
  std::tuple <Double_t,Double_t> results_fit_R_tau = fit_unbinned_Gauss(  R_tau, 
                                                                    {0.006, 0.005}, 
                                                                    {0.0001, 0.0001}, 
                                                                    {0., 0.}, 
                                                                    {0., 0.}, 
                                                                    {"mu", "sigma"},
                                                                    &hfR_tau, 
                                                                    outfile, 
                                                                    "R_tau", 
                                                                    -0.025, 
                                                                    0.035
                                                                  );   

  outfile->Close();
  dummy->Close();

}
