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


#define analysis_cxx
#include "fit_unbinned_TotalSpectrum_pseudo.h"
#include "fit_unbinned_InvariantMass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;
std::vector<double> fsigTrues;
std::vector<double> mins_band;
std::vector<double> maxs_band;
std::vector<double> x_band;
std::vector<double> y_band;


void toy() 
{

  TFile * outfile = TFile::Open("outfile.root","RECREATE");
  TFile * dummy = TFile::Open("dummy.root","RECREATE");

  gRandom ->SetSeed (12345);
  const double mMin = 1.60;
  const double mMax = 2.10;
  const int nBins = 100;
  const double binWidth = (mMax - mMin) / nBins;
  const int nToys = 500;

  const double nTotTrue = 10000.0;
  double fsigStart = 0.001;
  int bandpoints = 4;

  // Known model used to generate the pseudo - experiments .
  // TF1 :: GetRandom uses only the shape of the function.
  TF1 generatorModel("generatorModel", total_mass_spectrum_pdf , mMin , mMax, 4);

  for (int i=0; i<bandpoints; i++){
    fsigTrues.push_back(fsigStart+i*0.001);
  }
  for (int j = 0; j<bandpoints; j++)
  {
    std::vector<double> fsig;
    std::vector<double> fdiffsig;

  double fSigTrue = fsigTrues[j];
  generatorModel.SetParameters(f_D_PhiPi, f_Ds_PhiMuNu, f_Ds_PhiPi, fSigTrue);
  TH1D hfSigFit("hfSigFit", "Fitted signal yield;#hat{f}_{sig};Pseudo -experiments", 200, 0.0, 5*fSigTrue);
  TH1D hdifffSigFit("hdifffSigFit", "Fitted signal yield - Simulated signal yield;#hat{f}_{sig}-f_{sig};Pseudo -experiments", 200, -0.02, 0.02);

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

    Double_t fSigfromFIT = fit_unbinned_TotalSpectrum(
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

    hfSigFit.Fill(fSigfromFIT);
    fsig.push_back(fSigfromFIT);


  }
    std::tuple<Double_t, Double_t> fit_pars =  fit_unbinned_InvariantMass_results(fsig, 
                            {1, fSigTrue, fSigTrue, 0.005, 0}, 
                            {0, 0.001, 0, 0.0001, 0}, 
                            {0., 0., 0., 0., 0.}, 
                            {1., 1., 1., 1., 1.}, 
                            {"fraction","mu1","mu2","sigma1","sigma2"}, 
                            &hfSigFit, 
                            outfile, 
                            "fsig", 
                            0, 
                            1
                          );
  Double_t mean  = std::get<0>(fit_pars);
  Double_t sigma = std::get<1>(fit_pars);
  Double_t min_band = mean-2*sigma;
  Double_t max_band = mean+2*sigma;
  mins_band.push_back(min_band);
  maxs_band.push_back(max_band);
  //push back twice because there are two points for each fSigTrue
  x_band.push_back(fSigTrue);
  x_band.push_back(fSigTrue);

  y_band.push_back(min_band);
  y_band.push_back(max_band);
  }

  TGraph *gr = new TGraph(y_band.size(), x_band.data(), y_band.data());

    TCanvas *c = new TCanvas("c", "Band points", 800, 600);

    gr->SetTitle("Band points;fsig;min -/- max");
    gr->SetMarkerStyle(20);
    gr->Draw("AP");



    outfile -> cd();

    gr->Write();
  // hfSigFit.Write();
  outfile->Close();
  dummy->Close();

}
