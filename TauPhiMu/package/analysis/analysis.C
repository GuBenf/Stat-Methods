bool USE_MC = false;
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
Double_t TOT_MASS_HIGH = 2.1;

int TOT_MASS_BINS = 200;

Double_t TOT_MASS_BIN_WIDTH = (TOT_MASS_HIGH - TOT_MASS_LOW)/TOT_MASS_BINS;

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

// --> UNBINNED Double_t TAU_FIT_COMBINATORIAL = 1.03132e+00;
// --> UNBINNED 
// --> UNBINNED Double_t f_D_PhiPi = 0.02;
// --> UNBINNED Double_t f_Ds_PhiNuMu = 0.65;
// --> UNBINNED Double_t f_Ds_PhiPi = 0.04;
// --> UNBINNED 
// --> UNBINNED Double_t f_D_PhiPi_FIT = 0.0215021;
// --> UNBINNED Double_t f_Ds_PhiMuNu_FIT = 0.653769;
// --> UNBINNED Double_t f_Ds_PhiPi_FIT = 0.0431105;
// --> UNBINNED 
// --> UNBINNED Double_t sigma_f_D_PhiPi_FIT = 0.000960993;
// --> UNBINNED Double_t sigma_f_Ds_PhiMuNu_FIT = 0.0134033;
// --> UNBINNED Double_t sigma_f_Ds_PhiPi_FIT = 0.000912701;

// --> BINNED 100 BIN Double_t TAU_FIT_COMBINATORIAL = 0.966229;
// --> BINNED 100 BIN 
// --> BINNED 100 BIN Double_t f_D_PhiPi = 0.02;
// --> BINNED 100 BIN Double_t f_Ds_PhiNuMu = 0.65;
// --> BINNED 100 BIN Double_t f_Ds_PhiPi = 0.04;
// --> BINNED 100 BIN 
// --> BINNED 100 BIN Double_t f_D_PhiPi_FIT = 0.0212275;
// --> BINNED 100 BIN Double_t f_Ds_PhiMuNu_FIT = 0.649015;
// --> BINNED 100 BIN Double_t f_Ds_PhiPi_FIT = 0.0426402;
// --> BINNED 100 BIN 
// --> BINNED 100 BIN Double_t sigma_f_D_PhiPi_FIT = 0.000968707;
// --> BINNED 100 BIN Double_t sigma_f_Ds_PhiMuNu_FIT = 0.0150055;
// --> BINNED 100 BIN Double_t sigma_f_Ds_PhiPi_FIT = 0.000907028;

Double_t TAU_FIT_COMBINATORIAL = 1.02535;

Double_t f_D_PhiPi = 0.02;
Double_t f_Ds_PhiNuMu = 0.65;
Double_t f_Ds_PhiPi = 0.04;

Double_t f_D_PhiPi_FIT = 0.0214183;
Double_t f_Ds_PhiMuNu_FIT = 0.65332;
Double_t f_Ds_PhiPi_FIT = 0.0429899;

Double_t sigma_f_D_PhiPi_FIT = 0.000961335;
Double_t sigma_f_Ds_PhiMuNu_FIT = 0.0152777;
Double_t sigma_f_Ds_PhiPi_FIT = 0.00091669;

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

inline Double_t sample_val(Double_t mu, Double_t sigma)
{
    return gRandom->Gaus(mu,sigma);
}


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

void analysis::Loop()
{

  TFile * outfile = new TFile("outfile_n1.root","RECREATE");

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
     
    if(id == 0 && D_M >= TOT_MASS_LOW && D_M <= TOT_MASS_HIGH)
    {
      vec_tot_mass_DATA.push_back(D_M);
      tot_mass_DATA->Fill(D_M);
    }

    //if(id == 0 && TMath::Abs(D_M-m_tau_mc)>3*m_sigma_tau_mc && D_M >= TOT_MASS_LOW && D_M <= TOT_MASS_HIGH) //DON'T LOOK IN THE HISTO!!
    //{
    //  tot_mass_DATA->Fill(D_M);
    //}; //Mass blinding

    //if(id == 0 && TMath::Abs(D_M-m_tau_mc)>3*m_sigma_tau_mc && D_M >= MASS_COMBINATORIA_LOW && D_M <= MASS_COMBINATORIAL_HIGH) 
    //{
    //  vec_mass_combinatorial.push_back(D_M);
    //  mass_combinatorial_DATA->Fill(D_M);
    //}; //Mass blinding
     
  }; //End Loop Events


  if(USE_MC)
  {
    cout << endl << "FIT INVARIANT MASS D_PLUS -> PHI PI" << endl << endl;
    fit_unbinned_InvariantMass(vec_mass_Dplus_PhiPi_MC, 
                              {0.5, 1.86, 1.88, 0.005, 0.007}, 
                              {0.001, 0.001, 0.001, 0.001, 0.001}, 
                              {0., 0., 0., 0., 0.}, 
                              {1., 3., 3., 1., 1.}, 
                              {"fraction","mu1","mu2","sigma1","sigma2"}, 
                              mass_Dplus_PhiPi_MC, 
                              outfile, 
                              "D_Plus_Phi_Pi", 
                              MASS_D_PLUS_PHI_PI_LOW, 
                              MASS_D_PLUS_PHI_PI_HIGH
                            );
  }

  if(USE_MC)
  {
    cout << endl << "FIT INVARIANT MASS DS_PLUS -> PHI PI" << endl << endl;
    fit_unbinned_InvariantMass(vec_mass_DSplus_PhiPi_MC, 
                              {0.5, 1.966, 1.968, 0.005, 0.007}, 
                              {0.001, 0.001, 0.001, 0.001, 0.001}, 
                              {0., 0., 0., 0., 0.}, 
                              {1., 3., 3., 1., 1.}, 
                              {"fraction","mu1","mu2","sigma1","sigma2"}, 
                              mass_DSplus_PhiPi_MC, 
                              outfile, 
                              "D_Plus_Phi_Mu_Nu", 
                              MASS_DS_PLUS_PHI_PI_LOW, 
                              MASS_DS_PLUS_PHI_PI_HIGH
                            );
  }

  if(USE_MC)
  {
    cout << endl << "FIT INVARIANT MASS DS_PLUS -> TAU NU" << endl << endl;
    fit_unbinned_InvariantMass(vec_mass_DSplus_TauNu_MC, 
                              {0.5, 1.76, 1.78, 0.005, 0.007}, 
                              {0.001, 0.001, 0.001, 0.001, 0.001}, 
                              {0., 0., 0., 0., 0.}, 
                              {1., 3., 3., 1., 1.}, 
                              {"fraction","mu1","mu2","sigma1","sigma2"}, 
                              mass_DSplus_TauNu_MC, 
                              outfile, 
                              "D_Plus_Phi_Tau_Nu", 
                              MASS_DS_PLUS_TAU_NU_LOW, 
                              MASS_DS_PLUS_TAU_NU_HIGH
                            );
  }

  if(USE_MC)
  {
   cout << endl << "FIT INVARIANT MASS DS_PLUS -> PHI MU NU" << endl << endl;
   fit_unbinned_Argus(vec_mass_DSplus_PhiMuNu_MC, 
                     //{1.98, -4.0,  1.5}, 
                     //{0.001,  0.001,  0.001}, 
                     //{1.95, -10., 0.1}, 
                     //{2.05, -1., 2.5}, 
                     //{"m0", "c", "p"}, 
                     {-4.0,  1.5}, 
                     {0.001,  0.001}, 
                     {-10., 0.1}, 
                     {-1., 2.5}, 
                     {"c", "p"}, 
                     mass_DSplus_PhiMuNu_MC, 
                     outfile, 
                     "DS_Plus_Phi_Mu_Nu",
                     MASS_DS_PLUS_PHI_MU_NU_LOW,
                     MASS_DS_PLUS_PHI_MU_NU_HIGH
                   );
  }

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

  if(!USE_MC)
  {
    cout << endl << "FIT INVARIANT MASS TOTAL SPECTRUM" /*<< rnd_offset*/ << endl << endl;
    // fit_unbinned_TotalSpectrum(vec_tot_mass_DATA, 
    //                        //{0.15, 0.3, 0.25, 0.01 - rnd_offset, 2.}, 
    //                        {f_D_PhiPi, f_Ds_PhiNuMu, f_Ds_PhiPi, 0.01, 1.03},
    //                        {0.01, 0.01, 0.01, 0.0001, 0.01}, 
    //                        //{0., 0., 0., 0. - rnd_offset, 0.01}, 
    //                        //{1., 1., 1., 1. - rnd_offset, 10.}, 
    //                        {0., 0., 0., 0., 0.}, 
    //                        {0., 0., 0., 0., 0.},
    //                        {"f_D_PhiPi","f_DS_PhiMuNu","f_DS_PhiPi","f_DS_TauNu","tau"}, 
    //                        tot_mass_DATA, 
    //                        outfile, 
    //                        "Total_Invariant_Mass_Spectrum", 
    //                        TOT_MASS_LOW, 
    //                        TOT_MASS_HIGH
    //                      );

    TMatrixDSym cov;

    cov.Print();
    std::pair<std::vector<double>,std::vector<double>> FITres = fit_binned_TotalSpectrum(
                            {0.021, 0.654, 0.043, 0.001, 1., tot_mass_DATA->GetEntries()}, 
                            tot_mass_DATA, 
                            outfile, 
                            "DATA_Invariant_Mass_Spectrum", 
                            TOT_MASS_LOW, 
                            TOT_MASS_HIGH,
                            cov
                            );

    std::vector<double> pars_fit_res = FITres.first;
    std::vector<double> pars_err_fit_res = FITres.second;

    Double_t f_D_PhiPi_fromFIT = pars_fit_res[0];
    Double_t err_f_D_PhiPi_fromFIT = pars_err_fit_res[0];

    Double_t f_Ds_PhiMuNu_fromFIT = pars_fit_res[1];
    Double_t err_f_Ds_PhiMuNu_fromFIT = pars_err_fit_res[1];

    Double_t f_Ds_PhiPi_fromFIT = pars_fit_res[2];
    Double_t err_f_Ds_PhiPi_fromFIT = pars_err_fit_res[2];

    Double_t f_Sig_fromFIT = pars_fit_res[3];
    Double_t err_f_Sig_fromFIT = pars_err_fit_res[3];

    Double_t f_comb_fromFIT = pars_fit_res[4];
    Double_t err_f_comb_fromFIT = pars_err_fit_res[4];


    cout << "f_D_PhiPi_fromFIT : " << f_D_PhiPi_fromFIT << " pm " << err_f_D_PhiPi_fromFIT << endl;
    cout << "f_Ds_PhiMuNu_fromFIT : " << f_Ds_PhiMuNu_fromFIT << " pm " << err_f_Ds_PhiMuNu_fromFIT << endl;
    cout << "f_Ds_PhiPi_fromFIT : " << f_Ds_PhiPi_fromFIT << " pm " << err_f_Ds_PhiPi_fromFIT << endl;
    cout << "f_comb_fromFIT : " << f_comb_fromFIT << " pm " << err_f_comb_fromFIT << endl;

    cout << "chi2 / ndof = " <<  pars_fit_res[6] << " / " << pars_fit_res[7] << " = " <<  pars_fit_res[6] / pars_fit_res[7] << " | p value = " << pars_fit_res[8] << endl; 

    double r_phipi = (f_D_PhiPi_fromFIT / f_Ds_PhiPi_fromFIT) * (EFF_Ds_PhiPi / EFF_D_PhiPi);
    double r_phipi_theo = 834/353 * 2.69e-3/2.25e-2;

    //double BR = ..;

    //double rtau = ..;

    cout << "R_PhiPi " << r_phipi  << "  th:  " << r_phipi_theo << endl;


    double s  = f_Sig_fromFIT;
    double m  = f_Ds_PhiMuNu_fromFIT;

    double BR = ( s / ((1.0-s)*m) )
          * ( BR_Ds_PhiMuNu * EFF_Ds_PhiMuNu )
          / ( BR_Ds_TauNu   * EFF_Sig );
    double dBR_ds = BR * (1.0/s + 1.0/(1.0-s));
    double dBR_dm = -BR / m;

    double Vss = cov(3,3);
    double Vmm = cov(1,1);
    double Vsm = cov(3,1);

    double var_fit =
          dBR_ds*dBR_ds * Vss
        + dBR_dm*dBR_dm * Vmm
        + 2.0*dBR_ds*dBR_dm * Vsm;
    double rel2_ext =
          pow(SIGMA_BR_Ds_PhiMuNu  / BR_Ds_PhiMuNu,  2)
        + pow(SIGMA_EFF_Ds_PhiMuNu / EFF_Ds_PhiMuNu, 2)
        + pow(SIGMA_BR_Ds_TauNu    / BR_Ds_TauNu,    2)
        + pow(SIGMA_EFF_Sig        / EFF_Sig,        2);

    double var_ext = BR*BR*rel2_ext;
    double errBR = std::sqrt(var_fit + var_ext);

    cout << "BR = " << BR
     << " +/- " << errBR
     << endl;
    // double BR = ( f_Sig_fromFIT / ((1 - f_Sig_fromFIT) * f_Ds_PhiMuNu_fromFIT) ) * ( ( BR_Ds_PhiMuNu * EFF_Ds_PhiMuNu ) / ( BR_Ds_TauNu * EFF_Sig ) );

    double r_tau =
        (s * EFF_Ds_PhiMuNu) /
        (m * EFF_Sig);

    double dr_ds =  r_tau / s;
    double dr_dm = -r_tau / m;

    double var_fit_r =
      dr_ds*dr_ds * Vss
    + dr_dm*dr_dm * Vmm
    + 2.0*dr_ds*dr_dm * Vsm;

    double rel2_ext_r =
      pow(SIGMA_EFF_Ds_PhiMuNu / EFF_Ds_PhiMuNu, 2)
    + pow(SIGMA_EFF_Sig        / EFF_Sig,        2);

    double var_ext_r = r_tau*r_tau * rel2_ext_r;
    double err_r_tau = std::sqrt(var_fit_r + var_ext_r);

    cout << "Rtau = " << r_tau
     << " +/- " << err_r_tau
     << endl;

    Double_t LR;

    /* // --> BLINDING

    double BR = ( f_Sig_fromFIT / ((1 - f_Sig_fromFIT) * f_Ds_PhiMuNu_fromFIT) ) * ( ( sample_val(BR_Ds_PhiMuNu,SIGMA_BR_Ds_PhiMuNu) * sample_val(EFF_Ds_PhiMuNu,SIGMA_EFF_Ds_PhiMuNu) ) / ( sample_val(BR_Ds_TauNu,SIGMA_BR_Ds_TauNu) * sample_val(EFF_Sig,SIGMA_EFF_Sig) ) );
    Double_t logL = TMath::Log(TMath::Gaus(BR,BR,0.000890334,1));

    Double_t logL_null = TMath::Log(TMath::Gaus(BR,0.,0.000890334,1));

    
    if(BR >= 0)
    { 
      LR = -2*(logL_null - logL);
      cout << "LR for sensitivity is: " << LR << endl;
    }
    else
    {
      cout << "fSig is negative, there is not sensitivity" << endl;
    }
    */
  

  }

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

