bool USE_MC = true;

#define analysis_cxx
#include "analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

Double_t MASS_D_PLUS_PHI_PI_LOW = 1.84;
Double_t MASS_D_PLUS_PHI_PI_HIGH = 1.9;
int MASS_D_PLUS_PHI_PI_BINS = 120;

Double_t MASS_DS_PLUS_PHI_MU_NU_LOW = 1.6;
Double_t MASS_DS_PLUS_PHI_MU_NU_HIGH = 2.;
int MASS_DS_PLUS_PHI_MU_NU_BINS = 800;

Double_t MASS_DS_PLUS_PHI_PI_LOW = 1.94;
Double_t MASS_DS_PLUS_PHI_PI_HIGH = 2.;
int MASS_DS_PLUS_PHI_PI_BINS = 120;

Double_t MASS_DS_PLUS_TAU_NU_LOW = 1.75;
Double_t MASS_DS_PLUS_TAU_NU_HIGH = 1.81;
int MASS_DS_PLUS_TAU_NU_BINS = 120;

Double_t TOT_MASS_LOW = 1.6;
Double_t TOT_MASS_HIGH = 2;
int TOT_MASS_BINS = 800;

int ID_D_PLUS_PHI_PI = 34;
int ID_DS_PLUS_PHI_MU_NU = 42;
int ID_DS_PLUS_PHI_PI = 41;
int ID_DS_PLUS_TAU_NU = 44;

std::vector<double> vec_mass_Dplus_PhiPi_MC;
std::vector<double> vec_mass_DSplus_PhiMuNu_MC;
std::vector<double> vec_mass_DSplus_PhiPi_MC;
std::vector<double> vec_mass_DSplus_TauNu_MC;
std::vector<double> tot_mass_MC;

std::vector<double> xvar;

Double_t MASS_FIT_D_PLUS_PHI_PI;
Double_t SIGMA_FIT_D_PLUS_PHI_PI;

Doublt_t M0_FIT_DS_PLUS_PHI_MU_NU;
Double_t C_FIT_DS_PLUS_PHI_MU_NU;
Double_t P_FIT_DS_PLUS_PHI_MU_NU;

Double_t MASS_FIT_DS_PLUS_PHI_PI;
Double_t SIGMA_FIT_DS_PLUS_PHI_PI;

Double_t MASS_FIT_MASS_DS_PLUS_TAU_NU;
Double_t SIGMA_FIT_MASS_DS_PLUS_TAU_NU;

Double_t TAU_FIT_COMBINATORIAL;

Double_t invariant_mass_pdf(Double_t x, Double_t *par)
{
  Double_t mu = par[0];
  Double_t sigma = par[1];

  return TMath::Gauss(x,mu,sigma,1);
  
}

Double_t pdf_proj_invariant_mass(Double_t *x, Double_t *par, Int_t max, Double_t bin_width)
{
  return (max*bin_width)*invariant_mass_pdf(x[0],par);
}

void fcn_invariant_mass(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t )
{
  Double_t Like=0;
  Int_t i=0;
  for(size_t i=0; i<x_var.size(); ++i)
    {
      Double_t p = invariant_mass_pdf(x_var[i],par);
      Like += TMath::Log(p);
  };
   f= - 2. * Like;
}

double pdf_Argus ( double *x , double * par )
{
  const double m = x [0];
  const double m0 = par [0];
  const double c = par [1];
  const double p = par [2];

  if ( m >= m0 ) {
    return 0.0;
  }

  const double xL = 1.0 - ( L_taumass_range / m0 ) *( L_taumass_range / m0 ) ;
  const double xH = 1.0 - ( H_taumass_range / m0 ) *( H_taumass_range / m0 ) ;

  const double gammaA = ROOT :: Math :: tgamma (1.0 + p ) ;
  const double dL = gammaA * ROOT :: Math :: inc_gamma_c (1.0 + p , -c * xL ) ;
  const double dH = gammaA * ROOT :: Math :: inc_gamma_c (1.0 + p , -c * xH ) ;
  const double norm = ( m0 * m0 ) /(2.0* c * std :: pow ( -c , p ) ) * ( dL - dH ) ;

  const double u = 1.0 - ( m / m0 ) *( m / m0 ) ;
  return m * std :: pow (u , p ) * std :: exp ( c * u ) / norm ;

}

Double_t combinatorial_pdf(Double_t x, Double_t *par)
{
  Double_t tau = par[0];

  return 1./tau * TMath::Exp(-1 * x * (1./tau)); 
}


Double_t total_mass_spectrum_pdf(Double_t x, Double_t *par)
{
  Double_t f_D_PhiPi = par[0]; 
  Double_t f_DS_PhiMuNu = par[1];
  Double_t f_DS_PhiPi = par[2];
  Double_t f_DS_TauNu = par[3];

  Double_t par_D_PhiPi[2] = {MASS_FIT_D_PLUS_PHI_PI, SIGMA_FIT_D_PLUS_PHI_PI};
  Double_t par_DS_PhiMuNu[3] = {M0_FIT_DS_PLUS_PHI_MU_NU, C_FIT_DS_PLUS_PHI_MU_NU, P_FIT_DS_PLUS_PHI_MU_NU};
  Double_t par_DS_PhiPi[2] = {MASS_FIT_DS_PLUS_PHI_PI, SIGMA_FIT_DS_PLUS_PHI_PI};
  Double_t par_DS_TauNu[2] = {MASS_FIT_MASS_DS_PLUS_TAU_NU, SIGMA_FIT_MASS_DS_PLUS_TAU_NU};

  Double_t val = f_D_PhiPi * invariant_mass_pdf(x,&par_D_PhiPi) + 
                 f_DS_PhiMuNu * pdf_Argus(x,&par_DS_PhiMuNu) + 
                 f_DS_PhiPi * invariant_mass_pdf(x,&par_DS_PhiPi) +
                 f_DS_TauNu * invariant_mass_pdf(x,&par_DS_TauNu) + 
                 (1 - f_D_PhiPi - f_DS_PhiMuNu - f_DS_PhiPi - f_DS_TauNu) * combinatorial_pdf(x,&TAU_FIT_COMBINATORIAL);
  return val;

}

Double_t pdf_proj_total_mass_spectrum(Double_t *x, Double_t *par, Int_t max, Double_t bin_width)
{
  return (max*bin_width)*total_mass_spectrum_pdf(x[0],par);
}

void fcn_total_mass_spectrum(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t )
{
  Double_t Like=0;
  Int_t i=0;
  for(size_t i=0; i<x_var.size(); ++i)
    {
      Double_t p = total_mass_spectrum_pdf(x_var[i],par);
      Like += TMath::Log(p);
  };
   f= - 2. * Like;
}

void analysis::Loop()
{

  TFile * outfile = new TFile("outfile.root","RECREATE");

  TH1D * mass_Dplus_PhiPi_MC = new TH1D("mass_Dplus_PhiPi_MC","",MASS_D_PLUS_PHI_PI_BINS,MASS_D_PLUS_PHI_PI_LOW,MASS_D_PLUS_PHI_PI_HIGH);
  TH1D * mass_DSplus_PhiMuNu_MC = new TH1D("mass_DSplus_PhiMuNu_MC","",MASS_DS_PLUS_PHI_MU_NU_BINS,MASS_DS_PLUS_PHI_MU_NU_LOW,MASS_DS_PLUS_PHI_MU_NU_HIGH);
  TH1D * mass_DSplus_PhiPi_MC = new TH1D("mass_DSplus_PhiPi_MC","",MASS_DS_PLUS_PHI_PI_BINS,MASS_DS_PLUS_PHI_PI_LOW,MASS_DS_PLUS_PHI_PI_HIGH);
  TH1D * mass_DSplus_TauNu_MC = new TH1D("mass_DSplus_TauNu_MC","",MASS_DS_PLUS_TAU_NU_BINS,MASS_DS_PLUS_TAU_NU_LOW,MASS_DS_PLUS_TAU_NU_HIGH);
  TH1D * tot_mass_MC = new TH1D("tot_mass_MC","",TOT_MASS_BINS,TOT_MASS_LOW,TOT_MASS_HIGH);
  
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

    if(id == ID_D_PLUS_PHI_PI)
    {
      mass_Dplus_PhiPi_MC->Fill(D_M);
      vec_mass_Dplus_PhiPi_MC.push_back(D_M);
    }
    if(id == ID_DS_PLUS_PHI_MU_NU)
    {
      mass_DSplus_PhiMuNu_MC->Fill(D_M);
      vec_mass_DSplus_PhiMuNu_MC.push_back(D_M);
    }
    if(id == ID_DS_PLUS_PHI_PI)
    {
      mass_DSplus_PhiPi_MC->Fill(D_M);
      vec_mass_DSplus_PhiPi_MC.push_back(D_M);
    }
    if(id == ID_DS_PLUS_TAU_NU)
    {
      mass_DSplus_TauNu_MC->Fill(D_M);
      mass_DSplus_TauNu_MC.push_back(D_M);
    }
    if(id != 0)
    {
      tot_mass_MC->Fill(D_M);
      tot_mass_MC.push_back(D_M);
    }
     
    if(id == 0 && TMath::Abs(D_M-m_tau_mc)>3*m_sigma_tau_mc) //DON'T LOOK IN THE HISTO!!
    {
	 
    }; //Mass blinding
     
  }; //End Loop Events

  

  mass_Dplus_PhiPi_MC->Write();
  mass_DSplus_PhiMuNu_MC->Write();
  mass_DSplus_PhiPi_MC->Write();
  mass_DSplus_TauNu_MC->Write();
  tot_mass_MC->Write();
   
}


