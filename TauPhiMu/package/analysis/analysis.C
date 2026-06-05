bool USE_MC = true;

#define analysis_cxx
#include "analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

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


using namespace std;

void analysis::Loop()
{



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
    }
    if(id == ID_DS_PLUS_PHI_MU_NU)
    {
      mass_DSplus_PhiMuNu_MC->Fill(D_M);
    }
    if(id == ID_DS_PLUS_PHI_PI)
    {
      mass_DSplus_PhiPi_MC->Fill(D_M);
    }
    if(id == ID_DS_PLUS_TAU_NU)
    {
      mass_DSplus_TauNu_MC->Fill(D_M);
    }
    if(id != 0)
    {
      tot_mass_MC->Fill(D_M);
    }
     
    if(id == 0 && TMath::Abs(D_M-m_tau_mc)>3*m_sigma_tau_mc) //DON'T LOOK IN THE HISTO!!
    {
	 
    }; //Mass blinding
     
  }; //End Loop Events
   
}