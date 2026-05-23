#define lifetimeANA_cxx
#include "lifetimeANA.h"

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <TRandom3.h>  

Double_t MU_A_FIXED = -0.000513713; //--> NOT USED SINCE WE ARE RESTRICTED TO A BIAS=0 REGION (fairly accurate)
Double_t MU_B_FIXED = 0.00536732;   //--> NOT
Double_t S_A_FIXED = 0.0825203;
Double_t S_B_FIXED = 0.130739;
Double_t F_FIXED = 0.495402;

Double_t MAX_ACCEPTANCE = 0.115673;

Double_t ACCEPTANCE_PAR_FIXED[7] = {0.58397, -1.77281, 1.75128, 1.41352, -0.825079,  0.969095, -1.88183};

using namespace std;

Double_t acceptance_func_mod(Double_t *x, Double_t *par)
{
    Double_t xx = x[0];

    Double_t term1 = par[0] * TMath::Erf((xx - par[1]) / par[2]);
    Double_t term2 = par[3] * TMath::Erf((xx - par[4]) / par[5]);

    Double_t returned_func;
    if (x[0]<0.5) returned_func = 0;
    else returned_func = (term1 + term2 + par[6])/MAX_ACCEPTANCE; //--> TO NORMALIZE AT 1 THE MAXIMUM OF THE ACCEPTANCE

    return returned_func;
}

Double_t pdf_exp_x_gauss(Double_t x, Double_t *par)
{
   Double_t tau;
   Double_t fA;
   Double_t muA;
   Double_t sigmaA;
   Double_t fB;
   Double_t muB;
   Double_t sigmaB;

   tau = par[0];
   fA = F_FIXED;
   //muA = MU_A_FIXED;
   sigmaA = S_A_FIXED;
   fB = 1-F_FIXED;
   //muB = MU_B_FIXED;
   sigmaB = S_B_FIXED;

   Double_t value = fA*( 1./(2*tau)*TMath::Exp((sigmaA*sigmaA - 2*tau*x)/(2*tau*tau))*( TMath::Erf((tau*x - sigmaA*sigmaA)/(TMath::Sqrt(2)*tau*sigmaA)) +1 ) );
   
   value = value + fB*( 1./(2*tau)*TMath::Exp((sigmaB*sigmaB - 2*tau*x)/(2*tau*tau))*( TMath::Erf((tau*x - sigmaB*sigmaB)/(TMath::Sqrt(2)*tau*sigmaB)) +1 ) );

   value = value * acceptance_func_mod(x[0],ACCEPTANCE_PAR_FIXED);
   
   return value;
}

void fcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t )
{

  Double_t Like=0;
  Int_t i=0;
  for(size_t i=0; i<xvar.size(); ++i)
    {
      Double_t p = pdf_exp_x_gauss(xvar[i],par);
      Like += TMath::Log(p);
  };
   f= - 2. * Like;
}

Double_t pdf_proj(Double_t *x, Double_t *par, Int_t max, Double_t bin_width)
{
  return (max*bin_width)*pdf_exp_x_gauss(x[0],par);
}

void lifetimeANA::Loop()
{  
  
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // #### My code #########
      
      if (id==13) //MC
      {
	
      }; //end MC
      
   };// #### end loop over jentry


}
