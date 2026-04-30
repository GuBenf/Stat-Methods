#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "TFile.h"
#include <TCanvas.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include "TH1.h"
#include "TF1.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TStyle.h"
#include <TH2.h>
#include <TH2F.h>
#include <TH3.h>
#include <TH3F.h>
#include "TF2.h"
#include "TGraph2D.h"
#include "TH1F.h"
#include "TMath.h"
#include <string>
#include <iomanip>

#include "labstyle.C"

using namespace std;

#define MAX 10000 //<---Numero di eventi
Double_t xvar[MAX]; 

TRandom3 randGen(0);
Double_t mA= -1.;
Double_t mB= 1.;
Double_t sA= 0.5;
Double_t sB= 0.5;

Double_t pdf(Double_t x, Double_t *par)
{

 Double_t _fA;
 Double_t _mA;
 Double_t _mB;
 Double_t _sA;
 Double_t _sB;
 
 _fA  = par[0];				
 _mA  = par[1];
 _mB  = par[2];
 _sA  = par[3];
 _sB  = par[4];

 Double_t value = _fA*(TMath::Gaus(x,_mA,_sA,1)) + (1-_fA)*(TMath::Gaus(x,_mB,_sB,1));
 return value;
}

Double_t pdf_proj(Double_t *x, Double_t *par)
{
  Double_t bin_width;
  bin_width =0.08;
  return (MAX*bin_width)*pdf(x[0],par);

}


void fcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t )
{

  Double_t Like=0;
  Int_t i=0;
  for(i=0;i<MAX;++i){
   Like = Like + TMath::Log(pdf(xvar[i],par));
  };
   f= - 2. * Like;
}



void fitterLike() {
  
  //Setting style
   labstyle();
   
   // Create the canvas
   TCanvas *c_histo = new TCanvas("c_histo","canvas histo",600,600);
   c_histo->cd();
   
   //Create the histogram
   TH1F *histo = new TH1F("histo","",100,-4.,4.);
   
   // Event generation (Double Gaussian)
   Int_t i=0;
   for(i=0;i<7000;++i){
     xvar[i]= randGen.Gaus(mA,sA);
     histo->Fill(xvar[i]);
  };
   
   for(i=7000;i<10000;++i){
     xvar[i]= randGen.Gaus(mB,sB);
     histo->Fill(xvar[i]);
   };
   
   // Draw the histogram
   histo->GetXaxis()->SetTitle("lenght [cm] "); // Let us assume a lenght on x axis
   histo->GetYaxis()->SetTitle("Events per 0.08 cm "); 
   histo->Draw();
   
   TMinuit *my_gMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
   my_gMinuit->SetFCN(fcn);      // set the FCN
   
   Double_t arglist[2];
   Int_t ierflg = 0;  // Error return code: 0 if the command was correctly executed, >0 otherwise. 
   
   // arglist[0] = 1;                           
   // my_gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
   
   // Set starting values and step sizes for parameters
   Double_t vstart[5] = {0.3,  -1.  , 1.  , 0.5 , 0.5 };
   //Double_t vstart[5] = {0.53,  -0.34  , -0.2  , 0.014 , 0.017 };
   Double_t step[5]   = {0.0001, 0.0001, 0.0001, 0.0001, 0.0001};   //step 0 li rende costanti
   my_gMinuit->mnparm(0, "fA", vstart[0], step[0], 0., 1., ierflg);
   my_gMinuit->mnparm(1, "mA", vstart[1], step[1], 0., 0., ierflg);
   my_gMinuit->mnparm(2, "mB", vstart[2], step[2], 0., 0., ierflg);
   my_gMinuit->mnparm(3, "sA", vstart[3], step[3], 0., 0., ierflg);
   my_gMinuit->mnparm(4, "sB", vstart[4], step[4], 0., 0., ierflg);
   
   // Now ready for minimization step
   // Start minimization with arglist[0] max iterations and arglist[1] tolerance:
   // the default tolerance is 0.1, and the minimization will stop when
   // the estimated vertical distance to the minimum (EDM) is less than 
   // 0.001*[tolerance]*ERRDEF
   
   arglist[0] = 0.;//500;
   arglist[1] = 0.1;
   my_gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
      
   // Print results
   Double_t amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
   my_gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);     
   my_gMinuit->mnprin(3,amin);
   
   Double_t val0, val1, val2, val3, val4, error0,error1,error2, error3, error4, bnd1, bnd2;
   Int_t ivar=0;
   TString chnam;
   
   my_gMinuit->mnpout(0, chnam, val0, error0, bnd1, bnd2, ivar);
   my_gMinuit->mnpout(1, chnam, val1, error1, bnd1, bnd2, ivar);
   my_gMinuit->mnpout(2, chnam, val2, error2, bnd1, bnd2, ivar);
   my_gMinuit->mnpout(3, chnam, val3, error3, bnd1, bnd2, ivar);
   my_gMinuit->mnpout(4, chnam, val4, error4, bnd1, bnd2, ivar);
   
   Double_t matrix[5][5];
   Int_t n=5;
   my_gMinuit->mnemat (&matrix[0][0],5); // Accessing Covariance Matrix

   // Write on .txt file covariance Matrix 
   FILE *file_cov = fopen ("_txt/covariance_matrix.txt", "w");
   Int_t icol=0;
   Int_t irow=0;
   if (n==1)
     {
       fprintf (file_cov, "%e \n", matrix[0][0]); 
     }
   else
     {
       for (irow=0; irow<n; irow++)
	 {
	   for (icol=0; icol<n; icol++)
	     { 
               fprintf (file_cov, "%e ", matrix[irow][icol]); 
	     };
	   fprintf (file_cov, "\n");
	 };
       fprintf (file_cov, "\n");   
     }
   fclose (file_cov);
   
   

   //Projection
   TF1 *f_proj = new TF1("f_proj",pdf_proj,-4.,4.,5);
   
   // Sets values 
   f_proj->SetParameter(0,val0);//fA
   f_proj->SetParameter(1,val1);//mA
   f_proj->SetParameter(2,val2);//mB
   f_proj->SetParameter(3,val3);//sA
   f_proj->SetParameter(4,val4);//sB
   
   f_proj->Draw("SAME");

   c_histo->Print("./_fig/histo.pdf");
   c_histo->Print("./_fig/histo.eps");
         
}


