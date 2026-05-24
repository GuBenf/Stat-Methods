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

const double XMIN = 0.5;
const double XMAX = 10.0;

Double_t MAX_ACCEPTANCE = 0.115673;

Double_t ACCEPTANCE_PAR_FIXED[7] = {0.58397, -1.77281, 1.75128, 1.41352, -0.825079,  0.969095, -1.88183};

using namespace std;

std::vector<double> xvar;
std::vector<double> tau_grid;
std::vector<double> norm_grid;

TH1D *h_time = new TH1D("h_time","",100,0,10);


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

   value = value * acceptance_func_mod(&x,ACCEPTANCE_PAR_FIXED);
   
   return value;
}

Double_t computeNorm(Double_t tau)
{
    Double_t sum = 0;
    Double_t dx = 0.001;  // fine enough

    for (double x = XMIN; x < XMAX; x += dx)
    {
        double p = pdf_exp_x_gauss(x, &tau);
        sum += p * dx;
    }
    return sum;
}

for (double tau = 0.1; tau < 2.0; tau += 0.01)
{
    tau_grid.push_back(tau);
    norm_grid.push_back(computeNorm(tau));
}

double getNorm(double tau)
{
    // linear interpolation (fast)
    int i = int((tau - tau_grid[0]) / 0.01);
    if (i < 0) i = 0;
    if (i >= tau_grid.size()-1) i = tau_grid.size()-2;

    double t1 = tau_grid[i], t2 = tau_grid[i+1];
    double n1 = norm_grid[i], n2 = norm_grid[i+1];
    return n1 + (tau - t1) * (n2 - n1) / (t2 - t1);
}

double pdf_norm(double x, double *par)
{
    double raw = pdf_exp_x_gauss(x, par);
    return raw / getNorm(par[0]);
}

void fcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t )
{

  Double_t Like=0;
  Int_t i=0;
  for(size_t i=0; i<xvar.size(); ++i)
    {
      double p = pdf_exp_x_gauss(xvar[i], par);
      p /= getNorm(par[0]);
      Like += log(p);
  };
   f= - 2. * Like;
}

Double_t pdf_proj(Double_t *x, Double_t *par, Double_t N, Double_t bin_width)
{
    return N * pdf_norm(x[0], par) * bin_width;
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
            double t = M0_time / 410.3e-15;

            h_time->Fill(t);

            if (t >= XMIN && t <= XMAX)
            {
            xvar.push_back(t);
            }

      }; //end MC

      };// #### end loop over jentry

      // ***** FIT ***** 

      const int nparam = 1;
      int n_fit = xvar.size();

      TMinuit *my_gMinuit = new TMinuit(nparam);  //initialize TMinuit with a maximum of 5 params
      //gMinuit->SetPrintLevel(-1);
      my_gMinuit->SetFCN(fcn);      // set the FCN
        
      Double_t arglist[2];
      Int_t ierflg = 0;  // Error return code: 0 if the command was correctly executed, >0 otherwise. 
        
      // Set starting values and step sizes for parameters
      Double_t vstart[nparam] = {1.};
      Double_t step[nparam]   = {0.01};   //step 0 li rende costanti
     
      my_gMinuit->mnparm(0, "tau", vstart[0], step[0], 0.1, 2., ierflg);

      arglist[0] = 500.;//500;
      arglist[1] = 0.1;

      my_gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
      
      // Print results
      Double_t amin,edm,errdef;
      Int_t nvpar,nparx,icstat;
      my_gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);     
      my_gMinuit->mnprin(3,amin);
   
      Double_t tau_fit, etau_fit, bnd1, bnd2;
      Int_t ivar=0;
      TString chnam;
        
      my_gMinuit->mnpout(0, chnam, tau_fit, etau_fit, bnd1, bnd2, ivar);
      std::cout << "\n========== FIT RESULT ==========\n";
      std::cout << "Parameter : " << chnam << std::endl;
      std::cout << "tau_fit   = " << tau_fit*410.3 << " fs" << std::endl;
      std::cout << "error     = " << etau_fit*410.3 << " fs" << std::endl;
      std::cout << "================================\n";

      std::vector<double> x_points;
      std::vector<double> y_points;

      Double_t pars[nparam] = {tau_fit};

      for(int xi=0; xi<10000; xi++)
      {
            double xmin = 0;
            double xmax = 10;
        
            double x = xmin + xi * (xmax-xmin)/10000;

            x_points.push_back(x);
            y_points.push_back(pdf_proj(&x,pars,n_fit,h_time->GetBinWidth(1)));
      }

      TGraph * fit_function_plot = new TGraph(10000,x_points.data(),y_points.data());

      TCanvas *c_fit = new TCanvas("c_fit", "fit", 800, 600);
      c_fit->cd();

      const float split = 0.3;

      // -------------------- TOP PAD --------------------
      TPad *pad_fit = new TPad("pad_fit", "", 0., split, 1., 1.);
      pad_fit->SetBottomMargin(0.02);
      pad_fit->Draw();
      pad_fit->cd();

      gStyle->SetOptStat(0);
      gPad->SetTicks(1,1);

      h_time->GetXaxis()->SetRangeUser(0., 10.0);
      h_time->Draw("E");

      // overlay fit
      fit_function_plot->Draw("L SAME");

      // IMPORTANT: go back to canvas
      c_fit->cd();

      // -------------------- BOTTOM PAD --------------------
      TPad *pad_res = new TPad("pad_res", "", 0., 0., 1., split);
      pad_res->SetTopMargin(0.02);
      pad_res->SetBottomMargin(0.25);
      pad_res->Draw();
      pad_res->cd();

      // residuals
      int n = h_time->GetNbinsX();

      TGraphErrors *residuals = new TGraphErrors(n);
      residuals -> SetMarkerStyle(7);
      residuals -> GetXaxis()-> SetTitle("t / 410.3 fs");
      residuals -> GetYaxis()-> SetTitle("residuals");

      for (int i = 0; i < n; i++) 
      {
            Double_t x, y;
            x = h_time->GetBinCenter(i+1);
            y = h_time->GetBinContent(i+1);
            Double_t y_fit  = pdf_proj(&x,pars,n_fit,h_time->GetBinWidth(1));
            Double_t ey     = h_time->GetBinError(i+1);
            Double_t res    = y - y_fit;

            residuals->SetPoint(i, x, res);
            residuals->SetPointError(i, 0., ey);
      }
      residuals->Draw("P same");

      // baseline
      residuals->GetXaxis()->SetLimits(0.0, 10.0);

      TLine *line_zero = new TLine(0.0, 0., 10., 0.);
      line_zero->SetLineStyle(2);
      fit_function_plot->SetLineColor(kRed);

      line_zero->Draw("same");

      // final update
      c_fit->Update();
      c_fit->Draw();

      // // ***** DRAW *****
      // TCanvas *c_fit = new TCanvas("c_fit", "");
      // c_fit->cd();
      // const float split = 0.3;

      // // --- FIT ---
      // TPad *pad_fit = new TPad("pad_fit", "", 0., split, 1., 1.);

      // pad_fit->cd();

      // gStyle->SetOptStat(0);
      // gPad->SetTicks(1,1);

      // h_time->SetMarkerStyle(20);
      // h_time->SetMarkerSize(0.8);
      // h_time->SetLineColor(kBlack);
      // h_time->SetTitle("");
      // h_time->GetXaxis()->SetTitle("t / 410.3 fs");
      // h_time->GetYaxis()->SetTitle("Events");
      // h_time->GetXaxis()->SetRangeUser(0.5, 10.0);

      // h_time->Draw("E");

      // fit_function_plot->SetLineColor(kRed);
      // fit_function_plot->SetLineWidth(2);
      // fit_function_plot->Draw("L SAME");
      // c_fit->cd();


      // // --- RESIDUALS ---
      // c_fit->cd();

      // TPad *pad_res = new TPad("pad_res", "", 0., 0., 1., split);
      // pad_res->SetTopMargin(0.);
      // pad_res->Draw();
      // pad_res->cd();

      // int n = h_time->GetNbinsX();

      // TGraphErrors *residuals = new TGraphErrors(n);
      // residuals -> SetMarkerStyle(7);
      // residuals -> GetXaxis()-> SetTitle("t / 410.3 fs");
      // residuals -> GetYaxis()-> SetTitle("residuals");

      // for (int i = 0; i < n; i++) 
      // {
      //       Double_t x, y;
      //       x = h_time->GetBinCenter(i+1);
      //       y = h_time->GetBinContent(i+1);
      //       Double_t y_fit  = pdf_proj(&x,pars,h_time->GetEntries(),h_time->GetBinWidth(1));
      //       Double_t ey     = h_time->GetBinError(i+1);
      //       Double_t res    = y - y_fit;

      //       residuals->SetPoint(i, x, res);
      //       residuals->SetPointError(i, 0., ey);
      // }
      // residuals->Draw("P same");

      // TLine *line_zero = new TLine(0., 0., 10., 0.);
      // line_zero->SetLineStyle(2);
      // line_zero->SetLineColor(kGray);
      // line_zero->Draw("same");

      // c_fit->Update();
      // c_fit->Draw();
      // c_fit->Update();      
      TFile * f = new TFile("outfile_fit_MC.root","RECREATE");
      c_fit -> Write();
      f->Close();

}
