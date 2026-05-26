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

for (double tau = 0.01; tau < 10.0; tau += 0.001)
{
    tau_grid.push_back(tau);
    norm_grid.push_back(computeNorm(tau));
}

double getNorm(double tau)
{
    // linear interpolation (fast)
    int i = int((tau - tau_grid[0]) / 0.001);
    if (i < 0) i = 0;
    if (i >= tau_grid.size()-1) i = tau_grid.size()-2;

    double t1 = tau_grid[i], t2 = tau_grid[i+1];
    double n1 = norm_grid[i], n2 = norm_grid[i+1];
    return n1 + (tau - t1) * (n2 - n1) / (t2 - t1);
}


Double_t pdf_background(Double_t x, Double_t *par)
{
    Double_t first_exp = pdf_exp_x_gauss(x,&par[0])/getNorm(par[0]);

    Double_t second_exp = pdf_exp_x_gauss(x,&par[1])/getNorm(par[1]);

    return par[2] * first_exp + (1 - par[2]) * second_exp;
}

void fcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t )
{

  Double_t Like=0;
  Int_t i=0;
  for(size_t i=0; i<xvar.size(); ++i)
    {
      double p = pdf_background(xvar[i], par);
      //p /= getNorm(par[0]);
      Like += log(p);
  };
   f= - 2. * Like;
}

Double_t pdf_proj(Double_t *x, Double_t *par, Double_t N, Double_t bin_width)
{
    return N * pdf_background(x[0], par) * bin_width;
}

// 


//  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
//   1  tau          3.40304e+00   1.22615e-02   3.83827e-04   4.36550e-02

//========== FIT RESULT ==========
//Parameter : tau
//tau_fit   = 1396.27 fs
//error     = 5.03091 fs
//================================

// ******* FUNCTIONS TO FIT THE INVARIANT MASS ********* 

std::vector<double> x_invariant_mass;
double X_MIN_INVARIANT_MASS = 1.8;
double X_MAX_INVARIANT_MASS = 1.95;

//p0                        =      54150.5   +/-   578.437     
//p1                        =      1.86482   +/-   2.98666e-05 
//p2                        =   0.00521984   +/-   4.64803e-05 
//p3                        =      13.5657   +/-   0.460617    
//chi2 = 55.3241
//ndof = 56
//prob = 0.500404

//--> fittando tutto il fondo
//Chi2                      =      332.502
//NDf                       =          296
//Edm                       =  5.74985e-09
//NCalls                    =          126
//p0                        =      56509.3   +/-   289.488     
//p1                        =      1.86482   +/-   2.97709e-05 
//p2                        =   0.00540974   +/-   2.72867e-05 
//p3                        =      11.6043   +/-   0.0808608 

double MAX_HISTO = 1.95;
double MIN_HISTO = 1.8;
double BIN_HISTO = 300;
double BIN_WIDTH_INVARIANT_MASS = (MAX_HISTO-MIN_HISTO)/BIN_HISTO;

TH1D * histo_data_MKpi = new TH1D("histo_data_MKpi","",BIN_HISTO,MIN_HISTO,MAX_HISTO);

Double_t pdf_double_gauss(Double_t x, Double_t *par)
{

 Double_t _fA;
 Double_t _mA;
 Double_t _mB;
 Double_t _sA;
 Double_t _sB;
 
 _fA  = par[0];	
//  _fA  = 1;			
 _mA  = par[1];
 _mB  = par[2];
 _sA  = par[3];
 _sB  = par[4];

 Double_t value = _fA*(TMath::Gaus(x,_mA,_sA,1)) + (1-_fA)*(TMath::Gaus(x,_mB,_sB,1));
 return value;
}

Double_t pdf_proj_double_gauss(Double_t *x, Double_t *par, Int_t max, Double_t bin_width)
{
  return (max*bin_width)*pdf_double_gauss(x[0],par);
}


void fcn_double_gauss(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t )
{

  Double_t Like=0;
  Int_t i=0;
  for(size_t i=0; i<x_invariant_mass.size(); ++i)
    {
      Double_t p = pdf_double_gauss(x_invariant_mass[i],par);
      Like += TMath::Log(p);
  };
   f= - 2. * Like;
}


Double_t invariant_mass_Gauss_func(Double_t *x, Double_t *par)
{
    //return par[0]*(par[6] + par[3] * TMath::Gaus(x[0],par[1],par[2],1) + (1-par[3]) * TMath::Gaus(x[0],par[4],par[5],1));
    return par[0]*(par[3] + TMath::Gaus(x[0],par[1],par[2],1));
    //return par[0]*( TMath::Gaus(x[0],par[1],par[2],1) );
}

/*
Double_t invariant_mass_Gauss_func(Double_t *x, Double_t *par)
{
  Double_t _N ; // Total normalization
  Double_t _fs; // Signal fraction; 
  Double_t _fA; // Relative fraction of Gaussian A wrt Gaussian B
  Double_t _mA; // Mean Gaussian A
  Double_t _sA; // Sigma Gaussian A
  Double_t _mB; // Mean Gaussian A
  Double_t _sB; // Sigma Gaussian B
  Double_t _slope; // Exponential slope
 
  _N  = par[0];
  _fs = par[1];
  _fA = par[2];
  _mA = par[3];
  _sA = par[4];
  _mB = par[5];
  _sB = par[6];
  _slope = par[7];

  Double_t _fb;
  Double_t _fB;
  
  _fb = 1. - _fs;
  _fB = 1. - _fA;

  Double_t value = _N*( _fs * (_fA*(TMath::Gaus(x[0],_mA,_sA,1)) + _fB*(TMath::Gaus(x[0],_mB,_sB,1))) + 
			_fb * (_slope/(TMath::Exp(_slope * X_MAX_INVARIANT_MASS)-TMath::Exp(_slope * X_MIN_INVARIANT_MASS))) * TMath::Exp(_slope * x[0])  );
 return value;
}
*/

Double_t invariant_mass_Gauss_func_fit(Double_t *x, Double_t *par)
{
   return BIN_WIDTH_INVARIANT_MASS*invariant_mass_Gauss_func(x,par);
}


TH1D * left_bkg = new TH1D("left_bkg","",100,0,10);

TH1D * right_bkg = new TH1D("rigth_bkg","",100,0,10);

TH1D * average_bkg = new TH1D("average_bkg", "", 100,0,10);

TH1D * tot_bkg = new TH1D("tot_bkg", "", 100,0,10);

left_bkg -> Sumw2();
right_bkg -> Sumw2();
average_bkg -> Sumw2();
tot_bkg -> Sumw2();

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

      if (id==1) //DATA
      { 
	      histo_data_MKpi ->Fill(M0_MKpi);

          if(M0_MKpi < 1.86482 - 3 * 0.00540974)
          {
            left_bkg -> Fill(M0_time / 410.3e-15);

            if (M0_time / 410.3e-15 >= XMIN && M0_time / 410.3e-15 <= XMAX)
            {
                xvar.push_back(M0_time / 410.3e-15);
            }

            tot_bkg -> Fill(M0_time / 410.3e-15);
          }
          if(M0_MKpi > 1.86482 + 3 * 0.00540974)
          {
            right_bkg -> Fill(M0_time / 410.3e-15);

            if (M0_time / 410.3e-15 >= XMIN && M0_time / 410.3e-15 <= XMAX)
            {
                xvar.push_back(M0_time / 410.3e-15);
            }

            tot_bkg -> Fill(M0_time / 410.3e-15);
            
          }

      }; //end DATA

      if (id==13) //MC
      {
            //double t = M0_time / 410.3e-15;

            //h_time->Fill(t);

            //if (t >= XMIN && t <= XMAX)
            //{
            //xvar.push_back(t);
            //}

      }; //end MC

      };// #### end loop over jentry

      h_time = (TH1D*)tot_bkg -> Clone("h_time");

      for(int bin = 1; bin < right_bkg->GetNbinsX(); bin++)
      {
        double left_y = left_bkg -> GetBinContent(bin);
        double right_y = right_bkg -> GetBinContent(bin);
        double sigma_left_y = left_bkg -> GetBinError(bin);
        double sigma_rigth_y = right_bkg -> GetBinError(bin);
        average_bkg -> SetBinContent(bin,(left_y+right_y)/2);
        average_bkg -> SetBinError(bin,1./2*std::sqrt(std::pow(sigma_left_y,2)+std::pow(sigma_rigth_y,2)));
      }

      // ****** FIT INVARIANT MASS ********

      double bin_min = histo_data_MKpi -> FindBin(X_MIN_INVARIANT_MASS);
      double bin_max = histo_data_MKpi -> FindBin(X_MAX_INVARIANT_MASS);
      TF1 * invariant_mass_Gauss = new TF1("invariant_mass_Gauss",invariant_mass_Gauss_func_fit,1.8,1.95,4);
      invariant_mass_Gauss -> SetParameter(0,histo_data_MKpi->Integral());
      //invariant_mass_Gauss -> SetParameter(0,histo_data_MKpi->Integral());

      invariant_mass_Gauss -> SetParameter(1,1.864);
      invariant_mass_Gauss -> SetParameter(2,0.007);

      invariant_mass_Gauss -> SetParameter(3,11.6043);

      //invariant_mass_Gauss -> SetParameter(3,0.5);
      //invariant_mass_Gauss -> SetParameter(4,1.864);
      //invariant_mass_Gauss -> SetParameter(5,0.007);
      //invariant_mass_Gauss -> SetParameter(6,4);

      //invariant_mass_Gauss -> SetParameter(1,0.5);
      //invariant_mass_Gauss -> SetParameter(2,0.5);
      //invariant_mass_Gauss -> SetParameter(3,1.864);
      //invariant_mass_Gauss -> SetParameter(4,0.007);
      //invariant_mass_Gauss -> SetParameter(5,1.864);
      //invariant_mass_Gauss -> SetParameter(6,0.007);
      //invariant_mass_Gauss -> SetParameter(7,1);
      //invariant_mass_Gauss -> SetParameter(8,1);


      histo_data_MKpi -> Fit(invariant_mass_Gauss,"L0","",X_MIN_INVARIANT_MASS,X_MAX_INVARIANT_MASS);
 
      // -------------------- PULL PLOT (INVARIANT MASS) --------------------

      const int n_bins = histo_data_MKpi->GetNbinsX();
        
      TCanvas *c_mass_pull = new TCanvas("c_mass_pull","mass pull",800,600);
        
      TPad *pad_mass = new TPad("pad_mass","",0,0.3,1,1);
      TPad *pad_pull_mass = new TPad("pad_pull_mass","",0,0,1,0.3);
        
      pad_mass->SetBottomMargin(0.02);
      pad_pull_mass->SetTopMargin(0.02);
      pad_pull_mass->SetBottomMargin(0.25);
        
      pad_mass->Draw();
      pad_pull_mass->Draw();
        
      // ---------------- TOP PAD: histogram + fit ----------------
      pad_mass->cd();
        
      histo_data_MKpi->Draw("E");
        
      invariant_mass_Gauss->SetRange(
          X_MIN_INVARIANT_MASS,
          X_MAX_INVARIANT_MASS
      );
      invariant_mass_Gauss->SetNpx(5000);
      invariant_mass_Gauss->SetLineColor(kRed);
      invariant_mass_Gauss->Draw("SAME");
      
      // ---------------- BOTTOM PAD: pulls ----------------
      pad_pull_mass->cd();
      
      TGraphErrors *mass_pulls = new TGraphErrors(n_bins);
      mass_pulls->SetMarkerStyle(7);
      mass_pulls->GetYaxis()->SetTitle("Pull");
      mass_pulls->GetXaxis()->SetTitle("m(K#pi) [GeV]");
      
      for (int i = 1; i <= n_bins; i++)
      {
      
          double x = histo_data_MKpi->GetBinCenter(i);
          double y = histo_data_MKpi->GetBinContent(i);
          double ey = histo_data_MKpi->GetBinError(i);
      
          if (x < X_MIN_INVARIANT_MASS || x > X_MAX_INVARIANT_MASS)
              continue;
      
          double y_fit = invariant_mass_Gauss->Eval(x);
      
          double pull = 0;
          if (ey > 0) pull = (y - y_fit) / ey;
      
          mass_pulls->SetPoint(i-1, x, pull);
          mass_pulls->SetPointError(i-1, 0, 1);
      }
      
      mass_pulls->SetMinimum(-5);
      mass_pulls->SetMaximum(5);
      mass_pulls->Draw("AP");
      mass_pulls->GetHistogram()->GetXaxis()->SetLimits(MIN_HISTO, MAX_HISTO);
      
      
      // zero line
      TLine *zero_line = new TLine(1.8,0,1.95,0);
      zero_line->SetLineStyle(2);
      zero_line->Draw("same");
      
      c_mass_pull->Update();

      
      // ***** FIT ***** 

      const int nparam = 3;
      int n_fit = xvar.size();
      cout << "DIMENSION " << n_fit << endl;

      TMinuit *my_gMinuit = new TMinuit(nparam);  //initialize TMinuit with a maximum of 5 params
      //gMinuit->SetPrintLevel(-1);
      my_gMinuit->SetFCN(fcn);      // set the FCN
        
      Double_t arglist[2];
      Int_t ierflg = 0;  // Error return code: 0 if the command was correctly executed, >0 otherwise. 
        
      // Set starting values and step sizes for parameters
      Double_t vstart[nparam] = {2.,3.5,0.2};
      Double_t step[nparam]   = {0.01,0.01,0.01};   //step 0 li rende costanti
     
      my_gMinuit->mnparm(0, "tau1", vstart[0], step[0], 0.01, 10., ierflg);
      my_gMinuit->mnparm(1, "tau2", vstart[1], step[1], 0.01, 10., ierflg);
      my_gMinuit->mnparm(2, "fraction", vstart[2], step[2], 0., 1., ierflg);

      arglist[0] = 500.;//500;
      arglist[1] = 0.1;

      my_gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
      
      // Print results
      Double_t amin,edm,errdef;
      Int_t nvpar,nparx,icstat;
      my_gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);     
      my_gMinuit->mnprin(3,amin);
   
      Double_t tau1_fit, etau1_fit, tau2_fit, etau2_fit, fraction, efraction,  bnd1, bnd2;
      Int_t ivar=0;
      TString chnam;
        
      my_gMinuit->mnpout(0, chnam, tau1_fit, etau1_fit, bnd1, bnd2, ivar);
      my_gMinuit->mnpout(1, chnam, tau2_fit, etau2_fit, bnd1, bnd2, ivar);
      my_gMinuit->mnpout(2, chnam, fraction, efraction, bnd1, bnd2, ivar);

      //std::cout << "\n========== FIT RESULT ==========\n";
      //std::cout << "Parameter : " << chnam << std::endl;
      //std::cout << "tau_fit   = " << tau_fit*410.3 << " fs" << std::endl;
      //std::cout << "error     = " << etau_fit*410.3 << " fs" << std::endl;
      //std::cout << "================================\n";

      std::vector<double> x_points;
      std::vector<double> y_points;

      Double_t pars[nparam] = {tau1_fit,tau2_fit,fraction};

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
      fit_function_plot->SetLineColor(kRed);
      fit_function_plot->SetMarkerColor(kRed);
      fit_function_plot->Draw("L SAME");

      // IMPORTANT: go back to canvas
      c_fit->cd();

      // -------------------- BOTTOM PAD (PULLS) --------------------
      TPad *pad_pull = new TPad("pad_pull", "", 0., 0., 1., split);
      pad_pull->SetTopMargin(0.02);
      pad_pull->SetBottomMargin(0.25);
      pad_pull->Draw();
      pad_pull->cd();
        
      int n = h_time->GetNbinsX();
        
      TGraphErrors *pulls = new TGraphErrors(n);
      pulls->SetMarkerStyle(7);
        
      pulls->GetXaxis()->SetTitle("t / 410.3 fs");
      pulls->GetYaxis()->SetTitle("Pull");
        
      // optional: axis range for visibility
      pulls->SetMinimum(-5);
      pulls->SetMaximum(5);
        
      for (int i = 0; i < n; i++)
      {
          double x  = h_time->GetBinCenter(i+1);
          double y  = h_time->GetBinContent(i+1);
          double ey = h_time->GetBinError(i+1);
      
          double y_fit = pdf_proj(&x, pars, n_fit, h_time->GetBinWidth(1));
      
          // avoid division by zero
          double pull = 0;
          if (ey > 0) pull = (y - y_fit) / ey;
      
          pulls->SetPoint(i, x, pull);
          pulls->SetPointError(i, 0., 1.0);
      }
      
      pulls->Draw("AP");
      
      // baseline at 0
      TLine *line_zero = new TLine(0.0, 0., 10., 0.);
      line_zero->SetLineStyle(2);
      line_zero->Draw("same");
      // final update
      c_fit->Update();
      c_fit->Draw();
      
    
      TFile * f = new TFile("outfile_fit_DATA_BKG.root","RECREATE");
      //c_fit -> Write();
      histo_data_MKpi -> Write();
      left_bkg->Scale(1./left_bkg->Integral());
      right_bkg->Scale(1./right_bkg->Integral());
      average_bkg->Scale(1./average_bkg->Integral());
      left_bkg -> Write();
      right_bkg -> Write();
      average_bkg -> Write();
      c_fit -> Write();
      c_mass_pull->Write();

      f->Close();

}
