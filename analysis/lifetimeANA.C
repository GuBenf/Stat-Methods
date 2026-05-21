#define lifetimeANA_cxx
#include "lifetimeANA.h"

// #include "fitterLike.C"

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

using namespace std;

std::vector<double> xvar;

// one vector for each proper-time bin
std::array<std::vector<double>,100> timeDiffData;

TRandom3 rnd(0);

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
   fA = par[1];
   muA = par[2];
   sigmaA = par[3];
   fB = par[4];
   muB = par[5];
   sigmaB = par[6];

   Double_t value = 1./(2*tau)*TMath::Exp((sigmaA*sigmaA - 2*tau*x)/(2*tau*tau))*( TMath::Erf((tau*x - sigmaA*sigmaA)/(TMath::Sqrt(2)*tau*sigmaA)) +1 );
   return value;
}


Double_t pdf(Double_t x, Double_t *par)
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

Double_t pdf_proj(Double_t *x, Double_t *par, Int_t max, Double_t bin_width)
{
  return (max*bin_width)*pdf(x[0],par);
}


void fcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t )
{

  Double_t Like=0;
  Int_t i=0;
  for(size_t i=0; i<xvar.size(); ++i)
    {
      Double_t p = pdf(xvar[i],par);
      Like += TMath::Log(p);
  };
   f= - 2. * Like;
}

Double_t acceptance_func(Double_t *x, Double_t *par)
{
   Double_t f;
   if(x[0] > 0.5) f = par[0]*(1 + TMath::Erf((TMath::Power(x[0],par[3]) - par[1])/par[2]));
   else f = 0;
   return f;
}

// p0   =    0.0579004   +/-   0.000312769 
// p1   =      1.15429   +/-   0.13716     
// p2   =    -0.289437   +/-   0.247952    
// p3   =    -0.489007   +/-   0.422909 

//Create a new file to store histograms
TFile *histo_file = new TFile("./_root/histo_file_new.root","RECREATE","put a title");

void lifetimeANA::Loop()
{
//   In a ROOT session, you can do:
//      root> .L lifetimeANA.C
//      root> lifetimeANA t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//



  // Create histograms
  TH1D *histo_data_MKpi    = new TH1D("histo_data_MKpi","",100,1.8,1.95);
  TH1D *histo_mc_MKpi      = new TH1D("histo_mc_MKpi","",100,1.8,1.95);
  TH1D *histo_mc_time      = new TH1D("histo_mc_time","",100,0,10);
  TH1D *histo_mc_time_true = new TH1D("histo_mc_time_true","",100,0,10);
  
  TH2D *xy = new TH2D("xy","",100,0,10,100,0,1);
  TH1D *extracted_time_tmp = new TH1D("extracted_time_tmp","",100,0,10); 
  TH1D *extracted_time = new TH1D("extracted_time","",100,0,10); 

  TH1D *histo_mc_time_diff_meas_true = new TH1D("histo_mc_time_diff_meas_true","",200,-1,1);

  std::array<TH1D*,100> histo_diff_mc_time_bins;
  std::array<double,100> histo_diff_mc_time_bincenters;



  for(int h = 0; h < histo_diff_mc_time_bins.size(); h++)
  {
      histo_diff_mc_time_bins[h] = new TH1D(Form("histo_diff_mc_time_bin_%d",h),"",200,-1,1);
  }

  histo_data_MKpi        ->Sumw2();
  histo_mc_MKpi          ->Sumw2();
  histo_mc_time          ->Sumw2();
  
  
  

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // #### My code #########
      if (id==1) //DATA
      { 
	      histo_data_MKpi ->Fill(M0_MKpi);
      }; //end DATA
      
      if (id==13) //MC
      { 
	      histo_mc_MKpi ->Fill(M0_MKpi);
	      histo_mc_time->Fill(M0_time/(410.3e-15));
         histo_mc_time_true->Fill(M0_time_true/(410.3e-15));

         //to check time resolution
         histo_mc_time_diff_meas_true->Fill( (M0_time-M0_time_true)/410.3e-15 );

         double time_diff = (M0_time-M0_time_true)/410.3e-15;

         for(int h = 0; h < histo_diff_mc_time_bins.size(); h++)
         {
            //we scan from proper decay time 0 to 10 (in unit of 410.3 fs)
            double bin_min = 10./histo_diff_mc_time_bins.size()*h; 
            double bin_max = 10./histo_diff_mc_time_bins.size()*(h+1);
            histo_diff_mc_time_bincenters[h] = (bin_min+bin_max)/2;
            if(M0_time/410.3e-15 < bin_max && M0_time/410.3e-15 >= bin_min)
            {
                  histo_diff_mc_time_bins[h]->Fill(time_diff);
                  timeDiffData[h].push_back(time_diff);
            }
         }

      }; //end MC
      
   };// #### end loop over jentry
   
   // TMinuit *my_gMinuit = new TMinuit(3);  //initialize TMinuit with a maximum of 5 params
   // my_gMinuit->SetFCN(fcn);      // set the FCN
   
   // Double_t arglist[2];
   // Int_t ierflg = 0;  // Error return code: 0 if the command was correctly executed, >0 otherwise. 
   
   // arglist[0] = 1;                           
   // my_gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
   
   // Set starting values and step sizes for parameters
   // Double_t vstart[5] = {1,  0.015  , 0.1 };
   // Double_t vstart[5] = {0.53,  -0.34  , -0.2  , 0.014 , 0.017 };
   // Double_t step[5]   = {0.001, 0.001, 0.001};   //step 0 li rende costanti

   // my_gMinuit->mnparm(0, "fA", vstart[0], step[0], 0., 1., ierflg);
   // my_gMinuit->mnparm(1, "mA", vstart[1], step[1], 0., 0., ierflg);
   // my_gMinuit->mnparm(2, "mB", vstart[2], step[2], 0., 0., ierflg);
   // my_gMinuit->mnparm(2, "sA", vstart[3], step[3], 0., 0., ierflg);
   // my_gMinuit->mnparm(4, "sB", vstart[4], step[4], 0., 0., ierflg);
   // arglist[0] = 500.;//500;
   // arglist[1] = 0.1;

   /////////////////////////////////
   // RESOLUTION
   /////////////////////////////////

   TGraphErrors *g_resolution = new TGraphErrors();
   g_resolution->SetMaximum(0.15);
   g_resolution->SetMinimum(0.);

   TGraphErrors *g_resolution2 = new TGraphErrors();
   g_resolution2->SetMaximum(0.3);
   g_resolution2->SetMinimum(0.);


   for(int h = 0; h < histo_diff_mc_time_bins.size(); h++)
   {
      // if(h>= 4 && h<60){
      if(h>= 0 && h<100)
      {
         cout << endl;
         cout << "=================================" << endl;
         cout << "FITTING TIME BIN " << h << endl;
         cout << "ENTRIES = "
              << timeDiffData[h].size()
              << endl;
         cout << "=================================" << endl;

         // current active dataset
         xvar = timeDiffData[h];
         TMinuit *my_gMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
         my_gMinuit->SetFCN(fcn);      // set the FCN
        
         Double_t arglist[2];
         Int_t ierflg = 0;  // Error return code: 0 if the command was correctly executed, >0 otherwise. 
        
         // arglist[0] = 1;                           
         // my_gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
        
         // Set starting values and step sizes for parameters
         Double_t vstart[5] = {1,  0.015  , 0.015,  0.09 , 0.13};
         //Double_t vstart[5] = {0.53,  -0.34  , -0.2  , 0.014 , 0.017 };
         Double_t step[5]   = {0.01, 0.001, 0.001, 0.001, 0.001};   //step 0 li rende costanti
     
         my_gMinuit->mnparm(0, "fA", vstart[0], step[0], 0., 1., ierflg);
         my_gMinuit->mnparm(1, "mA", vstart[1], step[1], 0., 0., ierflg);
         my_gMinuit->mnparm(2, "mB", vstart[2], step[2], 0., 0., ierflg);
         my_gMinuit->mnparm(3, "sA", vstart[3], step[3], 0., 0., ierflg);
         my_gMinuit->mnparm(4, "sB", vstart[4], step[4], 0., 0., ierflg);
         arglist[0] = 500.;//500;
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

         g_resolution->SetPoint(h, histo_diff_mc_time_bincenters[h], val3);
         g_resolution->SetPointError(h, 0.0, error3);

         g_resolution2->SetPoint(h, histo_diff_mc_time_bincenters[h], val4);
         g_resolution2->SetPointError(h, 0.0, error4);


         std::vector<double> x_points;
         std::vector<double> y_points;

         Double_t pars[5] = {val0,val1,val2,val3,val4};

         for(int xi=0; xi<10000; xi++)
         {
            double xmin = histo_diff_mc_time_bins[h]->GetXaxis()->GetXmin();
            double xmax = histo_diff_mc_time_bins[h]->GetXaxis()->GetXmax();
        
            double x = xmin + xi * (xmax-xmin)/10000;

            // cout << x << endl;
            x_points.push_back(x);
            y_points.push_back(pdf_proj(&x,pars,histo_diff_mc_time_bins[h]->GetEntries(),histo_diff_mc_time_bins[h]->GetBinWidth(1)));
         }

         TGraph * plot_f = new TGraph(10000,x_points.data(),y_points.data());
         plot_f->SetMaximum(histo_diff_mc_time_bins[h]->GetMaximum()*1.05);
         plot_f->SetLineColor(kRed);

         TCanvas *c_plot = new TCanvas(Form("plot_fit_%d",h));
         c_plot->cd();
         plot_f->Draw();
         histo_diff_mc_time_bins[h]->Draw("same");
         histo_file->cd();
         c_plot->Write();
      }
   };
   
   // Double_t matrix[5][5];
   // Int_t n=5;
   // my_gMinuit->mnemat (&matrix[0][0],5); // Accessing Covariance Matrix

   // Write on .txt file covariance Matrix 
   // FILE *file_cov = fopen ("_txt/covariance_matrix.txt", "w");
   // Int_t icol=0;
   // Int_t irow=0;
   // if (n==1)
   // {
   //    fprintf (file_cov, "%e \n", matrix[0][0]); 
   // }
   //  else
   //  {
   //    for (irow=0; irow<n; irow++)
   //    {
   //       for (icol=0; icol<n; icol++)
   //       { 
   //         fprintf (file_cov, "%e ", matrix[irow][icol]); 
   //       };
   //       fprintf (file_cov, "\n");
   //    };
   //    fprintf (file_cov, "\n");   
   //  }
   //  fclose (file_cov);
   

   /////////////////////////////
   // ACCEPTANCE
   /////////////////////////////


   std::array<TH1D*,10> histo_mean_lives;
   std::array<TH1D*,10> histo_mc_scaled;

   int estrazione = 0;
   int eventi_buoni = 0;
   Double_t fixed_par[7] = {1., 4.63990e-01, 6.71361e-03, 8.41433e-02, 1-4.63990e-01 , 1.56500e-02, 1.28286e-01};
   Double_t mean_lives[10] = {0.1, 0.5, 0.8, 1., 1.2, 1.5, 1.8, 2., 2.5, 3.};
   Double_t fixed_par_h[7];
   for (int k = 1; k < 7; k++) {
      fixed_par_h[k] = fixed_par[k];
   }

   for(int h = 0; h < histo_mean_lives.size(); h++)
   {
         histo_mean_lives[h] = new TH1D(Form("histo_mean_lives_%g",mean_lives[h]),"",200,0,10);
         histo_mc_scaled[h] = (TH1D*)histo_mc_time->Clone(Form("histo_mc_scaled%g",mean_lives[h]));
   }

   for(int hist=0; hist < histo_mean_lives.size(); hist++){
      eventi_buoni = 0;
   while(eventi_buoni < 1e7)
   {
      fixed_par_h[0] = mean_lives[hist];

      double sim_time = rnd.Uniform(0,10);
      double y = rnd.Uniform(0,1.0/fixed_par_h[0]);

      //cout << sim_time << " " << y << " " << pdf_exp_x_gauss(sim_time,fixed_par) << endl;

      if(y < pdf_exp_x_gauss(sim_time,fixed_par_h))
      {
         // xy->Fill(sim_time,y);
         // extracted_time_tmp->Fill(sim_time);
         histo_mean_lives[hist]->Fill(sim_time);
         eventi_buoni++;
      }

      estrazione++;
   }
   }

   // Int_t control_integral = extracted_time_tmp->Integral(17,70);
   // Double_t control_ratio = (Double_t)control_integral/1e7;
   // Int_t mc_integral = histo_mc_time->Integral(17,70);
   // Double_t norm_ratio = (Double_t)mc_integral/control_integral;

   int estrazione2 = 0;
   int eventi_buoni2 = 0;
   // while(eventi_buoni2 < 1e7*norm_ratio)
   while(eventi_buoni2 < 1e7)
   {
      double sim_time = rnd.Uniform(0,10);
      double y = rnd.Uniform();

      //cout << sim_time << " " << y << " " << pdf_exp_x_gauss(sim_time,fixed_par) << endl;

      if(y < pdf_exp_x_gauss(sim_time,fixed_par))
      {
         xy->Fill(sim_time,y);
         extracted_time->Fill(sim_time);
         eventi_buoni2++;
      }

      estrazione2++;
   }
   
   std::array<TGraphErrors*, 10> g_ratios;
   std::array<TGraph*, 10> g_acceptances;


   for(int h=0; h< g_ratios.size(); h++){
      TGraphErrors *g_ratio = new TGraphErrors();
      TGraph *g_acceptance = new TGraphErrors();

      g_ratio->SetName(Form("g_ratio_lives_%g", mean_lives[h]));
      g_acceptance->SetName(Form("g_acceptance_%g", mean_lives[h]));

      for(int bin = 1; bin < histo_mean_lives[h]->GetNbinsX()+1; bin ++)
   {
      double ratio = histo_mean_lives[h]->GetBinContent(bin)/extracted_time->GetBinContent(bin);
      histo_mc_scaled[h]->SetBinContent(bin, histo_mc_scaled[h]->GetBinContent(bin)*ratio);
      double ratio_mc = histo_mean_lives[h]->GetBinContent(bin)/histo_mc_scaled->GetBinContent(bin);


      //is is ok to use the standard error propagation ? --> how about the errors on x ? 

      // ratio = B/A

      double A = extracted_time->GetBinContent(bin);
      double sigmaA = extracted_time->GetBinError(bin);
      double B = histo_mean_lives[h]->GetBinContent(bin);
      double sigmaB = histo_mean_lives[h]->GetBinError(bin);

      double ratio_error = ratio * TMath::Sqrt(std::pow(sigmaA/A,2) + std::pow(sigmaB/B,2)); //--> covariance ? 

      //Error on x: bin width/2, otherwise you get a doubled error

      g_ratio->SetPoint(bin,extracted_time->GetBinCenter(bin),ratio);
      g_ratio->SetPointError(bin,extracted_time->GetBinWidth(bin)/2.0,ratio_error);

      g_acceptance->SetPoint(bin,extracted_time->GetBinCenter(bin),ratio_mc);

   } 
   }

   std::array<Color_t,10> colors = {kOrange, kRed, kPink+1, kBlack, kMagenta, kViolet+1, kBlue,kCyan, kGreen, kAzure+7};
   std::array<std::string, 10> labels;

   for(int h = 0; h < labels.size(); h++)
   {   
      std::string s = std::format("{:.1f}", mean_lives[h]);

      labels[h] = ("Mean life: " + s);
   }


   TCanvas *c_ratio_lives = new TCanvas("ratio_lives","");
   c_ratio_lives -> cd();
   TLegend *l = new TLegend();
   g_ratios[0]-> Draw("A");
   g_ratios[0]-> SetLineColor(colors[0]);

   l->AddEntry(g_ratios[0], labels[0].c_str());

   for(int h=1; h< g_ratios.size(); h++){
      g_ratios[h]-> Draw("same");
      g_ratios[h]-> SetLineColor(colors[h]);

      l -> AddEntry(g_ratios[h], (labels[h]).c_str());

   }

   TCanvas *c_acceptances_lives = new TCanvas("acceptance_lives","");
   c_acceptances_lives -> cd();
   TLegend *l = new TLegend();
   g_acceptances[0]-> Draw("A");
   g_acceptances[0]-> SetLineColor(colors[0]);

   l->AddEntry(g_acceptances[0], labels[0].c_str());

   for(int h=1; h< g_acceptances.size(); h++){
      g_acceptances[h]-> Draw("same");
      g_acceptances[h]-> SetLineColor(colors[h]);

      l -> AddEntry(g_acceptances[h], (labels[h]).c_str());

   }


   TGraphErrors *acceptance = new TGraphErrors();

   for(int bin = 1; bin < histo_mc_time->GetNbinsX()+1; bin ++)
   {
      double ratio = histo_mc_time->GetBinContent(bin)/extracted_time->GetBinContent(bin);

      //is is ok to use the standard error propagation ? --> how about the errors on x ? 

      // ratio = B/A

      double A = extracted_time->GetBinContent(bin);
      double sigmaA = extracted_time->GetBinError(bin);
      double B = histo_mc_time->GetBinContent(bin);
      double sigmaB = histo_mc_time->GetBinError(bin);

      double ratio_error = ratio * TMath::Sqrt(std::pow(sigmaA/A,2) + std::pow(sigmaB/B,2)); //--> covariance ? 

      //Error on x: bin width/2, otherwise you get a doubled error

      cout << extracted_time->GetBinCenter(bin) << " " << ratio << endl;

      acceptance->SetPoint(bin,extracted_time->GetBinCenter(bin),ratio);
      acceptance->SetPointError(bin,extracted_time->GetBinWidth(bin)/2.0,ratio_error);

   } 

   // ACCEPTANCE FIT

   cout << "/////////////////////////////" << endl;
   cout << "START ACCEPTANCE FIT" << endl;

   TF1 *acceptance_tf1 = new TF1("acceptance_tf1",acceptance_func,0,10,4);
   acceptance_tf1 -> SetParameter(0,0.05); //--> accettanza è definita a meno di una costante moltiplicativa
   acceptance_tf1 -> SetParameter(1,0.9);
   acceptance_tf1 -> SetParameter(2,0.25);
   acceptance_tf1 -> SetParameter(3,0.6);

   acceptance -> Fit(acceptance_tf1,"","",0.,10);

   Double_t fit_chi2 = acceptance_tf1->GetChisquare();
   Int_t    fit_ndof = acceptance_tf1->GetNDF();
   Double_t fit_prob = acceptance_tf1->GetProb();

   cout << "chi2 = " << fit_chi2 << endl;
   cout << "ndof = " << fit_ndof << endl;
   cout << "prob = " << fit_prob << endl;

   TCanvas *c_acceptance_fit = new TCanvas("acceptance_fit","");
   c_acceptance_fit -> cd();
   acceptance -> Draw("AP");
   acceptance_tf1 -> Draw("same");
   acceptance_tf1 -> SetLineColor(kRed);
   

   /*
      //Create Canvas
      TCanvas *c_histo_data_MKpi = new TCanvas("c_histo_data_MKpi","canvas histo",500,500);
      c_histo_data_MKpi->cd();
      histo_data_MKpi->Draw();
      c_histo_data_MKpi->Print("./_fig/c_histo_data_MKpi.pdf");
      c_histo_data_MKpi->Print("./_fig/c_histo_data_MKpi.eps");


      //Create Canvas
      TCanvas *c_histo_mc_MKpi = new TCanvas("c_histo_mc_MKpi","canvas histo",500,500);
      c_histo_mc_MKpi->cd();
      histo_mc_MKpi->Draw();
      c_histo_mc_MKpi->Print("./_fig/c_histo_mc_MKpi.pdf");
      c_histo_mc_MKpi->Print("./_fig/c_histo_mc_MKpi.eps");

      //Create Canvas
      TCanvas *c_histo_mc_time = new TCanvas("c_histo_mc_time","canvas histo",500,500);
      c_histo_mc_time->cd();
      histo_mc_time->Draw();
      c_histo_mc_time->Print("./_fig/c_histo_mc_time.pdf");
      c_histo_mc_time->Print("./_fig/c_histo_mc_time.eps");
   */

   TCanvas *c1 = new TCanvas("c1", "Acceptance", 800, 600);
   c1 -> cd();
   acceptance->SetTitle("Acceptance;Time;Ratio");
   acceptance->SetMarkerStyle(20);

   acceptance->Draw("AP E1");  // axes + points + error bars

   c1->Update();
    
  
   histo_file->cd();

   histo_data_MKpi->Write();
   histo_mc_MKpi->Write();
   histo_mc_time->Write();
   histo_mc_time_true->Write();

   histo_mc_time->Scale(1./histo_mc_time->Integral());
   histo_mc_time->Write("histo_mc_time_norm");
   histo_mc_time_true->Scale(1./histo_mc_time_true->Integral());
   histo_mc_time_true->Write("histo_mc_time_true_norm");

   histo_mc_time_diff_meas_true->Write();

   for(int h = 0; h < histo_diff_mc_time_bins.size(); h++)
   {   
      histo_diff_mc_time_bins[h] -> Write();
   }

   for(int h = 0; h < histo_mean_lives.size(); h++)
   {   
      std::string s = std::format("{:.1f}", mean_lives[h]);

      histo_mean_lives[h] -> Write(("Mean life: " + s).c_str());
   }

   acceptance-> Write("acceptance_graph");
       g_resolution->SetTitle(
        "Resolution vs Proper Time;"
        "Proper Time;"
        "Resolution"
    );

   g_resolution->SetMarkerStyle(20);

   g_resolution->Draw("AP E1");

   g_resolution->Write("resolution_graph");

   g_resolution2->SetTitle(
        "Resolution2 vs Proper Time;"
        "Proper Time;"
        "Resolution2"
    );

   g_resolution2->SetMarkerStyle(20);

   g_resolution2->Draw("AP E1");

   g_resolution2->Write("resolution2_graph");


   xy->Write();
   extracted_time->Write();
   extracted_time->Scale(1./extracted_time->Integral());
   extracted_time->Write("extracted_time_norm");

   c_acceptance_fit -> Write();

   c_ratio_lives->Write();
   c_acceptances_lives->Write();

   histo_file->Close();




}
