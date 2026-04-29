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

using namespace std;

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

  std::vector<double> v_M0_time_true;

  std::array<std::vector<double>,50> time_bins;

  // Create histograms
  TH1D *histo_data_MKpi    = new TH1D("histo_data_MKpi","",100,1.8,1.95);
  TH1D *histo_mc_MKpi      = new TH1D("histo_mc_MKpi","",100,1.8,1.95);

  TH1D *histo_mc_true_M0_time_true = new TH1D("histo_mc_true_M0_time_true","",50,0,4e-12);

  histo_data_MKpi        ->Sumw2();
  histo_mc_MKpi          ->Sumw2();
  
  

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // #### My code #########
      if ( id==1){ //DATA
	
	//histo_data_MKpi ->Fill(M0_MKpi);
	// Print out the Mpipi values
	//std::cout <<M0_MKpi << std::endl;
	
      }; //end DATA
      
      if (id==13){ //MC
	//histo_mc_MKpi ->Fill(M0_MKpi);

      histo_mc_true_M0_time_true -> Fill(M0_time_true);

      v_M0_time_true.push_back(M0_time_true);

      int i = 0;
      double bin_width = 4e-12/50;

      double theta1 = TMath::ACos(h1_pz/(std::sqrt(h1_px*h1_px + h1_py*h1_py + h1_pz*h1_pz)));
      double theta2 = TMath::ACos(h2_pz/(std::sqrt(h2_px*h2_px + h2_py*h2_py + h2_pz*h2_pz)));

      while(i < 50)
      {

            double min_bin = bin_width * i;
            double max_bin = bin_width * (i+1);

            if(M0_time_true < max_bin && M0_time_true >= min_bin)
            {
                  time_bins[i].push_back(theta1);
                  time_bins[i].push_back(theta2);
                  break;
            }

            i++;
      }

      //cout << M0_time_true << endl;

	//std::cout <<M0_MKpi << std::endl;
	
      }; //end MC
      
   };// #### end loop over jentry



      TH1D *theta_time_bin = new TH1D("theta_time_bin","",100,0,0.35);
      for(int i = 0; i<time_bins[30].size(); i++)
      {
            theta_time_bin->Fill(time_bins[30][i]);
      }

      ////////////////////////
      ///    FIT      ////////
      ////////////////////////
      RooRealVar time("time", "M0 time true", 0., 5e-12);

      RooDataSet data("data", "data", RooArgSet(time));

      for(int i =0; i<int(v_M0_time_true.size()); i++)
        {

            time.setVal(v_M0_time_true[i]);
            if(v_M0_time_true[i]<0){cout << v_M0_time_true[i] << endl;}
            data.add(RooArgSet(time));
        }

      RooRealVar meanG("meanG", "Mean of Gaussian", 4e-13, 1e-13, 1e-12);
      RooRealVar sigmaG("sigmaG", "Gaussian sigma", 5e-13, 1e-14, 1e-12);
      RooGaussian gauss("gauss", "Gaussian PDF", time, meanG, sigmaG);

      RooRealVar lambda("lambda", "lambda of the Exponential", -1./4.1e-13, -1./1e-13, -1./1e-11);
      RooExponential exp("exp", "Exponential PDF", time, lambda);

      RooFFTConvPdf gauss_exp("gauss_exp", "Exponential convoluted with Gaussian", time, exp, gauss);

      RooFitResult* fitRes = gauss_exp.fitTo(data);

      RooPlot* frame = time.frame();

      data.plotOn(frame, RooFit::Name("data"));
 
      gauss_exp.plotOn(frame,RooFit::Name("fit"));
      //gauss_exp.paramOn(frame);

      TCanvas *c_fit = new TCanvas("c_fit","c_fit");
      c_fit->cd();
      frame->Draw();

      TCanvas *c_res = new TCanvas("c_res","c_res");
      c_res->cd();

      RooHist *hresid = frame->residHist("data", "fit");
      hresid->GetYaxis()->SetTitle("Data - Fit");
      hresid->GetXaxis()->SetTitle("M0 time true");
      hresid->Draw("AP"); 

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

*/
    
  //Create a new file to store histograms
   TFile *histo_file = new TFile("./_root/histo_file_new.root","RECREATE","put a title");

   histo_mc_true_M0_time_true->Write();
   theta_time_bin->Write();

   //histo_file->cd();
   //histo_data_MKpi->Write();
   //histo_mc_MKpi->Write();
   //histo_file->Close();
   

}
