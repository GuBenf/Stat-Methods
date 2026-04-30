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



  // Create histograms
  TH1D *histo_data_MKpi    = new TH1D("histo_data_MKpi","",100,1.8,1.95);
  TH1D *histo_mc_MKpi      = new TH1D("histo_mc_MKpi","",100,1.8,1.95);
  TH1D *histo_mc_time      = new TH1D("histo_mc_time","",100,0,10);
  TH1D *histo_mc_time_true = new TH1D("histo_mc_time_true","",100,0,10);

  TH1D *histo_mc_time_diff_meas_true = new TH1D("histo_mc_time_diff_meas_true","",200,-1,1);

  std::array<TH1D*,10> histo_diff_mc_time_bins;

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
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // #### My code #########
      if ( id==1){ //DATA
	
	histo_data_MKpi ->Fill(M0_MKpi);
	// Print out the Mpipi values
	//std::cout <<M0_MKpi << std::endl;
	
      }; //end DATA
      
      if (id==13){ //MC
	histo_mc_MKpi ->Fill(M0_MKpi);
	histo_mc_time->Fill(M0_time/(410.3e-15));
      histo_mc_time_true->Fill(M0_time_true/(410.3e-15));

      //to check time resolution
      histo_mc_time_diff_meas_true->Fill( (M0_time-M0_time_true)/410.3e-15 );

      double time_diff = (M0_time-M0_time_true)/410.3e-15;

      for(int h = 0; h < histo_diff_mc_time_bins.size(); h++)
      {
            //we scan from proper decay time 0 to 10 (in unit of 410.3 fs)
            int bin_min = 10./histo_diff_mc_time_bins.size()*h; 
            int bin_max = 10./histo_diff_mc_time_bins.size()*(h+1);
            if(M0_time/410.3e-15 < bin_max && M0_time/410.3e-15 >= bin_min)
            {
                  histo_diff_mc_time_bins[h]->Fill(time_diff);
            }
      }

      //to check if the time resolution varies as a function of the proper deacy time
      //--> check difference between measured decay time and sim in bins of the proper decay time

	//std::cout <<M0_MKpi << std::endl;
	
      }; //end MC
      
   };// #### end loop over jentry


   //here we should fit the histogram with the time difference to get the resolution -->
   //Micheal said that we should not use ROOFit because it is easy to write the function in root and to use the Fitter examples :)

   //here we should fit each of the 10 histograms to check the resolution as a function of the proper time -->

   //check the acceptance as a function of the proper decay time 
   // --> the empirical function will be 100 points with bins of proper time as x and the ratio between the proper time mc sim and meas
   // --> almeno io ho capito così oggi in classe :)

   TGraphErrors *acceptance = new TGraphErrors();

   for(int bin = 0; bin < histo_mc_time->GetNbinsX(); bin ++)
   {
      double ratio = histo_mc_time_true->GetBinContent(bin)/histo_mc_time->GetBinContent(bin);

      //is is ok to use the standard error propagation ? --> how about the errors on x ? 

      // ratio = A/B

      double A = histo_mc_time_true->GetBinContent(bin);
      double sigmaA = histo_mc_time_true->GetBinError(bin);
      double B = histo_mc_time->GetBinContent(bin);
      double sigmaB = histo_mc_time->GetBinError(bin);

      double ratio_error = ratio * TMath::Sqrt(std::pow(sigmaA/A,2) + std::pow(sigmaB/B,2)); //--> covariance ? 

      //cout << ratio << " " << ratio_error << endl;

      acceptance->SetPoint(bin,histo_mc_time_true->GetBinCenter(bin),ratio);
      acceptance->SetPointError(bin,histo_mc_time_true->GetBinWidth(bin),ratio_error);

   } //non mi plotta bene il TGraph ma i valori ci sono, non mi ricordo come si faceva con SetPoint :)



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


    
  //Create a new file to store histograms
   TFile *histo_file = new TFile("./_root/histo_file_new.root","RECREATE","put a title");
   histo_file->cd();

   histo_data_MKpi->Write();
   histo_mc_MKpi->Write();
   histo_mc_time->Write();

   histo_mc_time_diff_meas_true->Write();

   for(int h = 0; h < histo_diff_mc_time_bins.size(); h++)
   {   
      histo_diff_mc_time_bins[h] -> Write();
   }

   acceptance -> Write();

   histo_file->Close();

}
