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
std::array<std::vector<double>,10> timeDiffData;

TRandom3 randGen(0);
Double_t mA= -1.;
// Double_t mB= 1.;
Double_t sA= 0.5;
// Double_t sB= 0.5;

Double_t pdf(Double_t x, Double_t *par)
{

 Double_t _fA;
 Double_t _mA;
//  Double_t _mB;
 Double_t _sA;
//  Double_t _sB;
 
 _fA  = par[0];				
 _mA  = par[1];
//  _mB  = par[2];
 _sA  = par[2];
//  _sB  = par[4];

 Double_t value = _fA*(TMath::Gaus(x,_mA,_sA,1)); // + (1-_fA)*(TMath::Gaus(x,_mB,_sB,1));
 return value;
}

Double_t pdf_proj(Double_t *x, Double_t *par, Int_t max)
{
  Double_t bin_width;
  bin_width =0.01;
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
  std::array<double,10> histo_diff_mc_time_bincenters;

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
            double bin_min = 10./histo_diff_mc_time_bins.size()*h; 
            double bin_max = 10./histo_diff_mc_time_bins.size()*(h+1);
            histo_diff_mc_time_bincenters[h] = (bin_min+bin_max)/2;
            if(M0_time/410.3e-15 < bin_max && M0_time/410.3e-15 >= bin_min)
            {
                  histo_diff_mc_time_bins[h]->Fill(time_diff);
                  timeDiffData[h].push_back(time_diff);
            }
      }

      //to check if the time resolution varies as a function of the proper deacy time
      //--> check difference between measured decay time and sim in bins of the proper decay time

	//std::cout <<M0_MKpi << std::endl;
	
      }; //end MC
      
   };// #### end loop over jentry
   
   //here we should fit the histogram with the time difference to get the resolution -->
   //Micheal said that we should not use ROOFit because it is easy to write the function in root and to use the Fitter examples :)
   
//    TMinuit *my_gMinuit = new TMinuit(3);  //initialize TMinuit with a maximum of 5 params
//    my_gMinuit->SetFCN(fcn);      // set the FCN
   
//    Double_t arglist[2];
//    Int_t ierflg = 0;  // Error return code: 0 if the command was correctly executed, >0 otherwise. 
   
//    // arglist[0] = 1;                           
//    // my_gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
   
//    // Set starting values and step sizes for parameters
//    Double_t vstart[5] = {1,  0.015  , 0.1 };
//    //Double_t vstart[5] = {0.53,  -0.34  , -0.2  , 0.014 , 0.017 };
//    Double_t step[5]   = {0.001, 0.001, 0.001};   //step 0 li rende costanti

//    my_gMinuit->mnparm(0, "fA", vstart[0], step[0], 0., 1., ierflg);
//    my_gMinuit->mnparm(1, "mA", vstart[1], step[1], 0., 0., ierflg);
// //    my_gMinuit->mnparm(2, "mB", vstart[2], step[2], 0., 0., ierflg);
//    my_gMinuit->mnparm(2, "sA", vstart[3], step[3], 0., 0., ierflg);
// //    my_gMinuit->mnparm(4, "sB", vstart[4], step[4], 0., 0., ierflg);
//    arglist[0] = 500.;//500;
//    arglist[1] = 0.1;

   TGraphErrors *g_resolution = new TGraphErrors();
   for(int h = 0; h < histo_diff_mc_time_bins.size(); h++)
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
        TMinuit *my_gMinuit = new TMinuit(3);  //initialize TMinuit with a maximum of 5 params
        my_gMinuit->SetFCN(fcn);      // set the FCN
        
        Double_t arglist[2];
        Int_t ierflg = 0;  // Error return code: 0 if the command was correctly executed, >0 otherwise. 
        
        // arglist[0] = 1;                           
        // my_gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
        
        // Set starting values and step sizes for parameters
        Double_t vstart[3] = {1,  0.015  , 0.1 };
        //Double_t vstart[5] = {0.53,  -0.34  , -0.2  , 0.014 , 0.017 };
        Double_t step[3]   = {0.001, 0.001, 0.001};   //step 0 li rende costanti
     
        my_gMinuit->mnparm(0, "fA", vstart[0], step[0], 0., 1., ierflg);
        my_gMinuit->mnparm(1, "mA", vstart[1], step[1], 0., 0., ierflg);
      //   my_gMinuit->mnparm(2, "mB", vstart[2], step[2], 0., 0., ierflg);
        my_gMinuit->mnparm(2, "sA", vstart[2], step[2], 0., 0., ierflg);
      //   my_gMinuit->mnparm(4, "sB", vstart[4], step[4], 0., 0., ierflg);
        arglist[0] = 500.;//500;
        arglist[1] = 0.1;


        my_gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
      
        // Print results
        Double_t amin,edm,errdef;
        Int_t nvpar,nparx,icstat;
        my_gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);     
        my_gMinuit->mnprin(3,amin);
   
        Double_t val0, val1, val2, error0,error1,error2, bnd1, bnd2;
        Int_t ivar=0;
        TString chnam;
        
        my_gMinuit->mnpout(0, chnam, val0, error0, bnd1, bnd2, ivar);
        my_gMinuit->mnpout(1, chnam, val1, error1, bnd1, bnd2, ivar);
        my_gMinuit->mnpout(2, chnam, val2, error2, bnd1, bnd2, ivar);
      //   my_gMinuit->mnpout(3, chnam, val3, error3, bnd1, bnd2, ivar);
      //   my_gMinuit->mnpout(4, chnam, val4, error4, bnd1, bnd2, ivar);

      g_resolution->SetPoint(h, histo_diff_mc_time_bincenters[h], val2);
      g_resolution->SetPointError(h, 0.0, error2);
    };
   
//    Double_t matrix[5][5];
//    Int_t n=5;
//    my_gMinuit->mnemat (&matrix[0][0],5); // Accessing Covariance Matrix

//    // Write on .txt file covariance Matrix 
//    FILE *file_cov = fopen ("_txt/covariance_matrix.txt", "w");
//    Int_t icol=0;
//    Int_t irow=0;
//    if (n==1)
//      {
//        fprintf (file_cov, "%e \n", matrix[0][0]); 
//      }
//    else
//      {
//        for (irow=0; irow<n; irow++)
// 	 {
// 	   for (icol=0; icol<n; icol++)
// 	     { 
//                fprintf (file_cov, "%e ", matrix[irow][icol]); 
// 	     };
// 	   fprintf (file_cov, "\n");
// 	 };
//        fprintf (file_cov, "\n");   
//      }
//    fclose (file_cov);
   

   //here we should fit each of the 10 histograms to check the resolution as a function of the proper time -->

   //check the acceptance as a function of the proper decay time 
   // --> the empirical function will be 100 points with bins of proper time as x and the ratio between the proper time mc sim and meas
   // --> almeno io ho capito così oggi in classe :)

   TGraphErrors *acceptance = new TGraphErrors();

   for(int bin = 1; bin < histo_mc_time->GetNbinsX()+1; bin ++)
   {
      // double ratio = histo_mc_time_true->GetBinContent(bin)/histo_mc_time->GetBinContent(bin);
      double ratio = histo_mc_time->GetBinContent(bin)/histo_mc_time_true->GetBinContent(bin);

      //is is ok to use the standard error propagation ? --> how about the errors on x ? 

      // ratio = A/B

      double A = histo_mc_time_true->GetBinContent(bin);
      double sigmaA = histo_mc_time_true->GetBinError(bin);
      double B = histo_mc_time->GetBinContent(bin);
      double sigmaB = histo_mc_time->GetBinError(bin);

      double ratio_error = ratio * TMath::Sqrt(std::pow(sigmaA/A,2) + std::pow(sigmaB/B,2)); //--> covariance ? 

      //cout << ratio << " " << ratio_error << endl;

      //Error on x: bin width/2, otherwise you get a doubled error

      acceptance->SetPoint(bin,histo_mc_time_true->GetBinCenter(bin),ratio);
      acceptance->SetPointError(bin,histo_mc_time_true->GetBinWidth(bin)/2.0,ratio_error);

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

   TCanvas *c1 = new TCanvas("c1", "Acceptance", 800, 600);

   acceptance->SetTitle("Acceptance;Time;Ratio");
   acceptance->SetMarkerStyle(20);

   acceptance->Draw("AP E1");  // axes + points + error bars

   c1->Update();
    
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

   acceptance-> Write("acceptance_graph");
       g_resolution->SetTitle(
        "Resolution vs Proper Time;"
        "Proper Time;"
        "Resolution"
    );

   g_resolution->SetMarkerStyle(20);

   g_resolution->Draw("AP E1");

   g_resolution->Write("resolution_graph");

   histo_file->Close();

}
