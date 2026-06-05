#define analysis_cxx
#include "analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

void analysis::Loop()
{

  
  TH1D *histo_D_M         = new TH1D("histo_D_M","",100,1.640,2.100);
  TH1D *histo_D_M_S       = new TH1D("histo_D_M_S","",100,1.640,2.100);
  TH1D *histo_D_time      = new TH1D("histo_D_time","",100,0,5e-12);
  TH1D *histo_D_time_S    = new TH1D("histo_D_time_S","",100,0,5e-12);
  TH1D *histo_D_pt        = new TH1D("histo_D_pt","",100,0,10);
  TH1D *histo_D_pt_S      = new TH1D("histo_D_pt_S","",100,0,10);
  TH1D *histo_h3_p_back   = new TH1D("histo_h3_p_back","",1000,0,200);
  TH1D *histo_h3_mu_p   = new TH1D("histo_h3_mu_p","",1000,0.5,0.65);

      Double_t count_bkg_34 = 0;
      Double_t count_bkg_34_tot = 0;
      Double_t count_bkg_41 = 0;
      Double_t count_bkg_41_tot = 0;
      Double_t count_bkg_42 = 0;
      Double_t count_bkg_42_tot = 0;
      Double_t count_bkg_43 = 0;
      Double_t count_bkg_43_tot = 0;


  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);
     nbytes += nb;
     // if (Cut(ientry) < 0) continue;
     //std::cout << M0_MKK << std::endl;
     histo_h3_p_back->Fill(h3_p,weight);

     // from fit to signal MC
     //Mean    =  1.77699   +/-   1.98392e-05 
     //Sigma   =   0.00581298   +/-   1.49607e-05  	 (limited)
     Double_t m_tau_mc;
     Double_t m_sigma_tau_mc;
     m_tau_mc= 1.77699;
     m_sigma_tau_mc=  0.00581298;

     Double_t m_mu_mc;
     m_mu_mc= 0.1075;

     
  //    if (id == 0){ // DATA
  //     if (TMath::Abs(D_M-m_tau_mc)>3*m_sigma_tau_mc)
  //      {
	 
	//  if (h3_MuonID==1){
	//    //if (1){     
	//    histo_D_M->Fill(D_M) ;
	//    histo_D_time ->Fill(D_time);
	//    histo_D_pt     ->  Fill(D_pt);
	//    //std::cout << D_M << std::endl;
	//  }; 
	//  if (h3_MuonID==1 && D_time>0.25e-12 && D_pt>2.5  && D_FDt > 1e-3){
	//    //if (h3_MuonID==1 && D_time>0.25e-12 && D_pt>2.5 ){
	//    //if ( D_time>2.5e-12){
	//    histo_D_M_S->Fill(D_M) ;
	//    histo_D_time_S ->Fill(D_time);
	//    histo_D_pt_S     ->  Fill(D_pt);
	//  };
	 
  //      }; //Mass blinding
  //    }

    if (id > 1){ // MC
      // if (TMath::Abs(D_M-m_tau_mc)>3*m_sigma_tau_mc)
      //  {
	 
	 if (h3_MuonID==1){
	   //if (1){     
	   histo_D_M->Fill(D_M) ;
	   histo_D_time ->Fill(D_time);
	   histo_D_pt     ->  Fill(D_pt);
	   //std::cout << D_M << std::endl;
	 }; 
	 if (h3_MuonID==1 && D_time>0.25e-12 && D_pt>2.5  && D_FDt > 1e-3){
	   //if (h3_MuonID==1 && D_time>0.25e-12 && D_pt>2.5 ){
	   //if ( D_time>2.5e-12){
	   histo_D_M_S->Fill(D_M) ;
	   histo_D_time_S ->Fill(D_time);
	   histo_D_pt_S     ->  Fill(D_pt);
	 };
	 
  //  if (h3_MuonID==1 && (id == 44 || id == 42))
  //  if (id == 44 || id == 42)

  //  {

    Double_t mother_p = TMath::Sqrt(D_px*D_px+D_py*D_py+D_pz*D_pz);
    Double_t daughter_p = TMath::Sqrt(h3_px*h3_px+h3_py*h3_py+h3_pz*h3_pz);

    ROOT::Math::PxPyPzEVector mother(D_px,D_py,D_pz,TMath::Sqrt(D_p*D_p+m_tau_mc*m_tau_mc));
    ROOT::Math::PxPyPzEVector daughter(h3_px,h3_py,h3_pz,TMath::Sqrt(h3_p*h3_p+m_mu_mc*m_mu_mc));

    ROOT::Math::Boost boost(mother.BoostToCM());

    auto motherRest   = boost(mother);
    auto daughterRest = boost(daughter);

    double pCMS = daughterRest.P();
    double motherMass = mother.M();
    double daughterMass = daughter.M();


    // std::cout << "daughter P = " << daughter.P()
    //       << " mother M CM = " << motherRest.M()
    //       << " daughter P CM = " << pCMS
    //       << std::endl;


      histo_h3_mu_p->Fill(pCMS);

      if (id == 34){
        count_bkg_34_tot += 1;
        if((pCMS >0.57 && pCMS < 0.61))
        {
          count_bkg_34 += 1;
        }
      }
         
      if (id == 41){
        count_bkg_41_tot += 1;
        if((pCMS >0.57 && pCMS < 0.61))
        {
          count_bkg_41 += 1;
        }
      }
      if (id == 42){
        count_bkg_42_tot += 1;
        if((pCMS >0.57 && pCMS < 0.61))
        {
          count_bkg_42 += 1;
        }
      }
      if (id == 43){
        count_bkg_43_tot += 1;
        if((pCMS >0.57 && pCMS < 0.61))
        {
          count_bkg_43 += 1;
        }
      }

  //  }

      //  }; //Mass blinding
     }
       
     
   }; //End Loop Events
  cout << "survival ratio 34: " << count_bkg_34/count_bkg_34_tot << endl;
  cout << "survival ratio 41: " << count_bkg_41/count_bkg_41_tot << endl;
  cout << "survival ratio 42: " << count_bkg_42/count_bkg_42_tot << endl;
  cout << "survival ratio 43: " << count_bkg_43/count_bkg_43_tot << endl;

   
   //Create Canvas
    TCanvas *c_histo_D_M = new TCanvas("c_histo_D_M","canvas histo",500,500);
    c_histo_D_M->cd();
    histo_D_M->SetMinimum(0);
    histo_D_M->Draw();
    histo_D_M_S->SetLineColor(2);
    histo_D_M_S->Draw("same");

    TCanvas *c_histo_D_time = new TCanvas("c_histo_D_time","canvas histo",500,500);
    c_histo_D_time->cd();
    histo_D_time->SetMinimum(0);
    histo_D_time->Draw();
    histo_D_time_S->SetLineColor(2);
    histo_D_time_S->Draw("same");

    TCanvas *c_histo_D_pt = new TCanvas("c_histo_D_pt","canvas histo",500,500);
    c_histo_D_pt->cd();
    histo_D_pt->SetMinimum(0);
    histo_D_pt->Draw();
    histo_D_pt_S->SetLineColor(2);
    histo_D_pt_S->Draw("same");
    
    // Save canvas in .pdf and .epr format
    c_histo_D_M->Print("./_fig/c_histo.pdf");
  
  //Create a new file to store histograms
   TFile *histo_file = new TFile("./_root/histo_file_MC.root","RECREATE","put a title");
  //TFile *histo_file = new TFile("./_root/histo_p_back_DPLUS_PhiPi.root","RECREATE","put a title");
   histo_file->cd();
   histo_h3_p_back->Write();
   //histo_file->Write();
   histo_D_M ->Write();
   histo_D_M_S ->Write();
   histo_D_time->Write();
   histo_D_time_S->Write();
   histo_D_pt     ->Write();
   histo_D_pt_S   ->Write(); 
   histo_h3_mu_p->Write();

   histo_file->Close();
     
}
