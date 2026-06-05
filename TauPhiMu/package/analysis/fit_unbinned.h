void fit_unbinned(int nparam, std::vector<double> xvar, Double_t vstart[nparam], Double_t step[nparam], Double_t low[nparam], Double_t high[nparam], std::string par_names[nparam])
{
      const int nparam = 3;
      int n_fit = xvar.size();
      //cout << "DIMENSION " << n_fit << endl;

      TMinuit *my_gMinuit = new TMinuit(nparam);  //initialize TMinuit with a maximum of 5 params
      //gMinuit->SetPrintLevel(-1);
      my_gMinuit->SetFCN(fcn);      // set the FCN
        
      Double_t arglist[2];
      Int_t ierflg = 0;  // Error return code: 0 if the command was correctly executed, >0 otherwise. 
     
      for(int par = 0; par < nparam; par ++)
      {
        my_gMinuit->mnparm(par, par_names[par].c_str(), vstart[par], step[par], low[par], high[par], ierflg);
      }

      arglist[0] = 5000.;//500;
      arglist[1] = 0.1;

      my_gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
      
      // Print results
      Double_t amin,edm,errdef;
      Int_t nvpar,nparx,icstat;
      my_gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);     
      my_gMinuit->mnprin(3,amin);

      Double_t pars[nparam];
      Double_t pars_errors[nparam];

      Int_t ivar=0;
      TString chnam;
      
      for(int par = 0; par < nparam; par ++)
      {
        my_gMinuit->mnpout(par, chnam, pars[par], pars_errors[par], bnd1, bnd2, ivar);
      }

      std::cout << "\n========== FIT RESULT ==========\n";
      std::cout << "Parameter : " << chnam << std::endl;
      for(int par = 0; par < nparam; par ++)
      {
        std::cout << Form("%s   = ",par_names[0]) << pars[par] << std::endl;
      }
      std::cout << "================================\n";

      std::vector<double> x_points;
      std::vector<double> y_points;

      

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
      pad_fit->SetBottomMargin(0.0);
      pad_fit->Draw();
      pad_fit->cd();

      gStyle->SetOptStat(0);
      gPad->SetTicks(1,1);

      h_time->GetXaxis()->SetRangeUser(0., 10.0);
      h_time->SetLineColor(kBlack);
      h_time->SetMarkerStyle(8);
      h_time->SetMarkerSize(0.5);
      h_time ->GetYaxis() -> SetTitleOffset(0.85);
      h_time->GetYaxis()->SetTitle("counts / 0.1 [fs / 410.3]");
      h_time->GetYaxis()->SetTitleSize(0.06);
      h_time->GetYaxis()->SetLabelSize(0.05);
      h_time->Draw("PE");

      // overlay fit
      fit_function_plot->SetLineColor(kRed);
      fit_function_plot->SetMarkerColor(kRed);
      fit_function_plot->Draw("L SAME");

      // IMPORTANT: go back to canvas
      c_fit->cd();

      // -------------------- BOTTOM PAD (PULLS) --------------------
      TPad *pad_pull = new TPad("pad_pull", "", 0., 0., 1., split);
      pad_pull->SetTopMargin(0.0);
      pad_pull->SetBottomMargin(0.3);
      pad_pull->Draw();
      pad_pull->cd();
        
      int n = h_time->GetNbinsX();

      TH1D *h_pull = new TH1D("h_pull", "", 20, -5, 5);      
        
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

          if( x >= XMIN)
          {
      
          double y_fit = pdf_proj(&x, pars, n_fit, h_time->GetBinWidth(1));
      
          // avoid division by zero
          double pull = 0;
          if (ey > 0) pull = (y - y_fit) / ey;
      
          pulls->SetPoint(i, x, pull);
          pulls->SetPointError(i, 0., 1.0);
          h_pull->Fill(pull);
          }

      }
      
      pulls -> SetMarkerStyle(8);
      pulls -> SetTitle("");
      pulls -> SetMarkerSize(0.4);
      pulls -> GetYaxis() -> SetTitleOffset(0.3);
      pulls -> GetXaxis()-> SetTitle("t [fs / 410.3]");
      pulls -> GetYaxis()-> SetTitle("pulls");
      pulls -> GetYaxis()-> SetTitleSize(0.13);
      pulls -> GetXaxis()-> SetTitleSize(0.13);
      pulls -> GetYaxis()-> SetLabelSize(0.11);
      pulls -> GetXaxis()-> SetLabelSize(0.11);
      pulls->Draw("APE");
      
      // baseline at 0
      TLine *line_zero = new TLine(0.0, 0., 10., 0.);
      line_zero->SetLineStyle(2);
      line_zero->Draw("same");
      // final update

      c_fit -> Update();
      c_fit -> SaveAs("PLOT_REPORT/FIT_BKG_TIME.pdf");
}