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
Double_t TAU_BACKGROUND = 3.40304;

const double XMIN = 0.5;
const double XMAX = 10.0;

Double_t MAX_ACCEPTANCE = 0.115673;

Double_t ACCEPTANCE_PAR_FIXED[7] = {0.58397, -1.77281, 1.75128, 1.41352, -0.825079,  0.969095, -1.88183};

using namespace std;

std::vector<double> xvar;

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

//====================================================
// SIGNAL + BACKGROUND NORMALIZATION (FAST VERSION)
//====================================================

// ---------------------------------------------------
// GLOBAL NORMALIZATION TABLES
// ---------------------------------------------------

std::vector<double> tau_grid;
std::vector<double> signal_norm_grid;

double background_norm = 0.0;

// ---------------------------------------------------
// PURE SIGNAL PDF
// (without fsig mixing)
// ---------------------------------------------------

Double_t signal_pdf(Double_t x, Double_t tau)
{
    Double_t fA = F_FIXED;
    Double_t fB = 1.0 - F_FIXED;

    Double_t sigmaA = S_A_FIXED;
    Double_t sigmaB = S_B_FIXED;

    Double_t signalA =
        fA *
        (1.0 / (2.0 * tau) *
         TMath::Exp((sigmaA * sigmaA - 2.0 * tau * x) /
                    (2.0 * tau * tau)) *
         (TMath::Erf((tau * x - sigmaA * sigmaA) /
                     (TMath::Sqrt(2.0) * tau * sigmaA)) +
          1.0));

    Double_t signalB =
        fB *
        (1.0 / (2.0 * tau) *
         TMath::Exp((sigmaB * sigmaB - 2.0 * tau * x) /
                    (2.0 * tau * tau)) *
         (TMath::Erf((tau * x - sigmaB * sigmaB) /
                     (TMath::Sqrt(2.0) * tau * sigmaB)) +
          1.0));

    Double_t val = (signalA + signalB);

    val *= acceptance_func_mod(&x, ACCEPTANCE_PAR_FIXED);

    return val;
}

// ---------------------------------------------------
// PURE BACKGROUND PDF
// ---------------------------------------------------

Double_t background_pdf(Double_t x)
{
    Double_t tauB = TAU_BACKGROUND;

    Double_t fA = F_FIXED;
    Double_t fB = 1.0 - F_FIXED;

    Double_t sigmaA = S_A_FIXED;
    Double_t sigmaB = S_B_FIXED;

    Double_t bkgA =
        fA *
        (1.0 / (2.0 * tauB) *
         TMath::Exp((sigmaA * sigmaA - 2.0 * tauB * x) /
                    (2.0 * tauB * tauB)) *
         (TMath::Erf((tauB * x - sigmaA * sigmaA) /
                     (TMath::Sqrt(2.0) * tauB * sigmaA)) +
          1.0));

    Double_t bkgB =
        fB *
        (1.0 / (2.0 * tauB) *
         TMath::Exp((sigmaB * sigmaB - 2.0 * tauB * x) /
                    (2.0 * tauB * tauB)) *
         (TMath::Erf((tauB * x - sigmaB * sigmaB) /
                     (TMath::Sqrt(2.0) * tauB * sigmaB)) +
          1.0));

    Double_t val = (bkgA + bkgB);

    val *= acceptance_func_mod(&x, ACCEPTANCE_PAR_FIXED);

    return val;
}

// ---------------------------------------------------
// COMPUTE SIGNAL NORMALIZATION
// ---------------------------------------------------

Double_t computeSignalNorm(Double_t tau)
{
    Double_t dx = 0.001;
    Double_t sum = 0.0;

    for (double x = XMIN; x < XMAX; x += dx)
    {
        sum += signal_pdf(x, tau) * dx;
    }

    return sum;
}

// ---------------------------------------------------
// COMPUTE BACKGROUND NORMALIZATION
// ---------------------------------------------------

Double_t computeBackgroundNorm()
{
    Double_t dx = 0.001;
    Double_t sum = 0.0;

    for (double x = XMIN; x < XMAX; x += dx)
    {
        sum += background_pdf(x) * dx;
    }

    return sum;
}

// ---------------------------------------------------
// PRECOMPUTE NORMALIZATION TABLE
// ---------------------------------------------------

void initNormGrid()
{
    tau_grid.clear();
    signal_norm_grid.clear();

    for (double tau = 0.1; tau < 10.0; tau += 0.01)
    {
        tau_grid.push_back(tau);
        signal_norm_grid.push_back(computeSignalNorm(tau));
    }

    background_norm = computeBackgroundNorm();

    std::cout << "Normalization grids initialized\n";
}

// ---------------------------------------------------
// INTERPOLATED SIGNAL NORMALIZATION
// ---------------------------------------------------

double getSignalNorm(double tau)
{
    int i = int((tau - tau_grid[0]) / 0.01);

    if (i < 0)
        i = 0;

    if (i >= (int)tau_grid.size() - 1)
        i = tau_grid.size() - 2;

    double t1 = tau_grid[i];
    double t2 = tau_grid[i + 1];

    double n1 = signal_norm_grid[i];
    double n2 = signal_norm_grid[i + 1];

    return n1 + (tau - t1) * (n2 - n1) / (t2 - t1);
}

// ---------------------------------------------------
// FULL NORMALIZED PDF
// ---------------------------------------------------

Double_t pdf_norm(Double_t x, Double_t *par)
{
    Double_t tau = par[0];
    Double_t fsig = par[1];

    Double_t sig =
        signal_pdf(x, tau) / getSignalNorm(tau);

    Double_t bkg =
        background_pdf(x) / background_norm;

    return fsig * sig + (1.0 - fsig) * bkg;
}

// ---------------------------------------------------
// FCN
// ---------------------------------------------------

void fcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
    Double_t Like = 0.0;

    for (size_t i = 0; i < xvar.size(); ++i)
    {
        double p = pdf_norm(xvar[i], par);

        if (p <= 0)
        {
            f = 1e30;
            return;
        }

        Like += log(p);
    }

    f = -2.0 * Like;
}

// ---------------------------------------------------
// PROJECTION FOR HISTOGRAM PLOT
// ---------------------------------------------------

Double_t pdf_proj(Double_t *x,
                  Double_t *par,
                  Double_t N,
                  Double_t bin_width)
{
    return N * pdf_norm(x[0], par) * bin_width;
}
// 


//  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
//   1  tau          3.40304e+00   1.22615e-02   3.83827e-04   4.36550e-02

//========== FIT RESULT ==========
//Parameter : tau
//tau_fit   = 1396.27 fs
//error     = 5.03091 fs
//================================

//p0                        =      54150.5   +/-   578.437     
//p1                        =      1.86482   +/-   2.98666e-05 
//p2                        =   0.00521984   +/-   4.64803e-05 
//p3                        =      13.5657   +/-   0.460617    
//chi2 = 55.3241
//ndof = 56
//prob = 0.500404



void lifetimeANA::Loop()
{
    // ---------------------------------------------
    // INITIALIZATION
    // ---------------------------------------------

    initNormGrid();

    xvar.clear();
    h_time->Reset();

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;

    // ---------------------------------------------
    // EVENT LOOP
    // ---------------------------------------------

    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);

        if (ientry < 0)
            break;

        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        // -------------------------
        // DATA
        // -------------------------

        if (id == 1)
        {
            if (M0_MKpi > 1.86482 - 2 * 0.00540974 &&
                M0_MKpi < 1.86482 + 2 * 0.00540974)
            {
                double t = M0_time / 410.3e-15;

                if (t >= XMIN && t <= XMAX)
                {
                    xvar.push_back(t);
                    h_time->Fill(t);
                }
            }
        }

        // -------------------------
        // MC
        // -------------------------

        if (id == 13)
        {
            // optional MC studies
        }

    } // end event loop

    // ---------------------------------------------
    // FIT
    // ---------------------------------------------

    const int nparam = 2;

    int n_fit = xvar.size();

    cout << "DIMENSION " << n_fit << endl;

    TMinuit *my_gMinuit = new TMinuit(nparam);

    my_gMinuit->SetFCN(fcn);

    Double_t arglist[2];
    Int_t ierflg = 0;

    Double_t vstart[nparam] = {1.0, 0.5};
    Double_t step[nparam]   = {0.01, 0.01};

    my_gMinuit->mnparm(
        0,
        "tau",
        vstart[0],
        step[0],
        0.1,
        10.0,
        ierflg);

    my_gMinuit->mnparm(
        1,
        "signal_fraction",
        vstart[1],
        step[1],
        0.0,
        1.0,
        ierflg);

    arglist[0] = 500;
    arglist[1] = 0.1;

    my_gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // ---------------------------------------------
    // PRINT RESULTS
    // ---------------------------------------------

    Double_t amin, edm, errdef;
    Int_t nvpar, nparx, icstat;

    my_gMinuit->mnstat(
        amin,
        edm,
        errdef,
        nvpar,
        nparx,
        icstat);

    my_gMinuit->mnprin(3, amin);

    Double_t tau_fit;
    Double_t etau_fit;

    Double_t signal_fraction_fit;
    Double_t e_signal_fraction_fit;

    Double_t bnd1, bnd2;

    Int_t ivar = 0;

    TString chnam;

    my_gMinuit->mnpout(
        0,
        chnam,
        tau_fit,
        etau_fit,
        bnd1,
        bnd2,
        ivar);

    my_gMinuit->mnpout(
        1,
        chnam,
        signal_fraction_fit,
        e_signal_fraction_fit,
        bnd1,
        bnd2,
        ivar);

    std::cout << "\n========== FIT RESULT ==========\n";

    std::cout << "tau_fit = "
              << tau_fit * 410.3
              << " fs\n";

    std::cout << "error = "
              << etau_fit * 410.3
              << " fs\n";

    std::cout << "signal fraction = "
              << signal_fraction_fit
              << "\n";

    std::cout << "================================\n";

    // ---------------------------------------------
    // BUILD FIT CURVES
    // ---------------------------------------------

    std::vector<double> x_points;

    std::vector<double> y_total;
    std::vector<double> y_signal;
    std::vector<double> y_background;

    Double_t pars[nparam] =
    {
        tau_fit,
        signal_fraction_fit
    };

    for (int xi = 0; xi < 10000; xi++)
    {
        double xmin = 0.0;
        double xmax = 10.0;

        double x =
            xmin +
            xi * (xmax - xmin) / 10000.0;

        x_points.push_back(x);

        // normalized signal/background PDFs

        double sig =
            signal_pdf(x, tau_fit)
            /
            getSignalNorm(tau_fit);

        double bkg =
            background_pdf(x)
            /
            background_norm;

        // expected histogram counts

        double ysig =
            n_fit *
            signal_fraction_fit *
            sig *
            h_time->GetBinWidth(1);

        double ybkg =
            n_fit *
            (1.0 - signal_fraction_fit) *
            bkg *
            h_time->GetBinWidth(1);

        double ytot = ysig + ybkg;

        y_signal.push_back(ysig);
        y_background.push_back(ybkg);
        y_total.push_back(ytot);
    }

    // ---------------------------------------------
    // CREATE GRAPHS
    // ---------------------------------------------

    TGraph *fit_total =
        new TGraph(
            10000,
            x_points.data(),
            y_total.data());

    TGraph *fit_signal =
        new TGraph(
            10000,
            x_points.data(),
            y_signal.data());

    TGraph *fit_background =
        new TGraph(
            10000,
            x_points.data(),
            y_background.data());

    // ---------------------------------------------
    // STYLE
    // ---------------------------------------------

    fit_total->SetLineColor(kRed);
    fit_total->SetLineWidth(2);

    fit_signal->SetLineColor(kBlue);
    fit_signal->SetLineWidth(2);
    fit_signal->SetLineStyle(2);

    fit_background->SetLineColor(kGreen + 2);
    fit_background->SetLineWidth(2);
    fit_background->SetLineStyle(3);

    // ---------------------------------------------
    // CANVAS
    // ---------------------------------------------

    TCanvas *c_fit =
        new TCanvas(
            "c_fit",
            "fit",
            800,
            600);

    c_fit->cd();

    const float split = 0.3;

    // ---------------------------------------------
    // TOP PAD
    // ---------------------------------------------

    TPad *pad_fit =
        new TPad(
            "pad_fit",
            "",
            0.,
            split,
            1.,
            1.);

    pad_fit->SetBottomMargin(0.02);

    pad_fit->Draw();

    pad_fit->cd();

    gStyle->SetOptStat(0);

    gPad->SetTicks(1,1);

    h_time->GetXaxis()->SetRangeUser(0., 10.0);

    h_time->Draw("E");

    fit_total->Draw("L SAME");
    fit_signal->Draw("L SAME");
    fit_background->Draw("L SAME");

    // ---------------------------------------------
    // LEGEND
    // ---------------------------------------------

    TLegend *leg =
        new TLegend(
            0.60,
            0.65,
            0.88,
            0.88);

    leg->AddEntry(h_time, "Data", "lep");
    leg->AddEntry(fit_total, "Total fit", "l");
    leg->AddEntry(fit_signal, "Signal", "l");
    leg->AddEntry(fit_background, "Background", "l");

    leg->SetBorderSize(0);

    leg->Draw();

    c_fit->cd();

    // ---------------------------------------------
    // BOTTOM PAD
    // ---------------------------------------------

    TPad *pad_res =
        new TPad(
            "pad_res",
            "",
            0.,
            0.,
            1.,
            split);

    pad_res->SetTopMargin(0.02);
    pad_res->SetBottomMargin(0.25);

    pad_res->Draw();

    pad_res->cd();

    int n = h_time->GetNbinsX();

    TGraphErrors *residuals =
        new TGraphErrors(n);

    residuals->SetMarkerStyle(7);

    residuals->GetXaxis()->SetTitle("t / 410.3 fs");
    residuals->GetYaxis()->SetTitle("pull");

    for (int i = 0; i < n; i++)
    {
        Double_t x =
            h_time->GetBinCenter(i + 1);

        Double_t y =
            h_time->GetBinContent(i + 1);

        Double_t y_fit =
            pdf_proj(
                &x,
                pars,
                n_fit,
                h_time->GetBinWidth(1));

        Double_t ey =
            h_time->GetBinError(i + 1);

        Double_t pull = 0.0;

        if (ey > 0)
            pull = (y - y_fit) / ey;

        residuals->SetPoint(i, x, pull);

        residuals->SetPointError(i, 0., 1.0);
    }

    residuals->Draw("AP");

    residuals->GetXaxis()->SetLimits(0.0, 10.0);

    TLine *line_zero =
        new TLine(0.0, 0., 10., 0.);

    line_zero->SetLineStyle(2);

    line_zero->Draw("same");

    // ---------------------------------------------
    // FINALIZE
    // ---------------------------------------------

    c_fit->Update();
    c_fit->Draw();

    TFile *f =
        new TFile(
            "outfile_fit_DATA_SIGNAL_SINGLE_NORM.root",
            "RECREATE");

    c_fit->Write();

    f->Close();
}