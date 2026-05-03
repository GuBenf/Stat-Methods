
#ifndef LAB_STYLE_H
#define LAB_STYLE_H
#include <Riostream.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH2.h>
#include <TNamed.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveStats.h>



/******************************************************************************************
 * Style                                                                                  *
 ******************************************************************************************/
void labstyle()
{
	TStyle *labstyle = new TStyle("labstyle","The real perfect style");

	// canvas
	labstyle->SetCanvasColor(0);
	labstyle->SetCanvasBorderSize(10);
	labstyle->SetCanvasBorderMode(0);
	labstyle->SetCanvasDefW(600);
	labstyle->SetCanvasDefH(600);

	// pads
	labstyle->SetPadColor(0);
	labstyle->SetPadBorderSize(10);
	labstyle->SetPadBorderMode(0);
	labstyle->SetPadLeftMargin(.18);
	labstyle->SetPadRightMargin(.05);
	labstyle->SetPadBottomMargin(.12);
	labstyle->SetPadTopMargin(.1);
	//labstyle->SetPadGridX(1);
	//labstyle->SetPadGridY(1);
	labstyle->SetPadTickX(0);
	labstyle->SetPadTickY(0);

	// frame
	labstyle->SetFrameBorderMode(0);
	labstyle->SetFrameBorderSize(10);
	labstyle->SetFrameFillStyle(0);
	labstyle->SetFrameFillColor(0);
	labstyle->SetFrameLineColor(1);
	labstyle->SetFrameLineStyle(0);
	labstyle->SetFrameLineWidth(1);

	// histogram
	labstyle->SetHistFillColor(0);
	labstyle->SetHistFillStyle(1001);// solid
	labstyle->SetHistLineColor(1);
	labstyle->SetHistLineStyle(0);
	labstyle->SetHistLineWidth(1);
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,16,0)
	labstyle->SetHistMinimumZero();
#endif
	labstyle->SetOptStat(0);
	labstyle->SetOptFit(0);
	labstyle->SetStatColor(0);
	labstyle->SetStatBorderSize(1);
	labstyle->SetStatFontSize(.05);

	// graph
	labstyle->SetEndErrorSize(0);
   labstyle->SetErrorX(0.5);

	// marker
	labstyle->SetMarkerStyle(20);
	labstyle->SetMarkerColor(kBlack);
	labstyle->SetMarkerSize(0.6);

	// title 
	labstyle->SetOptTitle(1);
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,16,0)
	labstyle->SetTitleAlign(33);
	labstyle->SetTitleX(.95);
#else
	labstyle->SetTitleX(.15);
#endif
	labstyle->SetTitleFillColor(0);
	labstyle->SetTitleBorderSize(0);
	labstyle->SetTitleStyle(0);

	// axes
	labstyle->SetNdivisions(505,"X");
	labstyle->SetNdivisions(505,"Y");
	
	labstyle->SetTitleSize(.05,"X");//.055
	labstyle->SetTitleOffset(1.,"X");//1.2,0.9
	labstyle->SetLabelOffset(0.003,"X");
	labstyle->SetLabelSize(.05,"X");
	labstyle->SetLabelFont(42,"X");

	labstyle->SetTitleSize(.05,"Y");//.055
	labstyle->SetTitleOffset(1.8,"Y");
	labstyle->SetLabelOffset(0.008,"Y");
	labstyle->SetLabelSize(.05,"Y");
	labstyle->SetLabelFont(42,"Y");

	labstyle->SetStripDecimals(kFALSE);

	labstyle->SetTitleSize(0.05,"Z");
	labstyle->SetTitleOffset(1.800,"Z");
	labstyle->SetLabelOffset(0.008,"Z");
	labstyle->SetLabelSize(0.05,"Z");
	labstyle->SetLabelFont(42,"Z");
	
	// fonts
	labstyle->SetTextSize(.05);//.055
	labstyle->SetTextFont(42);
	labstyle->SetStatFont(42);
	labstyle->SetTitleFont(42,"");
	labstyle->SetTitleFont(42,"Z");
	labstyle->SetTitleFont(42,"X");
	labstyle->SetTitleFont(42,"Y");

	// function
	labstyle->SetFuncColor(kBlue);
	labstyle->SetFuncStyle(0);
	labstyle->SetFuncWidth(1);

	// legend
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,16,0)
	labstyle->SetLegendBorderSize(1);
#endif

	// palette
	labstyle->SetPalette(1);
	labstyle->SetNumberContours(20);
		
	// set labstyle as current style
	gROOT->SetStyle("labstyle");

	gSystem->ProcessEvents();

}



#endif // LAB_STYLE_C

