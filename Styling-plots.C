#include <TVector3.h>
#include <iostream>
#include <Math/RootFinder.h>
#include <Math/Functor.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCut.h>
#include <TArrow.h>
#include <TApplication.h>
#include <TTimer.h>
#include <TSystem.h>
#include <TEllipse.h>
#include <TArrow.h>

void stylehist1D(TH1D *h0)
{
  h0->SetTitleOffset(0.9,"X");
  h0->SetTitleOffset(1.,"Y");
  gStyle->SetTitleFontSize(0.08);
  h0->SetLineWidth(3);
  h0->GetXaxis()->SetLabelSize(0.05);
  h0->GetXaxis()->SetTitleSize(0.05);
  h0->GetYaxis()->SetLabelSize(0.05);
  h0->GetYaxis()->SetTitleSize(0.05);
  return;
}

void stylelegend(TLegend *leg)
{
  leg->SetTextSize(0.04);
  leg->SetMargin(0.05);
  leg->SetTextAlign(12);
}

void styletf1(TF1 *f1)
{
  f1->SetLineWidth(4);
  //f1->GetHistogram()->SetTitleOffset(1.,"X");
  //f1->GetHistogram()->SetTitleOffset(1.,"Y");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetLabelSize (0.05, "X");
  gStyle->SetLabelSize (0.05, "Y");
  gStyle->SetTitleSize (0.05, "X");
  gStyle->SetTitleSize (0.05, "Y");
  gStyle->SetTitleOffset (1., "X");
  gStyle->SetTitleOffset (1., "Y");
  //f1->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
  //f1->GetHistogram()->GetXaxis()->SetTitleSize(0.05);
  //f1->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
  //f1->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
}


void stylegraph1D(TGraph *t0)
{
  // t0->GetHistogram()->SetTitleOffset(0.8,"X");
  t0->GetYaxis()->SetTitleOffset(1.0);
  t0->GetXaxis()->SetTitleOffset(0.97);
  t0->GetXaxis()->SetLabelSize(0.05);
  t0->GetXaxis()->SetTitleSize(0.05);
  t0->GetYaxis()->SetLabelSize(0.05);
  t0->GetYaxis()->SetTitleSize(0.05);
  t0->SetMarkerSize(1.8);
  return;
}

void stylehist2d(TH2F *h0)
{
  h0->SetTitleOffset(1.,"X");
  h0->SetTitleOffset(1.,"Y");
  gStyle->SetTitleFontSize(0.08);
  h0->SetLineWidth(3);
  h0->GetXaxis()->SetLabelSize(0.05);
  h0->GetXaxis()->SetTitleSize(0.05);
  h0->GetYaxis()->SetLabelSize(0.05);
  h0->GetYaxis()->SetTitleSize(0.05);
  h0->SetTitleOffset(0.8,"Z");
  h0->GetZaxis()->SetLabelSize(0.05);
  h0->GetZaxis()->SetTitleSize(0.05);
  h0->SetMarkerSize(1.8);
  return;
}

void stylegrapherr1D(TGraphErrors *t0)
{
  t0->GetYaxis()->SetTitleOffset(1.);
  t0->GetXaxis()->SetTitleOffset(1.);
  t0->GetXaxis()->SetLabelSize(0.05);
  t0->GetXaxis()->SetTitleSize(0.05);
  t0->GetYaxis()->SetLabelSize(0.05);
  t0->GetYaxis()->SetTitleSize(0.05);
  t0->SetMarkerSize(1.8);
  return;
}