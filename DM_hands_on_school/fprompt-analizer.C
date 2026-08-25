#include <TMath.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <string>
#include <TTree.h>


void FpromptCalc(int wshort, int wlong) {
    
    
    TH2D *hFpromptE = new TH2D("hFpromptE","2D Histogram",1000,0.,10000.,110,0.4,1.);
    std::ifstream f("analysis/selection.runs",std::ifstream::in);

  while (true)
    {
      int run;
      
      f>> run;
      if ( f.eof() ) break;

    TFile *file = new TFile(Form("trees_v2/mytreefile_%d.root",run), "READ");
    if ( file->IsZombie() ) { printf("RedFile does not exist \n"); continue; }
    TTree *trace = (TTree*)file->Get("treetraces");
    
    std::vector<float> *energy  = nullptr;
	  int sa_start,sa_end;
	  
    trace->SetBranchAddress("energy",&energy);
	  trace->SetBranchAddress("sa_start",&sa_start);
	  trace->SetBranchAddress("sa_end",&sa_end);

    

    for (int ientry=0;ientry<trace->GetEntries();ientry++)
    {
      trace->GetEntry(ientry);
      if ( ientry%1000==0 ) printf("Run %d - Entry %d \n",run,ientry);

	    double qshort = 0;
      double qlong = 0;
      double totalenergy = 0;
      double energytoplot = 0;
      double fprompt = 0;

      for (int sa=0;sa<wshort;sa++)
    {
      qshort+=energy->at(sa);
    } 
    for (int sa=0;sa<wlong;sa++)
    {
      qlong+=energy->at(sa); 
    } 
    for (int sa=0;sa<(int) energy->size();sa++)
    {
      totalenergy+=energy->at(sa);
    } 
      fprompt = qshort/qlong;
      energytoplot = -totalenergy;
		  hFpromptE->Fill(energytoplot,fprompt);
	}
    trace->ResetBranchAddresses();
  }

  TCanvas *c1 = new TCanvas ("c1","c1",1920,1080);
  hFpromptE-> Draw("COLZ");

  TH1D *hfprompt=(TH1D*) hFpromptE->ProjectionY("hfprompt");
  
  TCanvas *c2 = new TCanvas ("c2","c2",1920,1080);
  hfprompt->Draw();
}
