#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include "/home/barbaras/red-deconv/barb/Styling-plots.C"
#include "TMarker.h"
#include "TCanvas.h"
#include "TColor.h"
#include <vector>

using namespace std;

Double_t funcECf(Double_t *val, Double_t *par);
Double_t func(Double_t *val, Double_t *par);
Double_t ziegler(Double_t *val, Double_t *par);
Double_t fnfunc(Double_t *val, Double_t *par);
Double_t moliere(Double_t *val, Double_t *par);
Double_t lenzjensen(Double_t *val, Double_t *par);
Double_t pscienergy(Double_t *val, Double_t *par);
Double_t Fprompt_PSci_Cut(Double_t *val, Double_t *par);
Double_t Fprompt_mu(Double_t *val, Double_t *par);
Double_t Fprompt_sigma(Double_t *val, Double_t *par);
Double_t EJ276D_efficiency(Double_t *val, Double_t *par);
auto *leg2 = new TLegend(0.61,0.66,0.97,0.86);
auto *leg3 = new TLegend(0.61,0.66,0.85,0.86);
Double_t ESTAR_stoppingpower(Double_t *val, Double_t *par);
Double_t EJ276D_LightOutput(Double_t *val, Double_t *par);

Double_t fnfunc_mcmc(Double_t *x, Double_t *par);
double CalculateNie_LenzJensen(double Er_keV, double cbox, double beta);
double CalculateLogLikelihood2D(TH2F* h2_data, const std::vector<double>& Nie_cache, double g2_guess);
void MCMC_Fit_LenzJensen(TH2F* h2_to_fit);


TGraph *grst = nullptr;
TF1 *lightout = nullptr;


TLatex l;


int randomhelp()
{
    cout<<"int n = number of events"<<endl;
    cout<<"int opt = (0,1,2,3) = (no fluc, s2 quant fluc, er fluc + s2 quant, all fluc)"<<endl;
    cout<<"const char funcopt = (ziegler, moliere, lenz-jensen)"<<endl;
    cout<<"const char drawopt = (empty, same)"<<endl;
    cout<<"int coloropt = choose any color number"<<endl;
    return 0;
}

void random(int n, int opt, const char* funcopt, const char* drawopt, int coloropt)
{   
    gStyle->SetOptStat(0000000000);
    if (n <= 0)
    {
        cout << "Invalid number" << endl;
        return;
    }

    auto f = new TF1("f", funcECf, 0.0001, 15., 0);
    auto fer = new TF1("fer", func, 0.0001, 15., 1);

    auto cutfprompt = new TF1("cutfprompt", Fprompt_PSci_Cut, 0.0001, 3000., 0);
    auto mufprompt = new TF1("mufprompt", Fprompt_mu, 0.0001, 3000., 0);
    auto sigmafprompt = new TF1("sigmafprompt", Fprompt_sigma, 0.0001, 3000., 0);
    auto pscyenergy = new TF1("pscyenergy", pscienergy, 0.0001, 20., 1);
    auto pscieff = new TF1("pscieff", EJ276D_efficiency, 0.0001, 20., 0);
    
    pscieff->SetNpx(1000);

    //pscieff->Draw();



    auto s2func = new TF1("s2func", ziegler, 0.0001, 10., 3);
    auto s2moli = new TF1("s2moli", moliere, 0.0001, 10., 3);
    auto s2lenz = new TF1("s2lenz", lenzjensen, 0.0001, 10., 3);

    /*auto sigmaneutron = new TF1("sigman", "pol2", 0.0001, 3000.);
    sigmaneutron->SetParameters(0.0304238,-2.55435e-05,9.71383e-09);
    auto meanneutron = new TF1("meann", "[0]*(1.-exp(-[1]*x))+[2]", 0.0001, 3000.);
    meanneutron->SetParameters(5.08055e-02,1.37182e-03,7.32697e-01);
    auto c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
    meanneutron->Draw();
    */
    auto fn = new TF1("fn", fnfunc, 0., 20., 3);
    
    pscyenergy->SetParameter(0,1.);


    double theta4[3] = {810.,6.8*pow(10,3),17.88};
    //s2func->SetParameters(theta4);    
    s2func->FixParameter(0,810.);
    s2func->FixParameter(1,6.8*pow(10,3));
    s2func->SetParameter(2,17.88);

    s2moli->FixParameter(0,810.);
    s2moli->FixParameter(1,6.8*pow(10,3));
    s2moli->SetParameter(2,17.88);

    s2lenz->FixParameter(0,810.);
    s2lenz->FixParameter(1,6.8*pow(10,3));
    s2lenz->SetParameter(2,17.88);


    TH1D *h1 = new TH1D("h1", "fprompt hist", 100, 0., 1.);
    TH1D *henergy = new TH1D("henergy", "neutrons energies", 100, 0.000001, 14000.);
    auto fpromptenergy = new TH2F("fpromptenergy", " ; Energy (keV); fprompt", 1000,0.,3000.,100,0.,1.);
    auto fpromptenergy2 = new TH2F("fpromptenergy2", " ; Energy (keV); fprompt", 1000,0.,3000.,100,0.,1.);
    auto h2 = new TH2F("h2", " ; Energy (keV); S2 (PE)", 100, 0., 10., 1000, 0., 1600.);
    TH1D *henergy2 = new TH1D("henergy2", "todos neutrons energies", 100, 0.000001, 14000.);
    
    TH1D *henergysecond = new TH1D("henergysecond", "todos x function - neutrons energies", 100, 0.000001, 14000.);
    auto fpromptenergy22 = new TH2F("fpromptenergy22", " ; Energy (keV); fprompt", 1000,0.,3000.,100,0.,1.);
    auto h22 = new TH2F("h22", " ; Energy (keV); S2 (PE)", 100, 0., 10., 1000, 0., 1600.);
    double fano = 0.107; // fano factor - > see directionality paper
    TRandom3 randGen;  
    randGen.SetSeed(4361);
    double massn = 939.565*1000. ; //MeV/c^2 - > neutron mass
    double dist = 2.; // m ->dist psci baf
    double deltat = 0.7e-9; //s -> uncertaint in time of flight
    double convfactor = 1.6e-16; //kev to Joule
    double c = 299792458.; //light speed
    double mu = 40./pow(41.,2.);
    int counter = 0;
    int counter2 = 0;

    grst = new TGraph("/home/barbaras/red-deconv/barb/g2/stoppingpower.txt");
    grst->GetN(); 
    
    lightout = new TF1("lightout",ESTAR_stoppingpower,0.01,1.,1);
    lightout->SetParameter(0, 0.01441);
    lightout->SetNpx(1000);
 
    auto light = new TF1("light",EJ276D_LightOutput,0.0100001,1.,0);
    light->SetNpx(1000);
    

    for (int i = 0; i < n; i++)
    {   
        double er;
        double cosangle = randGen.Uniform(TMath::Cos(TMath::Pi()/15.),TMath::Cos(17.*TMath::Pi()/180.));

        fer->SetParameter(0,cosangle);
        er = fer->GetRandom(); //recoil energy (MeV)
        
        //double s2 = s2func->Eval(er);

        double epsi = 0.0135*er; //reduced energy
        double beta = 6.8*pow(10.,3.); //normalization constant
        double se = 0.145*sqrt(epsi); //electronic stopping power
        double fz = 0.953; //energy correction factor
        double sn = 0.;
        
        if (strcmp(funcopt, "ziegler") == 0)
            sn = log(1. + 1.1383*fz*epsi)/(2.*(fz*epsi + 0.01321*pow((fz*epsi),0.21226) + 0.19593*pow((fz*epsi),0.5))); // nuclear stopping power

        else if (strcmp(funcopt, "moliere") == 0)
        {
            double part[3] = {0.216, 0.570, 2.37};
            fn->SetParameters(part);
            sn = (1./epsi)*(fn->Integral(0., epsi, 1e-7));
        }

        else if (strcmp(funcopt, "lenz-jensen") == 0)
        {
            double part[3] = {0.191, 0.512, 2.92};
            fn->SetParameters(part);
            sn = (1./epsi)*(fn->Integral(0., epsi, 1e-6));
        }
        
        
        double ni = beta*epsi*se/(sn+se); //number of ionized electrons
        double n0 = ni*1.2; //total number of quanta generated
        double edrift = 20000.0; //[V/m] -> drift electric field
        double gamma = 810./edrift; //another norm constant -> check the meaning of this later
        double prob = (1./(gamma*ni))*log(1. + gamma*ni); //(1-r) -> Thomas-Imel -> prob of a electron to survive recombination and arrive gas pocket
        double ne = ni*prob; //number of electrons that survive recombination

        //cout<<"ni = "<<ni<<"   n0 = "<<n0<<"     ne = "<<ne<<"   raiz-n0*fano = "<<sqrt(fano*n0)<<endl;

        double n0fluc = TMath::Nint(randGen.Gaus(n0,sqrt(fano*n0)));
        int nefluc = randGen.Binomial(n0fluc, ne/n0);
        double s2fq = nefluc*17.88;
        double s2fluc = randGen.Gaus(nefluc*17.88,sqrt(nefluc*17.88 + pow(0.001*nefluc*17.88,2.)));

        
        double energy = er/(2.*mu*(1.-cosangle)) - er; //(MeV)
        double efficiency = pscieff->Eval(energy);

        //double factor = randGen.Uniform(0.,1.);

        double pscidep = (pscyenergy->Eval(energy)); //MeV
        double psciobs = pscidep*efficiency*light->Eval(pscidep/1000.);

        //cerr<<" ENERGIA DEPOSITADA "<<pscidep<<"      BIRKS AQUIIIII   "<<light->Eval(pscidep)<<endl;

        double fprompt = randGen.Gaus(mufprompt->Eval(psciobs),sigmafprompt->Eval(psciobs));
        h1->Fill(fprompt);
        fpromptenergy->Fill(psciobs,fprompt);

        double deltaenergy = 2.*energy*c*sqrt(2.*energy/massn)*deltat/dist;

        double deltaer = 2.*mu*(1.-cosangle)*deltaenergy;

        double erfluc = randGen.Gaus(er,deltaer);
        double fcut = cutfprompt->Eval(psciobs);


        if (er > 0. && s2fluc > 0. && opt == 0)
        {
            h2->Fill(er,ne*17.88);
        }
        else if (er > 0. && s2fluc > 0. && opt == 1)
        {
            h2->Fill(er,s2fq);
        }
        else if (er > 0. && s2fluc > 0. && opt == 2)
        {
            h2->Fill(erfluc,s2fq);
        }
        else if (er > 0. && s2fluc > 4.0*17.88 && opt == 3 && cosangle != 0. && psciobs>200. && fprompt<fcut) 
        {

            //if (factor < efficiency)
            //{
            h2->Fill(erfluc,s2fluc);
            fpromptenergy2->Fill(psciobs,fprompt);
            //henergy->Fill(energy);
            //}
            //h22->Fill(erfluc,s2fluc);
            //fpromptenergy22->Fill(psciobs,fprompt);
            henergy2->Fill(pscidep);
            henergysecond->Fill(energy,pscieff->Eval(energy));
            
        }
        else 
        {
        counter = counter + 1;
        }
        
    }
    cout<<"counter = "<<counter<<endl;


    TCanvas* canvado = new TCanvas("canvado","",1920,1080);
    h2->SetMarkerColor(kViolet);
    h2->SetMarkerStyle(20);
    h2->SetMarkerSize(1.1);
    stylehist2d(h2);
    gStyle->SetTitleFontSize(0.06);
    h2->Draw(drawopt);

    if (strcmp(funcopt, "ziegler") == 0){
    s2func->SetLineColor(coloropt+2);
    styletf1(s2func);
    h2->Fit(s2func,"","",0.,10.);
    s2func->Draw("same");
    leg2->AddEntry("h2", "Ziegler - fitted g2 = 17.04 (0.11)","p");
    stylelegend(leg2);
    leg2->Draw("SAME");
    }
    else if (strcmp(funcopt, "moliere") == 0){
    s2moli->SetLineColor(coloropt+2);
    styletf1(s2moli);
    h2->Fit(s2moli,"","",0.,10.);
    s2moli->Draw("same");
    leg2->AddEntry("h2", "Moliere","p");
    stylelegend(leg2);
    //leg2->Draw("SAME");
    }
    else if (strcmp(funcopt, "lenz-jensen") == 0){
    s2lenz->SetLineColor(coloropt+2);
    styletf1(s2lenz);
    //h2->Fit(s2lenz,"","",0.,10.);
    //s2lenz->Draw("same");
    leg2->AddEntry("h2", Form("Lenz-Jensen - fitted g2 = %.2f (%.2f)",s2lenz->GetParameter(2),s2lenz->GetParError(2)),"p");
    stylelegend(leg2);
    //leg2->Draw("SAME");
    }
    else cout<<"Invalid option"<<endl;
    l.SetTextColor(kBlack);
    l.SetTextSize(0.05);
    //l.DrawLatex(7.7,1370,"Ziegler");
    //l.DrawLatex(7.7,1270,"Moliere");
    //l.DrawLatex(7.7,1170,"Lenz-Jensen");

    //TMarker *marker1 = new TMarker(7.5, 1400, 20);
    //TMarker *marker2 = new TMarker(7.5, 1300, 20);
    //TMarker *marker3 = new TMarker(7.5, 1200, 20);


    //marker1->SetMarkerColor(kBlue); 
    //marker1->SetMarkerSize(2.2);   
    //marker1->Draw("same");
    
    //marker2->SetMarkerColor(kOrange); 
    //marker2->SetMarkerSize(2.2);   
   // marker2->Draw("same");

    //marker3->SetMarkerColor(kGreen); 
    //marker3->SetMarkerSize(2.2);   
    //marker3->Draw("same");

    TCanvas *coelho1 = new TCanvas("coelho1","coelho1",800, 800);
    pscieff->SetTitle("Eficiencia - EJ-276; Energia (MeV); Eficiencia");
    pscieff->Draw();


    
    double intstart[7] = {0.5, 2.5,3.5, 4.5, 5.5, 6.5, 7.5};
    double intend[7] = {2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 10.};
    double percent[7] = {0.15, 0.05,0.14,0.05,0.08,0.12,0.21};

    TH1D* h2proj[7];
    TH2F* h2copy;
    h2copy = (TH2F*) h2->Clone("h2_clone");
    TCanvas* cevice = new TCanvas("cevice","",1920,1080);
    stylehist1D(h1);
    h1->Draw();
    TCanvas* fhist = new TCanvas("fhist","",1920,1080);
    fpromptenergy->SetMarkerStyle(20);
    stylehist2d(fpromptenergy);
    fpromptenergy->SetMarkerSize(1.5);
    fpromptenergy->Draw();
    cutfprompt->SetLineColor(kRed);
    cutfprompt->Draw("SAME");
    TCanvas* fhist2 = new TCanvas("fhist2","",1920,1080);
    fpromptenergy2->SetMarkerStyle(20);
    stylehist2d(fpromptenergy2);
    fpromptenergy2->SetMarkerSize(1.5);
    fpromptenergy2->Draw();
    cutfprompt->SetLineColor(kRed);
    cutfprompt->Draw("SAME");
    //TCanvas* canvenergy = new TCanvas("canvenergy","",1920,1080);
    //henergy->Draw();
    
    //TCanvas* fhist22 = new TCanvas("fhist22","",1920,1080);
    //fpromptenergy22->SetMarkerStyle(20);
    //stylehist2d(fpromptenergy22);
    //fpromptenergy22->SetMarkerSize(1.5);
    //fpromptenergy22->Draw();
    TCanvas* canvenergy2 = new TCanvas("canvenergy2","",1920,1080);
    henergy2->Draw();
    henergysecond->SetLineColor(kRed);
    henergysecond->Draw("SAME");


    //TCanvas* canvas[7];
    TCanvas* canvas;
    double means[7] = {0.,0.,0.,0.,0.,0.,0.};
    double er2,s2ran;

    
    for (int i = 0; i<7; i++)
    {
    std::string hname=Form("h2proj_%d",i);
	h2proj[i]=new TH1D(hname.c_str(),"",25,0.,80.);
    h2proj[i] = h2->ProjectionY("",h2->GetXaxis()->FindBin(intstart[i]), h2->GetXaxis()->FindBin(intend[i]));
    //h2proj[i]->RebinX(25);
    //canvas[i]=new TCanvas(Form("cdraw_%d",i),Form("cdraw_%d",i),1920,1080);
    //canvas=new TCanvas("cdraw1","",1920,1080);
    h2proj[i]->SetYTitle("Counts/1.0 PE");
    h2proj[i]->SetXTitle("N_{i.e.}");
    stylehist1D(h2proj[i]);
    //h2proj[i]->DrawCopy();
    auto *leg1 = new TLegend(0.45,0.66,0.87,0.75);
    stylelegend(leg1);
    leg1->SetTextSize(0.065);
    leg1->AddEntry((TObject*)0, Form("Energy interval: %.1f to %.1f",intstart[i],intend[i]), "");
    //leg1->DrawClone("same");
    int entries = h2proj[i]->GetEntries();
    int bkg = percent[i]*entries;
    cout<<"bkg numb = "<<bkg<<endl;
    for (int j = 0; j<bkg; j++)
    {
        s2ran = randGen.Uniform(0.,1610.);
        er2 = randGen.Uniform(intstart[i],intend[i]);
        h2copy->Fill(er2,s2ran);
    }
    }
    canvas=new TCanvas("cdraw1","",1920,1080);
    //h2copy->Fit(s2lenz,"","",0.,10.);
    //leg3->AddEntry("h2", Form("Fitted g2 = %.2f (%.2f)",s2lenz->GetParameter(2),s2lenz->GetParError(2)),"p");
    //stylelegend(leg3);
    //h2copy->SetTitle("Copy here");
    h2copy->Draw();
    //leg3->Draw("SAME");

    MCMC_Fit_LenzJensen(h2copy);


    

    /*
    h22->SetMarkerColor(kViolet);
    h22->SetMarkerStyle(20);
    h22->SetMarkerSize(1.1);
    TCanvas* canven2 = new TCanvas("canven2","",1920,1080);
    stylehist2d(h22);
    h22->Draw();
    leg3->Draw("SAME");
    */


    //cout<<"counter = "<<counter<<endl;

}

void compcalib() 
{
   auto c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
//   c1->SetFillColor(42);
   c1->SetGrid();
//   c1->GetFrame()->SetFillColor(21);
//   c1->GetFrame()->SetBorderSize(12);
   const Int_t n = 7;
   Double_t x1[n]  = {1.5,3.0,4.0,5.0,6.0,7.0,8.75}; 
   Double_t x2[n]  = {1.5,3.0,4.0,5.0,6.0,7.0,8.75};
   Double_t y1[n]  = {197.56/(17.88*1.5),319.71/(17.88*1.5),416.81/(17.88*1.5),465.01/(17.88*1.5),558.79/(17.88*1.5),585.74/(17.88*1.5),648.13/(17.88*1.5)};
   Double_t ey1[n] = {3.568/(17.88*1.5),4.739/(17.88*1.5),13.38/(17.88*1.5),20.412/(17.88*1.5),28.38/(17.88*1.5),25.92/(17.88*1.5),74.91/(17.88*1.5)};
   Double_t y4[n]  = {248.89/(17.88*1.5),364.75/(17.88*1.5),451.98/(17.88*1.5),520.84/(17.88*1.5),608.35/(17.88*1.5),620.36/(17.88*1.5),772.79/(17.88*1.5)};
   Double_t ey4[n] = {3.658/(17.88*1.5),5.366/(17.88*1.5),13.83/(17.88*1.5),27.92/(17.88*1.5),53.88/(17.88*1.5),132.92/(17.88*1.5),71.14/(17.88*1.5)};
   Double_t ex[n] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
   auto gr1 = new TGraphErrors(n,x1,y1,ex,ey1);
   auto gr4 = new TGraphErrors(n,x2,y4,ex,ey4);
   gr1->SetMarkerColor(kBlue);
   gr1->SetMarkerStyle(21);
   gr1->SetMarkerSize(2);
   gr1->SetTitle("Ziegler");
   gr4->SetMarkerColor(kGreen);
   gr4->SetMarkerStyle(21);
   gr4->SetMarkerSize(2);
   gr4->SetTitle("Lenz-Jensen");
   TMultiGraph *mg = new TMultiGraph();
   mg->SetTitle(" ;Energy (keV); Qy");
   mg->Add(gr1);
   mg->Add(gr4);
   mg->Draw("AP");
   c1->BuildLegend(0.7,0.7,0.9,0.9);
}

void comp_s2()
{

    auto s2func = new TF1("s2func", ziegler, 0.0001, 10., 3);
    auto s2moli = new TF1("s2moli", moliere, 0.0001, 10., 3);
    auto s2lenz = new TF1("s2lenz", lenzjensen, 0.0001, 10., 3);
    auto fn = new TF1("fn", fnfunc, 0., 20., 3);

    double theta4[3] = {810.,6.8*pow(10,3),17.88};
    //s2func->SetParameters(theta4);    
    s2func->FixParameter(0,810.);
    s2func->FixParameter(1,6.8*pow(10,3));
    s2func->SetParameter(2,17.88);

    s2moli->FixParameter(0,810.);
    s2moli->FixParameter(1,6.8*pow(10,3));
    s2moli->SetParameter(2,17.88);

    s2lenz->FixParameter(0,810.);
    s2lenz->FixParameter(1,6.8*pow(10,3));
    s2lenz->SetParameter(2,17.88);

    s2func->Draw();
    s2moli->SetLineColor(kGreen);
    s2moli->Draw("same");
    s2lenz->SetLineColor(kBlue);
    s2lenz->Draw("same");

}



Double_t funcECf(Double_t *val, Double_t *par)
{
    Double_t x = val[0];
    if (x <= 0) return 0;  // Prevent integration issues
    return TMath::Exp(-0.88*x) * TMath::SinH(sqrt(2.0*(x + 1e-6)));
}


Double_t func(Double_t *val, Double_t *par)
{
   Double_t x = val[0];
   Double_t alfa = (1-par[0])*80/(pow(41,2));
   Double_t f = (1/alfa)*TMath::Exp(-0.88*x/(alfa*1000))*TMath::SinH(sqrt(2.0*x/(alfa*1000)));
   return f;
}


Double_t ziegler(Double_t *x, Double_t *par) //par = {cbox,beta,g2};
{

    Float_t X = x[0]; //x = Er (keV)
    if (X == 0) return 0;
    double epsi = 0.0135*X;
    double beta = par[1];
    double se = 0.145*sqrt(epsi);
    double fz = 0.953;
    double sn = log(1. + 1.1383*fz*epsi)/(2.*(fz*epsi + 0.01321*pow((fz*epsi),0.21226) + 0.19593*pow((fz*epsi),0.5)));
    double ni = beta*epsi*se/(sn+se);

    double edrift = 20000.0; //[V/m]
    double gamma = par[0]/edrift;
    double prob = (1./(gamma*ni))*log(1 + gamma*ni); //(1-r) -> Thomas-Imel

    double nie = ni*prob;
    double s2 = par[2]*nie;

    return s2;
}

Double_t fnfunc(Double_t *x, Double_t *par) //par = {m,q,lambda};
{
    Float_t eta = x[0];
    double lambda = par[2];
    double q = par[1];
    double m = par[0];

    double fn = lambda*pow(eta,1.-2.*m)/pow(1. + pow(2.*lambda*pow(eta,2.*(1.-m)),q),1./q);

    return fn;
}

Double_t moliere(Double_t *x, Double_t *par) //par = {cbox,beta,g2};
{

    Float_t X = x[0]; //x = Er (keV)
    double epsi = 0.0135*X;
    double beta = par[1];
    double se = 0.145*sqrt(epsi);
    double fz = 0.953;
    
    double part[3] = {0.216, 0.570, 2.37};
    auto fn = new TF1("fn", fnfunc, 0., 250, 3);
    fn->SetParameters(part);
    double sn = (1./epsi)*(fn->Integral(0., epsi,1e-8));

    double ni = beta*epsi*se/(sn+se);

    double edrift = 20000.0; //[V/m]
    double gamma = par[0]/edrift;
    double prob = (1./(gamma*ni))*log(1 + gamma*ni); //(1-r) -> Thomas-Imel

    double nie = ni*prob;
    double s2 = par[2]*nie;

    delete fn;

    return s2;
}

Double_t lenzjensen(Double_t *x, Double_t *par) //par = {cbox,beta,g2};
{

    Float_t X = x[0]; //x = Er (keV)
    double epsi = 0.0135*X;
    double beta = par[1];
    double se = 0.145*sqrt(epsi);
    double fz = 0.953;

    double part[3] = {0.191, 0.512, 2.92};
    auto fn = new TF1("fn", fnfunc, 0., 250, 3);
    fn->SetParameters(part);
    double sn = (1./epsi)*(fn->Integral(0., epsi,1e-7));
    
    double ni = beta*epsi*se/(sn+se);

    double edrift = 20000.0; //[V/m]
    double gamma = par[0]/edrift;
    double prob = (1./(gamma*ni))*log(1 + gamma*ni); //(1-r) -> Thomas-Imel

    double nie = ni*prob;
    double s2 = par[2]*nie;

    delete fn;

    return s2;
}


Double_t pscienergy(Double_t *x, Double_t *par) //par = opt ---- model for the energy deposited in the psci, given the energy of the neutron
{
     Float_t X = x[0]; //x = Er (MeV)
    double opt = par[0];
    double  a=6.476e-01, b=1.714e+00, c=3.496e-01;  //astro-ph: 2011.01907
    //double  a=7.300e-01, b=3.500e+00, c=2.000e-01;  //NIM A 1051 (2023) 168216
    auto tquenching_nr_psci=new TF1("tquenching_nr_psci","([0]*x-[1]*(1.-exp(-[2]*x)))*1.2",0.,20.);
    tquenching_nr_psci->SetParameter(0,a);
    tquenching_nr_psci->SetParameter(1,b);
    tquenching_nr_psci->SetParameter(2,c);
    double fq=tquenching_nr_psci->Eval(X)/(X);
    
    double Qene=X*fq; 
    double quenchfactor = 0.884;

    //return Qene*quenchfactor;
    
    double Qobs;

    static TRandom3 rnd(0);
    double Fano=2.5; // 20% for Am241
    double Sigma=sqrt(Qene*Fano);
    Qobs=rnd.Gaus(Qene,Sigma);
    if ( Qobs<0. ) Qobs=0.;

    return Qobs*quenchfactor;
}

Double_t Fprompt_PSci_Cut(Double_t *x, Double_t *par)
{
  Float_t E0 = x[0];
  if ( E0>2000. ) E0=2000.;
  double parex[4]={  6.82e-01 , 1.01e-05 , 6.35e+01 ,  5.00e-01  };
  double val=parex[0]+parex[1]*E0-exp(-E0/parex[2])*parex[3];
  if ( val<0. ) val=0.;
  return val;
}

Double_t Fprompt_mu(Double_t *x, Double_t *par)
{
  Float_t E0 = x[0];
  if ( E0>3000. ) E0=3000.;
  double parex[2]={   6.3461e-01  ,  4.3481e-02 };
  double val=parex[0]*pow(E0/1000.,parex[1]);

  return val;
}

Double_t Fprompt_sigma(Double_t *x, Double_t *par)
{
  Float_t E0 = x[0];
  if ( E0>3000. ) E0=3000.;
  double parex[3]={ 5.9943e-03 ,  5.6147e-01 ,  -5.7931e-01 };
  double val=parex[0]+parex[1]*pow(E0,parex[2]);
  
  return val;
}


Double_t EJ276D_efficiency(Double_t *x, Double_t *par)
{
    Float_t X = x[0];

    if (X>9.) X = 9.;

    const Int_t n = 33;
    Double_t x1[n]  = {0.0123, 0.2952, 0.5781, 0.8487, 1.1316, 1.4022, 1.6851, 1.9557, 2.2386, 2.5092, 2.7921, 3.0627, 3.3456, 3.6285, 3.8991, 4.1820, 
        4.4649, 4.7355, 5.0062, 5.2891, 5.5720, 5.8426, 6.1132, 6.4084, 6.6667, 6.9496, 7.2325, 7.5031, 7.7860, 8.0689, 8.3272, 8.6224, 8.8930}; //before
    Double_t y1[n]  = {0.0855, 0.2214, 0.3128, 0.3547, 0.3613, 0.3564, 0.3597, 0.3383, 0.3259, 0.3053, 0.3119, 0.2922, 0.2774, 0.2790, 0.2617, 0.2683, 
        0.2395, 0.2296, 0.2305, 0.2123, 0.2206, 0.1975, 0.2156, 0.1942, 0.1984, 0.1753, 0.1975, 0.1761, 0.1654, 0.1572, 0.1885, 0.2033, 0.1868}; 

    
    
    auto gr1 = new TGraph(n,x1,y1);
    gr1->SetMarkerColor(4);
    gr1->SetMarkerStyle(21);
    gr1->SetMarkerSize(2);

    return gr1->Eval(X,0,"S");

    //if (X>1.) return sqrt(1. + X/1.)/15.;
    //else return 0.;

}

void birks()
{

    grst = new TGraph("/home/barbaras/red-deconv/barb/g2/stoppingpower.txt");
    grst->GetN(); 
    

    lightout = new TF1("lightout",ESTAR_stoppingpower,0.01,1.,1);
    lightout->SetParameter(0, 0.01441);

    TCanvas *cc = new TCanvas();
    lightout->SetNpx(1000);
    lightout->Draw();
 
    TCanvas *c1 = new TCanvas();
    auto light = new TF1("light",EJ276D_LightOutput,0.0100001,1.,0);
    light->SetNpx(1000);
    light->Draw();

}


Double_t EJ276D_LightOutput(Double_t *x, Double_t *par)
{
    double X = x[0];
    double ltout = lightout->Integral(0.01, X, 1.e-6);

    if (X>0.47743) return 1.;
    
    else 
    return (0.47743*ltout/(X*lightout->Integral(0.01,0.47743,1.e-10)));

}


Double_t ESTAR_stoppingpower(Double_t *x, Double_t *par)
{

    Float_t X = x[0]; //x = E (MeV)
    
    double st = grst->Eval(X,0,"");
    
    return 1./(1. + par[0]*st);

}

// ==========================================================
// FUNÇÕES DO MCMC (Atualizadas com Modelagem de Fundo)
// ==========================================================
Double_t fnfunc_mcmc(Double_t *x, Double_t *par) {
    double eta = x[0];
    double m = par[0];
    double q = par[1];
    double lambda = par[2];
    return lambda * pow(eta, 1. - 2.*m) / pow(1. + pow(2.*lambda * pow(eta, 2.*(1.-m)), q), 1./q);
}

double CalculateNie_LenzJensen(double Er_keV, double cbox, double beta) {
    if (Er_keV <= 0) return 0.0;
    double epsi = 0.0135 * Er_keV;
    double se = 0.145 * sqrt(epsi);
    TF1 *fn = new TF1("fn", fnfunc_mcmc, 0., 250., 3);
    fn->SetParameters(0.191, 0.512, 2.92); 
    double sn = (1. / epsi) * (fn->Integral(0., epsi, 1e-7));
    delete fn; 
    double ni = beta * epsi * se / (sn + se);
    double edrift = 20000.0; 
    double gamma = cbox / edrift;
    double prob = (1. / (gamma * ni)) * log(1. + gamma * ni); 
    return ni * prob;
}

double CalculateLogLikelihood2D_WithBkg(TH2F* h2_data, const vector<double>& Nie_cache, double g2_guess) {
    double logL = 0.0;
    int nBinsX = h2_data->GetNbinsX(); 
    int nBinsY = h2_data->GetNbinsY(); 
    
    // Largura do bin no eixo Y (S2) para calcular a densidade do fundo uniforme
    double s2_bin_width = h2_data->GetYaxis()->GetBinWidth(1);

    for (int ix = 1; ix <= nBinsX; ix++) {
        double N_ie = Nie_cache[ix]; 
        if (N_ie <= 0) continue;
        
        // 1. Conta o total de eventos reais (sinal + fundo) nesta fatia de energia
        double total_events_in_slice = 0;
        for (int iy = 1; iy <= nBinsY; iy++) {
            total_events_in_slice += h2_data->GetBinContent(ix, iy);
        }
        if (total_events_in_slice == 0) continue; 

        // 2. Estimativa da fração de fundo com base no seu vetor percent[] do SimS2.C
        // Como o fundo é uniformemente distribuído de 0 a 1610 PE (amplitude total de S2):
        double S2_range_total = 1610.0; 
        double bkg_fraction = 0.10; // Valor médio estimado do seu array percent[]
        double total_bkg_events = total_events_in_slice * bkg_fraction;
        double total_signal_events = total_events_in_slice - total_bkg_events;
        
        // Densidade esperada de fundo por bin de S2 nesta fatia
        double bkg_lambda_per_bin = (total_bkg_events / S2_range_total) * s2_bin_width;

        // 3. Física do Sinal
        double mean_S2 = g2_guess * N_ie;
        double sigma_S2 = sqrt(mean_S2 + pow(0.001 * mean_S2, 2.0)); 
        
        // 4. Soma de Verossimilhança de Poisson por bin (Sinal + Fundo)
        for (int iy = 1; iy <= nBinsY; iy++) {
            double n_obs = h2_data->GetBinContent(ix, iy);
            double S2_val = h2_data->GetYaxis()->GetBinCenter(iy);
            
            // Sinal esperado normalizado pelos eventos de sinal reais
            double signal_lambda = total_signal_events * TMath::Gaus(S2_val, mean_S2, sigma_S2, true) * s2_bin_width;
            
            // Expectativa total (lambda) = Sinal + Fundo
            double lambda = signal_lambda + bkg_lambda_per_bin;
            
            if (lambda < 1e-6) lambda = 1e-6; 
            
            logL += n_obs * log(lambda) - lambda;
        }
    }
    return logL;
}

// void MCMC_Fit_LenzJensen(TH2F* h2_to_fit) { 
//     double cbox = 810.0;
//     double beta = 6800.0;
//     int nBinsX = h2_to_fit->GetNbinsX();
//     vector<double> Nie_cache(nBinsX + 1, 0.0);
    
//     cout << "Pre-calculando integrais de Lenz-Jensen para o Cache... Aguarde." << endl;
//     for (int ix = 1; ix <= nBinsX; ix++) {
//         double E_val = h2_to_fit->GetXaxis()->GetBinCenter(ix);
//         Nie_cache[ix] = CalculateNie_LenzJensen(E_val, cbox, beta);
//     }
//     cout << "Cache concluido! Iniciando o MCMC com correcao de Fundo..." << endl;
    
//     TRandom3 randGen(0); 
//     int n_steps = 20000;          
//     int burn_in = 3000;           
//     double step_size = 0.25;      
//     double current_g2 = 16.0;     
    
//     TH1D* h_posterior = new TH1D("h_posterior", "Posterior distribution de g2; g2 (PE/e-); Freq", 100, 16.5, 18.5);
//     TGraph* gr_trace = new TGraph();
//     gr_trace->SetTitle("Trace Plot; Passo; Valor de g2");

//     // Usa a verossimilhança atualizada com o fundo
//     double current_logL = CalculateLogLikelihood2D_WithBkg(h2_to_fit, Nie_cache, current_g2);
//     int accepted_steps = 0;

//     for (int step = 0; step < n_steps; step++) {
//         if (step % 2000 == 0) cout << "MCMC Passo: " << step << " / " << n_steps << endl;
//         double proposed_g2 = randGen.Gaus(current_g2, step_size);
//         if (proposed_g2 <= 0) continue; 
        
//         double proposed_logL = CalculateLogLikelihood2D_WithBkg(h2_to_fit, Nie_cache, proposed_g2);
        
//         double log_u = log(randGen.Uniform(0.0, 1.0));
//         if (log_u < (proposed_logL - current_logL)) {
//             current_g2 = proposed_g2;
//             current_logL = proposed_logL;
//             accepted_steps++;
//         }
        
//         if (step > burn_in) h_posterior->Fill(current_g2);
//         gr_trace->SetPoint(step, step, current_g2);
//     }
    
//     cout << "------------------------------------------" << endl;
//     cout << "Taxa de aceitacao: " << (double)accepted_steps / n_steps * 100.0 << " %" << endl;
//     cout << "g2 estimado (Lenz-Jensen com Fundo): " << h_posterior->GetMean() << " +/- " << h_posterior->GetStdDev() << endl;
    
//     TCanvas* c_mcmc = new TCanvas("c_mcmc", "MCMC 2D Result", 1200, 500);
//     c_mcmc->Divide(2, 1);
//     c_mcmc->cd(1);
//     gr_trace->SetMarkerStyle(1);
//     gr_trace->Draw("AP"); 
//     c_mcmc->cd(2);
//     h_posterior->SetFillColor(kTeal-9);
//     h_posterior->Draw();
//     auto* leg_post = new TLegend(0.15, 0.75, 0.45, 0.88);
//     leg_post->AddEntry("", Form("Mean = %.2f)", h_posterior->GetMean()), "");
//     leg_post->AddEntry("", Form("Std dev = %.2f)", h_posterior->GetStdDev()), "");
//     leg_post->Draw("SAME");

//     // Plotagem da Curva do Modelo sobre o Espaço 2D
//     TCanvas* c_fit2d = new TCanvas("c_fit2d", "Ajuste MCMC no Espaço 2D", 1000, 800);
//     h2_to_fit->SetTitle("Ajuste MCMC com Fundo no Espaço E_{NR} x S2; Recoil Energy (keV); S2 (PE)");
//     h2_to_fit->Draw();
    
//     double g2_best = h_posterior->GetMean();
    
//     TF1* f_mean_s2 = new TF1("f_mean_s2", [=](Double_t *x, Double_t *par) {
//         double E_keV = x[0];
//         int binX = h2_to_fit->GetXaxis()->FindBin(E_keV);
//         if (binX >= 1 && binX <= h2_to_fit->GetNbinsX()) {
//             double N_ie = Nie_cache[binX];
//             return g2_best * N_ie;
//         }
//         return 0.0;
//     }, h2_to_fit->GetXaxis()->GetXmin(), h2_to_fit->GetXaxis()->GetXmax(), 0);
    
//     f_mean_s2->SetLineColor(kRed);
//     f_mean_s2->SetLineWidth(3);
//     f_mean_s2->Draw("SAME");
    
//     auto* leg_mcmc2d = new TLegend(0.15, 0.75, 0.45, 0.88);
//     leg_mcmc2d->AddEntry(h2_to_fit, "Amostra com Background (Monte Carlo)", "lep");
//     leg_mcmc2d->AddEntry(f_mean_s2, Form("Modelo Lenz-Jensen (g_{2} = %.2f)", g2_best), "l");
//     leg_mcmc2d->Draw("SAME");
// }


void MCMC_Fit_LenzJensen(TH2F* h2_to_fit) { 
    double cbox = 810.0;
    double beta = 6800.0;
    int nBinsX = h2_to_fit->GetNbinsX();
    std::vector<double> Nie_cache(nBinsX + 1, 0.0);
    
    cout << "Pre-calculando integrais de Lenz-Jensen para o Cache... Aguarde." << endl;
    for (int ix = 1; ix <= nBinsX; ix++) {
        double E_val = h2_to_fit->GetXaxis()->GetBinCenter(ix);
        Nie_cache[ix] = CalculateNie_LenzJensen(E_val, cbox, beta);
    }
    cout << "Cache concluido! Iniciando o MCMC com correcao de Fundo e Lag..." << endl;
    
    TRandom3 randGen(0); 
    int n_steps = 20000;          
    int burn_in = 0;           
    int lag = 10;                 // <--- NOVO: Fator de Thinning (Descorrelação)
    double step_size = 0.25;      
    double current_g2 = 16.0;     
    
    TH1D* h_posterior = new TH1D("h_posterior", "Posterior de g2 (Descorrelacionada); g2 (PE/e-); Freq", 100, 16.5, 18.5);
    h_posterior->SetDirectory(0); // <--- Evita que o ROOT delete o histograma
    
    TGraph* gr_trace = new TGraph();
    gr_trace->SetTitle("Trace Plot; Passo; Valor de g2");
    
    std::vector<double> full_chain; // <--- NOVO: Armazena a cadeia para calcular a ACF

    double current_logL = CalculateLogLikelihood2D_WithBkg(h2_to_fit, Nie_cache, current_g2);
    int accepted_steps = 0;

    for (int step = 0; step < n_steps; step++) {
        if (step % 2000 == 0) cout << "MCMC Passo: " << step << " / " << n_steps << endl;
        double proposed_g2 = randGen.Gaus(current_g2, step_size);
        if (proposed_g2 <= 0) continue; 
        
        double proposed_logL = CalculateLogLikelihood2D_WithBkg(h2_to_fit, Nie_cache, proposed_g2);
        
        double log_u = log(randGen.Uniform(0.0, 1.0));
        if (log_u < (proposed_logL - current_logL)) {
            current_g2 = proposed_g2;
            current_logL = proposed_logL;
            accepted_steps++;
        }
        
        if (step > burn_in) {
            full_chain.push_back(current_g2); // Guarda todos os passos pós burn-in para a ACF
            
            // Aplica o fator de lag: preenche a posterior apenas a cada 'lag' passos
            if ((step - burn_in) % lag == 0) {
                h_posterior->Fill(current_g2);
            }
        }
        gr_trace->SetPoint(step, step, current_g2);
    }
    
    // ==============================================================================
    // CÁLCULO DA AUTOCORRELAÇÃO (ACF)
    // ==============================================================================
    int n_chain = full_chain.size();
    double mean_chain = 0.0;
    for (double val : full_chain) mean_chain += val;
    mean_chain /= n_chain;
    
    double var_chain = 0.0;
    for (double val : full_chain) var_chain += (val - mean_chain) * (val - mean_chain);
    var_chain /= n_chain;
    
    TGraph* gr_acf = new TGraph();
    gr_acf->SetTitle("Funcao de Autocorrelacao (ACF); Lag k (Passos); ACF");
    
    int max_lag = 100; // Analisar a correlação em até 100 passos de distância
    for (int k = 0; k <= max_lag; k++) {
        double acf_k = 0.0;
        for (int i = 0; i < n_chain - k; i++) {
            acf_k += (full_chain[i] - mean_chain) * (full_chain[i+k] - mean_chain);
        }
        acf_k = acf_k / n_chain; 
        acf_k = acf_k / var_chain; // Normaliza para ter ACF = 1 no lag 0
        gr_acf->SetPoint(k, k, acf_k);
    }

    // ==============================================================================
    // RESULTADOS NO TERMINAL E PLOTAGEM DOS DIAGNÓSTICOS
    // ==============================================================================
    cout << "------------------------------------------" << endl;
    cout << "Taxa de aceitacao: " << (double)accepted_steps / n_steps * 100.0 << " %" << endl;
    cout << "g2 estimado (Lenz-Jensen com Fundo): " << h_posterior->GetMean() << " +/- " << h_posterior->GetStdDev() << endl;
    
    // Canvas dividido em 3 painéis para acomodar a ACF
    TCanvas* c_mcmc = new TCanvas("c_mcmc", "MCMC Diagnosticos", 1500, 500);
    c_mcmc->Divide(3, 1);
    
    c_mcmc->cd(1);
    gr_trace->SetMarkerStyle(1);
    gr_trace->Draw("AP"); 
    
    c_mcmc->cd(2);
    gr_acf->SetLineColor(kRed);
    gr_acf->SetMarkerStyle(20);
    gr_acf->SetMarkerSize(0.6);
    gr_acf->Draw("ALP"); 
    
    c_mcmc->cd(3);
    h_posterior->SetFillColor(kTeal-9);
    h_posterior->Draw();
    
    auto* leg_post = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg_post->AddEntry((TObject*)0, Form("Mean = %.2f", h_posterior->GetMean()), "");
    leg_post->AddEntry((TObject*)0, Form("Std dev = %.2f", h_posterior->GetStdDev()), "");
    leg_post->Draw("SAME");

    // ==============================================================================
    // PLOTAR A CURVA DO MODELO SOBRE O ESPAÇO E_NR x S2
    // ==============================================================================
    TCanvas* c_fit2d = new TCanvas("c_fit2d", "Ajuste MCMC no Espaço 2D", 1000, 800);
    h2_to_fit->SetTitle("Ajuste MCMC com Fundo no Espaço E_{NR} x S2; Recoil Energy (keV); S2 (PE)");
    h2_to_fit->Draw("COLZ");
    
    double g2_best = h_posterior->GetMean();
    
    TF1* f_mean_s2 = new TF1("f_mean_s2", [=](Double_t *x, Double_t *par) {
        double E_keV = x[0];
        int binX = h2_to_fit->GetXaxis()->FindBin(E_keV);
        if (binX >= 1 && binX <= h2_to_fit->GetNbinsX()) {
            double N_ie = Nie_cache[binX];
            return g2_best * N_ie;
        }
        return 0.0;
    }, h2_to_fit->GetXaxis()->GetXmin(), h2_to_fit->GetXaxis()->GetXmax(), 0);
    
    f_mean_s2->SetLineColor(kRed);
    f_mean_s2->SetLineWidth(3);
    f_mean_s2->Draw("SAME");
    
    auto* leg_mcmc2d = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg_mcmc2d->AddEntry(h2_to_fit, "Amostra com Background", "lep");
    leg_mcmc2d->AddEntry(f_mean_s2, Form("Modelo Lenz-Jensen (g_{2} = %.2f)", g2_best), "l");
    leg_mcmc2d->Draw("SAME");
}