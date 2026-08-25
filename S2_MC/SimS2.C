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

    auto f = new TF1("f", funcECf, 0.0001, 15., 0);
    auto fer = new TF1("fer", func, 0.0001, 15., 1);

    auto cutfprompt = new TF1("cutfprompt", Fprompt_PSci_Cut, 0.0001, 3000., 0);
    auto mufprompt = new TF1("mufprompt", Fprompt_mu, 0.0001, 3000., 0);
    auto sigmafprompt = new TF1("sigmafprompt", Fprompt_sigma, 0.0001, 3000., 0);
    auto pscyenergy = new TF1("pscyenergy", pscienergy, 0.0001, 20., 1);
    auto pscieff = new TF1("pscieff", EJ276D_efficiency, 0.0001, 20., 0);
    
    pscieff->SetNpx(1000);


    auto s2func = new TF1("s2func", ziegler, 0.0001, 10., 3);
    auto s2moli = new TF1("s2moli", moliere, 0.0001, 10., 3);
    auto s2lenz = new TF1("s2lenz", lenzjensen, 0.0001, 10., 3);
    auto fn = new TF1("fn", fnfunc, 0., 20., 3);
    
    pscyenergy->SetParameter(0,1.);

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
    auto fpromptenergy = new TH2F("fpromptenergy", " ; Energy (keV); fprompt", 1000,0.,3000.,100,0.,1.);
    auto fpromptenergy2 = new TH2F("fpromptenergy2", " ; Energy (keV); fprompt", 1000,0.,3000.,100,0.,1.);
    auto h2 = new TH2F("h2", " ; Energy (keV); S2 (PE)", 100, 0., 10., 1000, 0., 1600.);
    TH1D *henergy2 = new TH1D("henergy2", "todos neutrons energies", 100, 0.000001, 14000.);
    TH1D *henergysecond = new TH1D("henergysecond", "todos x function - neutrons energies", 100, 0.000001, 14000.);
    
	double fano = 0.107; // fano factor - > see directionality paper
    TRandom3 randGen;  
    randGen.SetSeed(4361); //Seed is fixed to enable direct comparison between samples
    double massn = 939.565*1000. ; //MeV/c^2 - > neutron mass
    double dist = 2.; // m ->dist psci baf
    double deltat = 0.7e-9; //s -> uncertaint in time of flight
    double convfactor = 1.6e-16; //kev to Joule
    double c = 299792458.; //light speed
    double mu = 40./pow(41.,2.); //reduced neutron-argon_nucleus mass
    int counter = 0;
    int counter2 = 0;

    grst = new TGraph("/home/barbaras/red-deconv/barb/g2/stoppingpower.txt");
    grst->GetN(); 
    
    auto light = new TF1("light",EJ276D_LightOutput,0.0100001,1.,0);
    light->SetNpx(1000);
    

    for (int i = 0; i < n; i++)
    {   
        double er;
        double cosangle = randGen.Uniform(TMath::Cos(TMath::Pi()/15.),TMath::Cos(17.*TMath::Pi()/180.));

        fer->SetParameter(0,cosangle);
        er = fer->GetRandom(); //recoil energy (MeV)

        double epsi = 0.0135*er; //reduced energy
        double beta = 6.8*pow(10.,3.); //normalization constant
        double se = 0.145*sqrt(epsi); //electronic stopping power
        double fz = 0.953; //energy correction factor
        double sn = 0.;
        
        if (strcmp(funcopt, "ziegler") == 0)
            sn = log(1. + 1.1383*fz*epsi)/(2.*(fz*epsi + 0.01321*pow((fz*epsi),0.21226) + 0.19593*pow((fz*epsi),0.5))); // ziegler nuclear stopping power

        else if (strcmp(funcopt, "moliere") == 0)
        {
            double part[3] = {0.216, 0.570, 2.37};
            fn->SetParameters(part);
            sn = (1./epsi)*(fn->Integral(0., epsi, 1e-7)); //moliere nuclear stopping power
        }

        else if (strcmp(funcopt, "lenz-jensen") == 0)
        {
            double part[3] = {0.191, 0.512, 2.92};
            fn->SetParameters(part);
            sn = (1./epsi)*(fn->Integral(0., epsi, 1e-6)); //lenz-jense nuclear stopping power
        }
        
        
        double ni = beta*epsi*se/(sn+se); //number of ionized electrons
        double n0 = ni*1.2; //total number of quanta generated
        double edrift = 20000.0; //[V/m] -> drift electric field
        double gamma = 810./edrift; //another norm constant -> check the meaning of this later
        double prob = (1./(gamma*ni))*log(1. + gamma*ni); //(1-r) -> Thomas-Imel -> prob of a electron to survive recombination and arrive gas pocket
        double ne = ni*prob; //number of electrons that survive recombination


        double n0fluc = TMath::Nint(randGen.Gaus(n0,sqrt(fano*n0))); //implementing fluctuation in the total quanta generated
        int nefluc = randGen.Binomial(n0fluc, ne/n0); //implementing fluctuation in the number of ionized electrons
        double s2fq = nefluc*17.88; //defining S2 with only quantum fluctuation
        double s2fluc = randGen.Gaus(nefluc*17.88,sqrt(nefluc*17.88 + pow(0.001*nefluc*17.88,2.))); //S2 with quantum and detection fluctuations

        
        double energy = er/(2.*mu*(1.-cosangle)) - er; //(MeV)
        double efficiency = pscieff->Eval(energy); //psci efficiency


        double pscidep = (pscyenergy->Eval(energy)); //energy that a neutron deposits in the psci (MeV)
        double psciobs = pscidep*efficiency*light->Eval(pscidep/1000.); //energy observed in the psci


        double fprompt = randGen.Gaus(mufprompt->Eval(psciobs),sigmafprompt->Eval(psciobs)); //generating a fprompt value for this neutron in the psci
        h1->Fill(fprompt); //filling histogram with all fprompt generated
        fpromptenergy->Fill(psciobs,fprompt); //filling histogram of fprompt vs. energy

        double deltaenergy = 2.*energy*c*sqrt(2.*energy/massn)*deltat/dist; //uncertainty on the reconstructed neutron energy

        double deltaer = 2.*mu*(1.-cosangle)*deltaenergy; //propagated uncertainty on the recoil energy measured on the tpc

        double erfluc = randGen.Gaus(er,deltaer); //recoil energy measured on the tpc with fluctuation given by the uncertainty
        double fcut = cutfprompt->Eval(psciobs); //fprompt selection cut calculation


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
            h2->Fill(erfluc,s2fluc);
            fpromptenergy2->Fill(psciobs,fprompt); 
        }
        else 
        {
        counter = counter + 1; //counting the number of discarded events
        }
        
    }
    cout<<"counter = "<<counter<<endl; //printing the number of discarded events


    TCanvas* canvando = new TCanvas("canvando","",1920,1080);
    h2->SetMarkerColor(kViolet);
    h2->SetMarkerStyle(20);
    h2->SetMarkerSize(1.1);
    stylehist2d(h2);
    gStyle->SetTitleFontSize(0.06);
    h2->Draw(drawopt); //drawing generated events, but yet without background

    if (strcmp(funcopt, "ziegler") == 0){
    s2func->SetLineColor(coloropt+2);
    styletf1(s2func);
    }
    else if (strcmp(funcopt, "moliere") == 0){
    s2moli->SetLineColor(coloropt+2);
    styletf1(s2moli);
    }
    else if (strcmp(funcopt, "lenz-jensen") == 0){
    s2lenz->SetLineColor(coloropt+2);
    styletf1(s2lenz);
    }

    l.SetTextColor(kBlack);
    l.SetTextSize(0.05);
    

    TCanvas *canvando1 = new TCanvas("canvando1","canvando1",800, 800);
    pscieff->SetTitle("Eficiencia - EJ-276; Energia (MeV); Eficiencia");
    pscieff->Draw(); //drawing the psci efficience interpolation being used

    
    double intstart[7] = {0.5, 2.5,3.5, 4.5, 5.5, 6.5, 7.5}; //energy interval beginings
    double intend[7] = {2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 10.}; //energy interval ends
    double percent[7] = {0.15, 0.05,0.14,0.05,0.08,0.12,0.21}; //background percentage for each energy interval

    TH1D* h2proj[7];
    TH2F* h2copy;
    h2copy = (TH2F*) h2->Clone("h2_clone");
    
	TCanvas* canvando2 = new TCanvas("canvando2","",1920,1080);
    stylehist1D(h1);
    h1->Draw(); //drawing fprompt distribution 

    TCanvas* fhist = new TCanvas("fhist","",1920,1080); // drawing fprompt vs. energy 2d histogram and fprompt-cut in function of energy for all events
    fpromptenergy->SetMarkerStyle(20);
    stylehist2d(fpromptenergy);
    fpromptenergy->SetMarkerSize(1.5);
    fpromptenergy->Draw();
    cutfprompt->SetLineColor(kRed);
    cutfprompt->Draw("SAME");

    TCanvas* canvando3 = new TCanvas("canvando3","",1920,1080); // drawing fprompt vs. energy 2d histogram and fprompt-cut in function of energy after selection cuts
    fpromptenergy2->SetMarkerStyle(20);
    stylehist2d(fpromptenergy2);
    fpromptenergy2->SetMarkerSize(1.5);
    fpromptenergy2->Draw();
    cutfprompt->SetLineColor(kRed);
    cutfprompt->Draw("SAME");
    

	//Now adding the background events

    TCanvas* canvas;
    double means[7] = {0.,0.,0.,0.,0.,0.,0.};
    double er2,s2ran;

    
    for (int i = 0; i<7; i++) //processing each energy interval individually 
    {
    std::string hname=Form("h2proj_%d",i);
	h2proj[i]=new TH1D(hname.c_str(),"",25,0.,80.);
    h2proj[i] = h2->ProjectionY("",h2->GetXaxis()->FindBin(intstart[i]), h2->GetXaxis()->FindBin(intend[i])); //Getting the projection in S2 for each slice of energy to calculate the number of background events
    h2proj[i]->SetYTitle("Counts/1.0 PE");
    h2proj[i]->SetXTitle("N_{i.e.}");
    stylehist1D(h2proj[i]);
    
    int entries = h2proj[i]->GetEntries();
    int bkg = percent[i]*entries; //calculating number of background events based on the percentage per slice. this percentage was extracted from real data estimatives
    
	for (int j = 0; j<bkg; j++) //adding bkg events
    {
        s2ran = randGen.Uniform(0.,1610.);
        er2 = randGen.Uniform(intstart[i],intend[i]);
        h2copy->Fill(er2,s2ran);
    }
    }

    canvas=new TCanvas("cdraw1","",1920,1080);
    h2copy->Draw(); //drawing the final sample


    MCMC_Fit_LenzJensen(h2copy); //applying the MCMC script to extract the fitted g2 value.

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

// ===================
//  MCMC FUNCTIONS 
// ===================

Double_t fnfunc_mcmc(Double_t *x, Double_t *par) //preparing for lenz-jensen stopping power calculation
{
    double eta = x[0];
    double m = par[0];
    double q = par[1];
    double lambda = par[2];
    return lambda * pow(eta, 1. - 2.*m) / pow(1. + pow(2.*lambda * pow(eta, 2.*(1.-m)), q), 1./q);
}

double CalculateNie_LenzJensen(double Er_keV, double cbox, double beta) //calculating the theoretical number of electrons extracted from LAr (thomas-imel + lenz-jensen)
{
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

double CalculateLogLikelihood2D_WithBkg(TH2F* h2_data, const vector<double>& Nie_cache, double g2_guess) //calculating loglikelihood 
{
    double logL = 0.0;
    int nBinsX = h2_data->GetNbinsX(); 
    int nBinsY = h2_data->GetNbinsY(); 
    
    // Bin width of Y axis (S2) to calculate the uniform background density
    double s2_bin_width = h2_data->GetYaxis()->GetBinWidth(1);

    for (int ix = 1; ix <= nBinsX; ix++) {
        double N_ie = Nie_cache[ix]; 
        if (N_ie <= 0) continue;
        
        // 1. Counts the total number of events (signal+background) in this energy interval
        double total_events_in_slice = 0;
        for (int iy = 1; iy <= nBinsY; iy++) {
            total_events_in_slice += h2_data->GetBinContent(ix, iy);
        }
        if (total_events_in_slice == 0) continue; 

        // 2. Estimative of the background fraction based on the uniform distrbution considering a background/signal ratio of 10%
        // The background has a uniform distribution from 0 to 1610 PE (total amplitude of S2):
        double S2_range_total = 1610.0; 
        double bkg_fraction = 0.10;
        double total_bkg_events = total_events_in_slice * bkg_fraction;
        double total_signal_events = total_events_in_slice - total_bkg_events;
        
        // Expected background density per S2 bin in this energy interval 
        double bkg_lambda_per_bin = (total_bkg_events / S2_range_total) * s2_bin_width;

        // 3. Signal 
        double mean_S2 = g2_guess * N_ie;
        double sigma_S2 = sqrt(mean_S2 + pow(0.001 * mean_S2, 2.0)); 
        
        // 4. Poisson likelihood sum per bin (signal + bakcground)
        for (int iy = 1; iy <= nBinsY; iy++) {
            double n_obs = h2_data->GetBinContent(ix, iy);
            double S2_val = h2_data->GetYaxis()->GetBinCenter(iy);
            
            // Expected singal normalized by the number of real signal events 
            double signal_lambda = total_signal_events * TMath::Gaus(S2_val, mean_S2, sigma_S2, true) * s2_bin_width;
            
            // Total number of events expected in this bin: (lambda) = signal + background
            double lambda = signal_lambda + bkg_lambda_per_bin;
            
            if (lambda < 1e-6) lambda = 1e-6; 
            
            logL += n_obs * log(lambda) - lambda;
        }
    }
    return logL;
}


void MCMC_Fit_LenzJensen(TH2F* h2_to_fit) // Metropolis-Hastings MCMC
{ 
    double cbox = 810.0;
    double beta = 6800.0;
    int nBinsX = h2_to_fit->GetNbinsX();
    std::vector<double> Nie_cache(nBinsX + 1, 0.0);
    
    cout << "Pre-calculating Lenz-Jensen model integrals... Wait a little bit." << endl;
    for (int ix = 1; ix <= nBinsX; ix++) {
        double E_val = h2_to_fit->GetXaxis()->GetBinCenter(ix);
        Nie_cache[ix] = CalculateNie_LenzJensen(E_val, cbox, beta);
    }
    cout << "Calculation completed! Now iniciating the MCMC..." << endl;
    
    TRandom3 randGen(0); 
    int n_steps = 20000;          
    int burn_in = 0;           
    int lag = 10;
    double step_size = 0.25;      
    double current_g2 = 16.0;     
    
    TH1D* h_posterior = new TH1D("h_posterior", "Posterior distribution of g2; g2 (PE/e-); Freq", 100, 16.5, 18.5);
    h_posterior->SetDirectory(0); // <--- Avoid ROOT from deleting our histogram :)
    
    TGraph* gr_trace = new TGraph();
    gr_trace->SetTitle("Trace Plot; Step; g2 value");
    
    std::vector<double> full_chain; // <--- Saves the g2 values to calculate the ACF

    double current_logL = CalculateLogLikelihood2D_WithBkg(h2_to_fit, Nie_cache, current_g2);
    int accepted_steps = 0;

    for (int step = 0; step < n_steps; step++) {
        if (step % 2000 == 0) cout << "MCMC step: " << step << " / " << n_steps << endl;
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
            full_chain.push_back(current_g2); // Saving values to calculate ACF
            
            // Applies the lag factor and fills the histogram only with the selected values
            if ((step - burn_in) % lag == 0) {
                h_posterior->Fill(current_g2);
            }
        }
        gr_trace->SetPoint(step, step, current_g2);
    }
    
    
	// AUTOCORRELATION FUNCTION CALCULATION	

    int n_chain = full_chain.size();
    double mean_chain = 0.0;
    for (double val : full_chain) mean_chain += val;
    mean_chain /= n_chain; //g2 mean
    
    double var_chain = 0.0;
    for (double val : full_chain) var_chain += (val - mean_chain) * (val - mean_chain);
    var_chain /= n_chain; //g2 variance
    
    TGraph* gr_acf = new TGraph();
    gr_acf->SetTitle("Autocorrelation function (ACF); Lag k (steps); ACF");
    
    int max_lag = 100; // Analysing the autocorrelation up to 100 steps
    for (int k = 0; k <= max_lag; k++) {
        double acf_k = 0.0;
        for (int i = 0; i < n_chain - k; i++) {
            acf_k += (full_chain[i] - mean_chain) * (full_chain[i+k] - mean_chain);
        }
        acf_k = acf_k / n_chain; 
        acf_k = acf_k / var_chain; // Normaliza para ter ACF = 1 no lag 0
        gr_acf->SetPoint(k, k, acf_k);
    }


    //  Ploting everything

    cout << "------------------------------------------" << endl;
    cout << "Acceptance ratio: " << (double)accepted_steps / n_steps * 100.0 << " %" << endl;
    cout << "estimated g2 (Lenz-Jensen + Thomas-Imel + Background ): " << h_posterior->GetMean() << " +/- " << h_posterior->GetStdDev() << endl;
    
    
    TCanvas* c_mcmc = new TCanvas("c_mcmc", "MCMC parameters", 1500, 500);
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
    TCanvas* c_fit2d = new TCanvas("c_fit2d", "MCMC fit result", 1000, 800);
    h2_to_fit->SetTitle(" MCMC fit result; Recoil Energy (keV); S2 (PE)");
    h2_to_fit->Draw("");
    
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
    leg_mcmc2d->AddEntry(h2_to_fit, "Sample with background", "lep");
    leg_mcmc2d->AddEntry(f_mean_s2, Form("Lenz-Jensen model (g_{2} = %.2f)", g2_best), "l");
    leg_mcmc2d->Draw("SAME");
}
