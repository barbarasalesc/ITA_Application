#include "BlueAnalyzer.hh"
#include "Root_Plot.cc"
#include "RDTrace.hh"
#include "RDChannel.hh"


double GetEnergy(int ipsci, double Q )
{
  if ( ipsci <0  || ipsci >17 ) return UNDEF;

  const double EnergyConstant_PSci_250sa[18]={  // Done by Maximo for wlong=250 with Am241 (ADC/energy)
    18.918,18.780,19.016,     20.144,18.200,17.296,       18.422,18.828,19.309,
    20.382,18.401,16.948,     17.909,22.833,15.586,       17.625,19.092,19.800  };

  return Q/EnergyConstant_PSci_250sa[ipsci];
  
}


void Process()
{  

	TFile *myRootFile = new TFile("/home/barbaras/red-deconv/barb/mytreefile_v2.root", "RECREATE");
    TTree *t1 = new TTree("treetraces", "treetraces");

	  std::vector<float> energy;
  	  int chanID;
  	  std::string chanType;
	  int sa_start,sa_end;

	  
	  t1->Branch("energy",&energy);
	  t1->Branch("chanID",&chanID);
	  t1->Branch("chanType",&chanType);
	  t1->Branch("sa_start",&sa_start);
	  t1->Branch("sa_end",&sa_end);

  std::ifstream file("analysis/selection.runs",std::ifstream::in);
  while (true)
    {
      int run;
      
      file>> run;
      if ( file.eof() ) break;

      //std::string BaseDir=RedHome;
      RDconfig::GetInstance()->Init(run,1,"",false);
   
      //-----------------------------------------------------------------------------------------------------------------

	  TTree *t;
	  
	  int type = 1; //PostFlash -- since we need gammas and neutrons

	  std::string filename;
	  filename=Form("results.red/fom/run_%d_1_%d.root",run,type);

	  TFile *f=new TFile(filename.c_str() );
	  if ( f->IsZombie() ) { printf("RedFile does not exist \n"); continue; }
	  f->GetObject("reco",t);

	  RedEvent* er=new RedEvent();
	  t->SetBranchAddress("recoevent",&er);

	  RDTrace* aTrace=new RDTrace();
	

	  printf("Processing run %d   | type=%d N=%d \n",run,type,(int) t->GetEntries() );
	
	  //-----------------------------------------------------------------------------------------------------------------
	  int nevt_sel=0;
	  for (int ientry=0;ientry<t->GetEntries();ientry++)
	    {
	      if ( ientry%1000==0 ) printf("Entry %d \n",ientry);
	      t->GetEntry(ientry);

	      //-----------------------------
	      //---FOM/Traces with traces Flash/PostFlash
	      //-----------------------------
	      if ( er->GetNTraces()==0 ) continue;
	      nevt_sel++;
	      
	      for (int ii=0;ii<er->GetNTraces();ii++)
		{
  	  	  chanID = 0;
  		  chanType.clear();
	  	  sa_start = 0;
		  sa_end = 0;
		  aTrace  =er->GetTrace(ii);
		  if (aTrace->chanType!="psci")  continue;
		  chanID = aTrace->chan;
  	  	  chanType = aTrace->chanType;
	  	  sa_start = aTrace-> sa_start;
		  sa_end = aTrace->sa_end;
		  energy.clear();
		  energy.resize(aTrace->adc.size());
		  for (int sa=0;sa<(int) aTrace->adc.size();sa++)
    	  {
		    energy.at(sa) = GetEnergy(chanID,aTrace->adc.at(sa));
		  }
		  t1->Fill();
		}// Loop over traces
	    }// Loop over entries
	  printf("Number of events with traces %d\n",nevt_sel);
		f->cd();
    	f->Close();
		delete f;
    }//Loop over event numbers
	myRootFile->cd();
	t1->Write();
	myRootFile->Close();
	delete myRootFile; 
}


