#ifndef __Run_C__
#define __Run_C__

#include "Run.h"

using namespace std;

/*---------------------------------------- 
   PID == 0 , Gamma Events
   PID == 1 , Beta like Events
   PID == 2 , Muon like Events
------------------------------------------*/

Run::Run(Int_t n,Int_t m,TSQLServer *sql)
{
  cout << "Attempting to open run " << n << endl;
  runnum      = n;    // Set Run Number
  runtype     = m;    // Set Run Type
  toverflow   = 0;    // Overflow counter set to zero.
  stupid      = 0;

  GammasEast  = 0.; // Set Gammas Counters to zero
  GammasWest  = 0.;
  GammasEastg = 0.;
  GammasWestg = 0.;

  passed1 = 0;      // Counter for Debugging Event selection cuts
  passed2 = 0;
  passed3 = 0;
  passed4 = 0;
  
  rotated = 3;
  octtype = new char[4];

  Hlist     = new TObjArray();
  HEastAn   = new TObjArray();
  HWestAn   = new TObjArray();
  RunDate   = new TDatime();
  SetReconEnergy(); // Set the Reconstruction parameters based on PENELOPE
  flipperOn = SetRunTime(sql);  // Connect to the database and get the run time
}

Run::~Run()
{
  
}

void Run::DeleteHistos()
{
  // Remove the histograms from the TObjArray's prior to deleting the pointers.
  HEastAn->RemoveRange(0,HEastAn->GetEntriesFast());
  HWestAn->RemoveRange(0,HWestAn->GetEntriesFast());
  // Delete the histogram pointers used by the class.
  delete hpe;
  delete hpw;
  delete heq;
  delete hwq;
//  delete heqB;
  delete heqC;
//  delete heqD;
  delete heqE;
  delete heqF;
  delete heqG;
//  delete heqH;
//  delete heqI;
//  delete hwqB;
  delete hwqC;
//  delete hwqD;
  delete hwqE;
  delete hwqF;
  delete hwqG;
//  delete hwqH;
//  delete hwqI;
  for(Int_t i = 0 ; i < 10 ; i++){
    delete hEHC[i];
    delete hWHC[i];
  }
  delete hte;
  delete htw;
  delete hApose;
  delete hAposw;
  delete hPsDfET23;
/*  delete hw1;
  delete hw2;
  delete hw3;
  delete hw4;
  delete he1;
  delete he2;
  delete he3;
  delete he4;
*/  delete hecy;
  delete hwcy;
  delete hEERef;
  delete hEERef1;
  delete hEERef2;
  delete hEWRef;
  delete hEWRef1;
  delete hEWRef2;
  delete hEmuon;
  delete hWmuon;
  delete hESigNos;
  delete hWSigNos;
  delete hENoMWPC;
  delete hWNoMWPC;
  delete hEMWPC;
  delete hWMWPC;
  delete hEType0_Multi;
  delete hEType1_Multi;
  delete hEType23_Multi;
  delete hWType0_Multi;
  delete hWType1_Multi;
  delete hWType23_Multi;
  delete hEtype_1;
  delete hEtype_23;
  delete hWtype_1;
  delete hWtype_23;
  for(Int_t i = 0 ; i < 12 ;i++){
    delete hERad[i];
    delete hWRad[i];
  }
  delete hTDCDiff;
  delete hRote;
  delete hRotw;
  delete hRoteI;
  delete hRotwI;
  delete hRote23;
  delete hRotw23;
  delete hRoteI23;
  delete hRotwI23;
  delete hPosDiffShifted;
  delete hPosDiffUnShifted;
  delete hEastType1XRot;
  delete hEastType1YRot;
  delete hEastType23XRot;
  delete hEastType23YRot;
  delete hWestType1XRot;
  delete hWestType1YRot;
  delete hWestType23XRot;
  delete hWestType23YRot;
  delete hETimeVsE;
  delete hETimeVsEP;
  delete hETimeVsES;
  delete hWTimeVsE;
  delete hWTimeVsEP;
  delete hWTimeVsES;
//  delete htdcE;
//  delete htdcW;
  delete hEAnode23;
  delete hWAnode23;
  delete hEType1_2d;
  delete hWType1_2d;
  delete hEType1_Primary;
  delete hEType1_Secondary;
  delete hWType1_Primary;
  delete hWType1_Secondary;
  delete hE23Anode2d;
  delete hW23Anode2d;
  delete hGammaCounts;
  delete hGammaCountsg;
  delete hETDC_cor;
  delete hETDC_cor_passed;
  delete hEvWTDC;
  delete hEvWTDC_cor;
  delete hEAnode;
  delete hECathS;
  delete hWAnode;
  delete hWCathS;
  delete c1;
  delete c2;
  delete hEKurie;
  delete hWKurie;
  delete xr; 
  delete hmr1;
  delete hCountTimeRecord;
}

Bool_t Run::OpenRun(Int_t n)
{
  // See the Bck_Run constructor for the relavent comments.
  const char* filepath = getenv("UCNAOUTPUTDIR"); 
  const char* author   = getenv("UCNA_ANA_AUTHOR");
  cout << "Attempting to open run spec_"<< n << ".root" << endl;
  f1     = new TFile(Form("%s/hists/spec_%d.root",filepath,n),"READ");
  //fHisto = new TFile(Form("%s/hists/hists_%d.root",filepath,n),"UPDATE");
  fHisto = new TFile(Form("/home/jwwexler/robbyWork/histOut/hists_%d.root",n),"UPDATE");
  
  AnalysisDirExist = kFALSE;
  if(!(f1->IsOpen())){
	cout << " Failed to open file "<<endl;
	return kFALSE;
  }
  // Get the tree..
  t1 = (TTree*)f1->Get("phys");
  
  if(!(fHisto->cd(Form("%s_Analysis_%d_%d",author,(Int_t)RAD,rotated)))){
    fHisto->mkdir(Form("%s_Analysis_%d_%d",author,(Int_t)RAD,rotated));
	cout << "Making Analysis Histogram" << endl;
  } else {
    fHisto->cd(Form("%s_Analysis_%d_%d",author,(Int_t)RAD,rotated));
    AnalysisDirExist = kTRUE;
	cout << "Reopening Analysis Histogram" << endl;
  }
  
  SetBranches();
  return kTRUE;
}
//--------------------------------------------------------------------------------
Int_t Run::Scale2Time(Int_t nex,Int_t nwx)
{
    cout << "Scaling to Run time " << endl;
    cout << "Counts East " << heqC->Integral(nlow,nhigh) << endl;
    cout << "Counts West " << hwqC->Integral(nlow,nhigh) << endl;
    
    ScaleList(HEastAn,rtime_e);
    ScaleList(HWestAn,rtime_w);
  
    //----------------------------------------------------------
    for(Int_t i = 1 ; i <= hEastType1XRot->GetNbinsX(); i+=2){
       for(Int_t j = 1 ; j <= hEastType1XRot->GetNbinsY() ; j+=2){
          Double_t eastevents = hRote->GetBinContent(i,j) + hRote->GetBinContent(i+1,j)+ hRote->GetBinContent(i,j+1)
                              + hRote->GetBinContent(i+1,j+1);
          Double_t westevents = hRotw->GetBinContent(i,j) + hRotw->GetBinContent(i+1,j)+ hRotw->GetBinContent(i,j+1)
                              + hRotw->GetBinContent(i+1,j+1);
          if(eastevents > 0){
            hEastType1XRot->SetBinContent(i,j,hEastType1XRot->GetBinContent(i,j)/eastevents);
            hEastType1YRot->SetBinContent(i,j,hEastType1YRot->GetBinContent(i,j)/eastevents);
          }
          if(westevents > 0){
            hWestType1XRot->SetBinContent(i,j,hWestType1XRot->GetBinContent(i,j)/westevents);
            hWestType1YRot->SetBinContent(i,j,hWestType1YRot->GetBinContent(i,j)/westevents);
          }

          eastevents = hRote23->GetBinContent(i,j) + hRote23->GetBinContent(i+1,j)+ hRote23->GetBinContent(i,j+1)
                     + hRote23->GetBinContent(i+1,j+1);
          westevents = hRotw23->GetBinContent(i,j) + hRotw23->GetBinContent(i+1,j)+ hRotw23->GetBinContent(i,j+1)
                     + hRotw23->GetBinContent(i+1,j+1);

          if(eastevents > 0){
            hEastType23XRot->SetBinContent(i,j,hEastType23XRot->GetBinContent(i,j)/eastevents);
            hEastType23YRot->SetBinContent(i,j,hEastType23YRot->GetBinContent(i,j)/eastevents);
          }
          if(westevents > 0){
            hWestType23XRot->SetBinContent(i,j,hWestType23XRot->GetBinContent(i,j)/westevents);
            hWestType23YRot->SetBinContent(i,j,hWestType23YRot->GetBinContent(i,j)/westevents);
          }
       }
    }

    T2_frac = hPsDfET23->Integral(40,60,50,70) / hPsDfET23->Integral(1,120,1,120);
  return 0;
}

void Run::CreateTH1F(TH1F *&h,const char* name,const char* title,Int_t nbins,Double_t xmin,Double_t xmax,Int_t nSide)
{

  h = new TH1F(Form("%s_%d",name,GetRunNumber()),title,nbins,xmin,xmax);
  Hlist->Add(h);

  if(nSide == 0)
     HEastAn->Add(h);
  else if(nSide==1)
     HWestAn->Add(h);
}

void Run::CreateTH2F(TH2F *&h,const char* name,const char* title,Int_t nbins,Double_t xmin,Double_t xmax,Int_t nybins,Double_t ymin,Double_t ymax,Int_t nSide)
{

  h = new TH2F(Form("%s_%d",name,GetRunNumber()),title,nbins,xmin,xmax,nybins,ymin,ymax);
  Hlist->Add(h);

  if(nSide == 0)
     HEastAn->Add(h);
  else if(nSide == 1)
     HWestAn->Add(h);

  return;
}

//---------------------------------------------------------------------
Int_t Run::Initialize_hist(Int_t n,Int_t nex,Int_t nwx)
{
  // This is a very long and boring function that just defines the histograms
  Int_t nEast = 0; 
  Int_t nWest = 1.;
  
  xr = new TRandom3();
  
  CreateTH1F(hEAnode23 ,"hEAnode23","East Trigger Side Anode 23;Energy(keV);Counts",100,0,20,nEast);
  CreateTH1F(hWAnode23 ,"hWAnode23","West Trigger Side Anode 23;Energy(keV);Counts",100,0,20,nWest);

  CreateTH2F(hpw,"hpw","West Position;X;Y",150,-75,75,150,-75,75,nWest);
  CreateTH2F(hpe,"hpe","East Position;X;Y",150,-75,75,150,-75,75,nEast);

  CreateTH2F(hEastType1XRot,"hEastType1XRot" ,"X Displacement vs. Position for East Type 1; X (mm) ; Y (mm)" ,75,-75,75,75,-75,75,nEast);
  CreateTH2F(hEastType1YRot,"hEastType1YRot" ,"Y Displacement vs. Position for East Type 1; X (mm) ; Y (mm)" ,75,-75,75,75,-75,75,nEast);
  CreateTH2F(hEastType23XRot,"hEastType23XRot","X Displacement vs. Position for East Type 23; X (mm) ; Y (mm)",75,-75,75,75,-75,75,nEast);
  CreateTH2F(hEastType23YRot,"hEastType23YRot","Y Displacement vs. Position for East Type 23; X (mm) ; Y (mm)",75,-75,75,75,-75,75,nEast);

  CreateTH2F(hWestType1XRot,"hWestType1XRot" ,"X Displacement vs. Position for West Type 1; X (mm) ; Y (mm)" ,75,-75,75,75,-75,75,nWest);
  CreateTH2F(hWestType1YRot,"hWestType1YRot" ,"Y Displacement vs. Position for West Type 1; X (mm) ; Y (mm)" ,75,-75,75,75,-75,75,nWest);
  CreateTH2F(hWestType23XRot,"hWestType23XRot","X Displacement vs. Position for West Type 23; X (mm) ; Y (mm)",75,-75,75,75,-75,75,nWest);
  CreateTH2F(hWestType23YRot,"hWestType23YRot","Y Displacement vs. Position for West Type 23; X (mm) ; Y (mm)",75,-75,75,75,-75,75,nWest);

  CreateTH2F(hPosDiffShifted,"hPosDiffShifted","Shifted Displacement; #Delta X (mm) ; #Delta Y (mm)"    ,150,-25,25,150,-25,25,-1);
  CreateTH2F(hPosDiffUnShifted,"hPosDiffUnShifted","UnShifted Displacement; #Delta X (mm) ; #Delta Y (mm)",150,-25,25,150,-25,25,-1);

 CreateTH1F(hEType0_Multi ,"hEType0_Multi","East Type 0 Multiplicity; Multiplicity;Cnts",20,0,20,nEast);
 CreateTH1F(hEType1_Multi ,"hEType1_Multi","East Type 1 Multiplicity; Multiplicity;Cnts",20,0,20,nEast);
 CreateTH1F(hEType23_Multi,"hEType23_Multi","East Type 23 Multiplicity; Multiplicity;Cnts",20,0,20,nEast);
 CreateTH1F(hWType0_Multi ,"hWType0_Multi","West Type 0 Multiplicity; Multiplicity;Cnts",20,0,20,nWest);
 CreateTH1F(hWType1_Multi ,"hWType1_Multi","West Type 1 Multiplicity; Multiplicity;Cnts",20,0,20,nWest);
 CreateTH1F(hWType23_Multi,"hWType23_Multi","West Type 23 Multiplicity; Multiplicity;Cnts",20,0,20,nWest);

 CreateTH1F( hGammaCounts  ,"hGammaCount","Gammas",2,0,2,-1);
 CreateTH1F( hGammaCountsg ,"hGammaCountg","Gammas",2,0,2,-1);	

 CreateTH1F(hENoMWPC ,"hENoMWPC","No MWPC cut ;Rate s^{-1} ; Energy (keV)",EBINS,0,ESCALE,nEast);
 CreateTH1F(hWNoMWPC ,"hWNoMWPC","No MWPC cut ;Rate s^{-1} ; Energy (keV)",EBINS,0,ESCALE,nWest);
 CreateTH1F(hEMWPC   ,"hEMWPC","MWPC cut ;Rate s^{-1} ; Energy (keV)",EBINS,0,ESCALE,nEast);
 CreateTH1F(hWMWPC   ,"hWMWPC","MWPC cut ;Rate s^{-1} ; Energy (keV)",EBINS,0,ESCALE,nWest);
 
  CreateTH2F(hETimeVsE,"hETimeVsE","East Timing vs. Energy",100,0.,1000.,100,0.,200.,nEast);
  CreateTH2F(hWTimeVsE,"hWTimeVsE","West Timing vs. Energy",100,0.,1000.,100,0.,200.,nWest);

 CreateTH1F(hESigNos ,"hESigNos","East Signal to Noise ; S / N ; Energy (keV)",EBINS,0,ESCALE,nEast);
 CreateTH1F(hWSigNos ,"hWSigNos","West Signal to Noise ; S / N ; Energy (keV)",EBINS,0,ESCALE,nWest);

 CreateTH1F(hEKurie ,"hEKurie","East Kurie Plot; E_{vis} (keV) ; #sqrt{N(E_{vis}) / F(Z,E) p^{2}(E_{vis})}",EBINS,0,ESCALE,nEast);
 CreateTH1F(hWKurie ,"hWKurie","West Kurie Plot; E_{vis} (keV) ; #sqrt{N(E_{vis}) / F(Z,E) p^{2}(E_{vis})}",EBINS,0,ESCALE,nWest);
 
 CreateTH2F(hPsDfET23,"hPsDfET23","Position Difference for Type 2/3",120,-60,60,120,-60,60,-2);

 CreateTH2F(hApose,"hApose","East Rate by Quadrant",2,-60,60,2,-60,60,nEast);
 CreateTH2F(hAposw,"hAposw","West Rate by Quadrant",2,-60,60,2,-60,60,nWest);
  
 CreateTH2F(hEvWTDC_cor,"hEvWTDC_cor","East vs. West TDC !BkhfGood; East (ns) ; West (ns)",200,0,200,200,0,200,-2);  
 CreateTH2F(hEvWTDC,"hEvWTDC","East vs. West TDC BkhfGood; East (ns) ; West (ns)",200,0,200,200,0,200,-2);

 CreateTH1F(hEAnode ,"hEAnode","East Anode Spectrum",100,0,20,nEast);
 CreateTH1F(hECathS ,"hECathS","East Cathode Sum Spectrum",5000,-1000,4000,nEast);
 CreateTH1F(hWAnode ,"hWAnode","West Anode Spectrum",100,0,20,nWest);
 CreateTH1F(hWCathS ,"hWCathS","West Cathode Sum Spectrum",5000,-1000,4000,nWest);
		     
 CreateTH2F(hE23Anode2d,"hE23Anode2d","East Type 23 E vs. W Anode; East Anode ; West Anode",100,0,3000,100,0,3000,nEast);
 CreateTH2F(hW23Anode2d,"hW23Anode2d","West Type 23 E vs. W Anode; East Anode ; West Anode",100,0,3000,100,0,3000,nWest);
  
  // Type 1 backscatter positions and rotations
  //------------------------------------------------------------------------------------------------------------------
  CreateTH1F(hTDCDiff ,"hTDCDiff","#Delta TDC ; TDCE - TDCW ; Counts",2000,-200,1800,-2);
  CreateTH2F(hRote,"hRote","East Type 1 Backscattering Rotation ; X_{east} - X_{west} (mm) ; Y_{east} - Y_{west} (mm)",150,-75,75,150,-75,75,nEast);
  CreateTH2F(hRotw,"hRotw","West Type 1 Backscattering Rotation ; X_{east} - X_{west} (mm) ; Y_{east} - Y_{west} (mm)",150,-75,75,150,-75,75,nWest);
  CreateTH2F(hRoteI,"hRoteI","East Type 1 Backscattering Rotation ; X_{east} - X_{west} (mm) ; Y_{east} - Y_{west} (mm)",150,-75,75,150,-75,75,nEast);
  CreateTH2F(hRotwI,"hRotwI","West Type 1 Backscattering Rotation ; X_{east} - X_{west} (mm) ; Y_{east} - Y_{west} (mm)",150,-75,75,150,-75,75,nWest);
  // Type 1 backscatter positions and rotations
  //------------------------------------------------------------------------------------------------------------------
  CreateTH2F(hRote23,"hRote23","East Type 23 Backscattering Rotation ; X_{east} - X_{west} (mm) ; Y_{east} - Y_{west} (mm)",
		     150,-75,75,150,-75,75,nEast);
  
  CreateTH2F(hRotw23,"hRotw23","West Type 23 Backscattering Rotation ; X_{east} - X_{west} (mm) ; Y_{east} - Y_{west} (mm)",
		     150,-75,75,150,-75,75,nWest);
  
  CreateTH2F(hRoteI23,"hRoteI23","East Type 32 Backscattering Rotation ; X_{east} - X_{west} (mm) ; Y_{east} - Y_{west} (mm)",
		      150,-75,75,150,-75,75,nEast);
  
  CreateTH2F(hRotwI23,"hRotwI23","West Type 23 Backscattering Rotation ; X_{east} - X_{west} (mm) ; Y_{east} - Y_{west} (mm)",
		      150,-75,75,150,-75,75,nWest);
  
  // West Scintillator energy
  //-------------------------------------------------------------------------------------------------------------------
 CreateTH1F(hw1 ,"hw1","Qadc4; Channel #;Cts",410,0,ESCALE,nWest);
 CreateTH1F(hw2 ,"hw2"," Qadc5; Channel #;Cts",410,0,ESCALE,nWest);
 CreateTH1F(hw3 ,"hw3","Qadc6; Channel #;Cts",410,0,ESCALE,nWest);
 CreateTH1F(hw4 ,"hw4"," Qadc7; Channel # ; Cts",410,0,ESCALE,nWest);
 CreateTH1F(hwq ,"hwq","West Energy ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nWest);
  
 CreateTH1F(hwqB ,"hwqB","West Energy B ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nWest);
 CreateTH1F(heqB ,"heqB","East Energy B ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nEast);
 CreateTH1F(hwqC ,"hwqC","West Energy C ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nWest); 
 CreateTH1F(heqC ,"heqC","East Energy C ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nEast);
 CreateTH1F(hwqD ,"hwqD","West Energy D ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nWest);
 CreateTH1F(heqD ,"heqD","East Energy D ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nEast);
 CreateTH1F(hwqE ,"hwqE","West Energy E ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nWest); 
 CreateTH1F(heqE ,"heqE","East Energy E ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nEast);
 CreateTH1F(hwqF ,"hwqF","West Energy F ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nWest);
 CreateTH1F(heqF ,"heqF","East Energy F ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nEast);
 CreateTH1F(hwqG ,"hwqG","West Energy G ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nWest); 
 CreateTH1F(heqG ,"heqG","East Energy G ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nEast);
 CreateTH1F(hwqH ,"hwqH","West Energy H ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nWest);
 CreateTH1F(heqH ,"heqH","East Energy H ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nEast);
 CreateTH1F(hwqI ,"hwqI","West Energy I ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nWest); 
 CreateTH1F(heqI ,"heqI","East Energy I ; E_{vis} (keV); Cts",EBINS,0,ESCALE,nEast);

 CreateTH1F(hEType1_Primary  ,"hEType1_Primary",";Primary Detector E_{vis} (keV); Cts",100,0,1,nEast);
 CreateTH1F(hEType1_Secondary,"hEType1_Secondary",";Primary Detector E_{vis} (keV); Cts",100,0,1,nEast);
 CreateTH1F(hWType1_Primary  ,"hWType1_Primary",";Primary Detector E_{vis} (keV); Cts",100,0,1,nWest);
 CreateTH1F(hWType1_Secondary,"hWType1_Secondary",";Primary Detector E_{vis} (keV); Cts",100,0,1,nWest);

 CreateTH2F(hETimeVsEP,"hETimeVsEP","East Timing vs. Energy",100,0.,1000.,100,0.,200.,nEast);
 CreateTH2F(hWTimeVsEP,"hWTimeVsEP","West Timing vs. Energy",100,0.,1000.,100,0.,200.,nWest);
 CreateTH2F(hETimeVsES,"hETimeVsES","East Timing vs. Energy",100,0.,1000.,100,0.,200.,nEast);
 CreateTH2F(hWTimeVsES,"hWTimeVsES","West Timing vs. Energy",100,0.,1000.,100,0.,200.,nWest);

 CreateTH1F(hWmuon ,"hWmuon","West #mu^{-} ; E_{vis} (keV); Cts",40,0,ESCALE,nWest);
 CreateTH1F(hEmuon ,"hEmuon","East #mu^{-} ; E_{vis} (keV); Cts",40,0,ESCALE,nEast);
  
 CreateTH1F(hEERef ,"hEERef","Reference ; Channel #;Cts",EBINS,0,ESCALE,-2);
 CreateTH1F(hEERef1 ,"hEERef1","Reference ; Channel #;Cts",EBINS,0,ESCALE,-2);
 CreateTH1F(hEERef2 ,"hEERef2","Reference ; Channel #;Cts",EBINS,0,ESCALE,-2);
 CreateTH1F(hEWRef ,"hEWRef","Reference ; Channel #;Cts",EBINS,0,ESCALE,-2);
 CreateTH1F(hEWRef1 ,"hEWRef1","Reference ; Channel #;Cts",EBINS,0,ESCALE,-2);
 CreateTH1F(hEWRef2 ,"hEWRef2","Reference ; Channel #;Cts",EBINS,0,ESCALE,-2);
  
 CreateTH1F(hETDC_cor ,"hETDC_cor","Corrupt TDC Energy; Channel #;Cts",EBINS,0,ESCALE,nEast);
 CreateTH1F(hWTDC_cor ,"hWTDC_cor","West Corrupt TDC Enegry; Channel #;Cts",EBINS,0,ESCALE,nWest);
 CreateTH1F(hETDC_cor_passed,"hETDC_cor_passed","Corrupt TDC Energy; Channel #;Cts",EBINS,0,ESCALE,nEast); 
 CreateTH1F(hWTDC_cor_passed,"hWTDC_cor_passed","West Corrupt TDC Enegry; Channel #;Cts",EBINS,0,ESCALE,nWest);
 CreateTH1F(htw,"htw","Timing West; Time ; Cts",2000,0,200,nWest);
 CreateTH1F(hWtype_1,"hWtype_1","Type 1 Events",EBINS,0,ESCALE,nWest);
 CreateTH1F(hWtype_23,"hWtype_23","Type 23 Events",EBINS,0,ESCALE,nWest);

 CreateTH1F(he1,"he1"," for Qadc0; Channel #;Cts",410,0,ESCALE,nEast);
 CreateTH1F(he2,"he2"," for Qadc1; Channel #;Cts",410,0,ESCALE,nEast);
 CreateTH1F(he3,"he3"," for Qadc2; Channel #;Cts",410,0,ESCALE,nEast);
 CreateTH1F(he4,"he4"," for Qadc3; Channel #;Cts",410,0,ESCALE,nEast);
 CreateTH1F(heq,"heq","East Energy ; E_{vis} (keV) ;Cts",EBINS,0,ESCALE,nEast);
  
 CreateTH1F(hte       ,"hte","Timing East for ; Time ; Cts",2000,0,200,nEast);
 CreateTH1F(hEtype_1  ,"hEtype_1","Type 1 Events for ",EBINS,0,ESCALE,nEast);
 CreateTH1F(hEtype_23 ,"hEtype_23" ,"Type 23 Events for ",EBINS,0,ESCALE,nEast);
 CreateTH2F(hEType1_2d ,"hEType1_2d","Type 1 event rates by quadrant",4,-60,60,4,-60,60,nEast);
 CreateTH2F(hWType1_2d ,"hWType1_2d","Type 1 event rates by quadrant",4,-60,60,4,-60,60,nWest);

  for(Int_t i = 0 ; i < 10 ; i++){
   CreateTH1F(hEHC[i] ,Form("hEHC[%d]",i),"Hard Cut",EBINS,0,ESCALE,-2);
   CreateTH1F(hWHC[i] ,Form("hWHC[%d]",i),"Hard Cut",EBINS,0,ESCALE,-2);
  }

  for(Int_t i = 0 ; i < 12 ; i++){
   CreateTH1F(hERad[i] ,Form("hERad[%d]",i),"Radius",EBINS,0,ESCALE,-2);
   CreateTH1F(hWRad[i] ,Form("hWRad[%d]",i),"Radius",EBINS,0,ESCALE,-2);
  }
  
  CreateTH1F(hmr1,"hmr1","UCN Monitor 2 Rate ; Time ; Rate",1000,0,4000,-2);
  CreateTH1F(hCountTimeRecord,"hCountTimeRecord","Count Times ; Index ; Counts",10,0,9,-2);

  return 0;
}

void Run::ScaleTDC()
{
  TDCE *= fEast_TDC_Scale;
  TDCW *= fWest_TDC_Scale;
};

Int_t Run::Count_Time_All()
{
  if(Side == 0 && CountTimeEAll==0){
    CountTimeEFirst=TDCE-3000;
    CountTimeEAll++;
  }
  else if(Side==0) {
    CountTimeEAll=TDCE-3000;
  }

  if(Side == 1 && CountTimeWBeta==0){
    CountTimeWFirst=TDCW-3000;
    CountTimeWAll++;
  }
  else if(Side==1) {
    CountTimeWAll=TDCW-3000;
  }

  hCountTimeRecord->SetBinContent(1,CountTimeEFirst);
  hCountTimeRecord->SetBinContent(2,CountTimeEAll);
  hCountTimeRecord->SetBinContent(5,CountTimeWFirst);
  hCountTimeRecord->SetBinContent(6,CountTimeWAll); 

  return -1;
}

Int_t Run::Count_Time_Beta()
{
  if(Side == 0 && CountTimeEBeta==0){
    CountTimeEFirstBeta=TDCE-3000;
    CountTimeEBeta++;
  }
  else if(Side==0) {
    CountTimeEBeta=TDCE-3000;
  }

  if(Side == 1 && CountTimeWBeta==0){
    CountTimeWFirstBeta=TDCW-3000;
    CountTimeWBeta++;
  }
  else if(Side==1) {
    CountTimeWBeta=TDCW-3000;
  }

  hCountTimeRecord->SetBinContent(3,CountTimeEFirstBeta);
  hCountTimeRecord->SetBinContent(4,CountTimeEBeta);
  hCountTimeRecord->SetBinContent(7,CountTimeWFirstBeta);
  hCountTimeRecord->SetBinContent(8,CountTimeWBeta); 

  return -1;
}

Int_t Run::Count_Gammas()
{
  if(Side == 0 && PID == 0) GammasEast++;
  if(Side == 1 && PID == 0) GammasWest++;
  return -1;
}

Int_t Run::Count_Gammas_Good()
{
  if(Side == 0 && PID == 0) GammasEastg++;
  if(Side == 1 && PID == 0) GammasWestg++;
  return -1;
}

Int_t Run::Set_Anode_Cut(Int_t n)
{
  /* ---------------------------------------------------------------
  /  The PADC spectrum from the MWPC is drawn with no cuts and the
  /  pedestal fit to a guassian.  The cut is set as the mean of 
  /  of the Gaussian plus nsig*sigma of the fit. 
  /  Get the Cathode Pedestal values.
  -----------------------------------------------------------------*/
  
  nsig = 3.;

  char *db = Form("mysql://%s/%s",getenv("UCNADBADDRESS"),getenv("UCNADB"));
  char *dbuser = Form("%s",getenv("UCNADBUSER"));
  char *dbpwd  = Form("%s",getenv("UCNADBPASS"));
  // Make a connection to the server at caltech
  //if(!(sql2->IsConnected()))ReConnect(sql2);
/*
  TSQLServer *sql2 = TSQLServer::Connect(db,dbuser,dbpwd);

  TSQLResult *res;
  TSQLRow *row;
  char buffer[400];

  // Mysql Statement for getting the runs...........................

  sprintf(buffer," select t1.run_number,t3.sensor_name,t2.sensors_sensor_id,t2.value from " 
                 " pedestal_set as t1, pedestals as t2,sensors as t3 where t1.run_number = %d "
                 " and t1.pedestal_set_id = t2.pedestal_set_id and t3.sensor_id = t2.sensors_sensor_id",runnum);

  res = (TSQLResult*)sql2->Query(buffer);

  if(res->GetRowCount() != 0){
     while((row = (TSQLRow*)res->Next())){ 
          if(atoi(row->GetField(2))> 50 && atoi(row->GetField(2)) < 67)
                  Exped[atoi(row->GetField(2)) - 51] = atof(row->GetField(3));
          else if(atoi(row->GetField(2))> 66 && atoi(row->GetField(2)) < 83)
	          Eyped[atoi(row->GetField(2)) - 67] = atof(row->GetField(3));
          else if(atoi(row->GetField(2))> 82 && atoi(row->GetField(2)) < 99)
                  Wxped[atoi(row->GetField(2)) - 83] = atof(row->GetField(3));
          else if(atoi(row->GetField(2))> 98 && atoi(row->GetField(2)) < 115)
                  Wyped[atoi(row->GetField(2)) - 67] = atof(row->GetField(3));
          delete row;
     }
  }

  delete res;
  delete sql2;*/

  return 0;
}

Int_t Run::Calculate_Backscatter(Int_t nex,Int_t nwx)
{
  // Count UP the backscatter fraction numbers.

/*  cout << "N_low, N_HIGH, (linebreak) heq entries " << endl;
  cout << nlow << " " << nhigh << " " << endl;
  cout << heq->GetEntries() << endl;
*/
  bscat.e_all     = heq->Integral(nlow,nhigh);
  bscat.etype_1   = heqF->Integral(nlow,nhigh);
  bscat.etype_23  = heqG->Integral(nlow,nhigh);
  
  bscat.e_alle    = sqrt(heq->Integral(nlow,nhigh)*rtime_e)/rtime_e;
  bscat.etype_1e  = sqrt(heqF->Integral(nlow,nhigh)*rtime_e)/rtime_e;
  bscat.etype_23e = sqrt(heqG->Integral(nlow,nhigh)*rtime_e)/rtime_e;
  
  bscat.w_all     = hwq->Integral(nlow,nhigh);
  bscat.wtype_1   = hwqF->Integral(nlow,nhigh);
  bscat.wtype_23  = hwqG->Integral(nlow,nhigh);
  
  bscat.w_alle    = sqrt(hwq->Integral(nlow,nhigh)*rtime_w)/rtime_w;
  bscat.wtype_1e  = sqrt(hwqF->Integral(nlow,nhigh)*rtime_w)/rtime_w;
  bscat.wtype_23e = sqrt(hwqG->Integral(nlow,nhigh)*rtime_w)/rtime_w;

  return 0;
}
//------------------------------------------------------------------------
Bool_t Run::IsType0(Int_t n)
{

  return kTRUE;
}
//-------------------------------------------------------------------------
Int_t Run::DefineCanvas(Int_t n)
{

  const char *name = Form(" ");
  
  if(GetRunType() == 1){
    name = Form("cr1");
  } else if(GetRunType() == 2) {
    name = Form("btr");
  } else if(GetRunType() == 3){
    name = Form("bckr");
  }
  
  Int_t ncols = 2;
  Int_t nrows = 2;
  
  c1 = new TCanvas(Form("%s[%i]->c1",name,n),
		   Form("Position Canvas for Run %i",runnum)
		       ,ncols*400,nrows*200);
  
  c2 = new TCanvas(Form("%s[%i]->c2",name,n),
		   Form("Backscattering Canvas for Run %i",runnum)
		       ,ncols*400,nrows*200);

  c1->Divide(ncols,nrows);
  c2->Divide(ncols+1,nrows);
  
  return 0;
}
//---------------------------------------------------------------------------
Int_t Run::Draw_Pos()
{

  c1->cd(1);
  hpw->Draw("colz");

  c1->cd(2);
  hpe->Draw("colz");

  return 0;
}
//----------------------------------------------------------------------------
void  Run::ReConnect(TSQLServer *sql)
{
  // Open file with the proper mysql databases  

  char *db = Form("mysql://%s/%s",getenv("UCNADBADDRESS"),getenv("UCNADB"));
  char *dbuser = Form("%s",getenv("UCNADBUSER"));
  char *dbpwd  = Form("%s",getenv("UCNADBPASS"));

  while(!(sql->IsConnected())){
    cout << "Having to retry connection (RunTime)  " << runnum <<  endl;
    sql->Connect(db,dbuser,dbpwd);
  }
  
  return;
}
//--------------------------------------------------------------------------
Int_t Run::SetRunTime(TSQLServer *sql)
{
  
  // read in the analysis database location. 
  // Check connection
  if(!(sql->IsConnected()))ReConnect(sql);
  
  TSQLResult *res;
  TSQLRow *row;
  char buffer[300];
  Int_t flipper = 0;
  // Mysql Statement for getting the runs...........................
  sprintf(buffer,"select t1.live_time_e,t1.live_time_w,t1.live_time,t2.flipper," 
	  "unix_timestamp(t2.end_time)-unix_timestamp(t2.start_time),t2.start_time "
	  "from analysis as t1, run as t2 where t1.run_number = t2.run_number and t1.run_number = %d ",runnum);
  // Send Query to Database;
  res = (TSQLResult*)sql->Query(buffer);
  // loop through results
  if(res->GetRowCount() != 0){
    while((row = (TSQLRow*)res->Next())){
      // Set east and west live times
      rtime_e = atof(row->GetField(0));
      rtime_w = atof(row->GetField(1));
      cout << "Run : " << runnum << " West live time " << rtime_w << " East live time " << rtime_e << endl;
      // Check flipper status
      if(!strncmp(row->GetField(3),"On",5))flipper = 1;
      // Look at log time
      log_time = atof(row->GetField(4));
      // calculate difference from the live time
      time_difference = log_time-rtime_e;
      // Set the Data of the run
      RunDate->Set(row->GetField(5));
      delete row;
    }
  }

  delete res;
  return flipper;
}
//-----------------------------------------------------------------------------------------------
Int_t Run::Book_Muons()
{
  
  // Fill some QADC spectra with backing veto tagged events.
  
  if(TaggedBackE && TaggedDriftW+TaggedDriftE+TaggedBackW+TaggedTopE == 0)hEmuon->Fill(EvisE);
  if(TaggedBackW && TaggedDriftW+TaggedDriftE+TaggedBackE+TaggedTopE == 0)hWmuon->Fill(EvisW);
  
  return -1;
}
//-----------------------------------------------------------------------------------------------
Int_t Run::SetOverFlow(Int_t i)
{
  toverflow = i;
  return -1;
}
//-----------------------------------------------------------------------------------------------
Int_t Run::SetBranches()
{
  
  // As the branch Address for the useful information in the tree
  // Things are commented out due to changes in the replay engine
  // over the years....
  t1->SetBranchAddress("ScintE"      ,&Escint);
//  t1->SetBranchAddress("TimeReal"    ,&TimeReal); //Seems to have been removed
					//  before I got here - JWW, 05132015
  t1->SetBranchAddress("ScintW"      ,&Wscint);
  t1->SetBranchAddress("DeltaT"      ,&DeltaT);
  t1->SetBranchAddress("PassedAnoE"  ,&PassedAnoE);
  t1->SetBranchAddress("PassedCathE" ,&PassedCathE);
  t1->SetBranchAddress("PassedAnoW"  ,&PassedAnoW);
  t1->SetBranchAddress("PassedCathW" ,&PassedCathW);
  t1->SetBranchAddress("TaggedTopE"  ,&TaggedTopE);
  t1->SetBranchAddress("TaggedBackE" ,&TaggedBackE);
  t1->SetBranchAddress("TaggedDriftE",&TaggedDriftE);
  t1->SetBranchAddress("TaggedBackW" ,&TaggedBackW);
  t1->SetBranchAddress("TaggedDriftW",&TaggedDriftW);
  t1->SetBranchAddress("TimeE"       ,&TimeE);
  t1->SetBranchAddress("TimeW"       ,&TimeW);
  t1->SetBranchAddress("TDCE"        ,&TDCE);
  t1->SetBranchAddress("TDCW"        ,&TDCW);
//  t1->SetBranchAddress("TofE"        ,&TofE);    //TofE/W deprecated, replaced							// by Tof
//  t1->SetBranchAddress("TofW"        ,&TofW);
  t1->SetBranchAddress("Tof"        ,&Tof);
  t1->SetBranchAddress("xEmpm"       ,&EastX);
  t1->SetBranchAddress("yEmpm"       ,&EastY);
  t1->SetBranchAddress("xWmpm"       ,&WestX);
  t1->SetBranchAddress("yWmpm"       ,&WestY);
  t1->SetBranchAddress("Sis00"       ,&Sis00);
  t1->SetBranchAddress("EvnbGood"    ,&EvnbGood);
  t1->SetBranchAddress("BkhfGood"    ,&BkhfGood);
  t1->SetBranchAddress("AnodeE"      ,&AnodeE);
  t1->SetBranchAddress("AnodeW"      ,&AnodeW);
  t1->SetBranchAddress("CathSumE"    ,&CathESum);
  t1->SetBranchAddress("CathSumW"    ,&CathWSum);
  t1->SetBranchAddress("EvisE"       ,&EvisE);
  t1->SetBranchAddress("EvisW"       ,&EvisW);
  t1->SetBranchAddress("PID"         ,&PID);
  t1->SetBranchAddress("Side"        ,&Side);
  t1->SetBranchAddress("Type"        ,&Type);
  //t1->SetBranchAddress("Etrue"       ,&Etrue);
  t1->SetBranchAddress("Erecon"       ,&Etrue);  // Change in 2011-2012 set
  t1->SetBranchAddress("EMWPC_E"     ,&EastMWPCEnergy);
  t1->SetBranchAddress("EMWPC_W"     ,&WestMWPCEnergy);

  return 0;
}

Int_t Run::Book_Raw(Int_t entry)
{
  // ------------------------------------------------------------
  // Analyze Each Event 
  //-----------------------------------------------------------

  // if the event counter in the header of the banks is misaligned then
  // the EvnbGood flag returns 1; 
  if(PID == 0) Count_Gammas();
  
  if(BkhfGood != 1 || EvnbGood != 1){
    // If BkhfGood flag is false then the header and or the footer of the TDC block
    // is missing in the event.   This routine is determine the affect of these events.
    // If the affect is small or zero we can use these events
    Handle_TDCCor_Evt();
//    return 0;
  }
  // Increment the passed counter if it was a good event
  passed1++;
  
  if( PID == 0 ) Count_Gammas_Good();
  // use the rotation matrix  
  //Rotate_Pos();
  // Define the Radius of the event in terms of the MPM unrotated variables.
  Float_t East_Rad = sqrt(EastX.center*EastX.center+EastY.center*EastY.center);
  Float_t West_Rad = sqrt(WestX.center*WestX.center+WestY.center*WestY.center);
  //----------------------------------------------------------------------------------------+
  // Fill Histograms with Evis based on a wire chamber or no wire chamber cut just          |
  // using the input register to determine trigger side.  this is for assessing neutron     |
  // generated backgrounds in the experiment.                                               |
  // Use the Raw ADC values with a stupid calibration parameter to get and an approximate   |
  // Evis. Since if the wire chamber cut fails we don't get a reliable energy calibration   |
  // in the replay.									    |
  //----------------------------------------------------------------------------------------+
  Float_t Esum = Etrue;  // Not really a difference between the scaling and using Evis so 
  Float_t Wsum = Etrue;  // will use Evis 
  //----------------------------------------------------------------------------------------
  if(PID != 2){  // Require DeltaT > 30 to remove after pulsing and PID != 2 
                                    // DeltaT cut removed for 2010 analysis
    if( Sis00 == 1 && East_Rad < RAD){ // removes the muons....
      hENoMWPC->Fill(Esum,1.);
      if(PID == 0 )hEMWPC->Fill(Esum,1.);
    }
    if( Sis00 == 2 && West_Rad < RAD){
      hWNoMWPC->Fill(Wsum,1.);
      if(PID == 1 )hWMWPC->Fill(Wsum,1.);
    }
  }
  //------------------------------------------------------------------------------------------
  // Make radial cut ....Require one side to be inside the radial cut.
  if((East_Rad > RAD && West_Rad > RAD ) )return 0;
  passed2++;
  //--------------------------------------------------------------------------------------
 // SetMultiplicity();
  // If a valid data stream event determine the event type.
  if((PID == 1) && Type == 0 && Side == 0 && East_Rad < RAD){
    // East Type 0 (no West side MWPC trigger)
    Handle_East_Type_0(entry);
  } else if((PID == 1) && Type == 0 && Side == 1 && West_Rad < RAD){
    // West Type 0 (no East side MWPC trigger)
    Handle_West_Type_0(entry);
  } else if((PID == 1) && Type > 0 ){
    // if both MWPC above threshold then must be a backscatter.
    Handle_Backscatters(entry);
  }
  
  if(PID!=2) Count_Time_All();
  if(PID==1) Count_Time_Beta();
  return -1;
}
//------------------------------------------------------------------------
void Run::SetMultiplicity()
{
   // Reset Multiplicities
   EastX.mult = 0.;
   EastY.mult = 0.; 
   WestX.mult = 0.;
   WestY.mult = 0.;
   // Recalculate 
   for(Int_t i = 0 ; i < 16 ; i++){
       if(eastX[i] > 100 + Exped[i]) EastX.mult++;
       if(eastY[i] > 100 + Eyped[i]) EastY.mult++;
       if(westX[i] > 100 + Wxped[i]) WestX.mult++;
       if(westY[i] > 100 + Wyped[i]) WestY.mult++;
   }

}
//---------------------------------------------------------------------------
Int_t Run::Handle_East_Type_0(Int_t entry)
{
  // Sum the two side together................
  EvisE += EvisW;
  // Fill the 1-d Energy historgrams
  ReconEnergy(&EvisE,0);
  // fill Radial 
  Double_t rad = sqrt(EastX.center*EastX.center + EastY.center*EastY.center);
  for(Int_t i = 0 ; i < 12 ; i++)
    if(rad > i*5 && rad < (i+1)*5)hERad[i]->Fill(Etrue,1.);
     
  // Fill Asymmetry Histograms...
  if(GetRunNumber() == 8300)cout << "Passed Choice D : " << entry << "   "  << EvisE << endl;
  heq->Fill(Etrue);
  heqB->Fill(Etrue);
  heqC->Fill(Etrue);
  heqD->Fill(Etrue);
  heqE->Fill(Etrue);
  // Fill Tube data
  he1->Fill(Escint.q1);
  he2->Fill(Escint.q2);
  he3->Fill(Escint.q3);
  he4->Fill(Escint.q4);
  // fill the position histogram for all events
  hpe->Fill(EastX.center,EastY.center,1);
  // fill a 2x2 grid histogram for the asymmetry calculation
  hApose->Fill(EastX.center,EastY.center);
  // Fill Multiplicity Plot
  hEType0_Multi->Fill((Float_t)(EastX.mult + EastY.mult),1.);
  hEAnode->Fill(EastMWPCEnergy,1.);
   
  return -1;
}
//--------------------------------------------------------------------------------
Int_t Run::Handle_West_Type_0(Int_t entry)
{

  EvisW += EvisE;
  
  ReconEnergy(&EvisW,0);
  
  Double_t rad = sqrt(WestX.center*WestX.center + WestY.center*WestY.center);

 // Hack for preNSAC analysis Etrue = 750./800. * Etrue;
  
  for(Int_t i = 0 ; i < 12 ; i++){
    if(rad > i*5 && rad < (i+1)*5){
      hWRad[i]->Fill(Etrue,1.);
    }
  }
  
   if(GetRunNumber() == 8300)cout << "Passed Choice D West : " << entry << "   "  << EvisW << endl;
  // Fill Asymmetry Histograms...
  hwq->Fill(Etrue);
  hwqB->Fill(Etrue);
  hwqC->Fill(Etrue);
  hwqD->Fill(Etrue);
  hwqE->Fill(Etrue);
  
  hw1->Fill(Wscint.q1);
  hw2->Fill(Wscint.q2);
  hw3->Fill(Wscint.q3);
  hw4->Fill(Wscint.q4);
  hpw->Fill(WestX.center,WestY.center,1);
  hAposw->Fill(WestX.center,WestY.center);
  hWType0_Multi->Fill((Float_t)(WestX.mult + WestY.mult),1.);
  hWAnode->Fill(WestMWPCEnergy,1.);
  
  return -1;
}

Int_t Run::Handle_Backscatters(Int_t entry)
{

  //--------------------------------------------------------------
  // Generate Random numbers for the likely-hood 2/3 separation
  
  Double_t type23guess = xr->Rndm();
  Double_t type2prob;
  Int_t mwpcbin;
  
  // Define the east and west radii  
  Double_t posrad  = sqrt(EastX.center*EastX.center + EastY.center*EastY.center);
  Double_t posradw = sqrt(WestX.center*WestX.center + WestY.center*WestY.center);

//if(Side==0) std::cout << "East: Radius " << posrad << ", x/y hit " << EastX.center << "/" << EastY.center << std::endl;
//if(Side==1) std::cout << "West: Radius " << posradw << ", x/y hit " << WestX.center << "/" << WestY.center << std::endl;

  // Define the difference in the east and west MWPC events.
  Double_t posdiff = TMath::Abs(sqrt(TMath::Power((WestX.center-EastX.center),2)
      + TMath::Power((WestY.center-EastY.center),2)));
  // There seems  to be a problem with the TDC jumping around
  // between runs so to line up the TDC histograms I've subtracted
  // off the cut time and added an arbitary scale of 150ns. This
  // doesn't resize, stretch or rebin the timing spectral, just
  // recenters the cut time to 150.
  Double_t NTDCE = TDCE - e_tdc_cut + TDC0;
  Double_t NTDCW = TDCW - w_tdc_cut + TDC0;

  NTDCE-=3000; NTDCW-=3000;
 
//  cout << "TEST THINGS FOREVER AND EVER " << NTDCE << " " << NTDCW << " Diff: " << NTDCE-NTDCW << endl;
 
  // Fill timing differenc histogram
  hTDCDiff->Fill(NTDCE-NTDCW,1);
  // Both MWPC above threshold, backscatter event
  htw->Fill(NTDCW);
  hte->Fill(NTDCE);
  // Fill this histogram with east vs west TDC values for events with valid header/footers
  hEvWTDC->Fill(NTDCE,NTDCW);
  // Sum total visible energy
  Float_t Etot = EvisW + EvisE;
  // the position different is greater than 25mm don't use this event
 // if(posdiff > 25)return -1; 
  passed3++;
			     // Based on the source studies there
			     // are accidentals in the backscattering
			     // by limiting the postion difference to 2.5cm
			     // we and put a tight cut on the accidental
			     // triggering.
  
  if( Side == 0  && Type == 1 && posrad < RAD){
    // East Type 1 backscatters
    ReconEnergy(&Etot,1);
    
    for(Int_t i = 0 ; i < 12 ; i++){
      if(posrad > i*5 && posrad < (i+1)*5){
	hERad[i]->Fill(Etrue,1);
      }
    }
 
    hEtype_1->Fill(Etot);
    // Fill the choice A analysis
    heq->Fill(Etrue);
    // fill the choise C analysis
    heqB->Fill(Etrue);
    heqC->Fill(Etrue);
    heqE->Fill(Etrue);
    heqF->Fill(Etrue);
    hEType1_Multi->Fill((Float_t)(EastX.mult + EastY.mult),1.);
    // fill the west position of east triggering events
    hRote->Fill(WestX.center,WestY.center,1);
    // look at the energy vs. normalized time
    hETimeVsE->Fill(Etrue,NTDCE,1);
    hRoteI->Fill(EastX.center,EastY.center,1);
    //-------------------------------------------------------------
    hEastType1XRot->Fill(EastX.center,EastY.center,EastX.center-WestX.center);
    hEastType1YRot->Fill(EastX.center,EastY.center,EastY.center-WestY.center);

    hPosDiffShifted->Fill(xErot-WestX.center,yErot-WestY.center);
    hPosDiffUnShifted->Fill(EastX.center-WestX.center,EastY.center-WestY.center);
  } else if(Side == 1 && Type == 1 && posradw < RAD){
    // West Type 1 backscatters
    ReconEnergy(&Etot,1);
    for(Int_t i = 0 ; i < 12 ; i++){
      if(posradw > i*5 && posradw < (i+1)*5){
	hWRad[i]->Fill(Etrue,1.);
      }
    }
    
    hWtype_1->Fill(Etrue);
    hRotw->Fill(EastX.center,EastY.center,1);
    hRotwI->Fill(WestX.center,WestY.center,1);

    hWestType1XRot->Fill(WestX.center,WestY.center,EastX.center-WestX.center);
    hWestType1YRot->Fill(WestX.center,WestY.center,EastY.center-WestY.center);
    hPosDiffShifted->Fill(xErot-WestX.center,yErot-WestY.center);
    hPosDiffUnShifted->Fill(EastX.center-WestX.center,EastY.center-WestY.center);
    hWTimeVsE->Fill(Etrue,NTDCW,1);
    // Fill Type 1 Analysis choices A,B,C,E,F
    hwq->Fill(Etrue);
    hwqC->Fill(Etrue);
    hwqB->Fill(Etrue);
    hwqE->Fill(Etrue);
    hwqF->Fill(Etrue);
    hWType1_Multi->Fill((Float_t)(WestX.mult + WestY.mult),1.);
  } else if(Type == 2 || Type == 3){
    
     ReconEnergy(&Etot,2);
     
     if(Side==0 && posrad < RAD){ 
       // East Type 2/3
      passed4++;
	for(Int_t i = 0 ; i < 12 ; i++){
	  if(posrad > i*5 && posrad < (i+1)*5){
	  hERad[i]->Fill(Etrue,1.);
	  }
	}
	hPsDfET23->Fill(EastX.center-WestX.center,EastY.center-WestY.center);
	hRoteI23->Fill(EastX.center,EastY.center);
	hRotw23->Fill(WestX.center,WestY.center);

        hEastType23XRot->Fill(EastX.center,EastY.center,EastX.center-WestX.center);
        hEastType23YRot->Fill(EastX.center,EastY.center,EastY.center-WestY.center);        

	hEtype_23->Fill(Etrue,1);
        hEType23_Multi->Fill((Float_t)(EastX.mult + EastY.mult),1.);
	// Fill Trigger side analysis
	heq->Fill(Etrue,1);	
	heqG->Fill(Etrue,1);
        //hwqG->Fill(Etrue,rtime_w/rtime_e);   //  PURELY FOR TESTING, NONPHYSICAL
        hEAnode23->Fill(EastMWPCEnergy,1);
	hE23Anode2d->Fill(AnodeE,AnodeW,1);
        // Fill hard cut analysis
	if(Type==3){//EastMWPCEnergy > Acut){
	  heqH->Fill(Etrue,1);
	  heqC->Fill(Etrue,1);
	} else {
	  hwqH->Fill(Etrue,1);
	  hwqC->Fill(Etrue,1);
	}
	for(Int_t i = 0 ; i < 10 ; i++){
          if(EastMWPCEnergy > (Double_t)i*0.1*Acut){
            hEHC[i]->Fill(Etrue,1);
          } else {
	    hWHC[i]->Fill(Etrue,1);
          }
        }
          mwpcbin = (Int_t)(EastMWPCEnergy/0.2);
	  if(mwpcbin < 0){
	    type2prob = 0.99;
	  } else if(mwpcbin > 99) {
	    type2prob = 0.20;
	  } else {
	    type2prob = seppar[mwpcbin];
	  }
	  if(type23guess < type2prob){
	    hwqE->Fill(Etrue,1);
	    hwqI->Fill(Etrue,1);
	  } else {
	    heqE->Fill(Etrue,1);
	    heqI->Fill(Etrue,1);
	  }
	 
     } else if(Side == 1 && posradw < RAD){
      // West Type 2/3
      for(Int_t i = 0 ; i < 12 ; i++){
	if(posradw > i*5 && posradw < (i+1)*5){
	  hWRad[i]->Fill(Etrue,1.);
	}
      }
      hWtype_23->Fill(Etrue,1);
      hWType23_Multi->Fill((Float_t)(WestX.mult + WestY.mult),1.);
      hRote23->Fill(EastX.center,EastY.center);
      hRotwI23->Fill(WestX.center,WestY.center);
      hPsDfET23->Fill(EastX.center-WestX.center,EastY.center-WestY.center);

      hWestType23XRot->Fill(WestX.center,WestY.center,EastX.center-WestX.center);
      hWestType23YRot->Fill(WestX.center,WestY.center,EastY.center-WestY.center);
      
      hwq->Fill(Etrue,1);
      //if(GetRunNumber() == 8296)cout << "Passed Choice on West G : " << entry << "  " << Etot <<  endl;
      hwqG->Fill(Etrue,1);
      hWAnode23->Fill(WestMWPCEnergy,1.);
      hW23Anode2d->Fill(AnodeE,AnodeW,1.);
      // fill the analysis choice c and h's based on trigger side anode cut
      if(Type==3){//WestMWPCEnergy > Acut){
	hwqH->Fill(Etrue,1);
	hwqC->Fill(Etrue,1);
      } else {
	heqH->Fill(Etrue,1);
	heqC->Fill(Etrue,1);
      }
      for(Int_t i = 0 ; i < 10 ; i++){
          if(WestMWPCEnergy > (Double_t)i*0.1*Acut){
            hWHC[i]->Fill(Etrue,1);
          } else {
            hEHC[i]->Fill(Etrue,1);
          }
      }

      mwpcbin = (Int_t)(WestMWPCEnergy/0.2);
	  if(mwpcbin < 0){
	    type2prob = 0.99;
	  } else if(mwpcbin > 99) {
	    type2prob = 0.20;
	  } else {
	    type2prob = seppar[mwpcbin];
	  }
	  if(type23guess < type2prob){
	    heqE->Fill(Etrue,1);
	    heqI->Fill(Etrue,1);
	  } else {
	    hwqE->Fill(Etrue,1);
	    hwqI->Fill(Etrue,1);
	  }
   }
  }
  return -1;
}

Int_t Run::Find_TDC_Cut(Int_t n,TSQLServer *sql2)
{
  
  //--------------------------------------------------------------------+
  //                                                                    |
  // Yawrrrrrrrrrrrrr                                                   |
  // This function finds the TDC cut, sensors_sensor_id = 143 and 144   |
  // are the id for the TDC cut chosen in the data quality inspections. |
  //                                                                    |
  // Because the cuts is not index by run number but analysis id the    |
  // analysis table must be cross linked by analysis id to determine    |
  // the run number.                                                    |
  //--------------------------------------------------------------------+

 // read in the analysis database location.
 //if(!(sql2->IsConnected()))ReConnect(sql2);

  e_tdc_cut = 130;
  w_tdc_cut = 130;
/*
  TH1F *hTDC_E_Scale = new TH1F("hTDC_E_Scale","Scaling",1000,100,4100);
  TH1F *hTDC_W_Scale = new TH1F("hTDC_W_Scale","Scaling",1000,100,4100);

  t1->Draw("TDCE >> hTDC_E_Scale","Type == 1","goff");
  t1->Draw("TDCW >> hTDC_W_Scale","Type == 1","goff");  

  Double_t xEast = hTDC_E_Scale->GetBinCenter(hTDC_E_Scale->GetMaximumBin());
  Double_t xWest = hTDC_W_Scale->GetBinCenter(hTDC_W_Scale->GetMaximumBin());
  */
  fEast_TDC_Scale = 1.;//fTDC_Max/xEast;
  fWest_TDC_Scale = 1.;//fTDC_Max/xWest;

//  delete hTDC_W_Scale;
 // delete hTDC_E_Scale;

  return 0;
}
//------------------------------------------------------------------------------
Int_t Run::Handle_TDCCor_Evt()
{
  
  Double_t NTDCE = (TDCE - e_tdc_cut)+TDC0;
  Double_t NTDCW = (TDCW - w_tdc_cut)+TDC0;
  
  // Plot total event energy for any event type
  if(BkhfGood == 0){
    hETDC_cor->Fill(EvisE,1.);
    hWTDC_cor->Fill(EvisW,1.);
  }
  
  hEvWTDC_cor->Fill(NTDCE,NTDCW,1.);
  
  if(EvnbGood == 0){
    hETDC_cor_passed->Fill(EvisE,1.);
    hWTDC_cor_passed->Fill(EvisW,1.);
  }
  
  return 0;
}
//-------------------------------------------------------------------------------
void Run::SetReconEnergy()
{
  
  // Set the Energy Reconstruction Parameters determined from
  // isotropic monoenergetic sources in PENELOPE
  
  // Also set the geometry type.
  
  if(runnum >= 7659 && runnum <= 9189){
    GeoType   = 0;
//  New parameters taken from Jianglai's code 
    epar0[1]= 71.1016;
    epar0[0]= 0.889869*1.1;
    epar0[2]=0.;

    epar1[1]= 184.128;
    epar1[0]= 0.775409*1.11;
    epar1[2]=  0.0;
    
    epar23[2]= 144.831;
    epar23[1]= 1.16586;
    epar23[0]= -0.000439037;
  } else if(runnum >= 9356 && runnum <= 10333){
    GeoType   = 1;

    epar0[0] = 0.87*1.1;
    epar0[1] = 88.6;
    epar0[2] = 0.;
     
    epar1[0] = 0.75*1.1;
    epar1[1] = 215.;
    epar1[2] = 0.;

    epar23[0] = -5.68e-4;
    epar23[1] = 1.20;
    epar23[2] = 171.7;
     
  } else if(runnum >= 10404){
    if(runnum<12000)
      GeoType = 2;
    else
      GeoType = 3;

    epar0[1]=41.2182;
    epar0[0]=0.929538*1.07;
    epar0[2]=0.;
     
    epar1[1]= 103.897;
    epar1[0]= 0.865697*1.07;
    epar1[2]= 0.;
    
    epar23[2]= 92.6774;
    epar23[1]= 1.19219;
    epar23[0]= -0.000400317;
  }
 
  
}
//------------------------------------------------------------------
void Run::ReconEnergy(Float_t *E, Int_t Type)
{
  // Reconstruct Energies 

  Float_t Etemp = *E;
  
   switch(Type){
     case 0 :
        Etemp = Etemp*epar0[0] + epar0[1];
       break;
     case 1 :
        Etemp = Etemp*epar1[0] + epar1[1];
       break;
     case 2 :
	Etemp = Etemp*Etemp*epar23[0] + Etemp*epar23[1] + epar23[2];
       break;
     default :
       break;
   }
    *E = Etemp;
}

//---------------------------------------------------------------------------------------
Int_t Run::SaveHistograms(Bool_t SAVE)
{

  // This routine should export the generated histograms to the root file.....
  if(SAVE){
    fHisto->cd(Form("%s_Analysis_%d_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated));
    Hlist->Write("",TObject::kOverwrite);
    fHisto->Write("",TObject::kOverwrite);
  }
  f1->Close("R");
  delete f1;
  return -1;

}
//---------------------------------------------------------------------------------------
Int_t Run::GetHistograms(Int_t dunce)
{
  fHisto->cd(Form("%s_Analysis_%d_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated));
  xr = new TRandom3();
  hpe          = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hpe_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hpw          = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hpw_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  heq          = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/heq_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hGammaCounts = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hGammaCount_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hGammaCountsg = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hGammaCountg_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hwq         = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hwq_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hENoMWPC    = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hENoMWPC_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWNoMWPC    = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWNoMWPC_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEMWPC      = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEMWPC_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWMWPC      = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWMWPC_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hETimeVsE   = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hETimeVsE_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWTimeVsE   = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWTimeVsE_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hETimeVsEP   = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hETimeVsEP_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWTimeVsEP   = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWTimeVsEP_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hETimeVsES   = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hETimeVsES_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWTimeVsES   = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWTimeVsES_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWAnode23   = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWAnode23_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEAnode23   = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEAnode23_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hESigNos    = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hESigNos_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWSigNos    = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWSigNos_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEKurie     = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEKurie_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWKurie     = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWKurie_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hPsDfET23   = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hPsDfET23_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hApose      = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hApose_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hAposw      = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hAposw_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEvWTDC_cor = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEvWTDC_cor_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEvWTDC     = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEvWTDC_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEastType1XRot  = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEastType1XRot_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEastType1YRot  = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEastType1YRot_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEastType23XRot = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEastType23XRot_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEastType23YRot = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEastType23YRot_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWestType1XRot  = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWestType1XRot_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWestType1YRot  = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWestType1YRot_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWestType23XRot = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWestType23XRot_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWestType23YRot = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWestType23YRot_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hPosDiffShifted = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hPosDiffShifted_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hPosDiffUnShifted = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hPosDiffUnShifted_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEAnode     = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEAnode_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWAnode     = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWAnode_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWCathS     = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWCathS_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hECathS     = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hECathS_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEType1_Primary   = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEType1_Primary_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEType1_Secondary = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEType1_Secondary_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWType1_Primary   = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWType1_Primary_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWType1_Secondary = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWType1_Secondary_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hE23Anode2d = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hE23Anode2d_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hW23Anode2d = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hW23Anode2d_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hTDCDiff    = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hTDCDiff_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hRote       = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hRote_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hRotw       = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hRotw_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hRoteI      = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hRoteI_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hRotwI      = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hRotwI_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hRote23     = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hRote23_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hRoteI23    = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hRoteI23_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hRotw23     = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hRotw23_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hRotwI23    = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hRotwI23_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
/*  hw1         = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hw1_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hw2         = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hw2_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hw3         = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hw3_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hw4         = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hw4_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
*/
//  hwqB        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hwqB_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
//  heqB        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/heqB_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hwqC        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hwqC_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  heqC        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/heqC_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
//  hwqD        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hwqD_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
//  heqD        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/heqD_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hwqE        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hwqE_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  heqE        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/heqE_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hwqF        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hwqF_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  heqF        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/heqF_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hwqG        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hwqG_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  heqG        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/heqG_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
//  hwqH        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hwqH_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
//  heqH        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/heqH_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
//  hwqI        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hwqI_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
//  heqI        = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/heqI_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWmuon      = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWmuon_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEmuon      = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEmuon_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hETDC_cor   = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hETDC_cor_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWTDC_cor   = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWTDC_cor_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hETDC_cor_passed = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hETDC_cor_passed_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWTDC_cor_passed = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWTDC_cor_passed_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  htw         = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/htw_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWtype_1    = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWtype_1_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWtype_23   = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWtype_23_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
/*  he1         = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/he1_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  he2         = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/he2_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  he3         = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/he3_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  he4         = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/he4_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
*/
  hte         = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hte_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEtype_1    = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEtype_1_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEtype_23   = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEtype_23_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hEType1_2d  = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEType1_2d_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hWType1_2d  = (TH2F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWType1_2d_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));

  for(Int_t i = 0 ; i < 12 ; i++){
    hERad[i]  = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hERad[%d]_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,i,runnum));
    hWRad[i]  = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWRad[%d]_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,i,runnum));
    if(i<10){
	hEHC[i] = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hEHC[%d]_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,i,runnum));
        hWHC[i] = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hWHC[%d]_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,i,runnum));
    }
  }

  hmr1   = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hmr1_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));
  hCountTimeRecord   = (TH1F*)fHisto->Get(Form("/%s_Analysis_%d_%d/hCountTimeRecord_%d",getenv("UCNA_ANA_AUTHOR"),(Int_t)RAD,rotated,runnum));

  hEERef  = new TH1F("hEERef","Reference ; Channel ;Cts"  ,EBINS,0,ESCALE);
  hEERef1 = new TH1F("hEERef1","Reference ; Channel ;Cts",EBINS,0,ESCALE);
  hEERef2 = new TH1F("hEERef2","Reference ; Channel ;Cts",EBINS,0,ESCALE);
  hEWRef  = new TH1F("hEWRef","Reference ; Channel ;Cts"  ,EBINS,0,ESCALE);
  hEWRef1 = new TH1F("hEWRef1","Reference ; Channel ;Cts",EBINS,0,ESCALE);
  hEWRef2 = new TH1F("hEWRef2","Reference ; Channel ;Cts",EBINS,0,ESCALE);
  
  Assign_to_Hist_List();
  
  return -1;
};

void Run::Assign_to_Hist_List()
{
  HEastAn->Add(heq);
//  HEastAn->Add(heqB);
  HEastAn->Add(heqC);
//  HEastAn->Add(heqD);
  HEastAn->Add(heqE);
  HEastAn->Add(heqF);
  HEastAn->Add(heqG);
//  HEastAn->Add(heqH);
//  HEastAn->Add(heqI);
  HEastAn->Add(hpe);
  HEastAn->Add(hte);
/*  HEastAn->Add(he4);
  HEastAn->Add(he1);
  HEastAn->Add(he2);
  HEastAn->Add(he3);             */
  HEastAn->Add(hEAnode23);
  HEastAn->Add(hEtype_1);
  HEastAn->Add(hETDC_cor_passed);
  HEastAn->Add(hEtype_23);
  HEastAn->Add(hEType1_Primary);
  HEastAn->Add(hEType1_Secondary);
  HEastAn->Add(hETDC_cor);
  HEastAn->Add(hApose);
  HEastAn->Add(hETimeVsE);
  HEastAn->Add(hETimeVsEP);
  HEastAn->Add(hETimeVsES);
  HEastAn->Add(hEMWPC);
  HEastAn->Add(hENoMWPC);
  HEastAn->Add(hRote);
  HEastAn->Add(hRoteI);
  HEastAn->Add(hRote23);
  HEastAn->Add(hRoteI23);
  
  HWestAn->Add(hwq);
//  HWestAn->Add(hwqB);
  HWestAn->Add(hwqC);
//  HWestAn->Add(hwqD);
  HWestAn->Add(hwqE);
  HWestAn->Add(hwqF);
  HWestAn->Add(hwqG);
//  HWestAn->Add(hwqH);
//  HWestAn->Add(hwqI);
  HWestAn->Add(hpw);
  HWestAn->Add(htw);
/*  HWestAn->Add(hw4);
  HWestAn->Add(hw1);
  HWestAn->Add(hw2);
  HWestAn->Add(hw3);        */
  HWestAn->Add(hWAnode23);
  HWestAn->Add(hWtype_1);
  HWestAn->Add(hWTDC_cor_passed);
  HWestAn->Add(hWtype_23);
  HWestAn->Add(hWType1_Primary);
  HWestAn->Add(hWType1_Secondary);
  HWestAn->Add(hWTDC_cor);
  HWestAn->Add(hAposw);
  HWestAn->Add(hWTimeVsE);
  HWestAn->Add(hWTimeVsEP);
  HWestAn->Add(hWTimeVsES);
  HWestAn->Add(hWMWPC);
  HWestAn->Add(hWNoMWPC);
  HWestAn->Add(hRotw);
  HWestAn->Add(hRotwI);
  HWestAn->Add(hRotw23);
  HWestAn->Add(hRotwI23);
  
   for(int i = 0; i < (int)(RAD/5.) ; i++){
     HEastAn->Add(hERad[i]);
     HWestAn->Add(hWRad[i]);
    }
  
};

void Run::Rotate_Pos()
{
   //================================================================================
   // This functions performs the rotation of the east and west MWPC positions.
   // A global variable "rotated" is set in the header file to choose between
   // two methods of applications : 
   //  1 - the east position is rotated and shifted by the full amount.
   //  2 - the east position is rotated and both sides are shift by 1/2 the 
   //      total amount in opposite directions.
   //
   // The TVector libraries are used to handle the vector, matrix operations
   //================================================================================
   // position variables are cast as Double_t's

   Double_t axE = EastX.center; 
   Double_t ayE = EastY.center;
   Double_t axW = WestX.center; 
   Double_t ayW = WestY.center;

   TVectorD posE(2); posE(0) = axE; posE(1) = ayE;
   TVectorD posW(2); posW(0) = axW; posW(1) = ayW;
   //
   switch (rotated)
   {
    case 1 :   
      posE *= matRotation;
      posE += vecShift;
      EastX.center = posE(0);
      EastY.center = posE(1);
      break;
    case 2 : // east/west split the shift
      posE *= matRotation;
      posE += 0.5*vecShift;
      EastX.center = posE(0);
      EastY.center = posE(1);
      posW -= 0.5*vecShift;
      WestX.center = posW(0);
      WestY.center = posW(1);
      break;
    case 3 : 
      posE -= vecShift;
      posE *= matRotation;
      posE += vecShift;
      xErot = posE(0);
      yErot = posE(1);
    default :
      break;
   }
   //
   return;
}

void Run::Load_Rotation_Matrix(Int_t nGeo)
{
  // the Following code is from Jianglai's analyzer.
  // It builds the 2x2 matrix for rotating and 
  // shifting the position vectors
 
  matRotation.Clear();
  vecShift.Clear();
  matRotation.ResizeTo(2,2);
  vecShift.ResizeTo(2);
  
  switch(nGeo)
  { 
    case 0 :
      matRotation(0,0) = 0.97;
      matRotation(0,1) = -0.077;
      matRotation(1,0) = 0.075;
      matRotation(1,1) = 0.97;
      vecShift(0) = 4.35;
      vecShift(1) = 1.79;
      break;
    case 1: 
      matRotation(0,0) = 0.98;
      matRotation(0,1) = -0.076;
      matRotation(1,0) = 0.072;
      matRotation(1,1) = 0.97;
      vecShift(0) = 5.11;
      vecShift(1) = 1.92;
      break;
    case 2 :
      matRotation(0,0) = 0.97;
      matRotation(0,1) = -0.077;
      matRotation(1,0) = 0.068;
      matRotation(1,1) = 0.98;
      vecShift(0) = 5.16;
      vecShift(1) = 1.2;
      break;
    case 3 :      
      matRotation(0,0) = 0.99901;
      matRotation(0,1) = -0.043084;
      matRotation(1,0) = 0.043084;
      matRotation(1,1) = 0.99901;
      vecShift(0)      = 7.11561;
      vecShift(1)      = 85.0963;
      break;
    default :
      matRotation(0,0) = 1;
      matRotation(0,1) = 0;
      matRotation(1,0) = 0;
      matRotation(1,1) = 1;
      vecShift(0) = 0;
      vecShift(1) = 0;
      break;
  }

  return;                                              
}

void Run::Load_Background(Int_t nRunNumber,Int_t nRunType,Int_t nGeo)
{
  // ----------------------------------------------
  const Int_t BetaRun    = 2;
  const Int_t BackGround = 3;

  fstream fparon;
  Int_t bin =0;
  
  if(nRunType == BackGround){
    if(nRunNumber > 9300 && nRunNumber < 9810){
      fparon.open("input_files/bkg_rad_gas.rate",fstream::in);
    } else {
      fparon.open("input_files/bkg_normal.rate",fstream::in);
    }
  } else if(nRunType == BetaRun){
    if(nGeo == 0){
      if(flipperOn == 0) fparon.open("input_files/beta_A_off.rate",fstream::in);
      if(flipperOn == 1) fparon.open("input_files/beta_A_on.rate",fstream::in);
    } else if(nGeo == 1){
      if(flipperOn == 0) fparon.open("input_files/beta_B_off.rate",fstream::in);
      if(flipperOn == 1) fparon.open("input_files/beta_B_on.rate",fstream::in);
    } else if(nGeo == 2){
      if(flipperOn == 0) fparon.open("input_files/beta_C_off.rate",fstream::in);
      if(flipperOn == 1) fparon.open("input_files/beta_C_on.rate",fstream::in);
    }   
  }
  
  for(Int_t i = 0 ; i< 53 ; i++)
    fparon >> bin >> parente[i][0] >> parente[i][1] >> parentw[i][0] >> parentw[i][1];

  fparon.close();
  
  return;
}
//======================================================================================
Float_t Run::CalcBinError(Int_t nBin,TH1F* h,Float_t time,Int_t Method)
{

  Float_t MINCNTS = 20.;
  Float_t MAXBIN  = 48.;
  Float_t BinCnt = h->GetBinContent(nBin)*time;
  
  if(Method == 0){
    if(BinCnt> MINCNTS || nBin > MAXBIN){ 
      if( nBin > MAXBIN ){
	h->SetBinError(nBin , (1. +sqrt(0.75 + BinCnt))/time);
      } else {
	h->SetBinError(nBin  , sqrt(BinCnt)/time);
      }
    } else {
      h->SetBinError(nBin  , sqrt(h->Integral(nlow,nhigh)/parente[49][0]*parente[nBin-1][0]*time)/time);
    }
  } else {
    if(BinCnt < MINCNTS ){
      h->SetBinError(nBin , (1. +sqrt(0.75 + BinCnt))/time);
    } else {
      h->SetBinError(nBin  , sqrt(BinCnt)/time);
    }
  }
  
  return -1.;
}
//====================================================================================================
void Run::Diagnosis_Run()
{
  using namespace std;

  cout << "  -----------------------------------------------------------" << endl;
  cout << "  Run : " << GetRunNumber() << " of type " << GetRunType() << " of Geometry " << GetGeo() << endl;
  cout << "  Octtype " << octtype << " Flipper State : " << flipperOn << endl;
  cout << "  Run Date : " << RunDate->AsString() << endl;
  cout << "  East Live Time : " << rtime_e << " \t West Live Time : " << rtime_w << endl;
  cout << "  Total Measured Rate : " << heq->Integral() + hwq->Integral() << " s^-1 " << endl;
  cout << "  -----------------------------------------------------------" << endl;
  cout << "  Backscatter Type \t East \t\t West " << endl;
//  cout << "  Type 0 \t \t"<< heqE->Integral() << "\t\t" << hwqE->Integral() << endl;
  cout << "  Type 1 \t \t"<< hEtype_1->Integral()  << "\t"<< hWtype_1->Integral() << endl;
  cout << "  Type 23\t \t"<< hEtype_23->Integral() << "\t"<< hWtype_23->Integral() << endl;
  cout << "  -----------------------------------------------------------" << endl;
  cout << "  TDC counts East   :  " << hte->Integral() << "\t West : " << htw->Integral() << endl;
  cout << "  Anode counts East :  " << hEAnode->Integral() << "\t West : " << hWAnode->Integral() << endl;
  cout << "  -----------------------------------------------------------" << endl;
  cout << "  Gamma East : " << GammasEast << " Gamma West : " << GammasWest << endl;
  cout << "  -----------------------------------------------------------" << endl;
  cout << "  Added : " << Hlist->GetEntries() << " objects to total histogram list. " << endl;
  cout << "  Added : " << HEastAn->GetEntries() << " objects to the East histogram list. " << endl;
  cout << "  Added : " << HWestAn->GetEntries() << " objects to the West histogram list. " << endl;
  cout << "************************** End of the Run Diagnostic *******************************" << endl;

};
//============================================================================================================
void Run::ScaleList(TObjArray *hL,Double_t time)
{
  
  TObjArrayIter *hLiter = new TObjArrayIter(hL);
  TH1* h;
  do{
    h = (TH1*)hLiter->Next();
    if(h) h->Scale(1./time);
  }while(h!=0);

  delete h; delete hLiter;

}

#endif
