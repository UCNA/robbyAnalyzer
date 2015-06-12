#ifndef runs_h
#define runs_h

#include <TMultiGraph.h>
#include <TH2.h>
#include <TH1.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraph.h> 
#include <TVectorD.h>
#include <TTree.h>
#include <TBranch.h>
#include <TDirectory.h>
#include <TLeaf.h>
#include <TAxis.h>
#include <TChain.h>
#include <TStyle.h>
#include <TMath.h>
#include <TApplication.h>
#include <TString.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TObjArray.h>
#include <TRandom3.h>
#include <TEllipse.h>
#include <TDatime.h>
#include <TMatrixD.h>
#include <TSQLServer.h>
#include <TSQLRow.h>
#include <TSQLResult.h>

#include "iostream"
#include "stdio.h"
#include "fstream"
#include "stdlib.h"
#include "string.h"

#include "typedefs.h"

#define TDC_MIN_CUT 20.  // TDC cut varaibles
#define TDC_MAX_CUT 135. // 
#define TDC_OFF_SET 6.   // 
#define ESCALE 2000.
#define EBINS  200
#define RAD  50.
#define nlow 23
#define Acut 4.14
#define nhigh 67
#define TDC0 150. // Define Renormalize cut position
#define fTDC_Max 150.

class Run {
  
  public:
    Int_t 	    PID;
    Int_t           Side;
    Int_t           Type;
    Int_t           PassedAnoE;
    Int_t           PassedCathE;
    Int_t           PassedAnoW;
    Int_t           PassedCathW;
    Int_t           TaggedTopE;
    Int_t           TaggedBackE;
    Int_t           TaggedDriftE;
    Int_t           TaggedBackW;
    Int_t           TaggedDriftW;
    Int_t           EvnbGood;
    Int_t           BkhfGood;
    
    Float_t         eastX[16];
    Float_t         eastY[16];
    Float_t         westX[16];
    Float_t         westY[16];
    
    Float_t         TimeReal;
    Float_t         TimeE;
    Float_t         TimeW;
    Int_t           Sis00;    // Micheal changed this to an Int_t after the june 2012 replay earlier replays will need a float
    Float_t         AnodeE;
    Float_t         AnodeW;
    Float_t         DeltaT;
    Float_t         CathESum;
    Float_t         CathWSum;
    Float_t         EvisE;
    Float_t         EvisW;
    Int_t           EastMultiX;
    Int_t           EastMultiY;
    Int_t           WestMultiX;
    Int_t           WestMultiY;

    Int_t           numevents;
    Int_t           numbetas;
    Float_t         TDCW;
    Float_t         TDCE;
    Float_t         TofE;
    Float_t         TofW;
    Float_t	    Tof;  // No TofE/W in 2011-2012 set, identical in 2013 set
				//  JWW 05/21
    Float_t         Etrue;
    Float_t         Erecon;  // True-Recon in 2011-2012 reanalysis (JWw 05/21)
    Float_t         WestMWPCEnergy;
    Float_t         EastMWPCEnergy;
    //------------------ 
    // Test position variables
    Float_t         xErot;
    Float_t         yErot;


  Scint_t Escint;
  Scint_t Wscint;

  
  mpmmwpc_t EastX;
  mpmmwpc_t EastY;
  mpmmwpc_t WestX;
  mpmmwpc_t WestY;

  public:
    Int_t toverflow;
    Int_t GeoType;

  private:
    Int_t runnum;
    Int_t runtype;
    
  public: 
    Run(Int_t n,Int_t m,TSQLServer *sql);
    virtual ~Run();
    char *octtype; 
    Int_t flipperOn;
    Double_t seppar[100];
    Double_t parente[53][2];
    Double_t parentw[53][2];
    TDatime *RunDate;  // Date
    TDirectory *dHists; // Directory for analysis
    Int_t runb;
    Bool_t AnalysisDirExist;
    Int_t passed1,passed2,passed3,passed4;
    Int_t rotated;
    Float_t epar0[3],epar1[3],epar23[3];
    Float_t Exped[16],Eyped[16],Wxped[16],Wyped[16];
    // Structure to hold backscatter numbers
    backs_t bscat;
    Int_t stupid;
    Double_t fEast_TDC_Scale;
    Double_t fWest_TDC_Scale;

    // Rad Cut rate
    Float_t erad[20],wrad[20];
    Float_t erader[20],wrader[20];
    Float_t eradu[20],wradu[20];
    Float_t eraduer[20],wraduer[20];
    
    // Anode Cut parameters
    Float_t nsig;
    Float_t e_acut;
    Float_t w_acut;
    Float_t e_tdc_cut;
    Float_t w_tdc_cut; 
    Float_t T2_frac;
    Float_t Mon_Rate;

    // Backscatter Counters
    Float_t E_tot_e,W_tot_e;
    Float_t E_type_I,E_type_23;
    Float_t W_type_I,W_type_23;
    TH1F *hEKurie,*hWKurie;
    
    //Simple Asymmerty
    Double_t fAsym_Run[40];
    Double_t fAsym_Run_Err[40];
    
    // Run time
    Float_t rtime_e,rtime_w;
    Float_t log_time,time_difference;
    Float_t CountTimeWAll,CountTimeEAll;
    Float_t CountTimeWFirst,CountTimeEFirst;
    Float_t CountTimeWBeta,CountTimeEBeta;
    Float_t CountTimeWFirstBeta,CountTimeEFirstBeta;
    Float_t EventsEast,EventsWest;
    Float_t GammasEast,GammasWest;
    Int_t BetasEast,BetasWest;
    Float_t GammasEastg,GammasWestg;
    Float_t EastClock[100000];
    Float_t WestClock[100000];

    // Array of Histograms
    TObjArray *Hlist,*HEastAn,*HWestAn;
    // Anode and Cathode histograms
    TH1F *hEAnode,*hECathS,*hWAnode,*hWCathS;
    TRandom3 *xr;
    
    // Tree and File for each run
    TFile *f1,*fHisto;
    TTree *t1;
    // Run Canvas 
    TCanvas *c1;
    TCanvas *c2;
    TCanvas *clocks;   
 
    // Rotation Parameters
    TMatrixD matRotation;
    TVectorD vecShift;

    // Histograms
    TH2F *hpe,*hpw;               // Over all position for all events
    // Asymmetry Analysis Spectra 
    TH1F *heq,*hwq;               // Qadc Spectrum for each peak
    TH1F *heqC,*hwqC;             // Analysis C option
    TH1F *heqB,*hwqB;             // Analysis B option
    TH1F *heqD,*hwqD;             // Analysis D option
    TH1F *heqE,*hwqE;             // Analysis E option
    TH1F *heqF,*hwqF;             // Analysis F option
    TH1F *heqG,*hwqG;             // Analysis G option
    TH1F *heqH,*hwqH;             // Analysis H option
    TH1F *heqI,*hwqI;             // Analysis I option
    TH1F *hEHC[10],*hWHC[10];
    //------------------------------------------------------------------------------------------------
    TH2F *hPsDfET23;                 // Difference in East-West Position for type 2/3
    TH2F *hApose,*hAposw;            // Quadrant based Asymmetry 
    TH1F *hte,*htw;                  // Tdc spectra for each peak
    TH1F *htwest,*hteast;            // East and West clocks for betas
    TH1F *hTimeE,*hTimeW;            // East and West clocks for all events
    TH1F *hTimeEBetas,*hTimeWBetas;  // East and West clocks for betas
    //------------------------------------------------------------------------------------------------
    TH1F *hw1,*hw2,*hw3,*hw4;     // Individual PMT 
    TH1F *he1,*he2,*he3,*he4;     // Qadc spectra
    TH1F *hecy,*hwcy;             // Projections
    TH1F *hEERef,*hEERef1,*hEERef2;// Monte Carlo Histogram
    TH1F *hEWRef,*hEWRef1,*hEWRef2;
    TH1F *hEmuon,*hWmuon;         // Backing Veto hits;
    TH1F *hESigNos,*hWSigNos;     // Signal to Noise Histograms
    TH1F *hENoMWPC,*hWNoMWPC;     // No MWPC cut used for looking at gamma backgrounds.
    TH1F *hEMWPC,*hWMWPC;         // MWPC cut used for looking at gamma backgrounds.
    TH1F *hEType0_Multi,*hWType0_Multi,*hEType1_Multi,*hWType1_Multi,*hEType23_Multi,*hWType23_Multi;
    TH1F *hERad[12],*hWRad[12];
    //--------------------------------------------------------------------------------------------------------
    // Backscatter Histograms
    TH1F *hEtype_1,*hEtype_23,*hWtype_1,*hWtype_23; // Energy spectra for backscatters of various types
    TH1F *hTDCDiff;
    TH2F *hRote,*hRotw,*hRoteI,*hRotwI;             // Initial and final positions of Type 1 backscatter events
    TH2F *hRote23,*hRotw23,*hRoteI23,*hRotwI23;     // Initial and final
    TH2F *hPosDiffShifted;
    TH2F *hPosDiffUnShifted;
    //--------------------------------------------------------------------------------------------------------
    //   Rotation directions............
    TH2F *hEastType1XRot,*hEastType1YRot;             
    TH2F *hWestType1XRot,*hWestType1YRot;  
    TH2F *hEastType23XRot,*hEastType23YRot;          
    TH2F *hWestType23XRot,*hWestType23YRot;
    //---------------------------------------------------------------------------------------------------------
    TH2F *hETimeVsE,*hWTimeVsE;
    TH2F *hETimeVsES,*hWTimeVsES;
    TH2F *hETimeVsEP,*hWTimeVsEP;

    TH1F *htdcE,*htdcW;           // 1-D TDC histos
    TH1F *hEAnode23,*hWAnode23;
    TH2F *hEType1_2d,*hWType1_2d; // 2-d backscattering 
    TH1F *hEType1_Primary,*hEType1_Secondary;
    TH1F *hWType1_Primary,*hWType1_Secondary;
    TH2F *hE23Anode2d,*hW23Anode2d;
    TH1F *hGammaCounts,*hGammaCountsg;
    //---------------------------------------------------------------------------------------------------------
    // TDC Corruption Histograms
    TH1F *hETDC_cor,*hWTDC_cor;
    TH2F *hETDC_cor_pos,*hWTDC_cor_pos;
    TH2F *hEvWTDC_cor,*hEvWTDC;
    TH1F *hETDC_cor_passed, *hWTDC_cor_passed;
    //----------------------------------------------------------------------------------------------------------
    // Functions
    Int_t Scale2Time(Int_t nex,Int_t nwx);
    void SetMultiplicity();
    Int_t Count_Time_Beta();
    Int_t Initialize_hist(Int_t n,Int_t nex,Int_t nwx);
    Int_t Set_Anode_Cut(Int_t n);//,TSQLServer *sql);
    Int_t Calculate_Backscatter(Int_t nex,Int_t nwx);
    Bool_t IsType0(Int_t n);
    Int_t DefineCanvas(Int_t n);
    Int_t Handle_TDCCor_Evt();
    Int_t SetRunTime(TSQLServer *sql);
    Int_t GetOverFlow() {return toverflow;}
    Int_t SetOverFlow(Int_t i); 
    Int_t GetRunNumber() {return runnum;}
    Int_t GetRunType() {return runtype;}
    Int_t GetBackgroundRun() {return runb;}
    Int_t GetGeo() {return GeoType;}
    Int_t Draw_Pos();
    Int_t Book_Raw(Int_t entry);
    void  SetReconEnergy();
    void ScaleList(TObjArray *hL,Double_t time);
    void  ReconEnergy(Float_t *E,Int_t Type);
    Int_t Handle_East_Type_0(Int_t entry);
    Int_t Handle_West_Type_0(Int_t entry);
    Int_t Handle_Backscatters(Int_t entry);
    Int_t Find_TDC_Cut(Int_t n,TSQLServer *sql);
    Bool_t OpenRun(Int_t n);
    //Bool_t Check_Vetos();
    Int_t Book_Muons();
    void  Load_Background(Int_t nRunNumber,Int_t nRunType,Int_t nGeo);
    Int_t SetBranches();
    Int_t SaveHistograms(Bool_t SAVE);
    Int_t GetHistograms();
    Int_t Count_Time_All();
    Int_t Count_Time_First();
    //Int_t Count_Time_Beta();
    Int_t Count_Events();
    Int_t Count_Gammas();
    Int_t Count_Betas();
    void  ReConnect(TSQLServer *sql);
    Int_t Count_Gammas_Good();
    void  Load_Rotation_Matrix(Int_t nGeo);
    Float_t CalcBinError(Int_t nbin,TH1F* h,Float_t time,Int_t Method);
    void  Rotate_Pos();
    void  ScaleTDC();
    void  DeleteHistos();
    void  Diagnosis_Run();
    void  Assign_to_Hist_List();

    void  CreateTH2F(TH2F *&h,const char* name,const char* title,Int_t nbins,Double_t xmin,Double_t xmax,Int_t nybins,Double_t ymin,Double_t ymax,Int_t nSide);
    void  CreateTH1F(TH1F *&h,const char* name,const char* title,Int_t nbins,Double_t xmin,Double_t xmax,Int_t nSide);
    
};

#endif
