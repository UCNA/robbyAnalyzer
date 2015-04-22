#ifndef __analysis_h__
#define __analysis_h__

#include "get_peds.h"
#include "Run.h"
#include "Octet.h"
#include "iomanip"
#include "Beta_Run.h"
#include "Bck_Run.h"
#include "analysisMathFunc.h"
#include "vector"
#include "typedefs.h"
#include <sys/resource.h>

// Root includes
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TPad.h>
#include <TLegend.h>
#include <TVector3.h>
#include <TGaxis.h>
#include <TGraphErrors.h> 
#include <TArrow.h>
#include <TRandom3.h>
#include <TSQLServer.h>
#include <TSQLRow.h>
#include <TSQLResult.h>

#define ECUTMIN  100
#define ECUTMAX  1600
#define ARRAYBIN 10000
#define radcut 50.
#define Rbins  20

using namespace std;
//-------------------------------------
// Analysis functions
Int_t analyze_background_runs(Int_t n,std::vector<Bck_Run *>bk,vector<Int_t> nrun,Int_t remake);
Int_t analyze_beta_runs(Int_t n, std::vector<Beta_Run *>bta,vector<Int_t> nrun,Int_t remake);
Int_t Subtract_Backgrounds(std::vector<Beta_Run *>bta,std::vector<Bck_Run *>bk,Int_t nb,Int_t nbk);
Int_t analyze_octets();
void  PlotRunTimes();
Int_t ParseMPMOctetList(vector<Int_t> &RunListMPM);
void  Fill_Asymmetry_Vector(Asym_t &A,Double_t oct,Double_t octer,Int_t &index);
void  Fill_Asymmetry_Vector_Oct(Asym_t &A,Double_t oct,Double_t octer,Int_t index);
void  SetROOTOpt();
void  CallAnalysisTasks();
void average_type1();
void AddRates(Double_t *Cter,TH1F *hOut,TH1F *hIn,Int_t ibin);
void  GetSingleRun();
bool  GetListofRuns(Int_t n1,Int_t n2);
void  Plot_RawAsymmetries(TGraphErrors *g1,TGraphErrors *g2,TGraphErrors *g3,Int_t noct,Int_t NoctC,Int_t NQuat);
//void  SetRate(Double_t *xR,Double_t *Rate,Double_t *RateEr,Double_t Run,Double_t mRate,Double_t mRateEr,Int_t &index);
void  DrawBackFractionPad(TCanvas *cCan,Int_t nPad,TGraphErrors *gEast,TGraphErrors *gWest,TF1 *fEast,
                          TF1 *fWest,TLegend *legBack,Double_t Ymax);
void  DrawMWPCPanel1(TH1F *h1,TH1F *h2,TH1F *h3,TH1F *h4,TLegend *leg1,Double_t ErrIntOff,Double_t ErrRatOff,Double_t ErrIntOn,Double_t ErrRatOn);
void  Collect_Pos();
void  ZeroThings();
void  CollectRates();
void  CollectTypeRot();
void  Average_A();
void  PrintPDF();
void  Plot_E_Chis();
void  Plot_ChiDis();
void  Collect_Stuff();
void  Collect_Energy_Spectra();
void  DefineRotationCollectionHistos(Int_t rbins,Float_t rad);
void  CollectTDCCor();
void  Collect_TvsE();
void  Plot_23_Diff();
bool  SetAnalysisType(Int_t *remake);
void  Collect_Rad();
void  Collect_Gammas();
void  CollectAsym();
void  Plot_Multi();
void  OutPutToMPMDB(Int_t noct);
void  Collect_23Anode();
void  Collect_TDCDiff();
void  Plot_MonRate();
void  Load23SepArray();
void  CalcSimplSuper();
void  Plot_Timing();
void  Fill_Timing();
void  Get_Base_Super(Int_t i, char oct1[4],char oct2[4]);
void  Collect_Octets();
Double_t  DrawEnergyPanel(TH1F *hEdraw,TH1F *hERef,TH1F *hTot,TH1F *hRefTot,Int_t npad,TCanvas *c,Double_t scaling = 1.);
void  Define_E_Spec();
void  GammaBack();
void  DrawRates();
void  DrawRadialCounts(Double_t *ARader,Double_t Aave);
Bool_t checkruns(Int_t currentrun,vector<Int_t> openruns,Int_t nruns);
void  TrackAnodeMPV();
void  TrackStats();
//---------------------------------------------------------------------------

void Fill_Foreground(Double_t *t1,Double_t *t2, Double_t *t3, Double_t *t4,Double_t *tl,Double_t *nmiss);
void Fill_Background(Double_t *t1,Double_t *t2, Double_t *t3, Double_t *t4);

void Set_Error_Bars(Double_t t1,Double_t t2, Double_t t3, Double_t t4,
		     Double_t tb1,Double_t tb2, Double_t tb3, Double_t tb4);
//--------------------------------------------------------

#endif
