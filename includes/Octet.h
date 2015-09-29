#ifndef Octet_h
#define Octet_h
 
#include <TMultiGraph.h>
#include <TH2.h>
#include <TH1.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TTree.h>
#include <TBranch.h>
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
#include <TMySQLServer.h>
#include <TMySQLRow.h>
#include <TMySQLResult.h>
#include <Beta_Run.h>
#include <Bck_Run.h>
#include <TEllipse.h>

#include "Enums.hh"
#include "AnalysisDB.hh"
#include "iostream"
#include "stdio.h"
#include "fstream"
#include "stdlib.h"
#include "string.h"

#define nchoices 20

class Octet {

  public :
    Octet(Int_t n,Int_t j, Int_t m,Double_t radius);
    virtual ~Octet();
  
  private :
    // Octet number, first and last runs, index, and radial cut
    Int_t Num;
    Int_t First,Last;
    Int_t QuartetAFirst,QuartetALast;
    Int_t QuartetBFirst,QuartetBLast;
    Int_t SuperA1First,SuperA1Last;
    Int_t SuperA2First,SuperA2Last;
    Int_t SuperB1First,SuperB1Last;
    Int_t SuperB2First,SuperB2Last;


    Int_t BetaIndex;
    Double_t radial_cut;
  
  public :
    // This set of vectors hold the run numbers for the octet.
    std::vector<Int_t> A2,A5,A7,A10,B2,B5,B7,B10;
    std::vector<Int_t> A1,A4,A9,A12,B1,B4,B9,B12;
    
    Double_t Asym;

    Double_t Asuper1[nchoices],Asuper2[nchoices],Bsuper1[nchoices],Bsuper2[nchoices];    
    Double_t Asuper1er[nchoices],Asuper2er[nchoices],Bsuper1er[nchoices],Bsuper2er[nchoices]; 
    // Holds the parent distributions for the backgrounds
    // used in calculating the uncertainty in low statistics bins
    Float_t  Asupertime1,Asupertime2,Bsupertime1,Bsupertime2;
    Float_t  Asupertime1a,Asupertime2a,Bsupertime1a,Bsupertime2a;
    Double_t parente[53][2],parentw[53][2];
    // the Product and Sum Octet Super-Ratios
    Double_t A_multi_A[nchoices],A_sum_A[nchoices];
    Double_t A_multi_B[nchoices],A_sum_B[nchoices];
    Double_t A_multi[nchoices]  ,A_sum[nchoices];
    // Uncertainties....
    Double_t A_multi_Aer[nchoices],A_sum_Aer[nchoices];
    Double_t A_multi_Ber[nchoices],A_sum_Ber[nchoices];
    Double_t A_multier[nchoices]  ,A_sumer[nchoices];
    // Integral Count rates for each run and analysis choice.
    Double_t nA2[nchoices],nA5[nchoices],nA7[nchoices],nA10[nchoices],nA2w[nchoices],nA5w[nchoices],nA7w[nchoices],nA10w[nchoices];
    Double_t nB2[nchoices],nB5[nchoices],nB7[nchoices],nB10[nchoices],nB2w[nchoices],nB5w[nchoices],nB7w[nchoices],nB10w[nchoices];
    Double_t nA2b[nchoices],nA5b[nchoices],nA7b[nchoices],nA10b[nchoices],nA2bw[nchoices],nA5bw[nchoices],nA7bw[nchoices],nA10bw[nchoices];
    Double_t nB2b[nchoices],nB5b[nchoices],nB7b[nchoices],nB10b[nchoices],nB2bw[nchoices],nB5bw[nchoices],nB7bw[nchoices],nB10bw[nchoices];
    // Run times....
    Float_t  tA2ef,tA2ei,tA2wf,tA2wi,tA5ef,tA5ei,tA5wf,tA5wi;
    Float_t  tA7ef,tA7ei,tA7wf,tA7wi,tA10ef,tA10ei,tA10wf,tA10wi;
    Float_t  tB2ef,tB2ei,tB2wf,tB2wi,tB5ef,tB5ei,tB5wf,tB5wi;
    Float_t  tB7ef,tB7ei,tB7wf,tB7wi,tB10ef,tB10ei,tB10wf,tB10wi;
    Float_t  tA2efa,tA2eia,tA2wfa,tA2wia,tA5efa,tA5eia,tA5wfa,tA5wia;
    Float_t  tA7efa,tA7eia,tA7wfa,tA7wia,tA10efa,tA10eia,tA10wfa,tA10wia;
    Float_t  tB2efa,tB2eia,tB2wfa,tB2wia,tB5efa,tB5eia,tB5wfa,tB5wia;
    Float_t  tB7efa,tB7eia,tB7wfa,tB7wia,tB10efa,tB10eia,tB10wfa,tB10wia;
    Double_t tA2,tA5,tA7,tA10,tB2,tB5,tB7,tB10;
    Double_t tA2w,tA5w,tA7w,tA10w,tB2w,tB5w,tB7w,tB10w;
    Double_t tA2b,tA5b,tA7b,tA10b,tB2b,tB5b,tB7b,tB10b;
    Double_t tA2bw,tA5bw,tA7bw,tA10bw,tB2bw,tB5bw,tB7bw,tB10bw;
    // Radial counts for each run type and analysis choice..
    Double_t A2rad[nchoices] ,A5rad[nchoices] ,A7rad[nchoices] ,A10rad[nchoices];
    Double_t A2wrad[nchoices],A5wrad[nchoices],A7wrad[nchoices],A10wrad[nchoices];
    Double_t B2rad[nchoices] ,B5rad[nchoices] ,B7rad[nchoices] ,B10rad[nchoices];
    Double_t B2wrad[nchoices],B5wrad[nchoices],B7wrad[nchoices],B10wrad[nchoices];
    
    Double_t A2rade[nchoices] ,A5rade[nchoices] ,A7rade[nchoices] ,A10rade[nchoices];
    Double_t A2wrade[nchoices],A5wrade[nchoices],A7wrade[nchoices],A10wrade[nchoices];
    Double_t B2rade[nchoices] ,B5rade[nchoices] ,B7rade[nchoices] ,B10rade[nchoices];
    Double_t B2wrade[nchoices],B5wrade[nchoices],B7wrade[nchoices],B10wrade[nchoices];
    // Averaged asymmetry and total counts....
    Double_t A_rad[nchoices],A_rader[nchoices];
    Double_t TotCounts,TotRadCounts[nchoices];
    //=======================================================================
    //----Analsis Choice Histograms .......
    TH1F *hA2[nchoices] ,*hA5[nchoices] ,*hA7[nchoices] ,*hA10[nchoices];
    TH1F *hB2[nchoices] ,*hB5[nchoices] ,*hB7[nchoices] ,*hB10[nchoices];    
    TH1F *hA2w[nchoices],*hA5w[nchoices],*hA7w[nchoices],*hA10w[nchoices];
    TH1F *hB2w[nchoices],*hB5w[nchoices],*hB7w[nchoices],*hB10w[nchoices];
    //======================================================================
    TH1F *hA1[nchoices] ,*hA4[nchoices] ,*hA9[nchoices] ,*hA12[nchoices];
    TH1F *hB1[nchoices] ,*hB4[nchoices] ,*hB9[nchoices] ,*hB12[nchoices];    
    TH1F *hA1w[nchoices],*hA4w[nchoices],*hA9w[nchoices],*hA12w[nchoices];
    TH1F *hB1w[nchoices],*hB4w[nchoices],*hB9w[nchoices],*hB12w[nchoices];
    // Asymmetry vs. E.
    TH1F *hAsyA[nchoices],*hAsyB[nchoices],*hAsyTot[nchoices];

    Int_t FirstRunIndex;

    // ---- Analysis Functions ----------------------------------------------------------------

    Double_t Calc_Super();
    Float_t Calc_Super_Time();
    Float_t Calc_Super_Time_All();
    Double_t Calc_A_multi();
    Double_t Calc_A_sum();
    void CalculateSuperTime(Float_t R1M, Float_t R1MI, Float_t R1P, Float_t R1PI, 
                            Float_t R2M, Float_t R2MI, Float_t R2P, Float_t R2PI, 
                            Float_t &R);
    void CalculateSuper(Double_t R1,Double_t R2,Double_t R3,Double_t R4,
			Double_t B1,Double_t B2,Double_t B3,Double_t B4,
			Double_t t1,Double_t t1b,Double_t t2,Double_t t2b,
			Double_t t3,Double_t t3b,Double_t t4,Double_t t4b,
			Double_t &A,Double_t &Aer);
    void Get_Rad_A();
    void Initialize_Histo(Int_t n);
    void Remove_Octet();  // From long-shot dream of fixing this
    void Find_Runs(std::vector<Beta_Run*>bta,std::vector<Bck_Run*>bck,Int_t nMax);

    Int_t GetFirst() {return First;}
    Int_t GetLast() {return Last;}
    
    Int_t GetQuartetBStart(){return QuartetBFirst;}
    Int_t GetQuartetBEnd(){return QuartetBLast;}
    Int_t GetQuartetAStart(){return QuartetAFirst;}
    Int_t GetQuartetAEnd(){return QuartetALast;}
    Int_t GetSuperA1Start(){return SuperA1First;}
    Int_t GetSuperA1End(){return SuperA1Last;}
    Int_t GetSuperA2Start(){return SuperA2First;}
    Int_t GetSuperA2End(){return SuperA2Last;}
    Int_t GetSuperB1Start(){return SuperB1First;}
    Int_t GetSuperB1End(){return SuperB1Last;}
    Int_t GetSuperB2Start(){return SuperB2First;}
    Int_t GetSuperB2End(){return SuperB2Last;}
    
    Int_t GetIndex() {return BetaIndex;}
    Double_t GetRadCut() {return radial_cut;}
    void Load_Background();
    void Fill(TH1F *h,TH1F *h2,Double_t t1,Double_t* t2,Int_t inc);
    
    Double_t Average_Quartet_Rate(Double_t r1,Double_t t1,Double_t b1,Double_t tb1,
				  Double_t r2,Double_t t2,Double_t b2,Double_t tb2,Int_t side);
				  
    Double_t Average_Octet_Rate(Double_t r1,Double_t t1,Double_t b1,Double_t tb1,
				Double_t r2,Double_t t2,Double_t b2,Double_t tb2,
				Double_t r3,Double_t t3,Double_t b3,Double_t tb3,
				Double_t r4,Double_t t4,Double_t b4,Double_t tb4,Int_t side);

    Double_t Average_Octet_Bin2(Double_t r1,Double_t t1,Double_t b1,Double_t bt1,
				Double_t r2,Double_t t2,Double_t b2,Double_t bt2, 
				Double_t r3,Double_t t3,Double_t b3,Double_t bt3,
				Double_t r4,Double_t t4,Double_t b4,Double_t bt4, 
				Double_t &ett);
    void Calc_A_sum_Bin();
    void Fill_Rad(Double_t *x1,Double_t *x1e,Float_t *x2,Float_t *x2e,Float_t t1);
    void GetIntegralCounts();
    
    Double_t Average_Octet_Cts_Bin(Double_t r1,Double_t t1,Double_t r2,Double_t t2,
				   Double_t r3,Double_t t3,Double_t r4,Double_t t4, 
				   Double_t b1,Double_t bt1,Double_t b2,Double_t bt2,
				   Double_t b3,Double_t bt3,Double_t b4,Double_t bt4, 
				   Double_t &ett);
    
    Double_t Average_Quartet_Bin(Double_t r1,Double_t e1,Double_t r2,Double_t e2,Double_t &ett);
    Double_t Average_Octet_Bin(Double_t r1,Double_t e1,Double_t r2,Double_t e2,
				  Double_t r3,Double_t e3,Double_t r4,Double_t e4, Double_t &ett);
    
    Double_t ErrorSum(Double_t n1,Double_t n2,Double_t t1,Double_t t2,
		      Double_t n3,Double_t n4,Double_t t3,Double_t t4,Int_t side);
		      
    Double_t ErrorSum(Double_t n1,Double_t n2,Double_t n3,Double_t n4,
		      Double_t t1,Double_t t2,Double_t t3,Double_t t4,
		      Double_t n5,Double_t n6,Double_t n7,Double_t n8,
		      Double_t t5,Double_t t6,Double_t t7,Double_t t8,Int_t side);
		      
    Double_t ErrorMultiRate(Double_t n1,Double_t n2,Double_t t1,Double_t t2);
    void OutPutRatesToMPMDB(std::vector<Int_t> nRuns,AFPState afp,EventType EType,TH1F *hEast
			    ,TH1F *hWest,Double_t tEast,Double_t tWest,AnalysisDB *adb);
    void Debugger(Int_t nOct);
    void OutPutToDB();
};

#endif
