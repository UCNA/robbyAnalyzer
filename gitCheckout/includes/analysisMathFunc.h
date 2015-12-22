#ifndef __MathFunc_H__
#define __MathFunc_H__

#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF2.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TGraphErrors.h> 
#include <TGraph.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TArrow.h>
#include "stdio.h"
#include <typedefs.h>
#include <vector>
#include <iostream>
#include <fstream>

Double_t ChiSquareDis(Double_t *x,Double_t *par);
void     DrawEnerPanel(TH1F *hData,TH1F *hMC,TCanvas *c,Int_t npad,TLegend *lleg,Double_t time);
Double_t Rate_Error(Float_t r1,Float_t r2,Float_t t1,Float_t b1,Float_t b2);
Double_t Bin_Errors(Double_t r1,Double_t t1,Double_t r2,Double_t t2);
void     Average_Array_Vector(Asym_t A,Int_t n, Double_t &y_Average,Double_t &y_Average_er);
void     Average_Array(Double_t *x,Double_t *xer,Int_t n,Double_t &y_Average,Double_t &y_Average_er);
void     CleanArray(Double_t x[][150],Int_t xbin,Int_t ybin);
void     CollectAllRot(Double_t x[][150],Int_t xbin,Int_t ybin,TH2F *h,TH2F *h2,Int_t last);
void     Average_All_Hists(TH1F *h1,std::vector<Double_t> &ave,std::vector<Double_t> &avee);
void     Return_Asymmetry(std::vector<Double_t> &ave,std::vector<Double_t> &avee,Int_t bins);
void     GetNorm(TH2F *hAveXRot,TH2F *hAveYRot,TLine *&lNorm,TArrow *&aDis,TCanvas *cCan,Int_t nline);
void     Get2DGaussianFit(TH2F *h2,Double_t &X,Double_t &Y);
Double_t Return_Rotation_Angle(TArrow *aDis,Double_t X,Double_t Y);
Int_t    GetInterSection(TLine *Line1,TLine *Line2,Double_t &Xin,Double_t &Yin);
Double_t PositionRotate(Double_t Angle,Double_t X,Double_t Y,TArrow *aDis);
Double_t BetaSpec(Double_t *x,Double_t *par);
Double_t PoverE(Double_t *x,Double_t *par);

void     ColorGraphic(TH1 *h, Int_t color, Int_t Style, Int_t Width,Double_t Size=1.);
void     ColorGraphic(TGraphErrors *h, Int_t color=1, Int_t Style=20, Int_t Width=1 ,Float_t Size = 1.,
                   const char *Title=" ",const char *XTitle=" ",const char *YTitle=" ");
void     ColorGraphic(TGraph *h, Int_t color=1, Int_t Style=20, Int_t Width=1,Float_t Size = 1.,
                   const char *Title=" ",const char *XTitle=" ",const char *YTitle=" ");
void     ColorGraphic(TH1 *h, Int_t color, Int_t Style, Int_t Width,Int_t nmarker,Double_t msize=1.);


Double_t GetSuperRatioAsymmetry(TH1F *R1M, TH1F *R1P,TH1F *R2M,TH1F *R2P,Int_t nl2,Int_t nh1);
Double_t GetSuperRatio(TH1F *R1M, TH1F *R1P,TH1F *R2M,TH1F *R2P,Int_t nl2,Int_t nh1);
Float_t GetSuperRatioEightAss(Float_t R1M, Float_t R1MI, Float_t R1P, Float_t R1PI, Float_t R2M, Float_t R2MI, Float_t R2P, Float_t R2PI);
Double_t GetSuperRatio(Double_t R1M,Double_t R1P,Double_t R2M,Double_t R2P);
Float_t GetSuperRatioEight(Float_t R1M, Float_t R1MI, Float_t R1P, Float_t R1PI, Float_t R2M, Float_t R2MI, Float_t R2P, Float_t R2PI);
Double_t GetSuperRatioAsymmetry(Double_t R1M,Double_t R1P,Double_t R2M,Double_t R2P);

Double_t GetSuperRatioError(TH1F *R1M, TH1F *R1P,TH1F *R2M, TH1F *R2P,Double_t time1,
                          Double_t time2,Int_t nl2,Int_t nh2);
Double_t GetSuperRatioError(Double_t R1M,Double_t R1P,Double_t R2M,Double_t R2P,
                          Double_t time1,Double_t time2);
Float_t GetDiff(Float_t R1M, Float_t R1P, Float_t R2M, Float_t R2P);
Float_t GetDifference(Float_t R1M, Float_t R1P);

void SetRate(Double_t *xR,Double_t *Rate,Double_t *RateEr,Double_t Run,
             Double_t mRate,Double_t mRateEr,Int_t &index);

void SetRateV(std::vector<Double_t> &xR,std::vector<Double_t> &Rate,std::vector<Double_t> &RateEr, Double_t Run,
	      Double_t mRate, Double_t mRateEr,Int_t nBck);

#endif
