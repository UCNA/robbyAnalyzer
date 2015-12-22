#ifndef __analysisPlotFunc_h__
#define __analysisPlotFunc_h__

#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include "stdio.h"
#include <iostream>
#include <fstream>
#include <algorithm>

Double_t DrawEnergyPanel(TH1F *hEdraw,TH1F *hERef,TH1F *hTot,TH1F *hRefTot,Int_t npad,TCanvas *c,Double_t scaling);
Double_t* Vector2Array(std::vector<Double_t> v,const Int_t sizev);
void DrawMWPCPanel1(TH1F *h1,TH1F *h2,TH1F *h3,TH1F *h4,TLegend *leg1,
                   Double_t ErrIntOff,Double_t ErrRatOff,Double_t ErrIntOn,Double_t ErrRatOn);



 
#endif 
