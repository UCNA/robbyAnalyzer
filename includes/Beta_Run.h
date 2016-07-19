#ifndef betarun_h
#define betarun_h

#include "Run.h"
#include "Bck_Run.h"
 
class Beta_Run : public Run {
   
  public :
    Beta_Run(Int_t n,Int_t m,TSQLServer *sql);
    virtual ~Beta_Run();
    // Variables to hold the mwpc ratio calculations.......................
    Double_t MWPC_RatioE,MWPC_RatioW,MWPC_RatioE_fit,MWPC_RatioE_fit_e;
    Double_t MWPC_RatioE_e,MWPC_RatioW_e,MWPC_RatioW_fit,MWPC_RatioW_fit_e;
    Double_t Diff[160];
    Double_t Chi;
    Double_t Residual_Bkg[2],Resid_Bkger[2];
    Double_t BkgRtE[10],BkgRtW[10];
    Double_t BkgTe,BkgTw;
    Double_t btime_e,btime_w;
    Int_t    Bkg_index;
    Float_t  E_Endpoint,W_Endpoint;
    Float_t  E_EndError,W_EndError;
    Float_t  E_Sig_Nos,W_Sig_Nos;
    Double_t emuon_rate,wmuon_rate;

    // some functions that do shit.............................
    
    Int_t    Fill(Int_t n,Int_t remake,Double_t *sep,Int_t nrun);
    void     FillKurie(TH1F *hESpec,TH1F *hKurie);
    Int_t    Draw_2d(Int_t nr,Int_t n);
    Int_t    Draw_Hists(Int_t n);
    Int_t    stuff;
    Int_t    SubBck(Bck_Run *br);
    Bool_t   Check_Vetos();
    Int_t    GetBackGround(TSQLServer *sql);
    Int_t    GetSimpleAsym();
    Double_t GetEnergyChi();
    Int_t    Make_Kurie();
    void     Load_Histograms(Bck_Run *br,Bool_t SUBBCK, Int_t run);
    void     Remove_Histograms(Bck_Run *br);
    
};

#endif
