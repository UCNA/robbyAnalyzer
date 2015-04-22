#ifndef calrun_h
#define calrun_h

#include "Run.h"
#include "Bck_Run.h"

//-----------------------------------------------------------------
class Cal_Run : public Run {
  
 public:
  Cal_Run(Int_t n,Int_t m,Int_t b);
  virtual ~Cal_Run();
  
  struct peak_t{
    Float_t x;
    Float_t sig;
  };

  struct calp_t{
    Float_t adc;
    Float_t sig;
    Float_t err;
  };
  
  
  
  // Peak counters
  Int_t nex,nwx;
  Int_t nep,nwp;
  // Run time
  Float_t rtime;
  // Structures
  peak_t *expeaks,*wypeaks,*eypeaks,*wxpeaks;
  calp_t *ecal,*wcal;

  // Class Functions
  Int_t get_peak_pos(peak_t *peaks,TH1D* h1,Int_t &npeaks);
 
  Int_t Draw_2d(Int_t nr,Int_t n);
  Int_t Fit_Response(Int_t n);
  Int_t FillGraph(Int_t n);
  Int_t Fill(Int_t n); 
  Int_t Locate_peaks(Int_t n,Int_t nrun);
  Int_t Draw_Hists(Int_t n);
  Bool_t Cut_on_peak(Float_t x,Float_t y,Float_t sig1,Float_t sig2,Float_t xp,Float_t yp);
  Int_t Draw_Pos(Int_t nex, Int_t nwx);
  Int_t SubBck(Bck_Run *br);
  // Ellipses
  TEllipse *epw[5];
  TEllipse *wpw[5];

  // Graphs
  TGraphErrors *gcale,*gcalw;


 // ClassDef(Cal_Run,1);   // Calibration Run 
};
#endif
