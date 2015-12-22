#ifndef bckrun_h
#define bckrun_h

#include "Run.h"

class Bck_Run : public Run {
   
  public :
    Bck_Run(Int_t n,Int_t m,TSQLServer *sql);
    virtual ~Bck_Run();
    
    //Functions
    Int_t Fill(Int_t n,Int_t remake,Double_t *sep,Int_t nrunl);
    Int_t Draw_2d(Int_t nr,Int_t n);
    Int_t Draw_Hists(Int_t n);
    Int_t stuff;
    Bool_t Check_Vetos();
    void Remove_Histograms();
    void Load_Histograms();
    
};
#endif
