#ifndef __variables_h__
#define __variables_h__

#include <fstream>
#include <Beta_Run.h>
#include <Run.h>
#include <Octet.h>

TCanvas *ct,*ciw,*cie;
 
Int_t runtot,nbeta,nbck,nbetatwo;
std::vector<Int_t>nbcks;
std::vector<Int_t>nbetas;
std::vector<Int_t>nrun;
std::vector<Int_t>ntype;
Int_t nrunb[2000];
Double_t sep23[100];
Int_t runopt,irun;
Int_t rty;
Int_t rb,runstart,runstop;
Int_t cnt2;
Int_t cnt3;
Int_t noct,noctstart,noctstop;
Int_t nlist;
Int_t NOctet;
Int_t numevents;
Int_t numbetas;
Int_t BetasWest,BetasEast;
Int_t BetasWestOne[2000],BetasEastOne[2000];
Int_t BkgdWestOne[2000],BkgdEastOne[2000];
Float_t CountTimeWAll,CountTimeEAll;
Float_t CountTimeWFirst,CountTimeEFirst;
Float_t CountTimeWBeta,CountTimeEBeta;
Float_t CountTimeWFirstBeta,CountTimeEFirstBeta;
Float_t FirstWTime[2000],FirstETime[2000];
Float_t FirstWTimeFive[2000],FirstETimeFive[2000];
Float_t FirstWTimeMore[2000],FirstETimeMore[2000];
Float_t FirstWTimeReal[2000],FirstETimeReal[2000];
Float_t FirstWTimePThree[2000],FirstETimePThree[2000];
Float_t LastWTime[2000],LastETime[2000];
Float_t LastWTimeFive[2000],LastETimeFive[2000];
Float_t LastWTimeMore[2000],LastETimeMore[2000];
Float_t LastWTimeReal[2000],LastETimeReal[2000];
Float_t LastWTimePThree[2000],LastETimePThree[2000];
Int_t lowbineast[2000],lowbinwest[2000];
Int_t lowbine[2000],lowbinw[2000];
Int_t highbineast[2000],highbinwest[2000];
Double_t SumLastPThreeE[2000],SumLastPThreeW[2000];
Double_t SumLastFiveE[2000],SumLastFiveW[2000];
Double_t RateLastFiveE[2000],RateLastFiveW[2000];
Double_t ProbLastFiveE[2000],ProbLastFiveW[2000];
Int_t SumLastFiveES[2000],SumLastFiveWS[2000];
Int_t lowbineastf[2000],lowbinwestf[2000];
Int_t highbineastf[2000],highbinwestf[2000];
Int_t highbinef[2000],highbinwf[2000];
Double_t SumFirstPThreeE[2000],SumFirstPThreeW[2000];
Double_t SumFirstFiveE[2000],SumFirstFiveW[2000];
Double_t RateFirstFiveE[2000],RateFirstFiveW[2000];
Double_t ProbFirstFiveE[2000],ProbFirstFiveW[2000];
Int_t SumFirstFiveES[2000],SumFirstFiveWS[2000];
Float_t FirstWBkgTime[2000],FirstEBkgTime[2000];
Float_t FirstWBkgTimeFive[2000],FirstEBkgTimeFive[2000];
Float_t FirstWBkgTimeMore[2000],FirstEBkgTimeMore[2000];
Float_t FirstWBkgTimeReal[2000],FirstEBkgTimeReal[2000];
Float_t FirstWBkgTimePThree[2000],FirstEBkgTimePThree[2000];
Float_t LastWBkgTime[2000],LastEBkgTime[2000];
Float_t LastWBkgTimeFive[2000],LastEBkgTimeFive[2000];
Float_t LastWBkgTimeMore[2000],LastEBkgTimeMore[2000];
Float_t LastWBkgTimeReal[2000],LastEBkgTimeReal[2000];
Float_t LastWBkgTimePThree[2000],LastEBkgTimePThree[2000];
Int_t lowbineastbkg[2000],lowbinwestbkg[2000];
Int_t lowbinebkg[2000],lowbinwbkg[2000];
Int_t highbineastbkg[2000],highbinwestbkg[2000];
Double_t SumLastPThreeEBkg[2000],SumLastPThreeWBkg[2000];
Double_t SumLastFiveEBkg[2000],SumLastFiveWBkg[2000];
Int_t lowbineastbkgf[2000],lowbinwestbkgf[2000];
Int_t highbineastbkgf[2000],highbinwestbkgf[2000];
Int_t highbinebkgf[2000],highbinwbkgf[2000];
Double_t SumFirstPThreeEBkg[2000],SumFirstPThreeWBkg[2000];
Double_t SumFirstFiveEBkg[2000],SumFirstFiveWBkg[2000];
Float_t EastClock[100000],WestClock[100000];
Float_t EastClockArray[2000][100000],WestClockArray[2000][100000];
Float_t EastClockBkgArray[2000][10000];
Float_t WestClockBkgArray[2000][10000];    
Double_t FFESIMPMin[2000], FFESIMLMin[2000];
Double_t FFWSIMPMin[2000], FFWSIMLMin[2000];
Double_t LFESIMPMax[2000], LFESIMLMax[2000];
Double_t LFWSIMPMax[2000], LFWSIMLMax[2000];
Double_t FirstFiveESIMP, FirstFiveESIML, FirstFiveWSIMP, FirstFiveWSIML;
Double_t LastFiveESIMP, LastFiveESIML, LastFiveWSIMP, LastFiveWSIML;
Double_t TSIMWP[2000], TSIMWL[2000], TSIMEP[2000], TSIMEL[2000];
Double_t SupSIM[2000];
Double_t FFESIMPMins[2000], FFESIMLMins[2000];
Double_t FFWSIMPMins[2000], FFWSIMLMins[2000];
Double_t LFESIMPMaxs[2000], LFESIMLMaxs[2000];
Double_t LFWSIMPMaxs[2000], LFWSIMLMaxs[2000];
Double_t FirstFiveESIMPS, FirstFiveESIMLS, FirstFiveWSIMPS, FirstFiveWSIMLS;
Double_t LastFiveESIMPS, LastFiveESIMLS, LastFiveWSIMPS, LastFiveWSIMLS;
Double_t TSIMWPS[2000], TSIMWLS[2000], TSIMEPS[2000], TSIMELS[2000];
Double_t SupSIMS[2000];
Double_t FirstFiveEAlP, FirstFiveEAlL, FEAlPSel, FEAlLSel;
Double_t FEAlP[2000], FEAlL[2000];
Double_t FirstFiveWAlP, FirstFiveWAlL, FWAlPSel, FWAlLSel;
Double_t FWAlP[2000], FWAlL[2000];
Double_t LastFiveEAlP, LastFiveEAlL, LEAlPSel, LEAlLSel;
Double_t LEAlP[2000], LEAlL[2000];
Double_t LastFiveWAlP, LastFiveWAlL, LWAlPSel, LWAlLSel;
Double_t LWAlP[2000], LWAlL[2000];
Double_t TAlWP[2000], TAlWL[2000], TAlEP[2000], TAlEL[2000];
Double_t SupAlSIM[2000];
Int_t mm[1000], p[1000], q[1000], r[1000];
Int_t ml[1000], pl[1000], ql[1000], rl[1000];
Double_t FirstFiveEAlPS, FirstFiveEAlLS, FEAlPSelS, FEAlLSelS;
Double_t FEAlPS[2000], FEAlLS[2000];
Double_t FirstFiveWAlPS, FirstFiveWAlLS, FWAlPSelS, FWAlLSelS;
Double_t FWAlPS[2000], FWAlLS[2000];
Double_t LastFiveEAlPS, LastFiveEAlLS, LEAlPSelS, LEAlLSelS;
Double_t LEAlPS[2000], LEAlLS[2000];
Double_t LastFiveWAlPS, LastFiveWAlLS, LWAlPSelS, LWAlLSelS;
Double_t LWAlPS[2000], LWAlLS[2000];
Double_t TAlWPS[2000], TAlWLS[2000], TAlEPS[2000], TAlELS[2000];
Double_t SupAlSIMS[2000];
Int_t ms[1000], ps[1000], qs[1000], rs[1000];
Int_t msl[1000], psl[1000], qsl[1000], rsl[1000];

std::vector<Double_t> erttt,ertet;
std::vector<Double_t> ertttf,ertetf;

fstream fsuprat1;

std::vector<Beta_Run*> btr;
std::vector<Bck_Run*> bckr;
std::vector<Octet*> octet;

Float_t WestDifference[2000],EastDifference[2000];
Float_t WestDiffFirst[2000],EastDiffFirst[2000];
/*Float_t LiveTimeBetaWest[2000], LiveTimeBetaEast[2000]; 
Float_t LiveTimeAllWest[2000], LiveTimeAllEast[2000];*/ 
Float_t DiffWestAllBeta[2000], DiffEastAllBeta[2000];
Float_t DiffBetaUp[2000],DiffBetaDown[2000];
Double_t SuperTiming[2000],SuperTimingA[2000];
Double_t Super[2000],SuperC[2000];
Double_t SuperE[2000],SuperCE[2000];
Double_t SuperPos1[2000],SuperPos2[2000],SuperPos3[2000],SuperPos4[2000];
Double_t SuperEPos1[2000],SuperEPos2[2000],SuperEPos3[2000],SuperEPos4[2000];
Double_t Xrun[2000];
Double_t Xer[2000];
Int_t    lookahead;
Int_t    Nsuper;

TH1F *hIn,*hOut;

TH1F *hA2TOTe,*hA5TOTe,*hA7TOTe,*hA10TOTe;
TH1F *hA2TOTw,*hA5TOTw,*hA7TOTw,*hA10TOTw;
TH1F *hB2TOTe,*hB5TOTe,*hB7TOTe,*hB10TOTe;
TH1F *hB2TOTw,*hB5TOTw,*hB7TOTw,*hB10TOTw;

TH1F *hEFlipperOn,*hEFlipperOff,*hWFlipperOn,*hWFlipperOff;
TH1F *hEFlipperOn_G,*hEFlipperOff_G,*hWFlipperOn_G,*hWFlipperOff_G;
TH1F *hEFlipperOn_S,*hEFlipperOff_S,*hWFlipperOn_S,*hWFlipperOff_S;
TH1F *hEFlipperOn_I,*hEFlipperOff_I,*hWFlipperOn_I,*hWFlipperOff_I;
TH1F *hEFlipperOn_2,*hEFlipperOff_2,*hWFlipperOn_2,*hWFlipperOff_2;

TH1F *hEFlipperOn_B,*hEFlipperOff_B,*hWFlipperOn_B,*hWFlipperOff_B;
TH1F *hEFlipperOn_B_I,*hEFlipperOff_B_I,*hWFlipperOn_B_I,*hWFlipperOff_B_I;
TH1F *hEFlipperOn_B_2,*hEFlipperOff_B_2,*hWFlipperOn_B_2,*hWFlipperOff_B_2;

TH2F *hTotRote,*hTotRotw,*hTotRoteI,*hTotRotwI;
TH2F *hTotRote23,*hTotRotw23,*hTotRoteI23,*hTotRotwI23;

TH2F *hEvWtot,*hEvWtotG,*hTotETimeVsE,*hTotWTimeVsE;
TH2F *hTotEPos,*hTotWPos,*hTE23Anode2d, *hTW23Anode2d;

// Collection Graphs
TGraphErrors *gre,*grw,*gre1,*grw1,*gre23,*grw23,*grtot,*grstot;
TGraphErrors *gref,*grwf,*gre1f,*grw1f,*gre23f,*grw23f,*grtotf,*grstotf;

// Mysql Server
TSQLServer *sql;

#endif
