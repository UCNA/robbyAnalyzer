#ifndef __analysis_c__
#define __analysis_c__

#include "analysis.h"
#include <algorithm>

void SetPlotOptions(TH1* h1,Double_t tsize,Double_t offset,Double_t lsize,Double_t lxsize);


int main()
{
  // Set up the ROOT environment------------------------------
  TApplication theApp("App",0,0);
  Int_t remake = 0;
  SetROOTOpt();
  ZeroThings();
  // Determine the type of analysis to perform

  //  cout << "Single Run or from list (1 of single, !=1 for list" << endl;
  //cin >> runopt;
  runopt = 2;
  runstop = 0;
  if(runopt == 1)
   GetSingleRun();
  else {
    do{
    cout << "Geometry A   : 7659  - 9189" << endl;
    cout << "Geometry B   : 9333  - 10333" << endl;
    cout << "Geometry C   : 10404 - 11095 "<< endl; 
    cout << "Geometry C09 : 12430 - 12700 " << endl;
    cout << "Adding List of Runs"  << endl;
    cout << "Enter first run : ";
    cin >> runstart;
    //cout << "Enter last run : ";
    //cin >> runstop;
    }while(runstart > 4);
    
    cout << "Remake Histograms ? ( 1 = Yes , 0 = No)" << endl;
    cin >> remake;
    //remake = 0;
    GetListofRuns(runstart,runstop);
  }
  
  fstream ref;
  Int_t j1,j2,j3,j4;
  
  // this get the file for the type23 separation fits from Monte Carlo....
  
  if(runstart == 1){
    nlist = 1;
    ref.open("type_23_output_a.txt",fstream::in);
  }
  else if(runstart == 2){
    nlist = 2;
    ref.open("type_23_output_b.txt",fstream::in);
  }
  else if(runstart == 3){
    nlist = 3;
    ref.open("type_23_output_c.txt",fstream::in);
  }
  else if(runstart == 4){
    ref.open("type_23_output_c.txt",fstream::in);
    nlist = 4;
  }
 
  for(Int_t i = 0 ; i < 100 ; i++){
    ref >> j1 >> sep23[i] >> j2 >> j3 >> j4;
  }
  ref.close();
  
  cout << "Finished Loading Run List from MySQL Database" << endl;
  cout << "Octet list is " << nlist << " There are " << nbck << " runs " << nbeta <<  endl;
  //--------------------------------------------------------------------------------------
  // Analyze  Run Data
  if(nbck  > 0) analyze_background_runs(nbck, bckr, nbcks,remake);
  if(nbeta > 0) analyze_beta_runs(nbeta, btr, nbetas,remake);
  //--------------------------------------------------------------------------------------
  cout << "Finished Analysis............" << endl;
  // Now Analyze generated histograms ....
  if(nbeta > 0){    
    Subtract_Backgrounds(btr,bckr,nbeta,nbck);
    cout << "done with backgrounds" << endl;
    analyze_octets();
  }
  
  cout << "done with the background subtraction" << endl;
  
  for(Int_t jrun = 0 ; jrun < nbeta ; jrun++){
     btr[jrun]->Calculate_Backscatter(1,1);
     btr[jrun]->GetSimpleAsym();
  }
  
  cout << "Collecting Rates " << endl;
  CollectRates();
  cout << "Collecting Type Rotation " << endl;
  CollectTypeRot();
  cout << "Collecting TDC Corruption Data " << endl;
  // CollectTDCCor();
  cout << "Lets Find the Octets " << endl;
  CalcSimplSuper();
  cout << "Plot Energy Chis" << endl;
  Plot_E_Chis();
  cout << "Plot Type 2/3 Difference" << endl;
  Plot_MonRate();
  cout << "Plotting Positions" << endl;
  Collect_Pos();  
  cout << "Clllect Octets" << endl;
  Collect_Octets();
  cout << "Finished Octets" << endl;
  Average_A();
  Collect_TvsE();
  Collect_23Anode();
  Collect_Stuff();
  Collect_Rad();
  Collect_Energy_Spectra();
  cout << "here " << endl;
  GammaBack();
  cout << " returned " << endl;
  Collect_Gammas();
  Collect_TDCDiff();
  Plot_ChiDis();
  average_type1();
  
  theApp.Run();
  return 0;
}
//---------------------------------------------------------------------------
Int_t analyze_background_runs(Int_t n, Bck_Run **bk, Int_t *nrun,Int_t remake)
{
// Analysis flow for background runs.
  for(Int_t i = 0 ; i < n ; i++){
    // bk[i]->Draw_2d(nrun[i],i);
    if(remake == 1)bk[i]->Find_TDC_Cut(i);
    bk[i]->Set_Anode_Cut(n);
    bk[i]->Fill(i,remake,sep23);
    bk[i]->Scale2Time(1,1);
    bk[i]->Calculate_Backscatter(1,1);
    if(bk[i]->GetRunNumber() == 8296){
      cout << "Passed cut 1 " << bk[i]->passed1 << endl;
      cout << "Passed cut 2 " << bk[i]->passed2 << endl;
      cout << "Passed cut 3 " << bk[i]->passed3 << endl;
      cout << "Passed cut 4 " << bk[i]->passed4 << endl;
    }
    cout << "At Run " << bk[i]->GetRunNumber() << endl;
 }
 return 0;

}
//---------------------------------------------------------------------------
Int_t analyze_beta_runs(Int_t n, Beta_Run **bta, Int_t *nrun,Int_t remake)
{
  
  for(Int_t i = 0 ; i < n ; i++){
    // bta[i]->Draw_2d(nrun[i],i);
    if(remake==1)bta[i]->Find_TDC_Cut(i);
    bta[i]->Set_Anode_Cut(n);
    bta[i]->Fill(i,remake,sep23);
    bta[i]->Scale2Time(1,1);
  }

 return 0;
}
//---------------------------------------------------------------------------
Int_t Subtract_Backgrounds(Beta_Run **bta,
			   Bck_Run **bk, Int_t nb,Int_t nbk)
{
  // Loop through the beta runs and subtract off the background spectra
  Int_t i,j,k;
  
  for(Int_t n = 0 ; n < nb ; n++){
    i = bta[n]->GetBackgroundRun();
    j = 0;
    k = -1;
    if(i!=0){
      do{
        k++;
	if(k == nbk) break;
        j = bk[k]->GetRunNumber();
      } while(i != j && k < nbk);
      if(i == j){
	bta[n]->Bkg_index = k;
	bta[n]->SubBck(bk[k]);
	bta[n]->GetEnergyChi();
	bta[n]->Make_Kurie();
      }
    }
  }
  return 0;
}
//---------------------------------------------------------------------------
Int_t analyze_octets()
{
  fstream runlist;
  cout << "getting octet list " << nlist << endl;
  runlist.open(Form("octet_list_%d.txt",nlist),fstream::in);
  // Here you would build the Octets ......
  Int_t noctets;
  noct = 0;
  
  if(runlist.is_open()){
    do{
      runlist >> noctets >> noctstart >> noctstop;
      cout << noctets << "\t" << noctstart << "\t" << noctstop << endl;      
      if(noctstart >= runstart && noctets != 1111){
        // Open the Octet object
        octet[noct] = new Octet(noct,noctstart,noctstop);
        // define histograms
        octet[noct]->Initialize_Histo(noct);
        // Use the MySQL database to determine the runs in each octet and fill the octet
        // defined runs.  If there are multiple DAQ runs for any such octet run they will
        // be summed such that the times are correctly counted and the rates are interms of
        // the total counts of all re
        octet[noct]->Find_Runs(btr,nbeta);
	// Need to adjust the octect calculation for the full octet analysis.
	octet[noct]->Calc_Super();
	octet[noct]->Calc_A_sum();
	octet[noct]->Calc_A_multi();
	octet[noct]->Calc_A_sum_Bin();
	octet[noct]->Get_Rad_A();
	octet[noct]->Debugger(btr[0]->GetGeo());
	noct++;
      }
    }while(!(runlist.eof()));
  }

  runlist.close(); 
  
  return -1;
}
//---------------------------------------------------------------------------------------
void GetSingleRun()
{
  cout << "Enter Run Number:";
  cin >> irun;
  cout << "Enter Run type  2 is beta decay ; 3 is background  : ";
  cin >> rty;  // 1 is calibration, 2 is beta decay ; 3 is background

  nrun[0]  = irun;
  ntype[0] = rty;

  if(rty == 2){
    nbeta++;
    nbck++;
    cout << "Enter Background :  " ;
    cin >> rb;
    nrun[1]   = rb;
    ntype[1]  = 3;
    nbetas[0] = irun;
    nrunb[0]  = rb;
    runtot    = 2;
    nbcks[0]  = rb;
  }
  if(rty == 3){
    nbck++;
    runtot   = 1;
    nbcks[0] = irun;
    nrunb[0] = irun;
  }

  return;
}
//---------------------------------------------------------------------------
void GetListofRuns(Int_t n1,Int_t n2)
{
  
  char db[100];
  fstream f;
  f.open("mysql_db.txt",fstream::in);
  f >> db;
  f.close();
  
  TMySQLServer *sql = new TMySQLServer(db,"ucnawrite","UCNBetAs!");  

  if(!(sql->IsConnected())){
    while(!(sql->IsConnected())){
      cout << "Having to retry connection (List) "<< n1 << "\t" << n2 <<  endl;
      sql->Connect(db,"ucn","UCNBetAs");
    }
  }

  TMySQLResult *res;
  TMySQLRow *row;
  char buffer[800];
  Int_t nst =0;
  Int_t nend=0;
  switch (n1){
  case 1 :
    nst = 7656;
    nend =9190;
    break;
  case 2 :
    nst = 9332;
    nend =10334;
    break;
  case 3 :
    nst = 10403;
    nend =11096;
    break;
  case 4 :
    nst = 12430;
    nend =12701;
    break;
  }
  /*
    cout << "Geometry A   : 7659  - 9189" << endl;
    cout << "Geometry B   : 9333  - 10333" << endl;
    cout << "Geometry C   : 10404 - 11095 "<< endl; 
    cout << "Geometry C09 : 12430 - 12700 " << endl;
  

  old query from initial analysis doesn't match the structure of the new DB

  sprintf(buffer,"select t1.run_number,t1.geometry,t1.asym_oct,t1.gate_valve"
	  ",t2.live_time_e from run as t1, analysis as t2 where t1.run_number between %d and %d and t1.run_type "
	  " = 'Asymmetry' and (t1.run_number < 9497 or t1.run_number > 9500) and "
	  " t1.run_number = t2.run_number"
	  ,n1,n2);
*/
  cout << "Nlist " << n1 << endl;
  /*  if(n1 == 1){
    sprintf(buffer,"select runnum,time from run_info where geotype = %d and octnumb >= 0 and (runnum < 9497 or"
	  " runnum > 9500) and (runnum < 9000 or runnum >9200) and (octtype = 'a2' || octtype = 'a5' || octtype = 'a7' "
	  "|| octtype = 'a10' ||  octtype = 'b2' || octtype = 'b5' || (octtype = 'b7' and octnumb != 4)  "
	  "|| octtype = 'b10') and time > 700",n1);
  } else if(n1 < 4 && n1 > 1){
    sprintf(buffer,"select t1.runnum,t1.time from run_info as t1, single_run as t2 where t1.geotype = %d and t1.octnumb >= 0 and (t1.runnum < 9497 or"
	  " t1.runnum > 9500) and (t1.runnum < 9000 or t1.runnum >9200) and (t1.octtype = 'a2' || t1.octtype = 'a5' || t1.octtype = 'a7' "
	  "|| t1.octtype = 'a10' ||  t1.octtype = 'b2' || t1.octtype = 'b5' || t1.octtype = 'b7' "
	  "|| t1.octtype = 'b10') and t1.time > 700 && t1.runnum = t2.run_number",n1);
	  } else {*/
  sprintf(buffer,"select run_number,time_w from single_run where "
	  "(tag = 'A2' || tag = 'A5' || tag = 'A7' "
	  "|| tag = 'A10' ||  tag = 'B2' || tag = 'B5' || tag = 'B7' "
	  "|| tag = 'B10') and time_w > 500 and run_number between %d and %d",nst,nend);
  // }
  cout << buffer << endl;
  res = (TMySQLResult*)sql->Query(buffer);
  
  int rn = 0;
  int chc = 0;
   if(res->GetRowCount() != 0){
     while((row = (TMySQLRow*)res->Next())){
       rn = atoi(row->GetField(0));
       chc=0;
       for(int i =0 ; i < nbeta;i++){
	 if(nbetas[i] == rn)chc++;
       }
       if(rn!= 7945 && rn!= 7946&& rn != 9812 && rn !=9007 && rn != 9808 && rn != 8539 && rn !=9010 && chc == 0) {
	if( atof(row->GetField(1)) > 0){
	  //----------------------------------------------
	  // Open a beta decay run
          cout << "At run " << rn << endl;
	  nbetas[nbeta] = rn;// atoi(row->GetField(0));
	  nrun[runtot]  = nbetas[nbeta];
	  ntype[runtot] = 2;
	  btr[nbeta]    = new Beta_Run(nbetas[nbeta],2);
	
	  nbeta++;
	  runtot++;

	  //-----------------------------------------------
	  // Open a background run for the above beta run

	  nbcks[nbck]   = btr[nbeta-1]->GetBackgroundRun();
	  cout << "found background run " << nbcks[nbck] <<  endl;
	  nrun[runtot]  = nbcks[nbck];
	  nrunb[nbck]   = nbcks[nbck];
	  ntype[runtot] = 3;
	  bckr[nbck]    = new Bck_Run(nbcks[nbck],3);
	  nbck++;
	  runtot++;
	}
      }
      delete row; 
     }
   } 
   
  delete res;
  delete sql;
  
}
//---------------------------------------------------------------------------
void SetROOTOpt()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(kFALSE); 
  gStyle->SetPadLeftMargin(0.2);
  gROOT->ForceStyle();
  gStyle->SetTextSize(1.01);
  gStyle->SetLabelSize(0.055,"xy");
  gStyle->SetTitleSize(0.06,"xy");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.16);
  
}
//---------------------------------------------------------------------------
void ZeroThings()
{
  
  nbck   = 0;
  nbeta  = 0;
  runtot = 0;
  cnt2   = 0;
  cnt3   = 0;
}
//---------------------------------------------------------------------------
void CollectRates()
{
  using namespace TMath;
  
  Double_t xrune[10000],xrunee[10000];
  Double_t xrunw[10000],xrunwe[10000];
  Double_t xrunw1[10000],xrunw1e[10000];
  Double_t xrunw23[10000],xrune1e[10000];
  Double_t xrune1[10000],xrune23e[10000];
  Double_t xrune23[10000],xrunw23e[10000];
  
  Double_t ertt[10000],erte[10000],wrtt[10000],wrte[10000];
  Double_t ert1[10000],er1e[10000],wrt1[10000],wr1e[10000];
  Double_t ert23[10000],er23e[10000],wrt23[10000],wr23e[10000];
  
  Double_t erttf[10000],ertef[10000],wrttf[10000],wrtef[10000];
  Double_t ert1f[10000],er1ef[10000],wrt1f[10000],wr1ef[10000];
  Double_t ert23f[10000],er23ef[10000],wrt23f[10000],wr23ef[10000];
  
  Double_t xrunef[10000],xruneef[10000];
  Double_t xrunwf[10000],xrunwef[10000];
  Double_t xrunw1f[10000],xrunw1ef[10000];
  Double_t xrunw23f[10000],xrune1ef[10000];
  Double_t xrune1f[10000],xrune23ef[10000];
  Double_t xrune23f[10000],xrunw23ef[10000];
  
  Int_t v1 = 0;
  Int_t v2 = 0; 
  Int_t v3 = 0;
  Int_t v4 = 0;
  Int_t v5 = 0; 
  Int_t v6 = 0;
  
  Int_t v1f = 0;
  Int_t v2f = 0; 
  Int_t v3f = 0;
  Int_t v4f = 0;
  Int_t v5f = 0; 
  Int_t v6f = 0;
   
  // Loop through background subtracted beta runs for 
  
  for(Int_t jrun = 0 ; jrun < nbeta ; jrun++){
    if(btr[jrun]->flipperOn == 1){
      xrune[v1]  = btr[jrun]->GetRunNumber();
      xrunee[v1] = 0.01;
      
      ertt[v1]   = btr[jrun]->bscat.e_all;
      erte[v1]   = TMath::Abs(btr[jrun]->bscat.e_alle);
      
      v1++;
      
      xrunw[v2] = btr[jrun]->GetRunNumber();
      xrunwe[v2] = 0.01;
      
      wrtt[v2] = btr[jrun]->bscat.w_all;
      wrte[v2] = TMath::Abs(btr[jrun]->bscat.w_alle);
      
      v2++;
      
      erttt[v1-1] = ertt[v1-1] + wrtt[v2-1];
      ertet[v1-1] = sqrt( wrte[v2-1]*wrte[v2-1] + erte[v1-1]*erte[v1-1]);
      
      xrunw1[v3] = btr[jrun]->GetRunNumber();
      xrunw1e[v3] = 0.01;
      
      wrt1[v3] = (btr[jrun]->bscat.wtype_1 / 
	   (1e-6+btr[jrun]->bscat.w_all)) * 100.;
      
      wr1e[v3] = wrt1[v3]*Rate_Error(btr[jrun]->bscat.wtype_1,
					 btr[jrun]->bscat.w_all,
                                         btr[jrun]->rtime_w,
                                         btr[jrun]->bscat.wtype_1_bck,
                                         btr[jrun]->bscat.w_all_bck);
      
      if(wr1e[v3] > 0.05)v3++;
      
      xrune1[v4] = btr[jrun]->GetRunNumber();
      xrune1e[v4] = 0.01;
      
      ert1[v4] = (btr[jrun]->bscat.etype_1 /(1e-6+ btr[jrun]->bscat.e_all)) * 100.;
      
      er1e[v4] = ert1[v4]*Rate_Error(btr[jrun]->bscat.etype_1,
				     btr[jrun]->bscat.e_all,
	                             btr[jrun]->rtime_e,
	                             btr[jrun]->bscat.etype_1_bck,
                                     btr[jrun]->bscat.e_all_bck);
      
      if(er1e[v4] > 0.05)v4++;
      
      
      xrunw23[v5] = btr[jrun]->GetRunNumber();
      xrunw23e[v5] = 0.01;
      
      wrt23[v5] = (btr[jrun]->bscat.wtype_23 / (1e-6+btr[jrun]->bscat.w_all)) * 100.;
      
      wr23e[v5] = wrt23[v5]*Rate_Error(btr[jrun]->bscat.wtype_23,
					   btr[jrun]->bscat.w_all,
	                                   btr[jrun]->rtime_w, 
                                           btr[jrun]->bscat.wtype_23_bck,
                                           btr[jrun]->bscat.w_all_bck);
      
      if(wr23e[v5] > 0.005) v5++;
      
      xrune23[v6] = btr[jrun]->GetRunNumber();
      xrune23e[v6] = 0.01;
      
      
      ert23[v6] = (btr[jrun]->bscat.etype_23 / (1e-6+btr[jrun]->bscat.e_all)) * 100.;
      
      er23e[v6] = ert23[v6]*Rate_Error(btr[jrun]->bscat.etype_23,
					   btr[jrun]->bscat.e_all,
	                                   btr[jrun]->rtime_e, 
				           btr[jrun]->bscat.etype_23_bck,
                                           btr[jrun]->bscat.e_all_bck);
      
      if(er23e[v6] > 0.005)v6++;
      
    } else if(btr[jrun]->flipperOn == 0){
      xrunef[v1f]  = btr[jrun]->GetRunNumber();
      xruneef[v1f] = 0.01;
      
      erttf[v1f]   = btr[jrun]->bscat.e_all;
      ertef[v1f]   = TMath::Abs(btr[jrun]->bscat.e_alle);
      
      v1f++;
      
      xrunwf[v2f] = btr[jrun]->GetRunNumber();
      xrunwef[v2f] = 0.01;
      
      wrtt[v2f] = btr[jrun]->bscat.w_all;
      wrte[v2f] = TMath::Abs(btr[jrun]->bscat.w_alle);
      
      v2f++;
      
      ertttf[v1f-1] = ertt[v1f-1] + wrtt[v2-1];
      ertetf[v1f-1] = sqrt( wrte[v2f-1]*wrte[v2f-1] + erte[v1f-1]*erte[v1f-1]);
      
      xrunw1f[v3f] = btr[jrun]->GetRunNumber();
      xrunw1ef[v3f] = 0.01;
      
      wrt1f[v3f] = (btr[jrun]->bscat.wtype_1 / 
		  (1e-6+btr[jrun]->bscat.w_all)) * 100.;
      
      wr1ef[v3f] = wrt1f[v3f]*Rate_Error(btr[jrun]->bscat.wtype_1,
					 btr[jrun]->bscat.w_all,
                                         btr[jrun]->rtime_w,
                                         btr[jrun]->bscat.wtype_1_bck,
                                         btr[jrun]->bscat.w_all_bck);
      
      if(wr1e[v3f] > 0.005)v3f++;
      
      xrune1f[v4f] = btr[jrun]->GetRunNumber();
      xrune1ef[v4f] = 0.01;
      
      ert1f[v4f] = (btr[jrun]->bscat.etype_1 /(1e-6 + btr[jrun]->bscat.e_all)) * 100.;
      
      er1ef[v4f] = ert1f[v4f]*Rate_Error(btr[jrun]->bscat.etype_1,
				     btr[jrun]->bscat.e_all,
	                             btr[jrun]->rtime_e,
	                             btr[jrun]->bscat.etype_1_bck,
                                     btr[jrun]->bscat.e_all_bck);
      
      if(er1ef[v4f] > 0.005)v4f++;
      
      
      xrunw23f[v5f] = btr[jrun]->GetRunNumber();
      xrunw23ef[v5f] = 0.01;
      
      wrt23f[v5f] = (btr[jrun]->bscat.wtype_23 / (1e-6+btr[jrun]->bscat.w_all)) * 100.;
      
      wr23ef[v5f] = wrt23f[v5f]*Rate_Error(btr[jrun]->bscat.wtype_23,
					   btr[jrun]->bscat.w_all,
	                                   btr[jrun]->rtime_w, 
                                           btr[jrun]->bscat.wtype_23_bck,
                                           btr[jrun]->bscat.w_all_bck);
      
      if(wr23ef[v5f] > 0.005) v5f++;
      
      xrune23f[v6f] = btr[jrun]->GetRunNumber();
      xrune23ef[v6f] = 0.01;
      
      
      ert23f[v6f] = (btr[jrun]->bscat.etype_23 / (1e-6+btr[jrun]->bscat.e_all)) * 100.;
      
      er23ef[v6f] = ert23f[v6f]*Rate_Error(btr[jrun]->bscat.etype_23,
					   btr[jrun]->bscat.e_all,
	                                   btr[jrun]->rtime_e, 
				           btr[jrun]->bscat.etype_23_bck,
                                           btr[jrun]->bscat.e_all_bck);
      
      if(er23ef[v6f] > 0.005)v6f++;
    }
      

  }

  grtot = new TGraphErrors(v1,xrune,erttt,xrunee,ertet); 
  gre   = new TGraphErrors(v1,xrune,ertt,xrunee,erte);
  grw   = new TGraphErrors(v2,xrunw,wrtt,xrunwe,wrte);
  gre1  = new TGraphErrors(v4,xrune1,ert1,xrune1e,er1e);
  grw1  = new TGraphErrors(v3,xrunw1,wrt1,xrunw1e,wr1e);
  gre23 = new TGraphErrors(v6,xrune23,ert23,xrune23e,er23e);
  grw23 = new TGraphErrors(v5,xrunw23,wrt23,xrunw23e,wr23e);
 
  cout << v3f << "  " << v4f << endl;

  grtotf = new TGraphErrors(v1f,xrunef,ertttf,xruneef,ertetf); 
  gref   = new TGraphErrors(v1f,xrunef,erttf,xruneef,ertef);
  grwf   = new TGraphErrors(v2f,xrunwf,wrttf,xrunwef,wrtef);
  gre1f  = new TGraphErrors(v4f,xrune1f,ert1f,xrune1ef,er1ef);
  grw1f  = new TGraphErrors(v3f,xrunw1f,wrt1f,xrunw1ef,wr1ef);
  gre23f = new TGraphErrors(v6f,xrune23f,ert23f,xrune23ef,er23ef);
  grw23f = new TGraphErrors(v5f,xrunw23f,wrt23f,xrunw23ef,wr23ef);
  
  fe1 = new TF1("fe1","[0]",runstart,runstop);
  fe2 = new TF1("fe2","[0]",runstart,runstop);
  fw1 = new TF1("fw1","[0]",runstart,runstop);
  fw2 = new TF1("fw2","[0]",runstart,runstop);
  
  fe1->SetLineColor(3);
  fe2->SetLineColor(3);
  
  fw1->SetLineColor(4);
  fw2->SetLineColor(4);
  
  TCanvas *cgr = new TCanvas("cgr","Total Rates");
  //cgr->Divide(1,3);
  cgr->cd(1);
  
  gre->SetTitle("Total Event Rates");
  gre->GetYaxis()->SetTitle("Total Rate s^{-1}");
  gre->GetXaxis()->SetTitle("Run Number");
  
  ColorGraphic(gre,2,20,2);
  ColorGraphic(grw,1,20,2);
  ColorGraphic(gre1,3,20,2);
  ColorGraphic(grw1,4,20,2);
  ColorGraphic(gre23,3,20,2);
  ColorGraphic(grw23,4,20,2);

  ColorGraphic(gref,2,4,2);
  ColorGraphic(grwf,1,4,2);
  ColorGraphic(gre1f,2,4,2);
  ColorGraphic(grw1f,1,4,2);
  ColorGraphic(gre23f,2,4,2);
  ColorGraphic(grw23f,1,4,2);
  
  gre->Draw("AP");
  grw->Draw("P");
  
  TLegend *lRates = new TLegend(0.6,0.7,0.9,0.9);
  lRates->AddEntry(gre,"East Rate","lp");
  lRates->AddEntry(grw,"West Rate","lp");
  lRates->Draw();
 
  TCanvas *cgrb = new TCanvas("cgrb","Backscattering");
  
  cgrb->Divide(1,2);
  cgrb->cd(1);
  gre1->Draw("AP");
  gre1->SetTitle("Type 1 Event Fractions");
  gre1->GetYaxis()->SetTitle("Fraction");
  gre1->GetXaxis()->SetTitle("Run Number");
  if(v3>0)gre1->Fit("fe1","REMQ+");
  grw1->Draw("P");
  gre1f->Draw("P");

  if(v4>0)grw1->Fit("fw1","REMQ+");
  grw1f->Draw("P");
  TLegend *lb1 = new TLegend(0.5,0.6,0.9,0.9);
  lb1->AddEntry(gre1,Form("East %5.3f #pm %5.3f",
		fe1->GetParameter(0),fe1->GetParError(0)),"lp");
      
  lb1->AddEntry(grw1,Form("West %5.3f #pm %5.3f", 
		fw1->GetParameter(0),fw1->GetParError(0)),"lp");
  
  lb1->SetFillColor(kWhite);
  lb1->Draw();
      
  cgrb->cd(2);
  
  gre23->SetTitle("Type 2/3 Event Fractions");
  gre23->GetYaxis()->SetTitle("Fraction");
  gre23->GetXaxis()->SetTitle("Run Number");
  gre23->Draw("AP");
  
  if(v5 >0)
    gre23->Fit("fe2","REMQ+");
  grw23->Draw("P");
  if(v6>0)
    grw23->Fit("fw2","REMQ+");
  grw23f->Draw("P");
  gre23f->Draw("P");
  TLegend *lb2 = new TLegend(0.5,0.6,0.9,0.9);
  if(v6 > 0 && v5 > 0){
    lb2->AddEntry(gre1,Form("East %5.3f #pm %5.3f",
  		fe2->GetParameter(0),fe2->GetParError(0)),"lp");
      
    lb2->AddEntry(grw1,Form("West %5.3f #pm %5.3f", 
		fw2->GetParameter(0),fw2->GetParError(0)),"lp");
  
    lb2->SetFillColor(kWhite);
    lb2->Draw();
  }
  gre->GetYaxis()->SetRangeUser(0,10);
  gre1->GetYaxis()->SetRangeUser(0,3);
  
}
//---------------------------------------------------------------------------
Double_t Rate_Error(Float_t r1,Float_t r2,Float_t t1,Float_t b1,Float_t b2)
{
  using namespace TMath;
  
  Double_t er;
  
  er  = (r1*t1 + b1)/Power((r1*t1 - b1),2);
  er += (r2*t1 + b2)/Power((r2*t1 - b2),2);
     
  er = sqrt(Abs(er));
  
 return er;
}
//---------------------------------------------------------------------------
void CollectTypeRot()
{
  Int_t rbins = 150;
  Float_t rad = 75.;
  Float_t radcut = 60.;
  Double_t ecount[rbins][rbins],wcount[rbins][rbins];
  
  hTotRote = new TH2F("hTotRote",
		   "East Type 1 Backscattering Rotation ; X_{west} (mm) ; Y_{west} (mm)"
		      ,rbins,-rad,rad,rbins,-rad,rad);
  
  hTotRotw = new TH2F("hTotRotw",
		      "West Type 1 Backscattering Rotation ; X_{east} (mm) ; Y_{east} (mm)"
		      ,rbins,-rad,rad,rbins,-rad,rad);
  
  for(Int_t ibin = 0 ; ibin < rbins ; ibin++){
    for(Int_t jbin = 0 ; jbin < rbins ; jbin++){
      ecount[ibin][jbin] = 0.;
      wcount[ibin][jbin] = 0.;
    }
  }
  
  for(Int_t i = 0 ; i <nbeta ; i++){
    for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
	ecount[ibin-1][jbin-1] += btr[i]->hRote->GetBinContent(ibin,jbin);
	wcount[ibin-1][jbin-1] += btr[i]->hRotw->GetBinContent(ibin,jbin);
      }
    }
  }
  
  for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
    for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
      hTotRote->SetBinContent(ibin,jbin,ecount[ibin-1][jbin-1]);
      hTotRotw->SetBinContent(ibin,jbin,wcount[ibin-1][jbin-1]);
    }
  }
  
  hTotRoteI = new TH2F("hTotRoteI",
		      "East Type Trigger Side Position ; X_{east} (mm) ; Y_{east} (mm)",
		       rbins,-rad,rad,rbins,-rad,rad);
  
  hTotRotwI = new TH2F("hTotRotwI",
		      "West Type Trigger Side Position ; X_{west} (mm) ; Y_{west} (mm)",
		       rbins,-rad,rad,rbins,-rad,rad);
  
  for(Int_t ibin = 0 ; ibin < rbins ; ibin++){
    for(Int_t jbin = 0 ; jbin < rbins ; jbin++){
      ecount[ibin][jbin] = 0.;
      wcount[ibin][jbin] = 0.;
    }
  }
  
  for(Int_t i = 0 ; i <nbeta ; i++){
    for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
	ecount[ibin-1][jbin-1] += btr[i]->hRoteI->GetBinContent(ibin,jbin);
	wcount[ibin-1][jbin-1] += btr[i]->hRotwI->GetBinContent(ibin,jbin);
      }
    }
  }
  
  for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
    for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
      hTotRoteI->SetBinContent(ibin,jbin,ecount[ibin-1][jbin-1]);
      hTotRotwI->SetBinContent(ibin,jbin,wcount[ibin-1][jbin-1]);
    }
  }
  
  
  
  hTotRote23 = new TH2F("hTotRote23",
			"East Type 23 Backscattering Rotation ; X_{west} (mm) ; Y_{west} (mm)"
			,rbins,-rad,rad,rbins,-rad,rad);
  
  hTotRotw23 = new TH2F("hTotRotw23",
			"West Type 23 Backscattering Rotation ; X_{east} (mm) ; Y_{east} (mm)"
			,rbins,-rad,rad,rbins,-rad,rad);
  
  for(Int_t ibin = 0 ; ibin < rbins ; ibin++){
    for(Int_t jbin = 0 ; jbin < rbins ; jbin++){
      ecount[ibin][jbin] = 0.;
      wcount[ibin][jbin] = 0.;
    }
  }
  
  for(Int_t i = 0 ; i <nbeta ; i++){
    for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
	ecount[ibin-1][jbin-1] += btr[i]->hRote23->GetBinContent(ibin,jbin);
	wcount[ibin-1][jbin-1] += btr[i]->hRotw23->GetBinContent(ibin,jbin);
      }
    }
  }
  
  for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
    for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
      hTotRote23->SetBinContent(ibin,jbin,ecount[ibin-1][jbin-1]);
      hTotRotw23->SetBinContent(ibin,jbin,wcount[ibin-1][jbin-1]);
    }
  }

  
  hTotRoteI23 = new TH2F("hTotRoteI23",
			 "East Type 23 trigger side Position ; X_{east} (mm) ; Y_{east} (mm)"
			 ,rbins,-rad,rad,rbins,-rad,rad);
  
  hTotRotwI23 = new TH2F("hTotRotwI23",
			 "West Type 23 trigger side position ; X_{west} (mm) ; Y_{west} (mm)"
			 ,rbins,-rad,rad,rbins,-rad,rad);
  
  for(Int_t ibin = 0 ; ibin < rbins ; ibin++){
    for(Int_t jbin = 0 ; jbin < rbins ; jbin++){
      ecount[ibin][jbin] = 0.;
      wcount[ibin][jbin] = 0.;
    }
  }
  
  for(Int_t i = 0 ; i <nbeta ; i++){
    for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
	ecount[ibin-1][jbin-1] += btr[i]->hRoteI23->GetBinContent(ibin,jbin);
	wcount[ibin-1][jbin-1] += btr[i]->hRotwI23->GetBinContent(ibin,jbin);
      }
    }
  }
  
  for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
    for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
      hTotRoteI23->SetBinContent(ibin,jbin,ecount[ibin-1][jbin-1]);
      hTotRotwI23->SetBinContent(ibin,jbin,wcount[ibin-1][jbin-1]);
    }
  }
  
  TCanvas *crot = new TCanvas("crot","Rotation Cuts",600,600);
  TEllipse *el1 = new TEllipse(0,0,radcut,radcut);
  TEllipse *el1c = new TEllipse(0,0,50,50);
  TEllipse *el1cc = new TEllipse(0,0,45,45);
  el1->SetLineColor(2);
  el1->SetLineWidth(2);
  el1->SetFillStyle(0);
  el1c->SetLineColor(2);
  el1c->SetLineWidth(2);
  el1c->SetLineStyle(2);
  el1c->SetFillStyle(0);
  el1cc->SetLineColor(2);
  el1cc->SetLineWidth(2);
  el1cc->SetLineStyle(3);
  el1cc->SetFillStyle(0);
  crot->Divide(2,2);
  crot->cd(1);
  hTotRote->Draw("colz");
  hTotRote->GetXaxis()->CenterTitle();
  hTotRote->GetYaxis()->CenterTitle();
  el1->Draw();
  el1c->Draw();
  el1cc->Draw();
  crot->cd(2);
  hTotRoteI->Draw("colz");
  hTotRoteI->GetXaxis()->CenterTitle();
  hTotRoteI->GetYaxis()->CenterTitle();
  el1->Draw();
  el1c->Draw();
  el1cc->Draw();
  crot->cd(3);
  hTotRotw->Draw("colz");
  hTotRotw->GetXaxis()->CenterTitle();
  hTotRotw->GetYaxis()->CenterTitle();
  TEllipse *el2 = new TEllipse(0,0,radcut,radcut);
  el2->SetLineColor(2);
  el2->SetLineWidth(2);
  el2->SetFillStyle(0);
  el2->Draw();
  el1c->Draw();
  el1cc->Draw();
  crot->cd(4);
  hTotRotwI->Draw("colz");
  hTotRotwI->GetXaxis()->CenterTitle();
  hTotRotwI->GetYaxis()->CenterTitle();
  el2->Draw();
  el1c->Draw();
  el1cc->Draw();
  crot->Print(Form("pdf_out/geo%d/type_1_rotation_%d.pdf",btr[0]->GetGeo(),btr[0]->GetGeo()));
  TCanvas *crot23 = new TCanvas("crot23","crot23",600,600);
  TEllipse *el23 = new TEllipse(0,0,radcut,radcut);
  el23->SetLineColor(2);
  el23->SetLineWidth(2);
  el23->SetFillStyle(0);
  crot23->Divide(2,2);
  crot23->cd(1);
  hTotRote23->Draw("colz");
  el23->Draw();
  crot23->cd(2);
  hTotRoteI23->Draw("colz");
  el23->Draw();
  crot23->cd(3);
  hTotRotw23->Draw("colz");
  TEllipse *el233 = new TEllipse(0,0,radcut,radcut);
  el233->SetLineColor(2);
  el233->SetLineWidth(2);
  el233->SetFillStyle(0);
  el233->Draw();
  crot23->cd(4);
  hTotRotwI23->Draw("colz");
  el233->Draw();
  
}
//---------------------------------------------------------------------------
void CollectTDCCor()
{
  Int_t rbins = 200;
  Double_t ecount2[rbins][rbins];
  Double_t x1[10000],x1e[10000],srate[10000],srater[10000],brate[10000],brater[10000];
  Double_t sut[10000],sute[10000];
  
  
  // Set ecount to zero
  for(Int_t ibin = 0 ; ibin < rbins ; ibin++){
    for(Int_t jbin = 0 ; jbin < rbins ; jbin++){
      ecount2[ibin][jbin] = 0.;
    }
  }
  
  // Fill srate array with the rate of TDC header footer failures
  for(Int_t i = 0 ; i < nbeta ; i++){
    x1[i]    = btr[i]->GetRunNumber();
    x1e[i]   = 0.00001;
    
    srate[i] = btr[i]->hETDC_cor->Integral();
    srater[i]= sqrt(btr[i]->rtime_e*srate[i])/btr[i]->rtime_e;
    
    brate[i] = btr[i]->hETDC_cor_passed->Integral();
    brater[i]= sqrt(btr[i]->rtime_e*brate[i])/btr[i]->rtime_e;
    
    sut[i] = srate[i] + brate[i] + erttt[i];
    sute[i] = sqrt( srater[i]*srater[i] + brater[i]*brater[i] + ertet[i]*ertet[i]);
  }
 
  // Collect those failures into a Graph
  
  TGraphErrors *gcorr = new TGraphErrors(nbeta,x1,srate,x1e,srater);
  TGraphErrors *gcorb = new TGraphErrors(nbeta,x1,brate,x1e,brater);
  grstot = new TGraphErrors(nbeta,x1,sut,x1e,sute);
  
  ColorGraphic(gcorr,2,20,2);
  ColorGraphic(gcorb,2,20,2);
  ColorGraphic(grstot,4,20,2);
  ColorGraphic(grtot,3,20,2);
  
  cHBrate = new TCanvas("cHBrate","cHBrate");
  cHBrate->cd();
  // Draw the TDC corrupted events
  gcorr->Draw("AP");
  grtot->Draw("P");
  gcorr->GetXaxis()->SetTitle("Run Number");
  gcorr->SetTitle("Events with no Header/Footer");
  gcorr->GetYaxis()->SetTitle("Summed East + West Rate s^{-1}");
  gcorb->Draw("P");
  //  grstot->Draw("P");
  gcorr->GetYaxis()->SetRangeUser(0,20);
  
  TCanvas *cTDC = new TCanvas("cTDC","cTDC");
  hEvWtot  = new TH2F("hEvWtot","East vs West No Header/Footer;East(ns);West(ns)"
		      ,rbins,0,rbins,rbins,0,rbins);
  hEvWtotG = new TH2F("hEvWtotG","East vs West Good;East (ns);West(ns)"
		      ,rbins,0,rbins,rbins,0,rbins);
  // Add up the hEvWTDC histogram for all beta runs
  for(Int_t i = 0 ; i <nbeta ; i++){
    for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
	ecount2[ibin-1][jbin-1] += btr[i]->hEvWTDC_cor->GetBinContent(ibin,jbin);
      }
    }
  }
  // Put the results in the collection histogram and 0 the array
  for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
    for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
      hEvWtot->SetBinContent(ibin,jbin,ecount2[ibin-1][jbin-1]);
      ecount2[ibin-1][ibin-1] = 0;
    }
  }
  // Do this again for the good events
  for(Int_t i = 0 ; i <nbeta ; i++){
    for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
        ecount2[ibin-1][jbin-1] += btr[i]->hEvWTDC->GetBinContent(ibin,jbin);
      }
    }
  }
  
  for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
    for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
      hEvWtotG->SetBinContent(ibin,jbin,ecount2[ibin-1][jbin-1]);
    }
  }
  
  cTDC->Divide(2,1);
  cTDC->cd(1);
  hEvWtot->Draw("colz");
  cTDC->cd(2);
  hEvWtotG->Draw("colz");
}
//------------------------------------------------------------------------------
void CollectAsym()
{
  using namespace TMath;
  
  Double_t Asym_Tot[40],Asym_Tot_Er[40],X[40],Xe[40];
  Double_t W[1500],Temp[1500];
  Double_t X1[1500],Xe1[1500];

  for(Int_t j = 0 ; j < 40 ; j++){
    // Reset the statistical weight and the value of the asymmetry
    W[j]    = 0.;
    Temp[j] = 0.;
    // Loop through to calculate the error weighted average
    for(Int_t i = 0 ; i < nbeta ; i++){
        W[j] += 1./Power(btr[i]->fAsym_Run_Err[j],2);
	Temp[j] += W[j]*btr[i]->fAsym_Run[j];
    }
    Asym_Tot[j]    = Temp[j] / W[j]; 
    Asym_Tot_Er[j] = sqrt(1./W[j]);
    X[j]           = (Double_t)(j+1)*10.;
    Xe[j]          = 0.01; 
  }
  
  for(Int_t i = 0 ; i < nbeta; i++){
     X1[i]   = btr[i]->GetRunNumber();
     Xe1[i]  = 0.01;
     W[i]    = btr[i]->fAsym_Run_Err[20];
     Temp[i] = btr[i]->fAsym_Run[20];
  }
  
  TGraphErrors *gAsym = new TGraphErrors(40,X,Asym_Tot,Xe,Asym_Tot_Er);
  TGraphErrors *gChck = new TGraphErrors(nbeta,X1,Temp,Xe1,W);
  ColorGraphic(gAsym,2,20,2);
  ColorGraphic(gChck,4,20,4);
  
  TCanvas *cAsyms = new TCanvas("cAsyms","Asymmetry Calculations");
  cAsyms->Divide(1,2);
  
  cAsyms->cd(1);
  gAsym->Draw("AP");
  
  cAsyms->cd(2);
  gChck->Draw("Ap");
}

void Plot_E_Chis()
{
  
  // Plot a random runs background subtracted Beta energy spectrum against the 
  // appropiate PENELOPE prediction.
  
  TCanvas *cEchis = new TCanvas("cEchis","Energy and Backgrounds");
  TRandom3 *xrand = new TRandom3();
  xrand->SetSeed(0);
  
  Double_t x1[10000],y1[10000],x1er[10000],y1er[10000];
  Double_t x1w[10000],y1w[10000],x1erw[10000],y1erw[10000];
  
  Int_t octn = xrand->Integer(nbeta);
  octn = 1;
  cEchis->Divide(2,1);
  cEchis->cd(1);
  
  cout << "Print Histogram for beta run " << btr[octn]->GetRunNumber() << endl;
/*
  btr[octn]->heq->Draw("E");
  bckr[btr[octn]->Bkg_index]->heq->SetLineColor(4);
  bckr[btr[octn]->Bkg_index]->heq->Draw("same");
  btr[octn]->hEtype_1->SetLineColor(3);
  btr[octn]->hEtype_1->Draw("same");
  btr[octn]->hEtype_23->SetLineColor(6);
  btr[octn]->hEtype_23->Draw("same");
  btr[octn]->hERef->SetLineColor(2);
  btr[octn]->hERef->Draw("same");
  */
  Int_t va = 0;
  for(Int_t i = 0 ; i < nbeta ; i++){
    if(btr[i]->flipperOn == 0){
      x1[va]   = btr[i]->GetRunNumber();
      y1[va]   = btr[i]->Residual_Bkg[0]*4.5/btr[i]->heq->Integral(20,65);
      y1er[va] = btr[i]->Resid_Bkger[0]*4.5/btr[i]->heq->Integral(20,65);
      va++;
    }
  }
  
  TGraphErrors *gEChi = new TGraphErrors(va,x1,y1,0,y1er);
  ColorGraphic(gEChi,2,20,2);

  cEchis->cd(1);
  gPad->SetGrid();
  gEChi->SetTitle("East");
  gEChi->GetXaxis()->SetTitle("Run Number");
  gEChi->GetYaxis()->SetTitle("Extrapolated Background under #beta-Spectra");
  gEChi->Draw("AP");
  
  cEchis->cd(2);
  /*
  btr[octn]->hwq->Draw("E");
  bckr[btr[octn]->Bkg_index]->hwq->SetLineColor(4);
  bckr[btr[octn]->Bkg_index]->hwq->Draw("same");
  btr[octn]->hWtype_1->SetLineColor(3);
  btr[octn]->hWtype_1->Draw("same");
  btr[octn]->hWtype_23->SetLineColor(6);
  btr[octn]->hWtype_23->Draw("same");
  btr[octn]->hERef->SetLineColor(2);
  btr[octn]->hERef->Draw("same");
  */
  // Extrapolate the Residual background beneath the beta spectrum..
  
  Int_t vb = 0;
  for(Int_t i = 0 ; i < nbeta ; i++){
    if(btr[i]->flipperOn == 0){
      x1w[vb]   = btr[i]->GetRunNumber();
      y1w[vb]   = btr[i]->Residual_Bkg[1]*4.5/btr[i]->hwq->Integral(20,65);
      y1erw[vb] = btr[i]->Resid_Bkger[1]*4.5/btr[i]->hwq->Integral(20,65);
      vb++;
    }
  }
  
  TGraphErrors *gEChiw = new TGraphErrors(vb,x1w,y1w,0,y1erw);
  ColorGraphic(gEChiw,2,20,2);
  
  cEchis->cd(4);
  gPad->SetGrid();
  gEChiw->SetTitle("West");
  gEChiw->GetXaxis()->SetTitle("Run Number");
  gEChiw->GetYaxis()->SetTitle("Extrapolated Background under #beta-Spectra");
  gEChiw->Draw("AP");
  

  cEchis->Print("output_files/energy_chis.pdf");

  delete gEChiw;
  delete xrand;
  delete cEchis;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
void CalcSimplSuper()
{
  
  // Quickly Calculate the super ratios.
  lookahead = -4;
  Nsuper    = 1;

  // Find Regular SuperRatios
  fsuprat1.open("super_ratio.txt",fstream::out);
  
  for(Int_t i = 0 ; i < nbeta ; i++){
    Get_Base_Super(i,(char*)"a2",(char*)"a5");
    Get_Base_Super(i,(char*)"a7",(char*)"a10");
    Get_Base_Super(i,(char*)"b2",(char*)"b5");
    Get_Base_Super(i,(char*)"b7",(char*)"b10");
  }
  
  fsuprat1.close();
  
  TCanvas *cSuper = new TCanvas("cSuper","Super Ratios");
  cSuper->Divide(1,2);
  cSuper->cd(1);
  
  TGraphErrors *gSup = new TGraphErrors(Nsuper-1,Xrun,Super,0,SuperE);
  TF1 *fSup = new TF1("fSup","[0]",0.,Nsuper+1.);
  fSup->SetParameter(0,0.044);
  fSup->SetLineWidth(2);
  fSup->SetLineColor(2);
  
  ColorGraphic(gSup,4,20,2);
  
  gSup->Draw("AP");
  // gSup->GetYaxis()->SetRangeUser(0,0.1);
  gSup->Fit("fSup","RMEQ","");
  gPad->SetGrid();
  gSup->SetTitle("Raw Super Ratio 0-600 keV (analysis choice A)");
  gSup->GetXaxis()->SetTitle("Quartet Pair");
  gSup->GetYaxis()->SetTitle("(1 - #sqrt{S}) / (1 + #sqrt{S})");
  gSup->GetXaxis()->CenterTitle();
  gSup->GetYaxis()->CenterTitle();
 
  TGraphErrors *gSupC = new TGraphErrors(Nsuper-1,Xrun,SuperC,Xer,SuperCE);
  TF1 *fSupC = new TF1("fSupC","[0]",0.,Nsuper+1.);
  fSupC->SetParameter(0,0.044);
  fSupC->SetLineWidth(2);
  fSupC->SetLineColor(1);
  
  ColorGraphic(gSupC,2,20,2);
  gSupC->Draw("P");
  gSupC->Fit("fSupC","REMQ","");
  
  
  TLegend *cl = new TLegend(0.6,0.6,0.9,0.9);
  cl->AddEntry(gSup,Form("A = %6.4f #pm %6.4f",
	         fSup->GetParameter(0),fSup->GetParError(0)),"lp");
  cl->AddEntry(gSup,Form("#chi^{2} / #nu = %6.4f",fSup->GetChisquare() 
			 / fSup->GetNDF()),"");
  cl->AddEntry(gSupC,Form("A_{c} = %6.4f #pm %6.4f",
	       fSupC->GetParameter(0),fSupC->GetParError(0)),"lp");
  cl->AddEntry(gSupC,Form("#chi^{2} / #nu = %6.4f",fSupC->GetChisquare() 
			  / fSupC->GetNDF()),"");
  cl->SetFillStyle(0);
  cl->Draw();
  
  cSuper->cd(2);
  
  Double_t Redis[100000],Rediser[100000];
  
  for(int i = 0 ; i < Nsuper-1 ; i++){
    Redis[i] = (Super[i] - fSup->GetParameter(0))/SuperE[i];
    Rediser[i] = 1.;
  }
  
  TGraphErrors *gRes = new TGraphErrors(Nsuper-1,Xrun,Redis,Xer,Rediser);
  ColorGraphic(gRes,4,20,2);
  
  gPad->SetGrid();
  gRes->Draw("AP");
  gRes->SetTitle("Residuals");
  gRes->GetXaxis()->SetTitle("Quartet Pair");
  gRes->GetYaxis()->SetTitle("Residuals in Units of #sigma");
  cSuper->Print("output_files/superratio_vs_q.pdf");

  //----------------------------------------------------------------------------------------
  TCanvas *cSup2d = new TCanvas("cSup2d","2d Super Ratios");
  cSup2d->Divide(2,2);
  cSup2d->cd(1);
  
  TGraphErrors *gSup1 = new TGraphErrors(Nsuper-1,Xrun,SuperPos1,0,SuperEPos1);
  TGraphErrors *gSup2 = new TGraphErrors(Nsuper-1,Xrun,SuperPos2,0,SuperEPos2);
  TGraphErrors *gSup3 = new TGraphErrors(Nsuper-1,Xrun,SuperPos3,0,SuperEPos3);
  TGraphErrors *gSup4 = new TGraphErrors(Nsuper-1,Xrun,SuperPos4,0,SuperEPos4);
  
  ColorGraphic(gSup1,4,20,2);
  ColorGraphic(gSup2,4,20,2);
  ColorGraphic(gSup3,4,20,2);
  ColorGraphic(gSup4,4,20,2);
  
  
  TF1 *fLSup1 = new TF1("fLSup1","[0]",0,Nsuper+1);
  TF1 *fLSup2 = new TF1("fLSup2","[0]",0,Nsuper+1);
  TF1 *fLSup3 = new TF1("fLSup3","[0]",0,Nsuper+1);
  TF1 *fLSup4 = new TF1("fLSup4","[0]",0,Nsuper+1);
  
  gSup1->Draw("AP");
  gSup1->Fit(fLSup1,"REMQ","");
  cSup2d->cd(2);
  gSup2->Draw("AP");
  gSup2->Fit(fLSup2,"REMQ","");
  cSup2d->cd(3);
  gSup3->Draw("AP");
  gSup3->Fit(fLSup3,"REMQ","");
  cSup2d->cd(4);
  gSup4->Draw("AP");
  gSup4->Fit(fLSup4,"REMQ","");

  cSup2d->Print("output_files/2d_asymmetries.pdf");

  //-----------------------------------------------------------------
  TCanvas *cPave = new TCanvas("cPave","Position Averages");
  cPave->cd();
  Double_t SupPos[4],SupPosE[4],Xpos[4],Xposer[4];
  
  SupPos[0] = fLSup1->GetParameter(0);
  SupPos[1] = fLSup2->GetParameter(0);
  SupPos[2] = fLSup3->GetParameter(0);
  SupPos[3] = fLSup4->GetParameter(0);
  
  SupPosE[0] = fLSup1->GetParError(0);
  SupPosE[1] = fLSup2->GetParError(0);
  SupPosE[2] = fLSup3->GetParError(0);
  SupPosE[3] = fLSup4->GetParError(0);
  
  for(int i = 0; i < 4; i++){
    Xpos[i] = i+1;
    Xposer[i] = 0.01;
  }
  
  
  TGraphErrors *gPave = new TGraphErrors(4,Xpos,SupPos,Xposer,SupPosE);
  ColorGraphic(gPave,4,20,2);

  gPave->Draw("AP");
  
 return;
}

void Get_Base_Super(Int_t i, char oct1[2],char oct2[2])
{
  
  // Eliminates the potential double counting of Quartets
  
  Int_t jstar = 0;
  
  if(i>0)
    
     if(!strcmp(btr[i-1]->octtype,btr[i]->octtype) && btr[i-1]->rtime_e > btr[i]->rtime_e) return;
  if(i < nbeta-1)
    if(!strcmp(btr[i+1]->octtype,btr[i]->octtype) && btr[i+1]->rtime_e > btr[i]->rtime_e) return;
  // Check for the existance of the quartet pair
  if(!strcmp(btr[i]->octtype,oct1)){
    // prevent from looking a negative array indices
    if(i < 4) lookahead = 0;
    // Loop around to find the correct run
    do{
      
      lookahead++;
      
    }while(strcmp(btr[i+lookahead]->octtype,oct2) && lookahead < 6);
    // if 5 then its a bad data set
    
    if(lookahead == 6){
      cout << "Matching Run for Run " << btr[i]->GetRunNumber() << " Not Found!! " << Nsuper << endl;
      if(btr[i]->GetRunNumber() == 11090){
         do{
	   jstar++;
	 }while(btr[jstar-1]->GetRunNumber() != 11093);
	 lookahead = jstar - i -1;
      }
      if(btr[i]->GetRunNumber() == 11084){
	do{
	   jstar++;
	 }while(btr[jstar-1]->GetRunNumber() != 11088);
	 lookahead = jstar-1 - i;
      }
    }// else {
      fsuprat1 <<"Super ratios " <<  btr[i]->GetRunNumber() << "\t" << btr[i+lookahead]->GetRunNumber() << endl;
      // if the quartet is found calculate the Asymmetry from the
      // super ratio.
      // This needs to be adjusted to add in the cases where there are 
      // multiple files for the same type of octet run.
      
      // Get Super Ratio
      Super[Nsuper-1]  = GetSuperRatio(btr[i]->hwq,btr[i]->heq,
				     btr[i+lookahead]->hwq,btr[i+lookahead]->heq);
      // Get Super Ratio error
      SuperE[Nsuper-1] = GetSuperRatioError(btr[i]->hwq,btr[i]->heq,
					  btr[i+lookahead]->hwq,btr[i+lookahead]->heq,
                                          btr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      // Get Analysis C Super Ratios
      SuperC[Nsuper-1]  = GetSuperRatio(btr[i]->hwqC,btr[i]->heqC,
				       btr[i+lookahead]->hwqC,btr[i+lookahead]->heqC);
      // Get Super Ratio error
      SuperCE[Nsuper-1] = GetSuperRatioError(btr[i]->hwqC,btr[i]->heqC,
					    btr[i+lookahead]->hwqC,btr[i+lookahead]->heqC,
	                                    btr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      // Get Position Dependent Super Ratios
      
      SuperPos1[Nsuper-1] = GetSuperRatio(btr[i]->hAposw->GetBinContent(1,1),
					  btr[i]->hApose->GetBinContent(1,1),
					  btr[i+lookahead]->hAposw->GetBinContent(1,1),
					  btr[i+lookahead]->hApose->GetBinContent(1,1));
      
      SuperEPos1[Nsuper-1] = GetSuperRatioError(btr[i]->hAposw->GetBinContent(1,1),
					    btr[i]->hApose->GetBinContent(1,1),
					    btr[i+lookahead]->hAposw->GetBinContent(1,1),
					    btr[i+lookahead]->hApose->GetBinContent(1,1),
	                                    btr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      //----------------------------------------------------------------------------------------
      SuperPos2[Nsuper-1] = GetSuperRatio(btr[i]->hAposw->GetBinContent(2,1),
					  btr[i]->hApose->GetBinContent(2,1),
					  btr[i+lookahead]->hAposw->GetBinContent(2,1),
					  btr[i+lookahead]->hApose->GetBinContent(2,1));
      
      SuperEPos2[Nsuper-1] = GetSuperRatioError(btr[i]->hAposw->GetBinContent(2,1),
	                                        btr[i]->hApose->GetBinContent(2,1),
					        btr[i+lookahead]->hAposw->GetBinContent(2,1),
					        btr[i+lookahead]->hApose->GetBinContent(2,1),
						btr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      //----------------------------------------------------------------------------------------
      
      SuperPos3[Nsuper-1] = GetSuperRatio(btr[i]->hAposw->GetBinContent(1,2),
					  btr[i]->hApose->GetBinContent(1,2),
					  btr[i+lookahead]->hAposw->GetBinContent(1,2),
				          btr[i+lookahead]->hApose->GetBinContent(1,2));
      
      SuperEPos3[Nsuper-1] = GetSuperRatioError(btr[i]->hAposw->GetBinContent(1,2),
	                                        btr[i]->hApose->GetBinContent(1,2),
					        btr[i+lookahead]->hAposw->GetBinContent(1,2),
					        btr[i+lookahead]->hApose->GetBinContent(1,2),
						btr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      //--------------------------------------------------------------------------------------
      SuperPos4[Nsuper-1] = GetSuperRatio(btr[i]->hAposw->GetBinContent(2,2),
					  btr[i]->hApose->GetBinContent(2,2),
					  btr[i+lookahead]->hAposw->GetBinContent(2,2),
				          btr[i+lookahead]->hApose->GetBinContent(2,2));
      
      SuperEPos4[Nsuper-1] = GetSuperRatioError(btr[i]->hAposw->GetBinContent(2,2),
	                                        btr[i]->hApose->GetBinContent(2,2),
					        btr[i+lookahead]->hAposw->GetBinContent(2,2),
					        btr[i+lookahead]->hApose->GetBinContent(2,2),
						btr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      
      // Incriment countsers.
      Xrun[Nsuper-1] = Nsuper;
      Xer[Nsuper-1]  = 0.02;
      Nsuper++;
    //}
    // reset lookahead variable to -4 so that we can look ahead and behind the current run
    // just incase.
    lookahead = -4;
    jstar = 0 ;
  }
  
  return;
}

Double_t GetSuperRatio(TH1F *R1M, TH1F *R1P,TH1F *R2M, TH1F *R2P)
{
  
  // Get the Super Ratio straight forward..
  
  using namespace TMath;

//  Int_t nlow = 20;
//  Int_t nhigh = 65;

  Double_t R = R1M->Integral(nlow,nhigh)*R2P->Integral(nlow,nhigh) / 
    (R1P->Integral(nlow,nhigh)*R2M->Integral(nlow,nhigh));
  Double_t S = Abs((1. - sqrt(R)) / ( 1. + sqrt(R)));
  
  return S;
}

Double_t GetSuperRatio(Double_t R1M, Double_t R1P,Double_t R2M, Double_t R2P)
{
  
  // Get the Super Ratio straight forward.. with Double values instead of from 
  // histograms
  
  using namespace TMath;
  
  Double_t R = R1M*R2P / (R1P*R2M);
  Double_t S = Abs((1. - sqrt(R)) / ( 1. + sqrt(R))); 
  
  return S;
}

Double_t GetSuperRatioError(TH1F *R1M, TH1F *R1P,TH1F *R2M, TH1F *R2P,Double_t time1,Double_t time2)
{
  
  // ----------------------------------------------------------------------------+
  // Calculate the Error in the Super Ratio Asymmetry which will be inflated     |
  // once the background numbers are added in.  For now I'm just trying to get   |
  // a number.                                                                   |
  //-----------------------------------------------------------------------------+
  using namespace TMath;
  
  //Int_t nlow = 20;
  //Int_t nhigh = 60;
  
  Double_t R = R1M->Integral(nlow,nhigh)*R2P->Integral(nlow,nhigh) / 
    (R1P->Integral(nlow,nhigh)*R2M->Integral(nlow,nhigh));
  
  Double_t K1 = sqrt(R) / Power(1.+sqrt(R),2);
  
  Double_t T1 = sqrt(R1M->Integral(nlow,nhigh)/time1) / R1M->Integral(nlow,nhigh);
  Double_t T2 = sqrt(R1P->Integral(nlow,nhigh)/time1) / R1P->Integral(nlow,nhigh);
  
  Double_t T3 = sqrt(R2M->Integral(nlow,nhigh)/time2) / R2M->Integral(nlow,nhigh);
  Double_t T4 = sqrt(R2P->Integral(nlow,nhigh)/time2) / R2P->Integral(nlow,nhigh);
  
  Double_t Error = K1*sqrt(T1*T1+T2*T2+T3*T3+T4*T4);
  
  return Error;
}

Double_t GetSuperRatioError(Double_t R1M,Double_t  R1P,Double_t R2M,Double_t R2P,Double_t time1,Double_t time2)
{
  
  // ----------------------------------------------------------------------------+
  // Calculate the Error in the Super Ratio Asymmetry which will be inflated     |
  // once the background numbers are added in.  For now I'm just trying to get   |
  // a number.                                                                   |
  //-----------------------------------------------------------------------------+
  using namespace TMath;
  
  Double_t R = R1M *R2P  / (R1P *R2M );
  
  Double_t K1 = sqrt(R) / Power(1.+sqrt(R),2);
  
  Double_t T1 = sqrt(R1M /time1) / R1M ;
  Double_t T2 = sqrt(R1P /time1) / R1P ;
  
  Double_t T3 = sqrt(R2M /time2) / R2M ;
  Double_t T4 = sqrt(R2P /time2) / R2P ;
  
  Double_t Error = K1*sqrt(T1*T1+T2*T2+T3*T3+T4*T4);
  
  return Error;
}

void ColorGraphic(TH1 *h, Int_t color, Int_t Style, Int_t Width,Int_t nmarker,Double_t msize)
{
  h->SetMarkerColor(color);
  h->SetLineColor(color);
  h->SetLineWidth(Width);
  h->SetLineStyle(Style);
  h->SetMarkerStyle(nmarker);
  h->SetMarkerSize(msize);
}

void ColorGraphic(TGraphErrors *h, Int_t color, Int_t Style, Int_t Width)
{
  h->SetLineColor(color);
  //h->SetLineWidth(Width);
  h->SetMarkerStyle(Style);
  h->SetMarkerColor(color);
}

void ColorGraphic(TGraph *h, Int_t color, Int_t Style, Int_t Width)
{
  h->SetLineColor(color);
  //h->SetLineWidth(Width);
  h->SetMarkerStyle(Style);
  h->SetMarkerColor(color);
}

void Collect_Pos()
{
  Int_t rbins = 150;
  Float_t rad = 75.;
  Double_t ecount[rbins][rbins],wcount[rbins][rbins];
  
  hTotEPos = new TH2F("hTotEPos",
		      "East Type 1 ; X_{west} (mm) ; Y_{west} (mm)",
		      rbins,-rad,rad,rbins,-rad,rad);
  
  hTotWPos = new TH2F("hTotWPos",
		      "West Type  ; X_{east} (mm) ; Y_{east} (mm)",
		      rbins,-rad,rad,rbins,-rad,rad);
  
  for(Int_t ibin = 0 ; ibin < rbins ; ibin++){
    for(Int_t jbin = 0 ; jbin < rbins ; jbin++){
      ecount[ibin][jbin] = 0.;
      wcount[ibin][jbin] = 0.;
    }
  }
  cout << "here " << endl;
  cout << "Integral " << btr[0]->hpe->Integral();

  for(Int_t i = 0 ; i <nbeta ; i++){
    for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
	ecount[ibin-1][jbin-1] += btr[i]->hpe->GetBinContent(ibin,jbin);
	wcount[ibin-1][jbin-1] += btr[i]->hpw->GetBinContent(ibin,jbin);
      }
    }
  }
  
  for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
    for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
      hTotEPos->SetBinContent(ibin,jbin,ecount[ibin-1][jbin-1]);
      hTotWPos->SetBinContent(ibin,jbin,wcount[ibin-1][jbin-1]);
    }
  }
  
  TCanvas *cPosAll = new TCanvas("cPosAll","Position All",600,600);
  TEllipse *elp1 = new TEllipse(0,0,58.5,58.5);
  TEllipse *elp1c = new TEllipse(0,0,50.0,50.0);
  TEllipse *elp1x = new TEllipse(0,0,45.0,45.0);
  
  elp1->SetLineColor(2);
  elp1->SetLineWidth(2);
  elp1->SetFillStyle(0);
  elp1c->SetLineColor(2);
  elp1c->SetLineWidth(2);
  elp1c->SetLineStyle(2);
  elp1c->SetFillStyle(0);
  elp1x->SetLineColor(2);
  elp1x->SetLineWidth(2);
  elp1x->SetLineStyle(3);
  elp1x->SetFillStyle(0);
  cPosAll->Divide(2,2);
  cPosAll->cd(1);  
  hTotEPos->Draw("colz");
  hTotEPos->GetXaxis()->CenterTitle();
  hTotEPos->GetYaxis()->CenterTitle();
  elp1->Draw();
  elp1c->Draw();
  elp1x->Draw();
  
  TLegend *lPosPlot = new TLegend(0.12,0.7,0.5,0.88);
  lPosPlot->SetFillColor(0);
  lPosPlot->AddEntry(elp1 ,"r=58.5mm ,Collimator","l");
  lPosPlot->AddEntry(elp1c,"r=50.0mm ,R_{MC} Max","l");
  lPosPlot->AddEntry(elp1x,"r=45.0mm ,Radial Cut","l");
  // lPosPlot->Draw();
  
  cPosAll->cd(2);
  hTotWPos->Draw("colz");
  elp1->Draw();
  elp1c->Draw();
  elp1x->Draw();


  Double_t ey = 0;
  Double_t wy = 0;
  Double_t earea = 0;
  Double_t warea = 0;
  
  TH1F *hEyC = new TH1F("hEyC","East height difference",rbins,-rad,rad);
  TH1F *hWyC = new TH1F("hWyC","West height difference",rbins,-rad,rad);
  
  TH1F *hEyCp = new TH1F("hEyCp","East height difference",75,0,rad);
  TH1F *hWyCp = new TH1F("hWyCp","West height difference",75,0,rad);
  TH1F *hEyCm = new TH1F("hEyCm","East height difference",75,0,rad);
  TH1F *hWyCm = new TH1F("hWyCm","West height difference",75,0,rad);
  
  for(Int_t i = 1 ; i <= rbins ; i++){
    ey =0;
    wy =0;
    earea = 0;
    warea = 0;
    for(Int_t j = 1 ; j <= rbins ;j++){
      ey += hTotEPos->GetBinContent(i,j);
      wy += hTotWPos->GetBinContent(i,j);
      if(ey > 0)earea++;
      if(wy > 0)warea++;
    }
    if(ey > 0){
      ey = ey/earea;
      hEyC->SetBinContent(i,ey);
      hEyC->SetBinError(i,sqrt(ey*earea)/earea);
      if(i <= 75){
	hEyCm->SetBinContent(76-i,ey);
	hEyCm->SetBinError(76-i,sqrt(ey*earea)/earea);
      } else {
        hEyCp->SetBinContent(i-75,ey);
	hEyCp->SetBinError(i-75,sqrt(ey*earea)/earea);
      }
    }
    if(wy > 0){
      wy = wy/warea;
      hWyC->SetBinContent(i,wy);
      hWyC->SetBinError(i,sqrt(wy*warea)/warea);
      if(i <= 75){
	hWyCm->SetBinContent(76-i,wy);
	hWyCm->SetBinError(76-i,sqrt(wy*warea)/warea);
      } else {
        hWyCp->SetBinContent(i-75,wy);
	hWyCp->SetBinError(i-75,sqrt(wy*warea)/warea);
      }
    }
  }
  
  cPosAll->cd(3);
  hEyC->Draw();
  hWyC->Draw("same");
  cPosAll->cd(4);
  hWyCm->Add(hWyCp,-1);
  hWyCm->Draw();
  hEyCm->Add(hEyCp,-1);
  hEyCm->SetLineColor(2);
  hEyCm->Draw("same");
  
  cPosAll->Print(Form("pdf_out/geo%d/position_all_%d.pdf",
		      btr[0]->GetGeo(),btr[0]->GetGeo()));
}
//--------------------------------------------------------------------
void Plot_23_Diff()
{
  TCanvas *cboring = new TCanvas("cboring","");
  cboring->Divide(2,1);
  cboring->cd(1);
  btr[1]->hPsDfET23->Draw("colz");
  
  Double_t nxr[10000],nt2f[10000];
  
  for(int i = 0 ; i < nbeta ; i++){
    nxr[i]  = btr[i]->heq->Integral()+btr[i]->hwq->Integral();
    nt2f[i] = btr[i]->T2_frac;
  }
  cboring->cd(2);
  TGraph *gT2 = new TGraph(nbeta,nxr,nt2f);
  gT2->SetMarkerColor(2);
  gT2->SetMarkerStyle(20);
  gT2->Draw("AP");
}
//-------------------------------------------------------------------------  
void Plot_MonRate()
{
  using namespace std;

  vector<Double_t> xt,y;
  Double_t ek[10000],wk[10000];
  Double_t eker[10000],wker[10000];
  Double_t emean[10000],wmean[10000];
  
  TH1F *hEastEndPoint = new TH1F("hEastEndPoint","East End Point ; E_{vis} (keV) ; Counts",
				 100,200,1200);
  TH1F *hWestEndPoint = new TH1F("hWestEndPoint","West End Point ; E_{vis} (keV) ; Counts",
				 100,200,1200);				  
  
  for(int i = 0 ; i < nbeta ; i++){
    xt.push_back(btr[i]->GetRunNumber());
    if(btr[i]->Mon_Rate >0)
      y.push_back((btr[i]->heq->Integral() + btr[i]->hwq->Integral())
		  /btr[i]->Mon_Rate);
    ek[i]   = btr[i]->E_Endpoint;
    wk[i]   = btr[i]->W_Endpoint;
    eker[i] = btr[i]->E_EndError;
    wker[i] = btr[i]->W_EndError;
    hEastEndPoint->Fill(btr[i]->E_Endpoint);
    hWestEndPoint->Fill(btr[i]->W_Endpoint);
    emean[i] = btr[i]->heq->GetMean();
    wmean[i] = btr[i]->hwq->GetMean();
  }
  
  TGraph *gBdM = new TGraph((int)y.size(),&xt[0],&y[0]);
  ColorGraphic(gBdM,2,20,1);
 
  TCanvas *cMon = new TCanvas("cMon","Beta-Rate Divided by Monitor 2");
  cMon->cd(2);
  gBdM->SetTitle("Normalized #beta Rate");
  gBdM->GetXaxis()->SetTitle("Run Number");
  gBdM->GetYaxis()->SetTitle("#beta -Rate / Monitor 2 ");
  gBdM->Draw("AP");
  
  TCanvas *cHom = new TCanvas("cHom","Stuff");
  cHom->cd();
  TGraphErrors *gKurieE = new TGraphErrors(nbeta,&xt[0],ek,0,eker);
  TGraphErrors *gKurieW = new TGraphErrors(nbeta,&xt[0],wk,0,wker);

  TGraph *gMeanE = new TGraph(nbeta,&xt[0],emean);
  TGraph *gMeanW = new TGraph(nbeta,&xt[0],wmean);
  
  ColorGraphic(gKurieE,2,20,1);
  ColorGraphic(gKurieW,4,20,1);
  ColorGraphic(gMeanE,2,20,1);
  ColorGraphic(gMeanW,4,20,1);
  
  cHom->Divide(2,2);
  cHom->cd(1);
  gKurieE->Draw("AP");
  gKurieE->GetXaxis()->SetRangeUser(0,1200);
  gKurieW->Draw("P");
  
  cHom->cd(2);
  
  ColorGraphic(hEastEndPoint,2,20,2,1);
  ColorGraphic(hWestEndPoint,4,20,2,1);
  hEastEndPoint->Draw();
  hWestEndPoint->Draw("same");
 
  cHom->cd(3);

  gMeanE->Draw("AP");
  gMeanW->Draw("P");
  
  /*
  for(int i = 0 ; i < nbeta ; i++){
    if(btr[i]->E_Endpoint  < 100){
      if(cter == 0)
	btr[i]->heq->Draw();
      else
	btr[i]->heq->Draw("same");
      cter++;
    }
    if(btr[i]->W_Endpoint  < 100){
      if(cter == 0)
	btr[i]->hwq->Draw();
      else
	btr[i]->hwq->Draw("same");
      cter++;
    }
  }*/
  
}

void Collect_Octets()
{
   
  using namespace TMath;
  
  TF1 *flA = new TF1("flA","[0]",0,noct);
  TF1 *flB = new TF1("flB","[0]",0,noct);

  flA->SetLineColor(2);
  flB->SetLineColor(4);

  Double_t x[100],A_sum_A[100],A_sum_B[100],A_multi_A[100],A_multi_B[100],xq[100],xqer[100];
  Double_t A_quat[100],A_quat_er[100];
  Double_t A_aveB[100],A_aveB_er[100];
  Double_t A_aveC[100],A_aveC_er[100];
  Double_t A_aveD[100],A_aveD_er[100];
  Double_t A_aveE[100],A_aveE_er[100];
  Double_t A_aveF[100],A_aveF_er[100];
  Double_t A_aveG[100],A_aveG_er[100];
  Double_t A_aveH[100],A_aveH_er[100];
  Double_t A_aveI[100],A_aveI_er[100]; 
  Double_t A_aves[100],A_aves_er[100];
  Double_t A_avesB[100],A_avesB_er[100];
  Double_t A_avesC[100],A_avesC_er[100];
  Double_t A_avesD[100],A_avesD_er[100];
  Double_t A_avesE[100],A_avesE_er[100];
  Double_t A_avesF[100],A_avesF_er[100];
  Double_t A_avesG[100],A_avesG_er[100];
  Double_t A_avesH[100],A_avesH_er[100];
  Double_t A_avesI[100],A_avesI_er[100]; 
  Double_t xer[100],A_sum_Aer[100],A_sum_Ber[100],A_multi_Aer[100],A_multi_Ber[100];
  Double_t A_ave[100],A_ave_er[100];
  Double_t A_super[400],A_superer[400],xsuper[400],xsuperer[400];
  Double_t A_ave_tot[9],A_ave_toter[9],A_aves_tot[9],A_aves_toter[9];

  Int_t NoctC   = 0;
  Int_t FullOct = 0;
  Int_t NQuat   = 0;
  
  fstream foct,fsup,fquart;
  foct.open(Form("oct_list_%d.txt",btr[0]->GetGeo()),fstream::out);
  fsup.open(Form("super_list_%d.txt",btr[0]->GetGeo()),fstream::out);
  fquart.open(Form("quart_list_%d.txt",btr[0]->GetGeo()),fstream::out);
  foct <<   "Noct\t Runs \t A " << endl;
  fquart << "Noct\t Runs \t A quartet A \t\t A quartet B" << endl;
  fsup <<   "Noct\t Runs \t A sup A1 \t\t A super A2 \t\t Asuper B1 \t\t Asuper B2"<<endl;
  
  
  for(Int_t i = 0; i < noct ; i++){

    foct << i << "\t" << octet[i]->GetFirst() << "-" << octet[i]->GetLast() << "\t" << octet[i]->A_sum[2] << "\t" << octet[i]->A_sumer[2] <<  endl;
    fquart << i<< "\t" << octet[i]->GetFirst() << "-" << octet[i]->GetLast() << "\t" << octet[i]->A_sum_A[2];
    fquart << "\t" <<octet[i]->A_sum_Aer[2] << "\t" << octet[i]->A_sum_B[2] << "\t"<< octet[i]->A_sum_Ber[2] << endl;
    fsup << i<< "\t" << octet[i]->GetFirst() << "-" << octet[i]->GetLast() << "\t";
    fsup << octet[i]->Asuper1 << "\t" << octet[i]->Asuper1er << "\t" << octet[i]->Asuper2 << "\t" << octet[i]->Asuper2er << "\t";
    fsup << octet[i]->Bsuper1 << "\t" << octet[i]->Bsuper1er << "\t" <<octet[i]->Bsuper2 << "\t" <<octet[i]->Bsuper2er << endl;
    
    // Output the oct list used in the analysis.
    if(IsNaN(octet[i]->Asuper2) || IsNaN(octet[i]->Asuper1) || IsNaN(octet[i]->Bsuper1) || IsNaN(octet[i]->Bsuper2) || 
     octet[i]->Asuper2 == 0 ||octet[i]->Asuper1 == 0 || octet[i]->Bsuper2 == 0 || octet[i]->Bsuper1 == 0   ){
      if(!(IsNaN(octet[i]->Asuper1)) && octet[i]->Asuper1 != 0  ){
	  A_super[NoctC]   = octet[i]->Asuper1;
	  A_superer[NoctC] = octet[i]->Asuper1er;
	  xsuper[NoctC]    = NoctC + 1;
	  xsuperer[NoctC]  = 0.001;
	  NoctC++;
      }
      if(!(IsNaN(octet[i]->Asuper2)) && octet[i]->Asuper2 != 0  ){
	  A_super[NoctC]   = octet[i]->Asuper2;
	  A_superer[NoctC] = octet[i]->Asuper2er;
	  xsuper[NoctC]    = NoctC + 1;
	  xsuperer[NoctC]  = 0.001;
	  NoctC++;
      }
      
      if(!(TMath::IsNaN(octet[i]->Bsuper1)) && octet[i]->Bsuper1 != 0 ){
	A_super[NoctC]   = octet[i]->Bsuper1;
	A_superer[NoctC] = octet[i]->Bsuper1er;
	xsuper[NoctC]    = NoctC + 1;
	xsuperer[NoctC]  = 0.001;
	NoctC++;
      }
      
      if(!(TMath::IsNaN(octet[i]->Bsuper2)) && octet[i]->Bsuper2 != 0 ){
	A_super[NoctC]   = octet[i]->Bsuper2;
	A_superer[NoctC] = octet[i]->Bsuper2er;
	xsuper[NoctC]    = NoctC + 1;
	xsuperer[NoctC]  = 0.001;
	NoctC++;
      }
    } else {
      A_super[NoctC]   = octet[i]->Asuper1;
      A_super[NoctC+1] = octet[i]->Asuper2;
      A_super[NoctC+2] = octet[i]->Bsuper1;
      A_super[NoctC+3] = octet[i]->Bsuper2;
      A_superer[NoctC]   = octet[i]->Asuper1er;
      A_superer[NoctC+1] = octet[i]->Asuper2er;
      A_superer[NoctC+2] = octet[i]->Bsuper1er;
      A_superer[NoctC+3] = octet[i]->Bsuper2er;
      xsuper[NoctC]      = NoctC + 1;
      xsuper[NoctC+1]    = NoctC + 2;
      xsuper[NoctC+2]    = NoctC + 3;
      xsuper[NoctC+3]    = NoctC + 4;
      xsuperer[NoctC]    = 0.001;
      xsuperer[NoctC+1]  = 0.001;
      xsuperer[NoctC+2]  = 0.001;
      xsuperer[NoctC+3]  = 0.001;
      NoctC+=4;
    }
      x[FullOct]           = FullOct+1;
      xer[FullOct]         = 0.01;
      if(octet[i]->A_sum_A[2] != 0 ){
	xq[NQuat]          = NQuat+1;
	xqer[NQuat]        = 0.01;
	A_quat[NQuat]      = Abs(octet[i]->A_sum_A[2]);
	A_quat_er[NQuat]   = Abs(octet[i]->A_sum_Aer[2]);
	NQuat++;
      }
      if(octet[i]->A_sum_B[2] != 0 ){
	xq[NQuat]          = NQuat+1;
	xqer[NQuat]        = 0.01;
	A_quat[NQuat]       = Abs(octet[i]->A_sum_B[2]);
	A_quat_er[NQuat]    = Abs(octet[i]->A_sum_Ber[2]);
	NQuat++;
      }
      
      A_ave[FullOct]       = octet[i]->A_multi[0];
      A_aveB[FullOct]      = octet[i]->A_multi[1];
      A_aveC[FullOct]      = octet[i]->A_multi[2];
      A_aveD[FullOct]      = octet[i]->A_multi[3];
      A_aveE[FullOct]      = octet[i]->A_multi[4];
      A_aveF[FullOct]      = octet[i]->A_multi[5];
      A_aveG[FullOct]      = octet[i]->A_multi[6];
      A_aveH[FullOct]      = octet[i]->A_multi[7];
      A_aveI[FullOct]      = octet[i]->A_multi[8];
      
      A_ave_er[FullOct]    = Abs(octet[i]->A_multier[0]);
      A_aveB_er[FullOct]   = Abs(octet[i]->A_multier[1]);
      A_aveC_er[FullOct]   = Abs(octet[i]->A_multier[2]);
      A_aveD_er[FullOct]   = Abs(octet[i]->A_multier[3]); 
      A_aveE_er[FullOct]   = Abs(octet[i]->A_multier[4]);
      A_aveF_er[FullOct]   = Abs(octet[i]->A_multier[5]);
      A_aveG_er[FullOct]   = Abs(octet[i]->A_multier[6]);
      A_aveH_er[FullOct]   = Abs(octet[i]->A_multier[7]); 
      A_aveI_er[FullOct]   = Abs(octet[i]->A_multier[8]);
      
      A_aves[FullOct]       = octet[i]->A_sum[0];
      A_avesB[FullOct]      = octet[i]->A_sum[1];
      A_avesC[FullOct]      = octet[i]->A_sum[2];
      A_avesD[FullOct]      = octet[i]->A_sum[3];
      A_avesE[FullOct]      = octet[i]->A_sum[4];
      A_avesF[FullOct]      = octet[i]->A_sum[5];
      A_avesG[FullOct]      = octet[i]->A_sum[6];
      A_avesH[FullOct]      = octet[i]->A_sum[7];
      A_avesI[FullOct]      = octet[i]->A_sum[8];
      
      A_aves_er[FullOct]    = Abs(octet[i]->A_sumer[0]);
      A_avesB_er[FullOct]   = Abs(octet[i]->A_sumer[1]);
      A_avesC_er[FullOct]   = Abs(octet[i]->A_sumer[2]);
      A_avesD_er[FullOct]   = Abs(octet[i]->A_sumer[3]); 
      A_avesE_er[FullOct]   = Abs(octet[i]->A_sumer[4]);
      A_avesF_er[FullOct]   = Abs(octet[i]->A_sumer[5]);
      A_avesG_er[FullOct]   = Abs(octet[i]->A_sumer[6]);
      A_avesH_er[FullOct]   = Abs(octet[i]->A_sumer[7]); 
      A_avesI_er[FullOct]   = Abs(octet[i]->A_sumer[8]);
      
      A_sum_A[FullOct]     = Abs(octet[i]->A_sum_A[2]);
      A_sum_B[FullOct]     = Abs(octet[i]->A_sum_B[2]);
      A_multi_A[FullOct]   = Abs(octet[i]->A_multi_A[2]);
      A_multi_B[FullOct]   = Abs(octet[i]->A_multi_B[2]);
      A_sum_Aer[FullOct]   = Abs(octet[i]->A_sum_Aer[2]);
      A_sum_Ber[FullOct]   = Abs(octet[i]->A_sum_Ber[2]);
      A_multi_Aer[FullOct] = Abs(octet[i]->A_multi_Aer[2]);
      A_multi_Ber[FullOct] = Abs(octet[i]->A_multi_Ber[2]);
      FullOct++;
  } 

  foct.close();
  fquart.close();
  fsup.close();
  
  Average_Array(A_ave , A_ave_er,FullOct, A_ave_tot[0], A_ave_toter[0]);
  Average_Array(A_aveB,A_aveB_er,FullOct, A_ave_tot[1], A_ave_toter[1]);
  Average_Array(A_aveC,A_aveC_er,FullOct, A_ave_tot[2], A_ave_toter[2]);
  Average_Array(A_aveD,A_aveD_er,FullOct, A_ave_tot[3], A_ave_toter[3]);
  Average_Array(A_aveE,A_aveE_er,FullOct, A_ave_tot[4], A_ave_toter[4]);
  Average_Array(A_aveF,A_aveF_er,FullOct, A_ave_tot[5], A_ave_toter[5]);
  Average_Array(A_aveG,A_aveG_er,FullOct, A_ave_tot[6], A_ave_toter[6]);
  Average_Array(A_aveH,A_aveH_er,FullOct, A_ave_tot[7], A_ave_toter[7]);
  Average_Array(A_aveI,A_aveI_er,FullOct, A_ave_tot[8], A_ave_toter[8]);
  
  Average_Array(A_aves , A_aves_er,FullOct, A_aves_tot[0], A_aves_toter[0]);
  Average_Array(A_avesB,A_avesB_er,FullOct, A_aves_tot[1], A_aves_toter[1]);
  Average_Array(A_avesC,A_avesC_er,FullOct, A_aves_tot[2], A_aves_toter[2]);
  Average_Array(A_avesD,A_avesD_er,FullOct, A_aves_tot[3], A_aves_toter[3]);
  Average_Array(A_avesE,A_avesE_er,FullOct, A_aves_tot[4], A_aves_toter[4]);
  Average_Array(A_avesF,A_avesF_er,FullOct, A_aves_tot[5], A_aves_toter[5]);
  Average_Array(A_avesG,A_avesG_er,FullOct, A_aves_tot[6], A_aves_toter[6]);
  Average_Array(A_avesH,A_avesH_er,FullOct, A_aves_tot[7], A_aves_toter[7]);
  Average_Array(A_avesI,A_avesI_er,FullOct, A_aves_tot[8], A_aves_toter[8]);
  
  fstream aveGout;
  aveGout.open("check_on_g.txt",fstream::out);
  for(int i = 0; i < FullOct ; i++){
     aveGout << A_aveG[i] << " +/- " << A_aveG_er[i] << "\t\t" << octet[i]->A_sum_A[6] << " +/- " << octet[i]->A_sum_Aer[6] << "\t\t" << octet[i]->A_sum_B[6] << " +/- " << octet[i]->A_sum_Ber[6];
     aveGout << "\t\t" <<  octet[i]->A_multi_A[6] << " +/- " << octet[i]->A_multi_Aer[6] << "\t\t" << octet[i]->A_multi_B[6] << " +/- " << octet[i]->A_multi_Ber[6] << endl;
  }
  aveGout.close();
   
  fstream ach;
  ach.open(Form("analysis_choice_%d.txt",btr[0]->GetGeo()),fstream::out);
	
  for(Int_t i = 0; i < 9 ; i++){
    cout << "Average A for Analysis " << i << "  " << A_ave_tot[i] << " +/- " << A_ave_toter[i] << "  " << A_aves_tot[i] << " +/- " << A_aves_toter[i]<< endl;
    ach << i << "\t" << A_ave_tot[i] << "\t" << A_ave_toter[i] << "\t" << A_aves_tot[i] << "\t" << A_aves_toter[i]<< endl;
  }

  TF1 *fbeta = new TF1("fbeta","sqrt(x*(x+2.*(510.998910)))/(x+510.998910)",0,2000);

  Double_t avebeta = fbeta->Integral(225,675)/450.;
  
  cout << "The average value for beta from 225,675 is " << avebeta << endl;
  for(Int_t i = 0; i < 9 ; i++){
    cout << "Average corrected A for Analysis " << i << "  " << A_aves_tot[i]/(0.5*avebeta) << " +/- " << A_aves_toter[i]/(0.5*avebeta) <<  endl;
  }
  
  TGraphErrors *gsumA    = new TGraphErrors(FullOct,x,A_sum_A,xer,A_sum_Aer);
  TGraphErrors *gsumB    = new TGraphErrors(FullOct,x,A_sum_B,xer,A_sum_Ber);
  TGraphErrors *gmultiA  = new TGraphErrors(FullOct,x,A_multi_A,xer,A_multi_Aer);
  TGraphErrors *gmultiB  = new TGraphErrors(FullOct,x,A_multi_B,xer,A_multi_Ber);
  TGraphErrors *gOctAave = new TGraphErrors(FullOct,x,A_avesC,xer,A_avesC_er);
  TGraphErrors *gSuperA  = new TGraphErrors(NoctC,xsuper,A_super,xsuperer,A_superer);
  TGraphErrors *gQuartetA = new TGraphErrors(NQuat,xq,A_quat,xqer,A_quat_er);
  
  ColorGraphic(gsumA,1,24,1);
  ColorGraphic(gsumB,1,24,1);
  ColorGraphic(gmultiA,1,24,1);
  ColorGraphic(gmultiB,1,24,1);
  ColorGraphic(gOctAave,1,24,1);
  ColorGraphic(gSuperA,1,24,1);
  ColorGraphic(gQuartetA,1,24,1);
   
  TCanvas *cdiff = new TCanvas("cdiff","Sum Product  Difference");
  cdiff->Divide(1,2);
  cdiff->cd(1);
  gsumA->Draw("AP");
  gmultiA->Draw("P");
  cdiff->cd(2);
  gsumB->Draw("AP");
  gmultiB->Draw("P");
  
  TCanvas *cass = new TCanvas("cass","Octet",300,800);
  cass->Divide(1,3);
  
  TF1 *flO = new TF1("flO","[0]",0,noct+1.);
  TF1 *flS = new TF1("flS","[0]",0,NoctC);
  TF1 *flQ = new TF1("flQ","[0]",0,NQuat+1);
  flO->SetLineColor(2);
  //  flO->SetLineStyle(2);
  flS->SetLineColor(2);
  //flS->SetLineStyle(2);
  flQ->SetLineColor(2);
  //flQ->SetLineStyle(2);
   
  Double_t chiSup[20],chiQuart[20],chiOct[20],xchi[20];
  cass->cd(1);
  gOctAave->Fit("flO","RMEQ","goff");
  gQuartetA->Fit("flQ","RMEQ","goff");
  gSuperA->Fit("flS","RMEQ","goff");
  
  Double_t bSup  = flS->GetParameter(0);
  Double_t bQuat = flQ->GetParameter(0);
  Double_t bOct  = flO->GetParameter(0);
  
  Double_t xmin1,xmin2,xmin3;
  
  for(Int_t i = 0 ; i < 20 ; i++){
    
    xchi[i]     = 0.94 + (double)i*0.006;
    chiSup[i]   = 0.;
    chiQuart[i] = 0.;
    chiOct[i]   = 0.;
    
    for(Int_t j = 0 ; j < NoctC ; j++){
      chiSup[i]     += Power((A_super[j] - bSup * xchi[i]) / (A_superer[j]),2);      
      if(j < NQuat  ){
	chiQuart[i] += Power((A_quat[j]  - bQuat* xchi[i]) / (A_quat_er[j]),2);
      }
      if(j < FullOct){
	chiOct[i]   += Power((A_avesC[j] - bOct * xchi[i]) / (A_avesC_er[j]),2);
      }
    }
  }
  
  xmin1 = chiSup[10];
  xmin2 = chiQuart[10];
  xmin3 = chiOct[10];
  
  for(Int_t i = 0 ; i < 20 ; i++){
    chiSup[i] = chiSup[i] - xmin1;
    chiQuart[i] = chiQuart[i] - xmin2;
    chiOct[i] = chiOct[i] - xmin3;
  }
  
  
  flS->SetParameter(0,bSup);
  flO->SetParameter(0,bQuat);
  flQ->SetParameter(0,bOct);
  
  TGraph *gSupChi  = new TGraph(20,xchi,chiSup);
  TGraph *gQuatChi = new TGraph(20,xchi,chiQuart);
  TGraph *gOctChi  = new TGraph(20,xchi,chiOct);
 
  cass->cd(1);
  gOctAave->SetTitle("");//Form("Octet Integral Asymmetry %d - %d keV",(Int_t)(btr[0]->heq->GetBinCenter(nlow)-btr[0]->heq->GetBinWidth(3)/2),
  //(Int_t)(btr[0]->heq->GetBinCenter(nhigh)+btr[0]->heq->GetBinWidth(3)/2)));
  gOctAave->GetYaxis()->SetTitle("Asym_{raw}");
  gOctAave->GetXaxis()->SetTitle("Octet Number");
  gOctAave->GetYaxis()->CenterTitle();
  gOctAave->GetXaxis()->CenterTitle();  
  gOctAave->GetYaxis()->SetRangeUser(0,.1);
  gOctAave->Draw("AP");
  gOctAave->Fit("flO","RMEQ","");
  

  TLegend *legO = new TLegend(0.2,0.7,0.87,0.85);
  legO->AddEntry(gOctAave,Form("A = %6.5f #pm %6.5f",
			       flO->GetParameter(0),flO->GetParError(0)),"lp");
  legO->AddEntry(gOctAave,Form("#chi^{2}/#nu = %6.5f",flO->GetChisquare()/flO->GetNDF()),"");
  legO->SetLineColor(0);
  legO->SetFillColor(0);
  legO->Draw();
  
  cass->cd(2);
  gQuartetA->SetTitle("");//Form("Quartet Integral Asymmetry %d-%d keV",
			  // (Int_t)(btr[0]->heq->GetBinCenter(nlow) -btr[0]->heq->GetBinWidth(3)/2),
  //(Int_t)(btr[0]->heq->GetBinCenter(nhigh)+btr[0]->heq->GetBinWidth(3)/2)));
  gQuartetA->GetYaxis()->SetTitle("Asym_{raw}");
  gQuartetA->GetXaxis()->SetTitle("Quartet Super-Ratio Number");
  gQuartetA->GetYaxis()->CenterTitle();
  gQuartetA->GetXaxis()->CenterTitle();  
  gQuartetA->GetYaxis()->SetRangeUser(0,.1);
  gQuartetA->Draw("AP");
  gQuartetA->Fit("flQ","RMEQ","");
  
  TLegend *legQ = new TLegend(0.2,0.7,0.87,0.85);
  legQ->AddEntry(gQuartetA,Form("A = %6.5f #pm %6.5f",
				flQ->GetParameter(0),flQ->GetParError(0)),"lp");
  legQ->SetFillColor(0);
  legQ->AddEntry(gQuartetA,Form("#chi^{2}/#nu = %6.5f",
			       flQ->GetChisquare()/flQ->GetNDF()),"");
  legQ->SetLineColor(0);
  legQ->Draw();
  
  
  cass->cd(3);
  gSuperA->SetTitle("");//Form("Super-Ratio Integral Asymmetry %d-%d keV",
  //(Int_t)(btr[0]->heq->GetBinCenter(nlow)-btr[0]->heq->GetBinWidth(3)/2),
  //			 (Int_t)(btr[0]->heq->GetBinCenter(nhigh)+btr[0]->heq->GetBinWidth(3)/2)));
  gSuperA->GetYaxis()->SetTitle("Asym_{raw}");
  gSuperA->GetXaxis()->SetTitle("Super-Ratio Number");
  gSuperA->GetYaxis()->CenterTitle();
  gSuperA->GetXaxis()->CenterTitle();  
  gSuperA->GetYaxis()->SetRangeUser(0,.1);
  gSuperA->Draw("AP");
  gSuperA->Fit("flS","RMEQ","");
  
  TLegend *legS = new TLegend(0.2,0.7,0.87,0.85);
  legS->AddEntry(gSuperA,Form("A = %6.5f #pm %6.5f",
			       flS->GetParameter(0),flS->GetParError(0)),"lp");
  legS->AddEntry(gSuperA,Form("#chi^{2}/#nu = %6.5f",
			      flS->GetChisquare()/flS->GetNDF()),"");
  legS->SetLineColor(0);
  legS->SetFillColor(0);
  legS->Draw();
  
  cass->Print(Form("pdf_out/geo%d/asy_groupings_%d.pdf",btr[0]->GetGeo(),btr[0]->GetGeo()));
  
  TCanvas *cChiTest = new TCanvas("cChiTest","Chi Test");
  cChiTest->cd();
  
  ColorGraphic(gSupChi,2,20,2);
  ColorGraphic(gQuatChi,3,20,3);
  ColorGraphic(gOctChi,4,20,4);
  
  gSupChi->GetXaxis()->SetTitle("#Delta A / A");
  gSupChi->GetXaxis()->CenterTitle();
  gSupChi->GetYaxis()->SetTitle("#chi^{2}");
  gSupChi->GetYaxis()->CenterTitle();
  gSupChi->SetTitle("#chi^{2} Contour");
  gSupChi->Draw("APC");
  gQuatChi->Draw("PC");
  gOctChi->Draw("PC");
  
  TLegend *lchi2 = new TLegend(0.3,0.6,0.8,0.85);
  lchi2->AddEntry(gSupChi,"Pulse Pair","pl");
  lchi2->AddEntry(gQuatChi,"Quartet","pl");
  lchi2->AddEntry(gOctChi,"Octet","pl");
  lchi2->SetFillColor(0);
  lchi2->Draw();
 
  cChiTest->Print(Form("chitests_%d.pdf",btr[0]->GetGeo()));
}
//====================================================================================================
void Average_A()
{
  using namespace TMath;
   
  Double_t x[200],y[200],yer[200];
  Double_t xb[200],yb[200],yber[200];
  Double_t ytot[200],xtot[200],ytoter[200];
  
  TF1 *fbeta = new TF1("fbeta","sqrt(x*(x+2.*(510.998910)))/(x+510.998910)",0,2000);
  
  for(Int_t i = 0;i<200;i++){

    x[i]      = 0.;
    y[i]      = 0.;
    yer[i]    = 0.;
    xb[i]     = 0.;
    yb[i]     = 0.;
    yber[i]   = 0.;
    xtot[i]   = 0.;
    ytot[i]   = 0.;
    ytoter[i] = 0.;
  } 
 
  // Calculate the weighted average of the Asymmetry
  //================================================================================================
  for(Int_t j = 0 ; j < noct ; j++){
    for(Int_t i = 1; i <= octet[0]->hAsyA[2]->GetNbinsX(); i++ ){
      if(!(IsNaN(octet[j]->hAsyA[2]->GetBinContent(i))) && octet[j]->hAsyA[2]->GetBinError(i) != 0 ){
	y[i-1] += octet[j]->hAsyA[2]->GetBinContent(i) / Power(octet[j]->hAsyA[2]->GetBinError(i),2);
	yer[i-1] += 1./Power(octet[j]->hAsyA[2]->GetBinError(i),2);
      }
      if(!(IsNaN(octet[j]->hAsyB[2]->GetBinContent(i)))&& octet[j]->hAsyB[2]->GetBinError(i) !=0){
	yb[i-1] += octet[j]->hAsyB[2]->GetBinContent(i) / Power(octet[j]->hAsyB[2]->GetBinError(i),2);
	yber[i-1] += 1./Power(octet[j]->hAsyB[2]->GetBinError(i),2);
      }
      if( !(IsNaN(octet[j]->hAsyTot[2]->GetBinContent(i))) && octet[j]->hAsyTot[2]->GetBinError(i) != 0){
	ytot[i-1] += octet[j]->hAsyTot[2]->GetBinContent(i) /
	  Power(octet[j]->hAsyTot[2]->GetBinError(i),2);
	ytoter[i-1] += 1./Power(octet[j]->hAsyTot[2]->GetBinError(i),2);
      }
    }
  }
  //=================================================================================================
  for(Int_t i = 0 ; i <= octet[0]->hAsyA[2]->GetNbinsX()-1 ;i++){
    x[i]   = octet[0]->hAsyA[2]->GetBinCenter(i+1);
    if(yer[i] > 0 ){
      y[i]   = y[i] / yer[i];
      yer[i] = sqrt( 1./yer[i]);
    }
    if(IsNaN(y[i])){
      y[i] = 0.0;
      yer[i] = 0.1;
    }
  }
  
  for(Int_t i = 0 ; i <= octet[0]->hAsyB[2]->GetNbinsX()-1 ;i++){
    xb[i]   = octet[0]->hAsyB[2]->GetBinCenter(i+1);
    yb[i]   = yb[i] / yber[i];
    yber[i] = sqrt( 1./yber[i]);
    if(IsNaN(yb[i])){
        yb[i] = 0.0;
	yber[i] = 0.1;
    }
  }

  
  for(Int_t i = 0 ; i <= octet[0]->hAsyTot[2]->GetNbinsX()-1 ;i++){
    xtot[i]   = octet[0]->hAsyTot[2]->GetBinCenter(i+1);
    ytot[i]   = ytot[i] / ytoter[i];
    ytoter[i] = sqrt( 1./ytoter[i]);
    if(IsNaN(ytot[i])){
        ytot[i] = 0.0;
	ytoter[i] = 0.1;
    }
  }
  
  fstream fasye;
  
  fasye.open(Form("asymmetry_e_%d.txt",btr[0]->GetGeo()),fstream::out);
  for(Int_t i = 0 ; i< octet[0]->hAsyTot[2]->GetNbinsX() ; i++){
    fasye << xtot[i] << "\t" << ytot[i] << "\t" << ytoter[i] << "\t" << y[i];
    fasye << "\t" << yer[i] << "\t" << yb[i] << "\t" << yber[i] << endl;
  }
  fasye.close();
  
  TGraphErrors *gAvsE    = new TGraphErrors(80,x,y,0,yer);
  TGraphErrors *gAvsEb   = new TGraphErrors(80,xb,yb,0,yber);
  TGraphErrors *gAvsEtot = new TGraphErrors(80,xtot,ytot,0,ytoter);

  ColorGraphic(gAvsE,2,20,2);
  ColorGraphic(gAvsEb,4,20,2);
  ColorGraphic(gAvsEtot,2,20,2);
  
  TCanvas *cGreat = new TCanvas("cGreat","Energy vs. A",500,800);
  cGreat->Divide(1,2);
  cGreat->cd(1);
  
  TF1 *fAbeta = new TF1("fAbeta","[0]",275,625);
  TF1 *fAbeta2 = new TF1("fAbeta2","[0]",275,625);
  
  gAvsE->Draw("AP");
  gAvsE->SetMarkerColor(1);
  gAvsE->SetMarkerStyle(24);
  gAvsE->SetMarkerSize(0.6);
  gAvsE->SetLineColor(1);
  gAvsEb->SetMarkerColor(1);
  gAvsEb->SetMarkerStyle(25);
  gAvsEb->SetMarkerSize(0.6);
  gAvsEb->SetLineColor(1);
  fAbeta->SetLineColor(2);
  fAbeta2->SetLineColor(2);
  gAvsE->GetXaxis()->SetRangeUser(0,900);
  gAvsE->GetYaxis()->SetRangeUser(-0.3,0.3);  
  gAvsE->SetTitle("");
  gAvsE->Fit("fAbeta","REMQ");
  gAvsE->GetYaxis()->SetTitle("#LT A_{o} #GT");
  gAvsE->GetXaxis()->SetTitle("Energy (keV)");
  gAvsE->GetXaxis()->CenterTitle();
  gAvsE->GetYaxis()->CenterTitle(); 

  gAvsEb->Draw("P");
  gAvsEb->Fit("fAbeta2","RMEQ","SAME");
  fAbeta->SetLineColor(2);
  fAbeta2->SetLineColor(2);
  gAvsE->GetYaxis()->SetRangeUser(-0.4,0.4);
  
  Double_t Aave = ((fAbeta->GetParameter(0)/Power(fAbeta->GetParError(0),2))
      + Abs((fAbeta2->GetParameter(0)/Power(fAbeta2->GetParError(0),2)))) /
      (1./Power(fAbeta->GetParError(0),2) + 1./Power(fAbeta2->GetParError(0),2));
  
  TLegend *lAs = new TLegend(0.3,0.7,0.9,0.9);
  lAs->SetFillColor(0);
  lAs->AddEntry(gAvsE,Form("A - Quartet A = %6.4f #pm %6.4f",fAbeta->GetParameter(0),
		fAbeta->GetParError(0)),"lp");
  lAs->AddEntry(gAvsE,Form("#chi^{2}/n = %6.4f",fAbeta->GetChisquare()/fAbeta->GetNDF()),"");
  lAs->AddEntry(gAvsEb,Form("B - Quartet A = %6.4f #pm %6.4f",fAbeta2->GetParameter(0),
		fAbeta2->GetParError(0)),"lp");
  lAs->AddEntry(gAvsEb,Form("#chi^{2}/n = %6.4f",fAbeta2->GetChisquare()/fAbeta2->GetNDF()),"");
  lAs->Draw("same");

  cGreat->cd(2);

  TF1 *fAbetat = new TF1("fAbetat","[0]",275,625);
  TF1 *fAsyOvr = new TF1("fAsyOvr","[0]",1000,1500);

  fAsyOvr->SetLineStyle(2);

  gAvsEtot->Draw("AP");
  gAvsEtot->SetTitle("");
  gAvsEtot->SetMarkerStyle(24);
  gAvsEtot->SetMarkerColor(1);
  gAvsEtot->SetLineColor(1);
  gAvsEtot->GetXaxis()->SetRangeUser(0,900);
  gAvsEtot->GetYaxis()->SetRangeUser(0.,0.2);  
  gAvsEtot->GetYaxis()->SetTitle("|#LT A_{o} #GT|");
  gAvsEtot->GetXaxis()->SetTitle("Energy (keV)");
  gAvsEtot->GetXaxis()->CenterTitle();
  gAvsEtot->GetYaxis()->CenterTitle();  
  fAbetat->SetLineColor(2);
  gAvsEtot->Fit("fAbetat","RMEQ");

  fAbetat->Draw("same");
  
  TLegend *lAst = new TLegend(0.2,0.2,0.7,0.4);
  lAst->SetFillColor(0);
  lAst->AddEntry(gAvsEtot,Form(" A_{0} = -%6.4f #pm %6.4f",fAbetat->GetParameter(0),
		fAbetat->GetParError(0)),"lp");
  lAst->AddEntry(gAvsEtot,
		 Form("#chi^{2}/#nu = %6.4f",fAbetat->GetChisquare()/fAbetat->GetNDF()),"");

  lAst->Draw("same");

  cGreat->Print(Form("pdf_out/geo%d/final_raw_asymmetry_%d.pdf",btr[0]->GetGeo(),btr[0]->GetGeo()));
    
}

void Collect_TvsE()
{

  // Calculate the average Energy vs. dt for the type 1 backscatters 
  // then compare it with the results from Monte Claro.

  Int_t rbins = 100;
  Double_t ecount[rbins][rbins],wcount[rbins][rbins];
  
  hTotETimeVsE = new TH2F("hTotETimeVsE",
			  "East Type 1 Backscattering Energy vs. time;  X_{west} (mm) ;"
			  "Y_{west} (mm)",rbins,0,1000,rbins,0,200);
  
  hTotWTimeVsE = new TH2F("hTotWTimeVsE",
			  "West Type 1 Backscattering Rotation ; X_{east} (mm) ;"
			  "Y_{east} (mm)",rbins,0,1000,rbins,0,200);

  TH2F *hTotWP = new TH2F("hTotWP","West Type 1",rbins,0,1000,rbins,0,200);
  TH2F *hTotWS = new TH2F("hTotWS","West Type 1",rbins,0,1000,rbins,0,200);
  TH2F *hTotEP = new TH2F("hTotEP","West Type 1",rbins,0,1000,rbins,0,200);
  TH2F *hTotES = new TH2F("hTotES","West Type 1",rbins,0,1000,rbins,0,200);
  // Fill the Collection 2d histogram

  for(Int_t ibin = 0 ; ibin < rbins ; ibin++){
    for(Int_t jbin = 0 ; jbin < rbins ; jbin++){
      ecount[ibin][jbin] = 0.;
      wcount[ibin][jbin] = 0.;
    }
  }
  
  for(Int_t i = 0 ; i <nbeta ; i++){
    for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
        if(btr[i]->hETimeVsE->GetBinContent(ibin,jbin) > 0)
	  ecount[ibin-1][jbin-1] += btr[i]->hETimeVsE->GetBinContent(ibin,jbin)*btr[i]->rtime_e;
	if(btr[i]->hWTimeVsE->GetBinContent(ibin,jbin) > 0)
	  wcount[ibin-1][jbin-1] += btr[i]->hWTimeVsE->GetBinContent(ibin,jbin)*btr[i]->rtime_w;
	hTotWP->SetBinContent(ibin,jbin,hTotWP->GetBinContent(ibin,jbin) + btr[i]->hWTimeVsEP->GetBinContent(ibin,jbin)*btr[i]->rtime_w);
	hTotWS->SetBinContent(ibin,jbin,hTotWS->GetBinContent(ibin,jbin) + btr[i]->hWTimeVsES->GetBinContent(ibin,jbin)*btr[i]->rtime_w);

	hTotEP->SetBinContent(ibin,jbin,hTotEP->GetBinContent(ibin,jbin) + btr[i]->hETimeVsEP->GetBinContent(ibin,jbin)*btr[i]->rtime_e);
	hTotES->SetBinContent(ibin,jbin,hTotES->GetBinContent(ibin,jbin) + btr[i]->hETimeVsES->GetBinContent(ibin,jbin)*btr[i]->rtime_e);
      }
    }
  }
  
  for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
    for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
      hTotETimeVsE->SetBinContent(ibin,jbin,ecount[ibin-1][jbin-1]);
      hTotWTimeVsE->SetBinContent(ibin,jbin,wcount[ibin-1][jbin-1]);
    }
  }
  // Finished Filling now compare...

  TCanvas *check = new TCanvas("check","Timing check");
  check->Divide(2,2);
  check->cd(1);
  hTotETimeVsE->Draw("colz");
  check->cd(2);
  hTotWTimeVsE->Draw("colz");

  TProfile *hTvsE_E = hTotETimeVsE->ProfileY("e",1,80);
  TProfile *hTvsE_W = hTotWTimeVsE->ProfileY("w",1,80);
 
  TProfile *hTvsE_WP = hTotWP->ProfileY("w1",1,80);
  TProfile *hTvsE_WS = hTotWS->ProfileY("w2",1,80);

  TProfile *hTvsE_EP = hTotEP->ProfileY("e1",1,80);
  TProfile *hTvsE_ES = hTotES->ProfileY("e2",1,80);
  
  fstream fback;
  fback.open(Form("backscatter_e_vs_t_%d.txt",btr[0]->GetGeo()),fstream::out);
  
  Double_t xe[100],ye[100],xw[100],yw[100],yet,ywt;
  Double_t yer[100],xer[100],ywr[100];

  for(Int_t i = 0 ; i < 100; i++){
    xe[i] = 1.+2.*i;
    xw[i] = 1.+2.*i;
    ye[i] = 0.;
    yw[i] = 0.;
    yet = 0.;
    ywt = 0.;
    for(Int_t j = 0 ; j < 70 ; j++){
      ye[i] +=  ecount[j][i]*j*10.;
      yw[i] +=  wcount[j][i]*j*10.;
      yet   +=  ecount[j][i];
      ywt   +=  wcount[j][i];
    }
    if(yet > 0){
      ye[i] = ye[i] / yet;
      yer[i] = ye[i]/sqrt(yet);
    }
    if(ywt > 0){
      yw[i]  = yw[i] / ywt;
      ywr[i] = yw[i] /sqrt(ywt);
    }
    fback << xe[i] <<"\t" << ye[i] << "\t" << yer[i] << "\t" << yw[i] << "\t" << ywr[i] << endl;
  } 

  fback.close();
  fstream fbck2;
  Double_t Toffset  = 0.;
  Double_t Toffsete = 0.;

  // Open and fill the expected average energy from monte claro
  if(btr[0]->GetGeo() == 0){
    fbck2.open("type_1_e_vs_time.txt",fstream::in);
 //   fbck2.open("evt_ave_a.txt",fstream::in);
    Toffset  = -10.0;
    Toffsete =  -4.; 
  } else if(btr[0]->GetGeo() == 1){
    fbck2.open("type_1_e_vs_time_1.txt",fstream::in);
    Toffset = -6.;
    Toffsete=  10.;
  }
  else if(btr[0]->GetGeo() == 2 || btr[0]->GetGeo() == 3 ){
    fbck2.open("type_1_e_vs_time_2.txt",fstream::in);
    if(btr[0]->GetGeo() == 2){
       Toffset = -10.;
       Toffsete = 10;
    } else {
      Toffset = -5.;
      Toffsete = 5.;
   }
  }

  Int_t ntbins = 100;
  Double_t xmc[ntbins],emc[ntbins],emcer[ntbins],xmcw[ntbins];

  for(Int_t i = 0 ; i < ntbins ; i++){
    fbck2 >> xmc[i] >> emc[i] >> emcer[i];
    xmcw[i] = xmc[i] - Toffset;
    xmc[i]  = xmc[i] + Toffsete;
  } 
  fbck2.close();
  
  TGraphErrors *gMcTe  = new TGraphErrors(ntbins,xmc,emc,0,emcer);
  TGraphErrors *gMcTw  = new TGraphErrors(ntbins,xmcw,emc,0,emcer);

  TGraphErrors *gAveTe = new TGraphErrors(100,xe,ye,0,yer);
  TGraphErrors *gAveTw = new TGraphErrors(100,xw,yw,0,ywr);
  
  ColorGraphic(gAveTe,2,20,2);
  ColorGraphic(gAveTw,4,20,2);
  ColorGraphic(gMcTe,1,7,2);
  ColorGraphic(gMcTw,1,7,2); 

  ColorGraphic(hTvsE_E,1,1,1,24);
  ColorGraphic(hTvsE_W,1,1,1,24);

  ColorGraphic(hTvsE_WS,2,1,1,24);
  ColorGraphic(hTvsE_WP,4,1,1,24);

  ColorGraphic(hTvsE_ES,1,1,1,24);
  ColorGraphic(hTvsE_EP,1,1,1,24);



  check->cd(3);
  gMcTe->GetXaxis()->SetTitle("Common Stop Timing (ns)");
  gMcTw->GetXaxis()->SetTitle("Common Stop Timing (ns)");
  gMcTe->GetYaxis()->SetTitle("#LT E_{vis} #GT (keV)");
  gMcTw->GetYaxis()->SetTitle("#LT E_{vis} #GT (keV)");

  gMcTe->GetXaxis()->CenterTitle();
  gMcTe->GetYaxis()->CenterTitle();
  gMcTw->GetXaxis()->CenterTitle();
  gMcTw->GetYaxis()->CenterTitle();
  
  gMcTw->GetXaxis()->SetRangeUser(20,150);

  gMcTe->SetTitle("");
  gMcTw->SetTitle("");

  TF1 *ffit0 = new TF1("ffit0","pol1(0)+expo(2)",20,140);
  TF1 *ffit1 = new TF1("ffit1","pol1(0)+expo(2)",20,140);
 
  ffit0->SetParameter(0,250);
  ffit1->SetParameter(0,250);
  ffit0->SetParLimits(0,220,320);
  ffit1->SetParLimits(0,220,320);

  ffit1->SetParLimits(1,0,10);
  ffit0->SetParLimits(1,0,10);

  ffit1->SetParameter(1,.4);
  ffit0->SetParameter(1,.4);

  ffit0->SetParameter(2,-5.);
  ffit1->SetParameter(2,-5.);
  ffit0->SetParLimits(2,-20,0);
  ffit1->SetParLimits(2,-20,0);

  ffit0->SetParameter(3,0.05);
  ffit1->SetParameter(3,0.05); 
  ffit0->SetParLimits(3,0,1);
  ffit1->SetParLimits(3,0,1);

  ffit0->SetLineColor(2);
  ffit0->SetLineStyle(2);
  ffit1->SetLineColor(2);

  gMcTe->Draw("AP");
  //gAveTe->Draw("P");
  hTvsE_E->Draw("same");
  check->cd(4);

  gMcTw->GetYaxis()->SetRangeUser(0,800);

  gMcTw->Draw("AP");
  gMcTw->Fit(ffit0,"me","same",40,138);
  ffit0->Draw("same");
  //gAveTw->Draw("P");
  hTvsE_W->Draw("same");
  hTvsE_W->Fit(ffit1,"me","same",40,138);
  hTvsE_WP->Draw("same");
  hTvsE_WS->Draw("same");

  TLegend *levt = new TLegend(0.22,0.6,0.5,0.85);
  levt->SetFillColor(0);
  levt->SetLineColor(0);
  levt->AddEntry(gMcTw,"MC","lp");
  levt->AddEntry(ffit0,"Fit to MC","l");
  levt->AddEntry(hTvsE_W,"Data","lp");
  levt->AddEntry(ffit1,"Fit to Data","l");
  levt->Draw();
  
  TLegend *levt2 = new TLegend(0.45,0.6,0.7,0.85);
  levt2->SetFillColor(0);
  levt2->SetLineColor(0);
  levt2->AddEntry(hTvsE_WP,"1^{st} Detector","pl");
  levt2->AddEntry(hTvsE_WS,"2^{nd} Detector","pl");
  levt2->Draw();


 check->Print(Form("pdf_out/geo%d/e_vs_timing_%d.pdf",btr[0]->GetGeo(),btr[0]->GetGeo()));

}

void Collect_23Anode()
{
   hTE23Anode2d = new TH2F("hTE23Anode2d",
		     "East Type 23 E vs. W Anode; East Anode ; West Anode"
		     ,100,0,3000,100,0,3000);
		     
   hTW23Anode2d = new TH2F("hTW23Anode2d",
		     "West Type 23 E vs. W Anode; East Anode ; West Anode"
		     ,100,0,3000,100,0,3000);
  
  Int_t rbins = 100;     
   Double_t ecount[rbins][rbins],wcount[rbins][rbins];
  
  for(Int_t ibin = 0 ; ibin < rbins ; ibin++){
    for(Int_t jbin = 0 ; jbin < rbins ; jbin++){
      ecount[ibin][jbin] = 0.;
      wcount[ibin][jbin] = 0.;
    }
  }
  
  for(Int_t i = 0 ; i <nbeta ; i++){
    for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
        if(btr[i]->hE23Anode2d->GetBinContent(ibin,jbin) > 0)
	ecount[ibin-1][jbin-1] += btr[i]->hE23Anode2d->GetBinContent(ibin,jbin)*btr[i]->rtime_e;
	if(btr[i]->hW23Anode2d->GetBinContent(ibin,jbin) > 0)
	wcount[ibin-1][jbin-1] += btr[i]->hW23Anode2d->GetBinContent(ibin,jbin)*btr[i]->rtime_w;
      }
    }
  }
  
  for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
    for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
      hTE23Anode2d->SetBinContent(ibin,jbin,ecount[ibin-1][jbin-1]);
      hTW23Anode2d->SetBinContent(ibin,jbin,wcount[ibin-1][jbin-1]);
    }
  }
  
  TCanvas *cAnode = new TCanvas("cAnode","cAnode");
  cAnode->Divide(2,1);
  cAnode->cd(1);
  hTE23Anode2d->Draw("colz");
  cAnode->cd(2);
  hTW23Anode2d->Draw("colz");
  
}

void Collect_Stuff()
{
  TCanvas *cfart = new TCanvas("cfart","Signal To Noise");
  cfart->cd(0);
  Double_t xS2N[10000],yS2N[10000],yS2Nw[10000];
  
  
  for(int i = 0 ; i < nbeta ; i++){
       xS2N[i] = btr[i]->GetRunNumber();
       yS2N[i] = btr[i]->E_Sig_Nos;
       yS2Nw[i] = btr[i]->W_Sig_Nos;
  }
  
  gS2Ne = new TGraph(nbeta,xS2N,yS2N);
  ColorGraphic(gS2Ne,2,20,2);
  gS2Nw = new TGraph(nbeta,xS2N,yS2Nw);
  ColorGraphic(gS2Nw,4,20,2);
  
  gS2Ne->Draw("AP");
  gS2Nw->Draw("P");

}

void Collect_Rad()
{


using namespace TMath;
   
  Double_t x[20],xer[20],y[20],yer[20];
  
  for(Int_t i = 0; i < 12 ;i++){
    x[i]      = 0.;
    xer[i]    = 0.;
    y[i]      = 0.;
    yer[i]    = 0.;
  } 
  
  for(Int_t j = 0 ; j < noct ; j++){
    for(Int_t i = 0; i < 12 ; i++ ){
      if(octet[j]->A_rader[i]!= 0 && !IsNaN(octet[j]->A_rader[i])){
	//cout << j << "\t" << octet[j]->A_rad[i] << "\t" << octet[j]->A_rader[i] << endl;	
	y[i] += TMath::Abs(octet[j]->A_rad[i] /
			   Power(octet[j]->A_rader[i],2));
	yer[i] += 1./Power(octet[j]->A_rader[i],2);
      }
    }
  }
  
  for(Int_t i = 0 ; i < 12 ;i++){
    x[i]   = (i+1)*5;
    y[i]   = y[i] / yer[i];
    yer[i] = sqrt( 1./yer[i]);
    if(IsNaN(y[i])){
	y[i] = 0.0;
	yer[i] = 0.1;
    }
  }
  
  
  TGraphErrors *gAvsRad = new TGraphErrors(12,x,y,0,yer);
  ColorGraphic(gAvsRad,2,20,2);
  
  TF1 *frad1 = new TF1("frad1","[0]",0,60);
  frad1->SetLineColor(2);

  TCanvas *cRadA = new TCanvas("cRadA","Radi Dependent");
  //cRadA->Divide(2,1);
  cRadA->cd();
  gAvsRad->SetTitle("Raw Octet Asymmetry for R_{i} < R < R_{i+1}");
  gAvsRad->GetXaxis()->SetTitle("Radius (mm)");
  gAvsRad->GetYaxis()->SetTitle("A_{raw} 0 < E < 0.8MeV");
  gAvsRad->Draw("AP");
  gAvsRad->Fit("frad1","REMQ");
  
  TLegend *lrado = new TLegend(0.6,0.6,0.9,0.9);
  lrado->AddEntry(gAvsRad,Form("A = %6.4f #pm %6.4f",frad1->GetParameter(0),
						    frad1->GetParError(0)),"lp");
  lrado->AddEntry(gAvsRad,Form("#chi^{2}/#nu = %6.4f",frad1->GetChisquare()/frad1->GetNDF()),"");
  lrado->SetFillColor(0);
  lrado->Draw();
  

}

void Collect_Energy_Spectra()
{
  using namespace TMath;

  // Initialize Histogram !!
  Define_E_Spec();
  // Signal Times
  Double_t west_time_on  = 0.;
  Double_t east_time_on  = 0.;
  Double_t west_time_off = 0.;
  Double_t east_time_off = 0.;
  // Background times
  Double_t west_time_onb  = 0.;
  Double_t east_time_onb  = 0.;
  Double_t west_time_offb = 0.;
  Double_t east_time_offb = 0.;
  //
  // Loop through and fill histograms for the signal and back ground spectra
  //
  Fill_Foreground(&west_time_on,&east_time_on,&west_time_off,&east_time_off);
  Fill_Background(&west_time_onb,&east_time_onb,&west_time_offb,&east_time_offb);
  
  Double_t totaltime = west_time_off + west_time_on + east_time_off+ east_time_on;
  
  // Set error bars in the collected histograms
  fstream parbck;
  parbck.open(Form("debug/back_ground_%d_on.txt",btr[0]->GetGeo()),fstream::out);
  fstream parbckoff;
  parbckoff.open(Form("debug/back_ground_%d_off.txt",btr[0]->GetGeo()),fstream::out);
  parbck << "E bin\tRe\tRe_er\tRw\tRw_er\tRe_1\tRe_1er\tRw_1\tRw_1er\tRe_2\tRe_2er\tRw_2\tRw_2er\t" << endl;
  parbckoff << "E bin\tRe\tRe_er\tRw\tRw_er\tRe_1\tRe_1er\tRw_1\tRw_1er\tRe_2\tRe_2er\tRw_2\tRw_2er\t" << endl;
  fstream parsig;
  parsig.open(Form("debug/signal_%d_on.txt",btr[0]->GetGeo()),fstream::out);
  fstream parsigoff;
  parsigoff.open(Form("debug/signal_%d_off.txt",btr[0]->GetGeo()),fstream::out);
  parsig << "E bin\tRe\tRe_er\tRw\tRw_er\tRe_1\tRe_1er\tRw_1\tRw_1er\tRe_2\tRe_2er\tRw_2\tRw_2er\t" << endl;
  parsigoff << "E bin\tRe\tRe_er\tRw\tRw_er\tRe_1\tRe_1er\tRw_1\tRw_1er\tRe_2\tRe_2er\tRw_2\tRw_2er\t" << endl;
  
  Double_t Etype1on   =  hEFlipperOn_I->Integral(nlow,nhigh)  / hEFlipperOn->Integral(nlow,nhigh);
  Double_t Etype23on  =  hEFlipperOn_2->Integral(nlow,nhigh)  / hEFlipperOn->Integral(nlow,nhigh);
  Double_t Etype1off  =  hEFlipperOff_I->Integral(nlow,nhigh) / hEFlipperOff->Integral(nlow,nhigh);
  Double_t Etype23off =  hEFlipperOff_2->Integral(nlow,nhigh) / hEFlipperOff->Integral(nlow,nhigh);
  Double_t Wtype1on   =  hWFlipperOn_I->Integral(nlow,nhigh)  / hWFlipperOn->Integral(nlow,nhigh);
  Double_t Wtype23on  =  hWFlipperOn_2->Integral(nlow,nhigh)  / hWFlipperOn->Integral(nlow,nhigh);
  Double_t Wtype1off  =  hWFlipperOff_I->Integral(nlow,nhigh) / hWFlipperOff->Integral(nlow,nhigh);
  Double_t Wtype23off =  hWFlipperOff_2->Integral(nlow,nhigh) / hWFlipperOff->Integral(nlow,nhigh);
  
  Double_t erEt1on    =  Etype1on*sqrt( hEFlipperOn_I->Integral(nlow,nhigh)/ east_time_on);
  Double_t erEt23on   =  Etype23on*sqrt( hEFlipperOn_2->Integral(nlow,nhigh)/ east_time_on);
  Double_t erEt1off   =  Etype1off*sqrt( hEFlipperOff_I->Integral(nlow,nhigh)/ east_time_off);
  Double_t erEt23off  =  Etype23off*sqrt( hEFlipperOff_2->Integral(nlow,nhigh)/ east_time_off);
  Double_t erWt1on    =  Wtype1on*sqrt( hWFlipperOn_I->Integral(nlow,nhigh)/ west_time_on);
  Double_t erWt23on   =  Wtype23on*sqrt( hWFlipperOn_2->Integral(nlow,nhigh)/ west_time_on);
  Double_t erWt1off   =  Wtype1off*sqrt( hWFlipperOff_I->Integral(nlow,nhigh)/ west_time_off);
  Double_t erWt23off  =  Wtype23off*sqrt( hWFlipperOff_2->Integral(nlow,nhigh)/ west_time_off);

  Double_t Etype1ont   =  hEFlipperOn_I->Integral()  / hEFlipperOn->Integral();
  Double_t Etype23ont  =  hEFlipperOn_2->Integral()  / hEFlipperOn->Integral();
  Double_t Etype1offt  =  hEFlipperOff_I->Integral() / hEFlipperOff->Integral();
  Double_t Etype23offt =  hEFlipperOff_2->Integral() / hEFlipperOff->Integral();
  Double_t Wtype1ont   =  hWFlipperOn_I->Integral()  / hWFlipperOn->Integral();
  Double_t Wtype23ont  =  hWFlipperOn_2->Integral()  / hWFlipperOn->Integral();
  Double_t Wtype1offt  =  hWFlipperOff_I->Integral() / hWFlipperOff->Integral();
  Double_t Wtype23offt =  hWFlipperOff_2->Integral() / hWFlipperOff->Integral();
  
  Double_t erEt1ont    =  Etype1on*sqrt( hEFlipperOn_I->Integral()/ east_time_on);
  Double_t erEt23ont   =  Etype23on*sqrt( hEFlipperOn_2->Integral()/ east_time_on);
  Double_t erEt1offt   =  Etype1off*sqrt( hEFlipperOff_I->Integral()/ east_time_off);
  Double_t erEt23offt  =  Etype23off*sqrt( hEFlipperOff_2->Integral()/ east_time_off);
  Double_t erWt1ont    =  Wtype1on*sqrt( hWFlipperOn_I->Integral()/ west_time_on);
  Double_t erWt23ont   =  Wtype23on*sqrt( hWFlipperOn_2->Integral()/ west_time_on);
  Double_t erWt1offt   =  Wtype1off*sqrt( hWFlipperOff_I->Integral()/ west_time_off);
  Double_t erWt23offt  =  Wtype23off*sqrt( hWFlipperOff_2->Integral()/ west_time_off);
  
  //-----------------------------------------------------------------------------
  Double_t binwidths = hWFlipperOff_2->GetBinWidth(2);
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
  cout << " Backscattering fraction for " << hWFlipperOff_2->GetBinCenter(nlow)-binwidths/2 << " - ";
  cout <<  hWFlipperOff_2->GetBinCenter(nhigh)+binwidths/2 << endl;
  cout << "East On Type 1 "   << Etype1on  << " +/- " << erEt1on <<   endl;
  cout << "East On Type 2/3 " << Etype23on << " +/- " << erEt23on << endl;
			      
  cout << "East Off Type 1 "   << Etype1off << " +/- " << erEt1off << endl;
  cout << "East Off Type 2/3 " << Etype23off<< " +/- " << erEt23off << endl;
			      
  cout << "West On Type 1 "   << Wtype1on << " +/- " << erWt1on<< endl;
  cout << "West On Type 2/3 " << Wtype23on << " +/- " << erWt23on<< endl;
			      
  cout << "West Off Type 1 "   << Wtype1off  << " +/- " << erWt1off<< endl;
  cout << "West Off Type 2/3 " << Wtype23off << " +/- " << erWt23off<< endl;
  cout << "----------------------------------------------------------" << endl;
  cout << "Total energy range" << endl;
  cout << "East On Type 1 "   << Etype1ont  << " +/- " << erEt1ont <<   endl;
  cout << "East On Type 2/3 " << Etype23ont << " +/- " << erEt23ont << endl;
			      
  cout << "East Off Type 1 "   << Etype1offt << " +/- " << erEt1offt << endl;
  cout << "East Off Type 2/3 " << Etype23offt<< " +/- " << erEt23offt << endl;
			      
  cout << "West On Type 1 "   << Wtype1ont << " +/- " << erWt1ont<< endl;
  cout << "West On Type 2/3 " << Wtype23ont << " +/- " << erWt23ont<< endl;
			      
  cout << "West Off Type 1 "   << Wtype1offt  << " +/- " << erWt1offt<< endl;
  cout << "West Off Type 2/3 " << Wtype23offt << " +/- " << erWt23offt<< endl;

  cout << "-------------------------------------------------------------" << endl; 

  cout << "Total East Counts On  (0-800keV) " << hEFlipperOn->Integral();
  cout << " Bck (0-800keV) " << hEFlipperOn_B->Integral();
  cout << " Sig / Noise " << hEFlipperOn->Integral()/ hEFlipperOn_B->Integral() << endl;

  cout << "Total East Counts Off (0-800keV) " << hEFlipperOff->Integral();
  cout << " Bck (0-800keV) " << hEFlipperOff_B->Integral();
  cout << " Sig / Noise " << hEFlipperOff->Integral()/ hEFlipperOff_B->Integral() << endl;

  cout << "Total West Counts On  (0-800keV) " << hWFlipperOn->Integral();
  cout << " Bck (0-800keV) " << hWFlipperOn_B->Integral();
  cout << " Sig / Noise " << hWFlipperOn->Integral()/ hWFlipperOn_B->Integral()<< endl;

  cout << "Total West Counts Off (0-800keV) " << hWFlipperOff->Integral();
  cout << " Bck (0-800keV) " << hWFlipperOff_B->Integral();
  cout << " Sig / Noise " << hWFlipperOff->Integral()/ hWFlipperOff_B->Integral() << endl;

  cout << "Total Counts (225-675keV) " << (hEFlipperOn->Integral(nlow,nhigh)*east_time_on + 
    east_time_on*hEFlipperOff->Integral(nlow,nhigh) + hWFlipperOn->Integral(nlow,nhigh)*west_time_on + west_time_off*hWFlipperOff->Integral(nlow,nhigh)) << endl;
  cout << "Total Counts (0-800keV) " << (hEFlipperOn->Integral()*east_time_on + east_time_on*hEFlipperOff->Integral() 
					   + hWFlipperOn->Integral()*west_time_on + west_time_off*hWFlipperOff->Integral()) << endl;
  
  
  for(Int_t i = 1; i < hEFlipperOn->GetNbinsX()+1 ; i++){
    
    parbck << i << "\t" << setprecision(4) << hEFlipperOn_B->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOn_B->GetBinError(i) << "\t";
    parbck << setprecision(4) << hWFlipperOn_B->GetBinContent(i)   << "\t"<< setprecision(4) << hWFlipperOn_B->GetBinError(i)   << "\t";
    parbck << setprecision(4) << hEFlipperOn_B_I->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOn_B_I->GetBinError(i) << "\t";
    parbck << setprecision(4) << hWFlipperOn_B_I->GetBinContent(i) << "\t"<< setprecision(4) << hWFlipperOn_B_I->GetBinError(i) << "\t";
    parbck << setprecision(4) << hEFlipperOn_B_2->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOn_B_2->GetBinError(i) << "\t";
    parbck << setprecision(4) << hWFlipperOn_B_2->GetBinContent(i) << "\t"<< setprecision(4) << hWFlipperOn_B_2->GetBinError(i) << endl;
    
    parbckoff << i << "\t" << setprecision(4) << hEFlipperOff_B->GetBinContent(i) << "\t"<< hEFlipperOff_B->GetBinError(i) << "\t";
    parbckoff << setprecision(4) << hWFlipperOff_B->GetBinContent(i)   << "\t"<< setprecision(4) << hWFlipperOff_B->GetBinError(i) << "\t";
    parbckoff << setprecision(4) << hEFlipperOff_B_I->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOff_B_I->GetBinError(i) << "\t";
    parbckoff << setprecision(4) << hWFlipperOff_B_I->GetBinContent(i) << "\t"<< setprecision(4) << hWFlipperOff_B_I->GetBinError(i) << "\t";
    parbckoff << setprecision(4) << hEFlipperOff_B_2->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOff_B_2->GetBinError(i) << "\t";
    parbckoff << setprecision(4) << hWFlipperOff_B_2->GetBinContent(i) << "\t"<< setprecision(4) << hWFlipperOff_B_2->GetBinError(i) << endl;

    parsig << i << "\t" << setprecision(4) << hEFlipperOn_B->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOn_B->GetBinError(i) << "\t";
    parsig << setprecision(4) << hWFlipperOn->GetBinContent(i)   << "\t"<< setprecision(4) << hWFlipperOn->GetBinError(i)   << "\t";
    parsig << setprecision(4) << hEFlipperOn_I->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOn_I->GetBinError(i) << "\t";
    parsig << setprecision(4) << hWFlipperOn_I->GetBinContent(i) << "\t"<< setprecision(4) << hWFlipperOn_I->GetBinError(i) << "\t";
    parsig << setprecision(4) << hEFlipperOn_2->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOn_2->GetBinError(i) << "\t";
    parsig << setprecision(4) << hWFlipperOn_2->GetBinContent(i) << "\t"<< setprecision(4) << hWFlipperOn_2->GetBinError(i) << endl;
    
    parsigoff << i << "\t" << setprecision(4) << hEFlipperOff->GetBinContent(i) << "\t"<< hEFlipperOff->GetBinError(i) << "\t";
    parsigoff << setprecision(4) << hWFlipperOff->GetBinContent(i)   << "\t"<< setprecision(4) << hWFlipperOff->GetBinError(i) << "\t";
    parsigoff << setprecision(4) << hEFlipperOff_I->GetBinContent(i) << "\t"<< setprecision(4)<< hEFlipperOff_I->GetBinError(i) << "\t";
    parsigoff << setprecision(4) << hWFlipperOff_I->GetBinContent(i) << "\t"<< setprecision(4)<< hWFlipperOff_I->GetBinError(i) << "\t";
    parsigoff << setprecision(4) << hEFlipperOff_2->GetBinContent(i) << "\t"<< setprecision(4)<< hEFlipperOff_2->GetBinError(i) << "\t";
    parsigoff << setprecision(4) << hWFlipperOff_2->GetBinContent(i) << "\t"<< setprecision(4)<< hWFlipperOff_2->GetBinError(i) << endl;
  }
  
  parbckoff.close();
  parbck.close();
  parsigoff.close();
  parsig.close();
  
  TCanvas *c1bk = new TCanvas("c1bk","One Clock back?",700,800);
  c1bk->Divide(3,4);
  TLegend *legSpc1,*legSpc2,*legSpc3,*legSpc4,*legSpc5,*legSpc6,*legSpc7,*legSpc8;
  TLegend *legSpc9,*legSpc10,*legSpc11,*legSpc12;

  DrawPanel(hWFlipperOff  ,btr[3]->hERef , c1bk,1,legSpc1,west_time_off);
  DrawPanel(hWFlipperOff_I,btr[3]->hERef1, c1bk,2,legSpc2,west_time_off);
  DrawPanel(hWFlipperOff_2,btr[3]->hERef2, c1bk,3,legSpc3,west_time_off);
  DrawPanel(hWFlipperOn   ,btr[4]->hERef , c1bk,4,legSpc4,west_time_on);
  DrawPanel(hWFlipperOn_I ,btr[4]->hERef1, c1bk,5,legSpc5,west_time_on);
  DrawPanel(hWFlipperOn_2 ,btr[4]->hERef2, c1bk,6,legSpc6,west_time_on);
  DrawPanel(hEFlipperOff  ,btr[5]->hERef , c1bk,7,legSpc7,east_time_off);
  DrawPanel(hEFlipperOff_I,btr[5]->hERef1, c1bk,8,legSpc8,east_time_off);
  DrawPanel(hEFlipperOff_2,btr[5]->hERef2, c1bk,9,legSpc9,east_time_off);
  DrawPanel(hEFlipperOn   ,btr[6]->hERef , c1bk,10,legSpc10,east_time_on);
  DrawPanel(hEFlipperOn_I ,btr[6]->hERef1, c1bk,11,legSpc11,east_time_on);
  DrawPanel(hEFlipperOn_2 ,btr[6]->hERef2, c1bk,12,legSpc12,east_time_on);
 
  btr[1]->hERef->Scale(hEFlipperOn->Integral(1,80)/btr[1]->hERef->Integral(1,80));
  btr[1]->hERef1->Scale(1./btr[1]->hERef1->Integral(1,80));
  btr[1]->hERef2->Scale(1./btr[1]->hERef2->Integral(1,80));

  c1bk->Print(Form("pdf_out/geo%d/e_spec_compares_%d.pdf",btr[0]->GetGeo(),btr[0]->GetGeo()));
  /*
  TCanvas *cTot = new TCanvas("cTot","Total Rates");
  cTot->Divide(2,2);
  cTot->cd(1);
  btr[1]->hERef->SetLineColor(6);
  hEFlipperOff->SetLineColor(2);
  hWFlipperOn->SetLineColor(3);
  hWFlipperOff->SetLineColor(4);
  
  hWFlipperOff->Draw();
  hEFlipperOff->Draw("same");
  hEFlipperOn->Draw("same");
  hWFlipperOn->Draw("same");
  // btr[1]->hERef->DrawNormalized("same");
  
  cTot->cd(2);
  
  hEFlipperOff_I->SetLineColor(2);
  hWFlipperOn_I->SetLineColor(3);
  hWFlipperOff_I->SetLineColor(4);
  
  hWFlipperOff_I->Draw();
  hEFlipperOff_I->Draw("same");
  hEFlipperOn_I->Draw("same");
  hWFlipperOn_I->Draw("same");
  
  //btr[1]->hERef1->SetLineColor(6);
  //btr[1]->hERef1->DrawNormalized("same");
  
  cTot->cd(3);
  hEFlipperOff_2->SetLineColor(2);
  hWFlipperOn_2->SetLineColor(3);
  hWFlipperOff_2->SetLineColor(4);
  
  hWFlipperOff_2->Draw();
  hEFlipperOff_2->Draw("same");
  hEFlipperOn_2->Draw("same");
  hWFlipperOn_2->Draw("same");
  
  btr[1]->hERef2->SetLineColor(6);
  // btr[1]->hERef2->DrawNormalized("same");
  */
  cout << "*********************************************" << endl;
  cout << " Flipper Off East " << hEFlipperOff->Integral(nlow,nhigh) << endl;
  cout << " Flipper Off West " << hWFlipperOff->Integral(nlow,nhigh) << endl;
  cout << " Flipper On  East " << hEFlipperOn->Integral(nlow,nhigh) << endl;
  cout << " Flipper On  West " << hWFlipperOn->Integral(nlow,nhigh) << endl;
  cout << "**********************************************" << endl;
  
  // Start Calculating the Chi^2/ndf
  
  Double_t chiE0on[80],chiE0off[80],chiE1on[80],chiE1off[80];
  Double_t chiW0on[80],chiW0off[80],chiW1on[80],chiW1off[80];
  Double_t chiE2on[80],chiE2off[80],chiW2on[80],chiW2off[80];
  Double_t x1[80];
  Double_t sig1[80],sig2[80],sig3[80];
  Double_t sig4[80],sig5[80],sig6[80];
  Double_t sig7[80],sig8[80],sig9[80];
  Double_t sig10[80],sig11[80],sig12[80];

  for(Int_t i = 0 ; i < 80 ; i++){
    x1[i]   = hEFlipperOn->GetBinCenter(i+1);
    sig1[i] = 0.;
    sig2[i] = 0.;
    sig3[i] = 0.;
    sig4[i] = 0.;
    sig5[i] = 0.;
    sig6[i] = 0.;
    sig7[i] = 0.;
    sig8[i] = 0.;
    sig9[i] = 0.;
    sig10[i] = 0.;
    sig11[i] = 0.;
    sig12[i] = 0.;

    chiE0on[i] = 0.;
    chiE1on[i] = 0.;
    chiE2on[i] = 0.;
    chiW0on[i] = 0.;
    chiW1on[i] = 0.;
    chiW2on[i] = 0.;
    chiE0off[i] = 0.;
    chiE1off[i] = 0.;
    chiE2off[i] = 0.;
    chiW0off[i] = 0.;
    chiW1off[i] = 0.;
    chiW2off[i] = 0.;
  }

   for(Int_t i = 0 ; i < 80 ; i++){
     if(hEFlipperOn->GetBinContent(i+1) > 0){
       chiE0on[i] = Power((hEFlipperOn->GetBinContent(i+1) -
			   btr[5]->hERef->GetBinContent(i+1))/
			  hEFlipperOn->GetBinContent(i+1),1);
       sig1[i] = sqrt(Power(hEFlipperOn->GetBinError(i+1)/hEFlipperOn->GetBinContent(i+1),2));
     }
     if(hEFlipperOn_I->GetBinContent(i+1) > 0){
       chiE1on[i]  = Power((hEFlipperOn_I->GetBinContent(i+1) -
			    btr[5]->hERef1->GetBinContent(i+1))/
			   hEFlipperOn_I->GetBinContent(i+1),1);
       sig2[i]    = sqrt(Power(hEFlipperOn_I->GetBinError(i+1)/hEFlipperOn_I->GetBinContent(i+1),2));
     }
     if(hEFlipperOn_2->GetBinContent(i+1) > 0){
       chiE2on[i] += Power((hEFlipperOn_2->GetBinContent(i+1) -
			    btr[5]->hERef2->GetBinContent(i+1))/
			   hEFlipperOn_2->GetBinContent(i+1),1);
       sig3[i]    = sqrt(Power(hEFlipperOn_2->GetBinError(i+1)/hEFlipperOn_2->GetBinContent(i+1),2));
     }
     if(hEFlipperOff->GetBinContent(i+1) > 0){
       chiE0off[i] += Power((hEFlipperOff->GetBinContent(i+1) -
			     btr[6]->hERef->GetBinContent(i+1))/
			    hEFlipperOff->GetBinContent(i+1),1);
       sig4[i]    = sqrt(Power(hEFlipperOff->GetBinError(i+1)/hEFlipperOff->GetBinContent(i+1),2));
     }
     if(hEFlipperOff_I->GetBinContent(i+1) > 0){
       chiE1off[i] += Power((hEFlipperOff_I->GetBinContent(i+1) -
			     btr[6]->hERef1->GetBinContent(i+1))/
			    hEFlipperOff_I->GetBinContent(i+1),1);
       sig5[i]    = sqrt(Power(hEFlipperOff_I->GetBinError(i+1)/hEFlipperOff_I->GetBinContent(i+1),2));
     }
     if(hEFlipperOff_2->GetBinContent(i+1) > 0){
       chiE2off[i] += Power((hEFlipperOff_2->GetBinContent(i+1) -
			     btr[6]->hERef2->GetBinContent(i+1))/
			    hEFlipperOff_2->GetBinContent(i+1),1);
       sig6[i]    = sqrt(Power(hEFlipperOff_2->GetBinError(i+1)/hEFlipperOff_2->GetBinContent(i+1),2));
     }
     if(hWFlipperOn->GetBinContent(i+1) > 0){
       chiW0on[i] += Power((hWFlipperOn->GetBinContent(i+1) -
			    btr[4]->hERef->GetBinContent(i+1))/
			   hWFlipperOn->GetBinContent(i+1),1);
       sig7[i]    = sqrt(Power(hWFlipperOn->GetBinError(i+1)/hWFlipperOn->GetBinContent(i+1),2));
     }
     if(hWFlipperOn_I->GetBinContent(i+1) > 0){	
       chiW1on[i] += Power((hWFlipperOn_I->GetBinContent(i+1) -
			    btr[4]->hERef1->GetBinContent(i+1))/
			   hWFlipperOn_I->GetBinContent(i+1),1); 
       sig8[i]    = sqrt(Power(hWFlipperOn_I->GetBinError(i+1)/hWFlipperOn_I->GetBinContent(i+1),2));
     }
     if(hWFlipperOn_2->GetBinContent(i+1) > 0){
       chiW2on[i] += Power((hWFlipperOn_2->GetBinContent(i+1) -
			    btr[4]->hERef2->GetBinContent(i+1))/
			   hWFlipperOn_2->GetBinContent(i+1),1);
       sig9[i]    = sqrt(Power(hWFlipperOn_2->GetBinError(i+1)/hWFlipperOn_2->GetBinContent(i+1),2));
     }
     if(hWFlipperOff->GetBinContent(i+1) > 0){
       chiW0off[i] += Power((hWFlipperOff->GetBinContent(i+1) -
			     btr[3]->hERef->GetBinContent(i+1))/
			    hWFlipperOff->GetBinContent(i+1),1);
       sig10[i]    = sqrt(Power(hWFlipperOff->GetBinError(i+1)/hWFlipperOff->GetBinContent(i+1),2));	      
     }
     if(hWFlipperOff_I->GetBinContent(i+1) > 0){
       chiW1off[i] += Power((hWFlipperOff_I->GetBinContent(i+1) -
			     btr[3]->hERef1->GetBinContent(i+1))/
			    hWFlipperOff_I->GetBinContent(i+1),1);
       sig11[i]    = sqrt(Power(hWFlipperOff_I->GetBinError(i+1)/hWFlipperOff_I->GetBinContent(i+1),2));	      
     }
     if(hWFlipperOff_2->GetBinContent(i+1) > 0){
       chiW2off[i] += Power((hWFlipperOff_2->GetBinContent(i+1) -
			     btr[3]->hERef2->GetBinContent(i+1))/
			    hWFlipperOff_2->GetBinContent(i+1),1);
       sig12[i]    = sqrt(Power(hWFlipperOff_2->GetBinError(i+1)/hWFlipperOff_2->GetBinContent(i+1),2));
     }
  }
  
  TGraphErrors *gEO0 = new TGraphErrors(80,x1,chiE0on,0,sig1);
  ColorGraphic(gEO0,1,20,2);
  TGraphErrors *gEO1 = new TGraphErrors(80,x1,chiE1on,0,sig2);
  ColorGraphic(gEO1,1,20,2);
  TGraphErrors *gEO2 = new TGraphErrors(80,x1,chiE2on,0,sig3);
  ColorGraphic(gEO2,1,20,2);
  
  TGraphErrors *gEF0 = new TGraphErrors(80,x1,chiE0off,0,sig4);
  ColorGraphic(gEF0,1,24,2);
  TGraphErrors *gEF1 = new TGraphErrors(80,x1,chiE1off,0,sig5);
  ColorGraphic(gEF1,1,24,2);
  TGraphErrors *gEF2 = new TGraphErrors(80,x1,chiE2off,0,sig6);
  ColorGraphic(gEF2,1,24,2);
  
  TGraphErrors *gWO0 = new TGraphErrors(80,x1,chiW0on,0,sig7);
  ColorGraphic(gWO0,1,21,2);
  TGraphErrors *gWO1 = new TGraphErrors(80,x1,chiW1on,0,sig8);
  ColorGraphic(gWO1,1,21,2);
  TGraphErrors *gWO2 = new TGraphErrors(80,x1,chiW2on,0,sig9);
  ColorGraphic(gWO2,1,21,2);
  
  TGraphErrors *gWF0 = new TGraphErrors(80,x1,chiW0off,0,sig10);
  ColorGraphic(gWF0,1,25,2);
  TGraphErrors *gWF1 = new TGraphErrors(80,x1,chiW1off,0,sig11);
  ColorGraphic(gWF1,1,25,2);
  TGraphErrors *gWF2 = new TGraphErrors(80,x1,chiW2off,0,sig12);
  ColorGraphic(gWF2,1,25,2);

  vector<Double_t> ave0 (4,0), ave1(4,0),ave2 (4,0), ave0er (4,0);
  vector<Double_t> ave1er (4,0), ave2er (4,0);
  
  Average_Array(chiE0on,sig1,30,ave0[0],ave0er[0]);
  Average_Array(chiE1on,sig2,30,ave1[0],ave1er[0]);
  Average_Array(chiE2on,sig3,30,ave2[0],ave2er[0]);
  Average_Array(chiE0off,sig4,30,ave0[1],ave0er[1]);
  Average_Array(chiE1off,sig5,30,ave1[1],ave1er[1]);
  Average_Array(chiE2off,sig6,30,ave2[1],ave2er[1]); 
  Average_Array(chiW0on,sig7,30,ave0[2],ave0er[2]);
  Average_Array(chiW1on,sig8,30,ave1[2],ave1er[2]);
  Average_Array(chiW2on,sig9,30,ave2[2],ave2er[2]);
  Average_Array(chiW0off,sig10,30,ave0[3],ave0er[3]);
  Average_Array(chiW1off,sig11,30,ave1[3],ave1er[3]);
  Average_Array(chiW2off,sig12,30,ave2[3],ave2er[3]); 
  
  Double_t AveResType0,AveResType1,AveResType2;
  Double_t AveResType0er,AveResType1er,AveResType2er;

  Average_Array(&ave0[0],&ave0er[0],4,AveResType0,AveResType0er);
  Average_Array(&ave1[0],&ave1er[0],4,AveResType1,AveResType1er);
  Average_Array(&ave2[0],&ave2er[0],4,AveResType2,AveResType2er);

  cout << "Average Residuals 0 : " << AveResType0 << " +/- "  << AveResType0er << endl;
  cout << "Average Residuals 1 : " << AveResType1 << " +/- "  << AveResType1er << endl;
  cout << "Average Residuals 2 : " << AveResType2 << " +/- "  << AveResType2er << endl;

   
  TCanvas *cETChi = new TCanvas("cETChi","Total Energy #chi^{2}",350,800);
  cETChi->Divide(1,3);
  
  cETChi->cd(1);
  gEO0->SetTitle("Type 0 Events");
  gEO0->GetXaxis()->SetRangeUser(0,800);
  gEO0->GetYaxis()->SetRangeUser(-1,1);
  gEO0->Draw("E0 AP");
  gEF0->Draw("E0 P");
  gWO0->Draw("E0 P");
  gWF0->Draw("E0 P");
  
  cETChi->cd(2);
  gEF1->SetTitle("Type 1 Events");
  gEF1->GetXaxis()->SetRangeUser(0,800);
  gEF1->GetYaxis()->SetRangeUser(-1,1);
  gEF1->Draw("E0 AP");
  gEO1->Draw("E0 P");
  gWO1->Draw("E0 P");
  gWF1->Draw("E0 P");
  
  cETChi->cd(3);
  gWO2->SetTitle("Type 2/3 Events");
  gWO2->GetXaxis()->SetRangeUser(0,800);
  gWO2->GetYaxis()->SetRangeUser(-1,1);
  gWO2->Draw("AP");
  gWF2->Draw("P");
  gEO2->Draw("P");
  gEF2->Draw("P");

  
  TCanvas *cFinal = new TCanvas("cFinal","Final");
  cFinal->Divide(2,1);

  cFinal->cd(1);
  hEFlipperOn->Draw();
  cFinal->cd(2);
  TH1F *hEFlipperClone = (TH1F*)hEFlipperOn->Clone("hEFlipperClone");
  hEFlipperClone->Draw();
  
  //------------------------------------------------------
  // output the rough super ratio 
  
  Double_t R1p = hEFlipperOn_2->Integral(nlow,nhigh);
  Double_t R1m = hWFlipperOn_2->Integral(nlow,nhigh);
  Double_t R2p = hEFlipperOff_2->Integral(nlow,nhigh);
  Double_t R2m = hWFlipperOff_2->Integral(nlow,nhigh);
  
  Double_t SS = (R1p*R2m)/(R2p*R1m);
  
  Double_t As_quick = (1.-sqrt(SS))/(1.+sqrt(SS));
  cout<<"Simple Super ratio for Type 2/3's is " << As_quick << endl;
  
  
}

void DrawPanel(TH1F *hData,TH1F *hMC,TCanvas *c,Int_t npad,TLegend *lleg,Double_t time)
{
  Int_t nt = (npad % 3)-1;
  if(nt ==-1) nt = 2;
  c->cd(npad);
  hData->GetXaxis()->SetRangeUser(0,1000);
   hMC->Scale(hData->Integral(1,20)/hMC->Integral(1,20));
  //  Double_t chi = hData->Chi2Test(hMC,"WW P");
  vector<Double_t> diff;
  Double_t chi2  = 0;
  Int_t goodbins = 0;
  for(Int_t i = 1 ; i<= 32;i++){
    diff.push_back(TMath::Abs(hData->Integral(1,i) - hMC->Integral(1,i)));
    if(hData->GetBinContent(i)!=0){
      double n1 = hData->GetBinContent(i);
      double n2 = hMC->GetBinContent(i);
      double MM = hMC->Integral(1,32);
      double NN = hData->Integral(1,32);
      double s1 = hData->GetBinError(i);
      double s2 = hMC->GetBinError(i);
      chi2 += TMath::Power((n1-n2)/s1,2);
      goodbins++;
    }
  }
  chi2 *= 1./(double)goodbins;
  Double_t n = hData->Integral(1,32);
  Double_t np = hMC->Integral(1,32);
  Double_t dnnp = *max_element(diff.begin(),diff.end());
  dnnp *= sqrt((n*np)/(n+np));
  hData->Draw("X0 E0 P0"); 
  hMC->Draw("same hist");
  lleg = new TLegend(0.55,0.6,0.89,0.89);
  lleg->SetHeader(Form("Type %d",nt));
  lleg->SetFillColor(0);
  lleg->SetLineColor(0);
  lleg->AddEntry(hData,"Data","pl");
  lleg->AddEntry(hMC,"MC","l");
  lleg->AddEntry(hMC,Form("K-test   : %5.3f",1-dnnp),"");
  lleg->AddEntry(hMC,Form("#chi/#nu : %5.3f",chi2),"");
  lleg->Draw();
};

void SetPlotOptions(TH1* h1,Double_t tsize,Double_t offset,Double_t lsize,Double_t lxsize)
{
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->CenterTitle();
  /*h1->GetYaxis()->SetTitleSize(tsize);
  h1->GetXaxis()->SetTitleSize(tsize);
  h1->GetXaxis()->SetLabelSize(lsize);
  h1->GetYaxis()->SetLabelSize(lxsize);*/
  h1->GetYaxis()->SetTitleOffset(offset);

  
};
void Define_E_Spec()
{
  Int_t nxbins = 80;

  hEFlipperOn = new TH1F("hEFlipperOn","; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOn = new TH1F("hWFlipperOn","; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hEFlipperOff = new TH1F("hEFlipperOff","; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOff = new TH1F("hWFlipperOff","; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
			 
  hEFlipperOn_I = new TH1F("hEFlipperOn_I","; Energy(keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOn_I = new TH1F("hWFlipperOn_I","; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hEFlipperOff_I = new TH1F("hEFlipperOff_I","; Energy (keV) ;"
			    "Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOff_I = new TH1F("hWFlipperOff_I","; Energy (keV) ; Rate"
			 "(Hz)"
			 ,nxbins,0,2000);
			 
  hEFlipperOn_2 = new TH1F("hEFlipperOn_2","; Energy(keV) ; Rate (Hz)"
		         ,nxbins,0,2000);
  hWFlipperOn_2 = new TH1F("hWFlipperOn_2","; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hEFlipperOff_2 = new TH1F("hEFlipperOff_2","; Energy (keV) ;"
			    "Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOff_2 = new TH1F("hWFlipperOff_2","; Energy (keV) ; Rate"
			 "(Hz)"
			 ,nxbins,0,2000);

  // ------------------------------------------------------------------------------
  // Define Background Histograms 
  // - This is calculate the error in each bin based on the signal to background
  //   ratio over the course of a geometry
  // ------------------------------------------------------------------------------
  
  hEFlipperOn_B = new TH1F("hEFlipperOn_B","East On_B ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOn_B = new TH1F("hWFlipperOn_B","West On_B ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hEFlipperOff_B = new TH1F("hEFlipperOff_B","East Off_B ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOff_B = new TH1F("hWFlipperOff_B","West Off_B ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
			 
  hEFlipperOn_B_I = new TH1F("hEFlipperOn_B_I","East On_B ; Energy(keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOn_B_I = new TH1F("hWFlipperOn_B_I","West On_B ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hEFlipperOff_B_I = new TH1F("hEFlipperOff_B_I","East Off_B ; Energy (keV) ;"
			    "Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOff_B_I = new TH1F("hWFlipperOff_B_I","West Off_B ; Energy (keV) ; Rate"
			 "(Hz)"
			 ,nxbins,0,2000);
			 
  hEFlipperOn_B_2 = new TH1F("hEFlipperOn_B_2","East On_B ; Energy(keV) ; Rate (Hz)"
		         ,nxbins,0,2000);
  hWFlipperOn_B_2 = new TH1F("hWFlipperOn_B_2","West On_B ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hEFlipperOff_B_2 = new TH1F("hEFlipperOff_B_2","East Off_B ; Energy (keV) ;"
			    "Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOff_B_2 = new TH1F("hWFlipperOff_B_2","West Off_B ; Energy (keV) ; Rate"
			 "(Hz)"
			 ,nxbins,0,2000);

  //void SetPlotOptions(TH1* h1,Double_t tsize,Double_t offset,Double_t lsize)
  Double_t textsize = 0.07;
  Double_t toffset  = 1.5;
  Double_t tlsize   = 0.04;
  Double_t txlsize  = 0.06;
   
  
  SetPlotOptions(hEFlipperOn   ,textsize,toffset,tlsize,txlsize);
  SetPlotOptions(hEFlipperOn_I ,textsize,toffset,tlsize,txlsize);
  SetPlotOptions(hEFlipperOn_2 ,textsize,toffset,tlsize,txlsize);
  SetPlotOptions(hEFlipperOff  ,textsize,toffset,tlsize,txlsize);
  SetPlotOptions(hEFlipperOff_I,textsize,toffset,tlsize,txlsize);
  SetPlotOptions(hEFlipperOff_2,textsize,toffset,tlsize,txlsize);
  SetPlotOptions(hWFlipperOn   ,textsize,toffset,tlsize,txlsize);
  SetPlotOptions(hWFlipperOn_I ,textsize,toffset,tlsize,txlsize);
  SetPlotOptions(hWFlipperOn_2 ,textsize,toffset,tlsize,txlsize);
  SetPlotOptions(hWFlipperOff  ,textsize,toffset,tlsize,txlsize);
  SetPlotOptions(hWFlipperOff_I,textsize,toffset,tlsize,txlsize);
  SetPlotOptions(hWFlipperOff_2,textsize,toffset,tlsize,txlsize);
  
  Int_t nmark = 24;
  Double_t markerstyle =0.6;
  ColorGraphic(hEFlipperOn_B,1,1,1,nmark);
  ColorGraphic(hEFlipperOn_B_I,1,1,1,nmark);
  ColorGraphic(hEFlipperOn_B_2,1,1,1,nmark);
  ColorGraphic(hEFlipperOff_B,1,1,1,nmark);
  ColorGraphic(hEFlipperOff_B_I,1,1,1,nmark);
  ColorGraphic(hEFlipperOff_B_2,1,1,1,nmark);
  ColorGraphic(hWFlipperOn_B,1,1,1,nmark);
  ColorGraphic(hWFlipperOn_B_I,1,1,1,nmark);
  ColorGraphic(hWFlipperOn_B_2,1,1,1,nmark);
  ColorGraphic(hWFlipperOff_B,1,1,1,nmark);
  ColorGraphic(hWFlipperOff_B_I,1,1,1,nmark);
  ColorGraphic(hWFlipperOff_B_2,1,1,1,nmark);

  ColorGraphic(hEFlipperOn,1,1,1,nmark, markerstyle );
  ColorGraphic(hEFlipperOn_I,1,1,1,nmark, markerstyle );
  ColorGraphic(hEFlipperOn_2,1,1,1,nmark, markerstyle );
  ColorGraphic(hEFlipperOff,1,1,1,nmark, markerstyle );
  ColorGraphic(hEFlipperOff_I,1,1,1,nmark, markerstyle );
  ColorGraphic(hEFlipperOff_2,1,1,1,nmark, markerstyle );
  ColorGraphic(hWFlipperOn,1,1,1,nmark, markerstyle );
  ColorGraphic(hWFlipperOn_I,1,1,1,nmark, markerstyle );
  ColorGraphic(hWFlipperOn_2,1,1,1,nmark, markerstyle );
  ColorGraphic(hWFlipperOff,1,1,1,nmark, markerstyle );
  ColorGraphic(hWFlipperOff_I,1,1,1,nmark, markerstyle );
  ColorGraphic(hWFlipperOff_2,1,1,1,nmark, markerstyle );
}
//---------------------------------------------------------------------------
void Fill_Foreground(Double_t *t1,Double_t *t2, Double_t *t3, Double_t *t4)
{

  using namespace TMath;
  
  Double_t t1val = *t1;
  Double_t t2val = *t2;
  Double_t t3val = *t3;
  Double_t t4val = *t4;

  Double_t EastOn[100],EastOn_I[100],EastOn_2[100];
  Double_t EastOff[100],EastOff_I[100],EastOff_2[100];
  Double_t WestOn[100],WestOn_I[100],WestOn_2[100];
  Double_t WestOff[100],WestOff_I[100],WestOff_2[100];
  
  for(Int_t i = 0; i < 100 ; i++){
    EastOn[i]    = 0.;
    EastOn_I[i]  = 0.;
    EastOn_2[i]  = 0.;
    EastOff[i]   = 0.;
    EastOff_I[i] = 0.;
    EastOff_2[i] = 0.;
    WestOn[i]    = 0.;
    WestOn_I[i]  = 0.;
    WestOn_2[i]  = 0.;
    WestOff[i]   = 0.;
    WestOff_I[i] = 0.;
    WestOff_2[i] = 0.;
  }
  
  t1val = 0.;
  t2val = 0.;
  t3val = 0.;
  t4val = 0.;
  
  // Loop over all beta runs and sum up the background subtracted spectra
  fstream f;

  for(Int_t i = 0 ; i < nbeta ; i++){
    for(Int_t bin = 1 ; bin <= hEFlipperOn->GetNbinsX(); bin++){
      if(btr[i]->flipperOn == 1){
	if(btr[i]->heq->GetBinError(bin)>0){
	  EastOn[bin-1] += 1./Power(btr[i]->heq->GetBinError(bin),2);
	  hEFlipperOn->SetBinContent(bin,hEFlipperOn->GetBinContent(bin) +
				     btr[i]->heq->GetBinContent(bin)/Power(btr[i]->heq->GetBinError(bin),2));
	}
	if(btr[i]->hwq->GetBinError(bin)>0){ 
	  WestOn[bin-1] += 1./Power(btr[i]->hwq->GetBinError(bin),2);
	  hWFlipperOn->SetBinContent(bin,hWFlipperOn->GetBinContent(bin)+
				     btr[i]->hwq->GetBinContent(bin)/Power(btr[i]->hwq->GetBinError(bin),2));
	}
	if(btr[i]->heqF->GetBinError(bin)>0){
	  EastOn_I[bin-1] += 1./Power(btr[i]->heqF->GetBinError(bin),2);
	  hEFlipperOn_I->SetBinContent(bin,hEFlipperOn_I->GetBinContent(bin)+
				       btr[i]->heqF->GetBinContent(bin)/Power(btr[i]->heqF->GetBinError(bin),2));//  *btr[i]->rtime_e);
	}
	if(btr[i]->hwqF->GetBinError(bin)>0){
	  WestOn_I[bin-1] += 1./Power(btr[i]->hwqF->GetBinError(bin),2);
	  hWFlipperOn_I->SetBinContent(bin,hWFlipperOn_I->GetBinContent(bin)+
				       btr[i]->hwqF->GetBinContent(bin)/Power(btr[i]->hwqF->GetBinError(bin),2)); //*btr[i]->rtime_w);
	}	
	if(btr[i]->heqG->GetBinError(bin)>0){
     	  EastOn_2[bin-1] += 1./Power(btr[i]->heqG->GetBinError(bin),2);
	  hEFlipperOn_2->SetBinContent(bin,hEFlipperOn_2->GetBinContent(bin)+
				       btr[i]->heqG->GetBinContent(bin)/Power(btr[i]->heqG->GetBinError(bin),2)); //*btr[i]->rtime_e);
	}
	
	if(btr[i]->hwqG->GetBinError(bin)>0){
	  WestOn_2[bin-1] += 1./Power(btr[i]->hwqG->GetBinError(bin),2);
	  hWFlipperOn_2->SetBinContent(bin,hWFlipperOn_2->GetBinContent(bin)+
				       btr[i]->hwqG->GetBinContent(bin)/Power(btr[i]->hwqG->GetBinError(bin),2)); //*btr[i]->rtime_w);
	}
	// Incriment the run time
	if(bin == 1){
	  t2val += btr[i]->rtime_e;
	  t1val += btr[i]->rtime_w;
	}
      } else if(btr[i]->flipperOn == 0){
	if(btr[i]->heq->GetBinError(bin)>0){
		  
	  EastOff[bin-1] += 1./Power(btr[i]->heq->GetBinError(bin),2);
	  hEFlipperOff->SetBinContent(bin,hEFlipperOff->GetBinContent(bin)+
				      btr[i]->heq->GetBinContent(bin)/ Power(btr[i]->heq->GetBinError(bin),2));//*btr[i]->rtime_e);
	}
	if(btr[i]->hwq->GetBinError(bin)>0){
	  
	  WestOff[bin-1] += 1./Power(btr[i]->hwq->GetBinError(bin),2);    
	  hWFlipperOff->SetBinContent(bin,hWFlipperOff->GetBinContent(bin)+
				      btr[i]->hwq->GetBinContent(bin)/Power(btr[i]->hwq->GetBinError(bin),2));   // *btr[i]->rtime_w);
	}
	if(btr[i]->heqF->GetBinError(bin)>0){

	  EastOff_I[bin-1] += 1./Power(btr[i]->heqF->GetBinError(bin),2);
	  hEFlipperOff_I->SetBinContent(bin,hEFlipperOff_I->GetBinContent(bin)+
					btr[i]->heqF->GetBinContent(bin)/Power(btr[i]->heqF->GetBinError(bin),2));  //*btr[i]->rtime_e);
	}
	if(btr[i]->hwqF->GetBinError(bin)>0){

	  WestOff_I[bin-1] += 1./Power(btr[i]->hwqF->GetBinError(bin),2);
	  hWFlipperOff_I->SetBinContent(bin,hWFlipperOff_I->GetBinContent(bin)+
					btr[i]->hwqF->GetBinContent(bin)/Power(btr[i]->hwqF->GetBinError(bin),2));  // *btr[i]->rtime_w);
	}
	if(btr[i]->heqG->GetBinError(bin)>0){

	  EastOff_2[bin-1] += 1./Power(btr[i]->heqG->GetBinError(bin),2);
	  hEFlipperOff_2->SetBinContent(bin,hEFlipperOff_2->GetBinContent(bin)+
					btr[i]->heqG->GetBinContent(bin)/Power(btr[i]->heqG->GetBinError(bin),2));  // *btr[i]->rtime_e);
	}
	if(btr[i]->hwqG->GetBinError(bin)>0){
	  
	  WestOff_2[bin-1] += 1./Power(btr[i]->hwqG->GetBinError(bin),2);
	  hWFlipperOff_2->SetBinContent(bin,hWFlipperOff_2->GetBinContent(bin)+
				      btr[i]->hwqG->GetBinContent(bin)/Power(btr[i]->hwqG->GetBinError(bin),2));  //*btr[i]->rtime_w);
	}
	// Incriment the run time
	if(bin == 1){
	  t4val += btr[i]->rtime_e;
	  t3val += btr[i]->rtime_w;
	}
      }
    }
  }

  for(Int_t i = 1 ; i <= btr[0]->heq->GetNbinsX() ; i++){
    if(EastOn[i-1] > 0){
      hEFlipperOn->SetBinContent(i,   hEFlipperOn->GetBinContent(i)   / EastOn[i-1]);
      hEFlipperOn->SetBinError(i, 1./Sqrt(EastOn[i-1]));
    }
    if(EastOn_I[i-1] > 0){
      hEFlipperOn_I->SetBinContent(i, hEFlipperOn_I->GetBinContent(i) / EastOn_I[i-1]);
      hEFlipperOn_I->SetBinError(i, 1./Sqrt(EastOn_I[i-1]));	  
    }
    if(EastOn_2[i-1] > 0){
      hEFlipperOn_2->SetBinContent(i, hEFlipperOn_2->GetBinContent(i) / EastOn_2[i-1]);
      hEFlipperOn_2->SetBinError(i, 1./Sqrt(EastOn_2[i-1]));
    }
    if(EastOff[i-1] > 0){
      hEFlipperOff->SetBinContent(i,   hEFlipperOff->GetBinContent(i)   / EastOff[i-1]);
      hEFlipperOff->SetBinError(i, 1./Sqrt(EastOff[i-1]));
    }
    if(EastOff_I[i-1] > 0){
      hEFlipperOff_I->SetBinContent(i, hEFlipperOff_I->GetBinContent(i) / EastOff_I[i-1]);
      hEFlipperOff_I->SetBinError(i, 1./Sqrt(EastOff_I[i-1]));	  
    }
    if(EastOff_2[i-1] > 0){
      hEFlipperOff_2->SetBinContent(i, hEFlipperOff_2->GetBinContent(i) / EastOff_2[i-1]);
      hEFlipperOff_2->SetBinError(i, 1./Sqrt(EastOff_2[i-1]));
    }
    if(WestOn[i-1] > 0){
      hWFlipperOn->SetBinContent(i,   hWFlipperOn->GetBinContent(i)   / WestOn[i-1]);
      hWFlipperOn->SetBinError(i, 1./Sqrt(WestOn[i-1]));
    }
    if(WestOn_I[i-1] > 0){
      hWFlipperOn_I->SetBinContent(i, hWFlipperOn_I->GetBinContent(i) / WestOn_I[i-1]);
      hWFlipperOn_I->SetBinError(i, 1./Sqrt(WestOn_I[i-1]));
    }
    if(WestOn_2[i-1] > 0){
      hWFlipperOn_2->SetBinContent(i, hWFlipperOn_2->GetBinContent(i) / WestOn_2[i-1]);
      hWFlipperOn_2->SetBinError(i, 1./Sqrt(WestOn_2[i-1]));
    }
    if(WestOff[i-1] > 0){
      hWFlipperOff->SetBinContent(i,   hWFlipperOff->GetBinContent(i)   / WestOff[i-1]);
      hWFlipperOff->SetBinError(i, 1./Sqrt(WestOff[i-1]));
    }
    if(WestOff_I[i-1] > 0){
      hWFlipperOff_I->SetBinContent(i, hWFlipperOff_I->GetBinContent(i) / WestOff_I[i-1]);
      hWFlipperOff_I->SetBinError(i, 1./Sqrt(WestOff_I[i-1]));
    }
    if(WestOff_2[i-1] > 0){
      hWFlipperOff_2->SetBinContent(i, hWFlipperOff_2->GetBinContent(i) / WestOff_2[i-1]);
      hWFlipperOff_2->SetBinError(i, 1./Sqrt(WestOff_2[i-1]));
    }
   
  }

  *t1 = t1val;
  *t2 = t2val;
  *t3 = t3val;
  *t4 = t4val;
 
}

//---------------------------------------------------------------------------
void Fill_Background(Double_t* t1,Double_t* t2, Double_t* t3, Double_t* t4)
{

  using namespace TMath;
  
  Double_t t1val = *t1;
  Double_t t2val = *t2;
  Double_t t3val = *t3;
  Double_t t4val = *t4;

  Double_t EastOn[100],EastOn_I[100],EastOn_2[100];
  Double_t EastOff[100],EastOff_I[100],EastOff_2[100];
  Double_t WestOn[100],WestOn_I[100],WestOn_2[100];
  Double_t WestOff[100],WestOff_I[100],WestOff_2[100];
  
  for(Int_t i = 0; i < 100 ; i++){
    EastOn[i]    = 0.;
    EastOn_I[i]  = 0.;
    EastOn_2[i]  = 0.;
    EastOff[i]   = 0.;
    EastOff_I[i] = 0.;
    EastOff_2[i] = 0.;
    WestOn[i]    = 0.;
    WestOn_I[i]  = 0.;
    WestOn_2[i]  = 0.;
    WestOff[i]   = 0.;
    WestOff_I[i] = 0.;
    WestOff_2[i] = 0.;
  }
  
  t1val = 0.;
  t2val = 0.;
  t3val = 0.;
  t4val = 0.;
  Int_t last=0;
  // Loop over all beta runs and sum up the background subtracted spectra

  for(Int_t i = 0 ; i < nbck ; i++){
    if(bckr[i]->GetRunNumber() != last){
    for(Int_t bin = 1 ; bin <= hEFlipperOn_B->GetNbinsX(); bin++){
      if(bckr[i]->flipperOn == 1){
	if(bckr[i]->heq->GetBinError(bin)>0){
	  EastOn[bin-1] += 1./Power(bckr[i]->heq->GetBinError(bin),2);
	  hEFlipperOn_B->SetBinContent(bin,hEFlipperOn_B->GetBinContent(bin)+
				    bckr[i]->heq->GetBinContent(bin)/ Power(bckr[i]->heq->GetBinError(bin),2));//*bckr[i]->rtime_e);
	}
	if(bckr[i]->hwq->GetBinError(bin)>0){
	  WestOn[bin-1] += 1./Power(bckr[i]->hwq->GetBinError(bin),2);    
	  hWFlipperOn_B->SetBinContent(bin,hWFlipperOn_B->GetBinContent(bin)+
				       bckr[i]->hwq->GetBinContent(bin)/Power(bckr[i]->hwq->GetBinError(bin),2));   // *bckr[i]->rtime_w);
	}
	if(bckr[i]->heqF->GetBinError(bin)>0){
	  
	EastOn_I[bin-1] += 1./Power(bckr[i]->heqF->GetBinError(bin),2);
	hEFlipperOn_B_I->SetBinContent(bin,hEFlipperOn_B_I->GetBinContent(bin)+
				      bckr[i]->heqF->GetBinContent(bin)/Power(bckr[i]->heqF->GetBinError(bin),2));  //*bckr[i]->rtime_e);
	}
	if(bckr[i]->hwqF->GetBinError(bin)>0){
		  
	  WestOn_I[bin-1] += 1./Power(bckr[i]->hwqF->GetBinError(bin),2);
	  hWFlipperOn_B_I->SetBinContent(bin,hWFlipperOn_B_I->GetBinContent(bin)+
					 bckr[i]->hwqF->GetBinContent(bin)/Power(bckr[i]->hwqF->GetBinError(bin),2));  // *bckr[i]->rtime_w);
	}
	if(bckr[i]->heqG->GetBinError(bin)>0){
	  
	  EastOn_2[bin-1] += 1./Power(bckr[i]->heqG->GetBinError(bin),2);
	  hEFlipperOn_B_2->SetBinContent(bin,hEFlipperOn_B_2->GetBinContent(bin)+
					 bckr[i]->heqG->GetBinContent(bin)/Power(bckr[i]->heqG->GetBinError(bin),2));  // *bckr[i]->rtime_e);
	}
	if(bckr[i]->hwqG->GetBinError(bin)>0){
	  WestOn_2[bin-1] += 1./Power(bckr[i]->hwqG->GetBinError(bin),2);
	  hWFlipperOn_B_2->SetBinContent(bin,hWFlipperOn_B_2->GetBinContent(bin) +
					 bckr[i]->hwqG->GetBinContent(bin)/Power(bckr[i]->hwqG->GetBinError(bin),2));
	}
	// Incriment the run time
	if(bin == 1){
	  t2val += bckr[i]->rtime_e;
	  t1val += bckr[i]->rtime_w;
	}
	
      } else {
	if(bckr[i]->heq->GetBinError(bin)>0){
	  EastOff[bin-1] += 1./Power(bckr[i]->heq->GetBinError(bin),2);
	  hEFlipperOff_B->SetBinContent(bin,hEFlipperOff_B->GetBinContent(bin)+
					bckr[i]->heq->GetBinContent(bin)/ Power(bckr[i]->heq->GetBinError(bin),2));//*bckr[i]->rtime_e);
	}
	if(bckr[i]->hwq->GetBinError(bin)>0){

	  WestOff[bin-1] += 1./Power(bckr[i]->hwq->GetBinError(bin),2);    
	  hWFlipperOff_B->SetBinContent(bin,hWFlipperOff_B->GetBinContent(bin)+
					bckr[i]->hwq->GetBinContent(bin)/Power(bckr[i]->hwq->GetBinError(bin),2));   // *bckr[i]->rtime_w);
	}		    
	if(bckr[i]->heqF->GetBinError(bin)>0){
	  EastOff_I[bin-1] += 1./Power(bckr[i]->heqF->GetBinError(bin),2);
	  hEFlipperOff_B_I->SetBinContent(bin,hEFlipperOff_B_I->GetBinContent(bin)+
					  bckr[i]->heqF->GetBinContent(bin)/Power(bckr[i]->heqF->GetBinError(bin),2));  //*bckr[i]->rtime_e);
	}	
	if(bckr[i]->hwqF->GetBinError(bin)>0){

	  WestOff_I[bin-1] += 1./Power(bckr[i]->hwqF->GetBinError(bin),2);
	  hWFlipperOff_B_I->SetBinContent(bin,hWFlipperOff_B_I->GetBinContent(bin)+
					  bckr[i]->hwqF->GetBinContent(bin)/Power(bckr[i]->hwqF->GetBinError(bin),2));  // *bckr[i]->rtime_w);
	}
	if(bckr[i]->heqG->GetBinError(bin)>0){
	  
	  EastOff_2[bin-1] += 1./Power(bckr[i]->heqG->GetBinError(bin),2);
	  hEFlipperOff_B_2->SetBinContent(bin,hEFlipperOff_B_2->GetBinContent(bin) +
					  bckr[i]->heqG->GetBinContent(bin)/Power(bckr[i]->heqG->GetBinError(bin),2));  // *bckr[i]->rtime_e);
	}
	if(bckr[i]->hwqG->GetBinError(bin)>0){

	  WestOff_2[bin-1] += 1./Power(bckr[i]->hwqG->GetBinError(bin),2);
	  hWFlipperOff_B_2->SetBinContent(bin,hWFlipperOff_B_2->GetBinContent(bin) +
					  bckr[i]->hwqG->GetBinContent(bin)/Power(bckr[i]->hwqG->GetBinError(bin),2));
	}
	
	// Incriment the run time
	if(bin == 1) {
	  t4val += bckr[i]->rtime_e;
	  t3val += bckr[i]->rtime_w;
	}
      }
       last = bckr[i]->GetRunNumber();
      }
    }
  }
 
   for(Int_t i = 1 ; i <= btr[0]->heq->GetNbinsX() ; i++){
     if(EastOn[i-1] > 0){
      hEFlipperOn_B->SetBinContent(i,   hEFlipperOn_B->GetBinContent(i)   / EastOn[i-1]);
      hEFlipperOn_B->SetBinError(i, 1./Sqrt(EastOn[i-1]));
    }
    if(EastOn_I[i-1] > 0){
      hEFlipperOn_B_I->SetBinContent(i, hEFlipperOn_B_I->GetBinContent(i) / EastOn_I[i-1]);
      hEFlipperOn_B_I->SetBinError(i, 1./Sqrt(EastOn_I[i-1]));	  
    }
    if(EastOn_2[i-1] > 0){
      hEFlipperOn_B_2->SetBinContent(i, hEFlipperOn_B_2->GetBinContent(i) / EastOn_2[i-1]);
      hEFlipperOn_B_2->SetBinError(i, 1./Sqrt(EastOn_2[i-1]));
    }
    if(EastOff[i-1] > 0){
      hEFlipperOff_B->SetBinContent(i,   hEFlipperOff_B->GetBinContent(i)   / EastOff[i-1]);
      hEFlipperOff_B->SetBinError(i, 1./Sqrt(EastOff[i-1]));
    }
    if(EastOff_I[i-1] > 0){
      hEFlipperOff_B_I->SetBinContent(i, hEFlipperOff_B_I->GetBinContent(i) / EastOff_I[i-1]);
      hEFlipperOff_B_I->SetBinError(i, 1./Sqrt(EastOff_I[i-1]));	  
    }
    if(EastOff_2[i-1] > 0){
      hEFlipperOff_B_2->SetBinContent(i, hEFlipperOff_B_2->GetBinContent(i) / EastOff_2[i-1]);
      hEFlipperOff_B_2->SetBinError(i, 1./Sqrt(EastOff_2[i-1]));
    }
    if(WestOn[i-1] > 0){
      hWFlipperOn_B->SetBinContent(i,   hWFlipperOn_B->GetBinContent(i)   / WestOn[i-1]);
      hWFlipperOn_B->SetBinError(i, 1./Sqrt(WestOn[i-1]));
    }
    if(WestOn_I[i-1] > 0){
      hWFlipperOn_B_I->SetBinContent(i, hWFlipperOn_B_I->GetBinContent(i) / WestOn_I[i-1]);
      hWFlipperOn_B_I->SetBinError(i, 1./Sqrt(WestOn_I[i-1]));
    }
    if(WestOn_2[i-1] > 0){
      hWFlipperOn_B_2->SetBinContent(i, hWFlipperOn_B_2->GetBinContent(i) / WestOn_2[i-1]);
      hWFlipperOn_B_2->SetBinError(i, 1./Sqrt(WestOn_2[i-1]));
    }
    if(WestOff[i-1] > 0){
      hWFlipperOff_B->SetBinContent(i,   hWFlipperOff_B->GetBinContent(i)   / WestOff[i-1]);
      hWFlipperOff_B->SetBinError(i, 1./Sqrt(WestOff[i-1]));
    }
    if(WestOff_I[i-1] > 0){
      hWFlipperOff_B_I->SetBinContent(i, hWFlipperOff_B_I->GetBinContent(i) / WestOff_I[i-1]);
      hWFlipperOff_B_I->SetBinError(i, 1./Sqrt(WestOff_I[i-1]));
    }
    if(WestOff_2[i-1] > 0){
      hWFlipperOff_B_2->SetBinContent(i, hWFlipperOff_B_2->GetBinContent(i) / WestOff_2[i-1]);
      hWFlipperOff_B_2->SetBinError(i, 1./Sqrt(WestOff_2[i-1]));
    }
   }
 

  *t1 = t1val;
  *t2 = t2val;
  *t3 = t3val;
  *t4 = t4val;
}

void Set_Error_Bars(Double_t t1,Double_t t2, Double_t t3, Double_t t4,
		     Double_t tb1,Double_t tb2, Double_t tb3, Double_t tb4)
{

  for(Int_t i = 1 ; i <= hWFlipperOff->GetNbinsX() ; i++){
    hEFlipperOn->SetBinError(i,Bin_Errors(hEFlipperOn->GetBinContent(i),t2,
					  hEFlipperOn_B->GetBinContent(i),tb2));
    hEFlipperOff->SetBinError(i,Bin_Errors(hEFlipperOff->GetBinContent(i),t4,
					   hEFlipperOff_B->GetBinContent(i),tb4));

    hWFlipperOn->SetBinError(i,Bin_Errors(hWFlipperOn->GetBinContent(i),t1,
					  hWFlipperOn_B->GetBinContent(i),tb1));
    hWFlipperOff->SetBinError(i,Bin_Errors(hWFlipperOff->GetBinContent(i),t3,
					   hWFlipperOff_B->GetBinContent(i),tb3));
    //- Repeat of type 1 events
    hEFlipperOn_I->SetBinError(i,Bin_Errors(hEFlipperOn_I->GetBinContent(i),t2,
					  hEFlipperOn_B_I->GetBinContent(i),tb2));
    hEFlipperOff_I->SetBinError(i,Bin_Errors(hEFlipperOff_I->GetBinContent(i),t4,
					   hEFlipperOff_B_I->GetBinContent(i),tb4));

    hWFlipperOn_I->SetBinError(i,Bin_Errors(hWFlipperOn_I->GetBinContent(i),t1,
					  hWFlipperOn_B_I->GetBinContent(i),tb1));
    hWFlipperOff_I->SetBinError(i,Bin_Errors(hWFlipperOff_I->GetBinContent(i),t3,
					   hWFlipperOff_B_I->GetBinContent(i),tb3));

    //- Repeat of type 1 events
    hEFlipperOn_2->SetBinError(i,Bin_Errors(hEFlipperOn_2->GetBinContent(i),t2,
					  hEFlipperOn_B_2->GetBinContent(i),tb2));
    hEFlipperOff_2->SetBinError(i,Bin_Errors(hEFlipperOff_2->GetBinContent(i),t4,
					   hEFlipperOff_B_2->GetBinContent(i),tb4));

    hWFlipperOn_2->SetBinError(i,Bin_Errors(hWFlipperOn_2->GetBinContent(i),t1,
					  hWFlipperOn_B_2->GetBinContent(i),tb1));
    hWFlipperOff_2->SetBinError(i,Bin_Errors(hWFlipperOff_2->GetBinContent(i),t3,
					   hWFlipperOff_B_2->GetBinContent(i),tb3));

    
    
  }
}


Double_t Bin_Errors(Double_t r1,Double_t t1,Double_t r2,Double_t t2)
{
  using namespace TMath;

  Double_t er = 0.;
  if(IsNaN(r1) != 0 ) r1 = 0.;
  
  er = Sqrt(Abs( (r1 + r2) + r2));
  
  if(IsNaN(er) != 0){
    cout << "Error calculation is bad \t " << r1 << "\t" << r2 <<endl;
    er = 0.;
  }
  return er;
}
//----------------------------------------------------------------------------------
void GammaBack()
{
  
  using namespace TMath;
  
  Int_t nxbins = 80;
  Int_t fon=0;
  Int_t foff=0;
  Double_t aveEO[nxbins],aveWO[nxbins],aveE1[nxbins],aveW1[nxbins];
  Double_t aveEOe[nxbins],aveWOe[nxbins],aveE1e[nxbins],aveW1e[nxbins];
  Double_t avesEO[nxbins],avesWO[nxbins],avesE1[nxbins],avesW1[nxbins];
  Double_t avesEOe[nxbins],avesWOe[nxbins],avesE1e[nxbins],avesW1e[nxbins];
  Double_t x[1000],REON[1000],REOFF[1000],RWON[1000],RWOFF[1000];
  Double_t xe[1000],REONe[1000],REOFFe[1000],RWONe[1000],RWOFFe[1000];
  Double_t xoff[1000];
  Double_t xoffe[1000];
  Double_t reonr,reonre,reofr,reofre,rwonr,rwonre,rwofr,rwofre;

  for(Int_t i = 0; i < nxbins ;i++){
    aveEO[i] = 0.;
    aveE1[i] = 0.;
    aveWO[i] = 0.;
    aveW1[i] = 0.;
    aveEOe[i] = 0.;
    aveE1e[i] = 0.;
    aveWOe[i] = 0.;
    aveW1e[i] = 0.;

    avesEO[i] = 0.;
    avesE1[i] = 0.;
    avesWO[i] = 0.;
    avesW1[i] = 0.;
    avesEOe[i] = 0.;
    avesE1e[i] = 0.;
    avesWOe[i] = 0.;
    avesW1e[i] = 0.; 
    reonr  = 0.;
    reonre = 0.;
    reofr  = 0.;
    reofre = 0.;
    rwonr  = 0.;
    rwonre = 0.;
    rwofr  = 0.;
    rwofre = 0.; 
  }
  
  hEFlipperOn_G = new TH1F("hEFlipperOn_G","East On ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOn_G = new TH1F("hWFlipperOn_G","West On ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hEFlipperOff_G = new TH1F("hEFlipperOff_G","East Off ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOff_G = new TH1F("hWFlipperOff_G","West Off ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);

  hEFlipperOn_S = new TH1F("hEFlipperOn_S","East On ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOn_S = new TH1F("hWFlipperOn_S","West On ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hEFlipperOff_S = new TH1F("hEFlipperOff_S","East Off ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOff_S = new TH1F("hWFlipperOff_S","West Off ; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);

  for(Int_t i = 1 ; i <= hEFlipperOn_G->GetNbinsX() ; i++){
    hEFlipperOn_G->SetBinContent(i,0.);
    hEFlipperOff_G->SetBinContent(i,0.);
    hWFlipperOn_G->SetBinContent(i,0.);
    hWFlipperOff_G->SetBinContent(i,0.);
    hEFlipperOn_S->SetBinContent(i,0.);
    hEFlipperOff_S->SetBinContent(i,0.);
    hWFlipperOn_S->SetBinContent(i,0.);
    hWFlipperOff_S->SetBinContent(i,0.);
  }
    
  for(Int_t i = 0 ; i < nbeta ; i++){
    for(Int_t j = 1 ; j <= hEFlipperOn_G->GetNbinsX(); j++){
    //  if(Abs(btr[i]->hENoMWPC->Integral(1,5)) < .5 && Abs(btr[i]->hWNoMWPC->Integral(1,5)) < .5 ){
      if(btr[i]->flipperOn == 1){
	
	if(btr[i]->hENoMWPC->GetBinError(j) > 0 && Abs(btr[i]->MWPC_RatioE) < 0.2){
	  
	  aveEO[j-1] += btr[i]->hENoMWPC->GetBinContent(j)
	    /Power(btr[i]->hENoMWPC->GetBinError(j),2);
	      aveEOe[j-1] += 1./Power(btr[i]->hENoMWPC->GetBinError(j),2);
	}
	
	if(btr[i]->hEMWPC->GetBinError(j) > 0 && Abs(btr[i]->MWPC_RatioE) < 0.2){ 
	  avesEO[j-1] += btr[i]->hEMWPC->GetBinContent(j)
	    /Power(btr[i]->hEMWPC->GetBinError(j),2);
	  avesEOe[j-1] += 1./Power(btr[i]->hEMWPC->GetBinError(j),2);
	  
	}
	
	if(btr[i]->hWNoMWPC->GetBinError(j) > 0&& Abs(btr[i]->MWPC_RatioW) < 0.2){

	  aveWO[j-1] += btr[i]->hWNoMWPC->GetBinContent(j)
	    /Power(btr[i]->hWNoMWPC->GetBinError(j),2);
	  
	  aveWOe[j-1] += 1./Power(btr[i]->hWNoMWPC->GetBinError(j),2);
	}
	if(btr[i]->hWMWPC->GetBinError(j) > 0&& Abs(btr[i]->MWPC_RatioW) < 0.2){ 
	  avesWO[j-1] += btr[i]->hWMWPC->GetBinContent(j)
	    /Power(btr[i]->hWMWPC->GetBinError(j),2);
	  
	  avesWOe[j-1] += 1./Power(btr[i]->hWMWPC->GetBinError(j),2);
	}
	
      }
      if(btr[i]->flipperOn == 0){
	
	if(btr[i]->hENoMWPC->GetBinError(j) > 0 && Abs(btr[i]->MWPC_RatioE) < 0.2){
	  
	  aveE1[j-1] += btr[i]->hENoMWPC->GetBinContent(j)
	    / Power(btr[i]->hENoMWPC->GetBinError(j),2);
	  
	  aveE1e[j-1] += 1./Power(btr[i]->hENoMWPC->GetBinError(j),2);
	}
	if(btr[i]->hEMWPC->GetBinError(j) > 0 && Abs(btr[i]->MWPC_RatioE) < 0.2){ 
	  avesE1[j-1] += btr[i]->hEMWPC->GetBinContent(j)
	    / Power(btr[i]->hEMWPC->GetBinError(j),2);
	  
	  avesE1e[j-1] += 1./Power(btr[i]->hEMWPC->GetBinError(j),2);
	}
	if(btr[i]->hWNoMWPC->GetBinError(j) > 0 && Abs(btr[i]->MWPC_RatioW) < 0.2){
	  
	  aveW1[j-1]  += btr[i]->hWNoMWPC->GetBinContent(j)
	    / Power(btr[i]->hWNoMWPC->GetBinError(j),2);
	  
	  aveW1e[j-1] += 1./Power(btr[i]->hWNoMWPC->GetBinError(j),2);
	}
	if(btr[i]->hWMWPC->GetBinError(j) > 0 && Abs(btr[i]->MWPC_RatioW) < 0.2){ 
	
	  avesW1[j-1]  += btr[i]->hWMWPC->GetBinContent(j)
	    / Power(btr[i]->hWMWPC->GetBinError(j),2);
	    
	  avesW1e[j-1] += 1./Power(btr[i]->hWMWPC->GetBinError(j),2);
	  
	}
      }
   // }
    }
    //-----------------------------------------------------------------
    // Look at individual runs
    //-----------------------------------------------------------------

       
    if(btr[i]->flipperOn == 1 && Abs(btr[i]->MWPC_RatioE) < 0.2 && Abs(btr[i]->MWPC_RatioW) < 0.2){
       
      x[fon]  = btr[i]->GetRunNumber();
      xe[fon] = 0.01;
      
 //     cout << "Run " << x[fon] << "  " << fon << endl;
      
      reonr +=  btr[i]->MWPC_RatioE / Power(btr[i]->MWPC_RatioE_e,2);
      rwonr +=  btr[i]->MWPC_RatioW / Power(btr[i]->MWPC_RatioW_e,2);

      reonre +=  1. / Power(btr[i]->MWPC_RatioE_e,2);
      rwonre +=  1. / Power(btr[i]->MWPC_RatioW_e,2);
      
      REON[fon] = btr[i]->MWPC_RatioE_fit;
      RWON[fon] = btr[i]->MWPC_RatioW_fit;

      REONe[fon] = 2.2*btr[i]->MWPC_RatioE_fit_e;
      RWONe[fon] = 2.2*btr[i]->MWPC_RatioW_fit_e;
      
      fon++;

    } else if(btr[i]->flipperOn == 0&& Abs(btr[i]->MWPC_RatioE) < 0.2 && Abs(btr[i]->MWPC_RatioW) < 0.2) {

      xoff[foff]  = btr[i]->GetRunNumber();
      xoffe[foff] = 0.01;

      reofr +=  btr[i]->MWPC_RatioE / Power(btr[i]->MWPC_RatioE_e,2);
      rwofr +=  btr[i]->MWPC_RatioW / Power(btr[i]->MWPC_RatioW_e,2);

      reofre +=  1. / Power(btr[i]->MWPC_RatioE_e,2);
      rwofre +=  1. / Power(btr[i]->MWPC_RatioW_e,2);
      
      REOFF[foff] = btr[i]->MWPC_RatioE_fit;
      RWOFF[foff] = btr[i]->MWPC_RatioW_fit;
      
      REOFFe[foff] = 2.2*btr[i]->MWPC_RatioE_fit_e;
      RWOFFe[foff] = 2.2*btr[i]->MWPC_RatioW_fit_e;

      foff++;

    }
    

  }

  TGraphErrors *gRen = new TGraphErrors(fon,x,REON,xe,REONe);
  ColorGraphic(gRen,2,20,1);
  TGraphErrors *gRef = new TGraphErrors(foff,xoff,REOFF,xoffe,REOFFe);
  ColorGraphic(gRef,4,20,1);
  TGraphErrors *gRwn = new TGraphErrors(fon,x,RWON,xe,RWONe);
  ColorGraphic(gRwn,2,20,1);
  TGraphErrors *gRwf = new TGraphErrors(foff,xoff,RWOFF,xoffe,RWOFFe);
  ColorGraphic(gRwf,4,20,1);

  Double_t interEon = 0.;
  Double_t interEof = 0.;
  Double_t interWon = 0.;
  Double_t interWof = 0.;
  
  reonre = 1./sqrt(reonre);
  rwonre = 1./sqrt(rwonre);
  reofre = 1./sqrt(reofre);
  rwofre = 1./sqrt(rwofre);

  reonr  = reonr*reonre*reonre;
  reofr  = reofr*reofre*reofre;
  rwonr  = rwonr*rwonre*rwonre;
  rwofr  = rwofr*rwofre*rwofre;

 //  cout << "rwofre " << rwofre << "  rwonr "  << rwonr << endl;
 // cout << "rwonre " << rwonre << "  rwofr "  << rwofr << endl;

  for(Int_t j = 1 ; j <=  hEFlipperOn_G->GetNbinsX(); j++){

    if(aveEOe[j-1] >0)
      aveEOe[j-1] = 1./sqrt(Abs(aveEOe[j-1]));
    else 
      aveEOe[j-1] = 0;
    
    if(aveWOe[j-1] >0)
      aveWOe[j-1] = 1./sqrt(Abs(aveWOe[j-1]));
    else 
      aveWOe[j-1] = 0;
    
    if(aveE1e[j-1] >0)
      aveE1e[j-1] = 1./sqrt(Abs(aveE1e[j-1]));
    else 
      aveE1e[j-1] = 0;
    
    if(aveW1e[j-1] >0)
      aveW1e[j-1] = 1./sqrt(Abs(aveW1e[j-1]));
    else 
      aveW1e[j-1] = 0;
    
    aveEO[j-1]  = aveEO[j-1] * Power(aveEOe[j-1],2);
    //aveWOe[j-1] = 1./sqrt(Abs(aveWOe[j-1]));
    aveWO[j-1]  = aveWO[j-1] * Power(aveWOe[j-1],2);
    //aveE1e[j-1] = 1./sqrt(Abs(aveE1e[j-1]));
    aveE1[j-1]  = aveE1[j-1] * Power(aveE1e[j-1],2);
    //aveW1e[j-1] = 1./sqrt(Abs(aveW1e[j-1]));
    aveW1[j-1]  = aveW1[j-1] * Power(aveW1e[j-1],2); 

    
    if(avesEOe[j-1] >0)
      avesEOe[j-1] = 1./sqrt(Abs(avesEOe[j-1]));
    else 
      avesEOe[j-1] = 0;
    
    //avesEOe[j-1] = 1./sqrt(Abs(avesEOe[j-1]));
    avesEO[j-1]  = avesEO[j-1] * Power(avesEOe[j-1],2);
    //avesWOe[j-1] = 1./sqrt(Abs(avesWOe[j-1]));
    
    if(avesWOe[j-1] >0)
      avesWOe[j-1] = 1./sqrt(Abs(avesWOe[j-1]));
    else 
      avesWOe[j-1] = 0;
    
    avesWO[j-1]  = avesWO[j-1] * Power(avesWOe[j-1],2);
    
    //avesE1e[j-1] = 1./sqrt(Abs(avesE1e[j-1]));
    if(avesE1e[j-1] >0)
      avesE1e[j-1] = 1./sqrt(Abs(avesE1e[j-1]));
    else 
      avesE1e[j-1] = 0;
    
    avesE1[j-1]  = avesE1[j-1] * Power(avesE1e[j-1],2);
    //avesW1e[j-1] = 1./sqrt(Abs(avesW1e[j-1]));
    if(avesW1e[j-1] >0)
      avesW1e[j-1] = 1./sqrt(Abs(avesW1e[j-1]));
    else 
      avesW1e[j-1] = 0;
    avesW1[j-1]  = avesW1[j-1] * Power(avesW1e[j-1],2); 

    

    if( j >=nlow && j <=nhigh){
      interEon += aveEOe[j-1]*aveEOe[j-1];
      interEof += aveE1e[j-1]*aveE1e[j-1];
      interWon += aveWOe[j-1]*aveWOe[j-1];
      interWof += aveW1e[j-1]*aveW1e[j-1];
    }
    
    if(aveEO[j-1] !=0 && aveEOe[j-1] > 0){
      hEFlipperOn_G->SetBinContent(j,aveEO[j-1]);
      hEFlipperOn_G->SetBinError(j,aveEOe[j-1]);
      hEFlipperOn_S->SetBinContent(j,avesEO[j-1]);
      hEFlipperOn_S->SetBinError(j,avesEOe[j-1]);
    } else if(IsNaN(aveEO[j-1]) == 0 ){
      hEFlipperOn_G->SetBinContent(j,0);
      hEFlipperOn_G->SetBinError(j,0);
      hEFlipperOn_S->SetBinContent(j,0);
      hEFlipperOn_S->SetBinError(j,0);
    }
   
    if(aveWO[j-1] !=0 && aveWOe[j-1] >0){      
      hWFlipperOn_G->SetBinContent(j,aveWO[j-1]);
      hWFlipperOn_G->SetBinError(j,aveWOe[j-1]);
      hWFlipperOn_S->SetBinContent(j,avesWO[j-1]);
      hWFlipperOn_S->SetBinError(j,avesWOe[j-1]);
    } else if(IsNaN(aveWO[j-1]) == 0 ){
      hWFlipperOn_G->SetBinContent(j,0);
      hWFlipperOn_G->SetBinError(j,0);
      hWFlipperOn_S->SetBinContent(j,0);
      hWFlipperOn_S->SetBinError(j,0);
    }
    if(aveE1[j-1] !=0 && aveE1e[j-1] > 0){
      hEFlipperOff_G->SetBinContent(j,aveE1[j-1]);
      hEFlipperOff_G->SetBinError(j,aveE1e[j-1]);
      hEFlipperOff_S->SetBinContent(j,avesE1[j-1]);
      hEFlipperOff_S->SetBinError(j,avesE1e[j-1]);
    } else if(IsNaN(aveE1[j-1])==0){
      hEFlipperOff_G->SetBinContent(j,0);
      hEFlipperOff_G->SetBinError(j,0);
      hEFlipperOff_S->SetBinContent(j,0);
      hEFlipperOff_S->SetBinError(j,0);
    }
    if(aveW1[j-1] != 0 && aveW1e[j-1] > 0){
      hWFlipperOff_G->SetBinContent(j,aveW1[j-1]);
      hWFlipperOff_G->SetBinError(j,aveW1e[j-1]);
      hWFlipperOff_S->SetBinContent(j,avesW1[j-1]);
      hWFlipperOff_S->SetBinError(j,avesW1e[j-1]);
    } else if(IsNaN(aveW1[j-1]) == 0){
      hWFlipperOff_G->SetBinContent(j,0);
      hWFlipperOff_G->SetBinError(j,0);
      hWFlipperOff_S->SetBinContent(j,0);
      hWFlipperOff_S->SetBinError(j,0);
    }
    
    if(j>93){
      hWFlipperOff_G->SetBinContent(j,0);
      hWFlipperOff_G->SetBinError(j,0);
      hWFlipperOff_S->SetBinContent(j,0);
      hWFlipperOff_S->SetBinError(j,0);
      hWFlipperOn_G->SetBinContent(j,0);
      hWFlipperOn_G->SetBinError(j,0);
      hWFlipperOn_S->SetBinContent(j,0);
      hWFlipperOn_S->SetBinError(j,0);
      hEFlipperOff_G->SetBinContent(j,0);
      hEFlipperOff_G->SetBinError(j,0);
      hEFlipperOff_S->SetBinContent(j,0);
      hEFlipperOff_S->SetBinError(j,0);
      hEFlipperOn_G->SetBinContent(j,0);
      hEFlipperOn_G->SetBinError(j,0);
      hEFlipperOn_S->SetBinContent(j,0);
      hEFlipperOn_S->SetBinError(j,0);
    }

    //cout << "East On " << j << " " << hEFlipperOn_G->GetBinContent(j-1) << "  " << hEFlipperOn_G->GetBinError(j-1) << endl;
    //cout << "East Off" << j << " " << hEFlipperOff_G->GetBinContent(j-1) << "  " << hEFlipperOff_G->GetBinError(j-1) << endl;
    //cout << "West On " << j << " " << hWFlipperOn_G->GetBinContent(j-1) << "  " << hWFlipperOn_G->GetBinError(j-1) << endl;  
    //cout << "West Off "<< j << " " << hEFlipperOn_S->GetBinContent(j) << "  " << hEFlipperOn_S->GetBinError(j) << endl;
    
  }

  interEon = sqrt(Abs(interEon));
  interEof = sqrt(Abs(interEof));
  interWon = sqrt(Abs(interWon));
  interWof = sqrt(Abs(interWof));

  TCanvas *cNoMWPC = new TCanvas("cNoMWPC","No wire chamber");
  cNoMWPC->Divide(2,2);
  cNoMWPC->cd(1);
  ColorGraphic(hEFlipperOn_G,2,1,1,4);
  ColorGraphic(hEFlipperOn_S,2,1,1,20);
  ColorGraphic(hEFlipperOff_G,1,1,1,21);
  ColorGraphic(hEFlipperOff_S,1,1,1,25);
  hEFlipperOff_S->Draw("P");
  hEFlipperOff_G->Draw("P same");
  hEFlipperOn_G->Draw("same P");
  hEFlipperOn_S->Draw("same P");


  TLegend *lEst = new TLegend(0.6,0.7,0.9,0.9);
  lEst->AddEntry(hEFlipperOff_G,
		 Form("Rate Off (200-600 keV) %5.4f #pm %5.4f",
		      hEFlipperOff_G->Integral(nlow,nhigh)
		      ,interEof)
		 ,"lp");
  lEst->AddEntry(hEFlipperOff_G,
		 Form("Ratio (200-600 keV) %5.4f #pm %5.4f",
		      hEFlipperOff_G->Integral(nlow,nhigh)/hEFlipperOff_S->Integral(nlow,nhigh)
		      ,reofre)
		 ,"lp");
  
  lEst->AddEntry(hEFlipperOn_G,
		 Form("Rate On (200-600 keV) %5.4f #pm %5.4f",
		      hEFlipperOn_G->Integral(nlow,nhigh)
		      ,interEon)
		 ,"lp");
  lEst->AddEntry(hEFlipperOn_G,
		 Form("Ratio (200-600 keV) %5.4f #pm %5.4f",
		      hEFlipperOn_G->Integral(nlow,nhigh)/hEFlipperOn_S->Integral(nlow,nhigh)
		      ,reonre)
		 ,"lp");
  lEst->SetFillColor(0);
  lEst->Draw();

  cNoMWPC->cd(2);

  
  TF1 *fEnl = new TF1("fEnl","[0]",runstart,runstop);//min(x[0],x[fon-1]),max(x[0],x[fon-1]));
  TF1 *fEfl = new TF1("fEfl","[0]",runstart,runstop);//min(xoff[0],xoff[foff-1]),max(xoff[0],xoff[foff-1]));
  
  TMultiGraph *mgRe = new TMultiGraph();
  mgRe->Add(gRen);
  mgRe->Add(gRef);

  fEnl->SetLineColor(2);
  gRen->Fit(fEnl,"RMEQ","goff");

  fEfl->SetLineColor(4);
  gRef->Fit(fEfl,"RMEQ","goff");

  mgRe->Draw("AP");

  TLegend *lEn = new TLegend(0.6,0.7,0.9,0.9);
  lEn->AddEntry(gRen,Form("(On) Average Ratio = %5.3f #pm %5.3f",
			  fEnl->GetParameter(0),fEnl->GetParError(0))
		,"lp");
  lEn->AddEntry(gRen,Form("#chi^{2} / #nu = %5.3f",fEnl->GetChisquare()/fEnl->GetNDF()),"lp");
  
  lEn->AddEntry(gRef,Form("(Off) Average Ratio = %5.3f #pm %5.3f",
			  fEfl->GetParameter(0),fEfl->GetParError(0))
		,"lp");
  lEn->AddEntry(gRef,Form("#chi^{2} / #nu = %5.3f",fEfl->GetChisquare()/fEfl->GetNDF()),"lp");
  
  lEn->SetFillColor(0);
  lEn->Draw();

  
  cNoMWPC->cd(3);


  ColorGraphic(hWFlipperOn_G,2,1,1,4);
  ColorGraphic(hWFlipperOff_G,1,1,1,25);
  ColorGraphic(hWFlipperOn_S,2,1,1,20);
  ColorGraphic(hWFlipperOff_S,1,1,1,21);

  hWFlipperOff_S->Draw("P");
  hWFlipperOff_G->Draw("same P");
  hWFlipperOn_G->Draw("same P");
  hWFlipperOn_S->Draw("same P");
  
  TLegend *lWst = new TLegend(0.6,0.7,0.9,0.9);
  lWst->AddEntry(hWFlipperOff_G,
		 Form("Rate Off (200-600 keV) %5.4f #pm %5.4f",
		      hWFlipperOff_G->Integral(nlow,nhigh),
		      interWof),"lp");
  lWst->AddEntry(hWFlipperOff_G,
		 Form("Ratio (200-600 keV) %5.4f #pm %5.4f",
		      hWFlipperOff_G->Integral(nlow,nhigh) / hWFlipperOff_S->Integral(nlow,nhigh),rwofre)
		 ,"lp");
  lWst->AddEntry(hWFlipperOn_G,
		 Form("Rate On  (200-600 keV) %5.4f #pm %5.4f",
		      hWFlipperOn_G->Integral(nlow,nhigh),
		      interWon),"lp");
  lWst->AddEntry(hWFlipperOn_G,
		 Form("Ratio (200-600 keV) %5.4f #pm %5.4f",
		      hWFlipperOn_G->Integral(nlow,nhigh) / hWFlipperOn_S->Integral(nlow,nhigh),reonre)
		 ,"lp");
  lWst->SetFillColor(0);
  lWst->Draw();
  
  cNoMWPC->cd(4);

  TF1 *fWnl = new TF1("fWnl","[0]",min(x[0],x[fon-1]),max(x[0],x[fon-1]));
  TF1 *fWfl = new TF1("fWfl","[0]",min(xoff[0],xoff[foff-1]),max(xoff[0],xoff[foff-1]));

  TMultiGraph *mgRw = new TMultiGraph();
  mgRw->Add(gRwn);
  mgRw->Add(gRwf);
  
  fWnl->SetLineColor(2);
  gRwn->Fit(fWnl,"RMEQ","goff");

  fWfl->SetLineColor(4);
  gRwf->Fit(fWfl,"RMEQ","goff");

  mgRw->Draw("AP");

  TLegend *lWn = new TLegend(0.6,0.7,0.9,0.9);
  lWn->AddEntry(gRwn,Form("(On) Average Ratio = %5.3f #pm %5.3f",
			  fWnl->GetParameter(0),fWnl->GetParError(0))
		,"lp");
  lWn->AddEntry(gRwn,Form("#chi^{2} / #nu = %5.3f",fWnl->GetChisquare()/fWnl->GetNDF()),"lp");
  lWn->AddEntry(gRwf,Form("(Off) Average Ratio = %5.3f #pm %5.3f",
			  fWfl->GetParameter(0),fWfl->GetParError(0))
		,"lp");
  lWn->AddEntry(gRwf,Form("#chi^{2} / #nu = %5.3f",fWfl->GetChisquare()/fWfl->GetNDF()),"lp");
  lWn->SetFillColor(0);
  lWn->Draw();
  
  TCanvas *fuckup = new TCanvas("fuckup");
  fuckup->Divide(2,1);
  fuckup->cd(1);
  btr[7]->hWMWPC->SetTitle(Form("Run : %d",btr[7]->GetRunNumber()));
  btr[7]->hWNoMWPC->SetLineColor(2);
  btr[7]->hWMWPC->SetLineColor(4);
  btr[7]->hWMWPC->Draw();
  bckr[btr[7]->Bkg_index]->hWMWPC->Draw("same");
  bckr[btr[7]->Bkg_index]->hWNoMWPC->SetLineColor(3);
  bckr[btr[7]->Bkg_index]->hWNoMWPC->Draw("same");
  btr[7]->hWNoMWPC->Draw("same");

  fuckup->cd(2);
  btr[6]->hWMWPC->SetTitle(Form("Run : %d",btr[6]->GetRunNumber()));
  btr[6]->hWNoMWPC->SetLineColor(2);
  btr[6]->hWMWPC->SetLineColor(4);
  btr[6]->hWMWPC->Draw();
  bckr[btr[6]->Bkg_index]->hWMWPC->Draw("same");
  bckr[btr[6]->Bkg_index]->hWNoMWPC->SetLineColor(3);
  bckr[btr[6]->Bkg_index]->hWNoMWPC->Draw("same");
  btr[6]->hWNoMWPC->Draw("same");

// cout << "Background Run " << bckr[btr[7]->Bkg_index]->GetRunNumber() << endl;
  
  return;
  
}

//-------------------------------------------------------------------------------------------------
void Collect_TDCDiff()
{
  
  TH1F *hTDCDiffTot_On  = new TH1F("hTDCDiffTot_On","#Delta TDC Flipper On;Counts; TDCE-TDCW"
				   ,1000,-200,200);
  TH1F *hTDCDiffTot_Off = new TH1F("hTDCDiffTot_Off","#Delta TDC Flipper Off;Counts; TDCE-TDCW"
				   ,1000,-200,200);
  

  for(Int_t i = 0 ; i < nbeta ;i++){
    for(Int_t j = 1 ; j <= hTDCDiffTot_Off->GetNbinsX() ; j++){
      if(btr[i]->flipperOn == 0){
	hTDCDiffTot_Off->SetBinContent(j,hTDCDiffTot_Off->GetBinContent(j)+
				       btr[i]->hTDCDiff->GetBinContent(j));
      } else {
	hTDCDiffTot_On->SetBinContent(j,hTDCDiffTot_On->GetBinContent(j)+
				       btr[i]->hTDCDiff->GetBinContent(j));
      }
    }
  }
  
  TCanvas *cTD = new TCanvas("cTD","TDC Difference");
  cTD->Divide(1,2);
  cTD->cd(1);
  hTDCDiffTot_On->Draw();
  cTD->cd(2);
  hTDCDiffTot_Off->Draw();
  

  return;
}


void Collect_Gammas()
{
  
 
  Double_t x[1000],fge[1000],fgw[1000],fr[1000];
  
  for(Int_t i = 0 ; i < nbeta ; i++){
    x[i]   = btr[i]->GetRunNumber();
    fge[i] = btr[i]->hGammaCounts->GetBinContent(1) / btr[i]->hGammaCountsg->GetBinContent(1);
    fgw[i] = btr[i]->hGammaCounts->GetBinContent(2) / btr[i]->hGammaCountsg->GetBinContent(2);
    if(fgw[i] !=0)fr[i]  = fge[i] / fgw[i];
  }
  
  TGraph *gGaRat  = new TGraph(nbeta,x,fge);
  TGraph *gGaRatw = new TGraph(nbeta,x,fgw);
  TGraph *gRat    = new TGraph(nbeta,x,fr);
  
  ColorGraphic(gGaRat,2,20,2);
  ColorGraphic(gGaRatw,4,20,2);
  ColorGraphic(gRat,4,20,2);
  
  TCanvas *cGamRat = new TCanvas("cGamRat","Gamma Ratios");
  cGamRat->Divide(1,3);
  cGamRat->cd(1);
  gGaRat->Draw("AP");
  cGamRat->cd(2);
  gGaRatw->Draw("AP");
  cGamRat->cd(3);
  gRat->Draw("AP");
}


void Average_Array(Double_t *x,Double_t *xer,Int_t n,Double_t &y_Average,Double_t &y_Average_er)
{
  using namespace TMath;
  
  //------------------------------------------------------------------
  // this function calcutes the error weighted average : 
  // 
  //   y_Average    = Sum_0^n( x[i]i/ xer[i]^2) / Sum_0^n(1./xer[i]^2)
  //   y_Average_er = Sqrt(Sum(1./xer[i]^2))
  //
  
  
  Double_t ytemp   = 0.;
  Double_t ytemper = 0.;
  
  for(Int_t j = 0 ; j < n ; j++){
    if(!(IsNaN(x[j])) && xer[j] > 0 && xer[j] != 1){
	  ytemp    += x[j] / Power(xer[j],2);
	  ytemper  += 1./Power(xer[j],2);
    }
  }
  
  y_Average    = ytemp/ytemper;  // Error weighted Average of the Array passed
  y_Average_er = sqrt(1./ytemper);
  
  return;
}

void  Plot_ChiDis()
{
 
  TF1 *fbeta2 = new TF1("fbeta2","sqrt(x*(x+2.*(510.998910)))/(x+510.998910)",0,2000);  
  TF1 *fline  = new TF1("fline","[0]",0,1000);
  TCanvas *cTest = new TCanvas("cTest","Chi Square Test");
  cTest->cd();
  Double_t chiA = 0.;

  fstream fch;
  
  if(btr[0]->GetGeo()!=2){
    fch.open(Form("chi_squares_%d.txt",btr[0]->GetGeo()),fstream::out);
  } else {
    if(btr[0]->GetRunNumber() <12000){
      fch.open("chi_squares_2.txt",fstream::out);
    } else {
      fch.open("chi_squares_3.txt",fstream::out);
    }
  }
  
  TH1F *hChiOct = new TH1F("hChiOct","#chi^{2} Distribution",40,0.,40.);
  TH1F *hTest = new TH1F("hTest","testing",octet[0]->hAsyTot[2]->GetNbinsX()
					  ,octet[0]->hAsyTot[2]->GetBinCenter(1) - octet[0]->hAsyTot[2]->GetBinWidth(1)/2.
					  ,octet[0]->hAsyTot[2]->GetBinCenter(octet[0]->hAsyTot[2]->GetNbinsX()) + octet[0]->hAsyTot[2]->GetBinWidth(1)/2.);
  for(Int_t i = 0 ; i < noct ; i++){
    for(Int_t j = 1 ; j <= octet[i]->hAsyTot[2]->GetNbinsX() ; j++){
      hTest->SetBinContent(j,octet[i]->hAsyTot[2]->GetBinContent(j) / fbeta2->Eval(octet[i]->hAsyTot[2]->GetBinCenter(j)));
      hTest->SetBinError(j,octet[i]->hAsyTot[2]->GetBinError(j) / fbeta2->Eval(octet[i]->hAsyTot[2]->GetBinCenter(j)));
    }
    hTest->Fit(fline,"REMQ","",octet[i]->hAsyTot[2]->GetBinCenter(nlow),octet[i]->hAsyTot[2]->GetBinCenter(nhigh));
    hChiOct->Fill(fline->GetChisquare());
    fch << i << "\t" << fline->GetChisquare() << endl;
  }
  
  hChiOct->Draw();
  TF1 *fchi = new TF1("fchi",ChiSquareDis,0,40.,2);
  fchi->SetParameter(0,hChiOct->Integral());
  fchi->SetParameter(1,(nhigh-nlow-1));
  fchi->SetLineColor(2);
  fchi->Draw("same");
  
}

Double_t ChiSquareDis(Double_t *x,Double_t *par)
{
  
  using namespace TMath;
  
  Double_t xx  = x[0];   // chi^2
  Double_t nd  = par[1]; // degrees of freedom
  Double_t scl = par[0]; // scale
  
  Double_t chip = scl*(Power(xx,nd/2.-1)*exp(-xx/2.)/(Gamma(nd/2.)*Power(2,nd/2.)));
  return chip;
  
}

void average_type1()
{
  TH1F *hE1Type_Pr = new TH1F("hE1Type_Pr",";Energy;Cts",100,0,1);
  TH1F *hE1Type_Sc = new TH1F("hE1Type_Sc",";Energy;Cts",100,0,1);
  TH1F *hW1Type_Pr = new TH1F("hW1Type_Pr",";Energy;Cts",100,0,1);
  TH1F *hW1Type_Sc = new TH1F("hW1Type_Sc",";Eenrgy;Cts",100,0,1);

  for(Int_t i = 0 ; i < nbeta ; i++){
    for(Int_t j = 1 ; j <= 100 ; j++){
      double epk1 = btr[i]->hEType1_Primary->GetBinCenter(j);
      hE1Type_Pr->SetBinContent(j,hE1Type_Pr->GetBinContent(j) + btr[i]->hEType1_Primary->GetBinContent(j)*btr[i]->rtime_e);
      hE1Type_Sc->SetBinContent(j,hE1Type_Sc->GetBinContent(j) + btr[i]->hEType1_Secondary->GetBinContent(j)*btr[i]->rtime_e);
      hW1Type_Pr->SetBinContent(j,hW1Type_Pr->GetBinContent(j) + btr[i]->hWType1_Primary->GetBinContent(j)*btr[i]->rtime_w);
      hW1Type_Sc->SetBinContent(j,hW1Type_Sc->GetBinContent(j) + btr[i]->hWType1_Secondary->GetBinContent(j)*btr[i]->rtime_w);
    }
  }
  
  TCanvas *cType1 = new TCanvas("cType1","Type 1 Primary vs. Secondary",600,600);
  cType1->cd();
  fstream penin;
  penin.open("type_1_e_fraction_3.txt",fstream::in);
  TH1F *hSimT1p = new TH1F("hSimT1p","East",100,0,1);
  TH1F *hSimT2p = new TH1F("hSimT2P","West",100,0,1);
  for(Int_t i = 1 ; i <= 100 ; i++){
    double h1,h2,h3;
    penin >> h1 >> h2 >> h3;
    hSimT1p->SetBinContent(i,h2);
    hSimT2p->SetBinContent(i,h3);
  }
  penin.close();
  hSimT2p->Sumw2();
  hSimT1p->Sumw2();
  ColorGraphic(hE1Type_Pr,1,1,1,20,0.8);
  ColorGraphic(hE1Type_Sc,1,1,1,20,0.8);
  ColorGraphic(hW1Type_Pr,2,1,1,24,0.8);
  ColorGraphic(hW1Type_Sc,1,1,1,24,0.8);
  ColorGraphic(hSimT1p,1,1,2,1,1);
  ColorGraphic(hSimT2p,2,1,2,1,1);

  hSimT1p->Scale(hE1Type_Pr->Integral()/hSimT1p->Integral());
  hSimT2p->Scale(hW1Type_Pr->Integral()/hSimT2p->Integral());

  hE1Type_Pr->Draw("e0 x0");
  hE1Type_Pr->GetXaxis()->SetTitle("E_{vis}(primary) / E_{vis}(total)");
  hE1Type_Pr->GetXaxis()->CenterTitle();
  hE1Type_Pr->GetYaxis()->CenterTitle();  
  hE1Type_Pr->GetXaxis()->SetNdivisions(509);
  // hE1Type_Sc->Draw("HIST same");
  hW1Type_Pr->Draw("same e0 x0");
  hSimT1p->Draw("hist same");
  hSimT2p->Draw("hist same");
  //hW1Type_Sc->Draw("HIST same");
  TLegend *lT1p = new TLegend(0.41,0.2,0.7,0.49);
  lT1p->SetLineColor(0);
  lT1p->SetFillColor(0);
  lT1p->AddEntry(hE1Type_Pr,"Data East","lp");
  lT1p->AddEntry(hSimT1p,"MC East","l");
  lT1p->AddEntry(hW1Type_Pr,"Data West","lp");
  lT1p->AddEntry(hSimT2p,"MC West","l");
  
  lT1p->Draw();
  Int_t ggo = btr[0]->GetGeo();
  cType1->Print(Form("pdf_out/geo%d/type_1_comp_pr_%d.pdf",ggo,ggo));
		
  
};
