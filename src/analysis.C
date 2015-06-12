#ifndef __analysis_c__
#define __analysis_c__

#include "analysis.h"
#include "variables.h"
#include "AnalysisDB.hh"
#include "Run.h"

#include <algorithm>

#define DEBUG      1
//#define NRUNFIRST  14206
#define NRUNFIRST  17125
#define NRUNLAST   19966

#define MAXRUNS    NRUNLAST-NRUNFIRST
 
int main(int argc,char *argv[])
{
  // Set up the ROOT environment------------------------------
  Int_t remake = 1;
  SetROOTOpt();
  ZeroThings();
  if(!(SetAnalysisType(&remake)))return -1;
  cout << endl;
  if(DEBUG == 1) cout << "Finished with getting the run list " << endl;
  // Load the 2/3 separation probability from simulation
  Load23SepArray();
  //--------------------------------------------------------------------------------------
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
    if(DEBUG == 1)cout << "Done with Backgrounds" << endl;
    analyze_octets();
    cout << "Finished Octets"<< endl;
  } 
//  return -1;
  //--------------------------------------------------------------------------------------
//  for(Int_t jrun = 0 ; jrun < nbeta ; jrun++){ 
//    btr[jrun]->Calculate_Backscatter(1,1);
//    btr[jrun]->GetSimpleAsym();
//  }
  //Close the connection to the database;
  sql->Close();
// cout << "Analysis Refill ";
//  if(nbeta > 0)  AnalysisTaskRefill(nbeta, btr, nbetas);
// cout << "Succeeded!" << endl;
  CallAnalysisTasks();

  if(DEBUG == 1)cout << "Quiting" << endl;

  return 0;

};
//-----------------------------------------------------------------------------------------
void CallAnalysisTasks()
{
  //Begin outputting things.....
  cout << "Collecting Rates " << endl;
  CollectRates();
  cout << "Drawing Rates " << endl;
  DrawRates();
/*  cout << "Plotting Positions" << endl;
  Collect_Pos();
  cout << "Collecting Type Rotation " << endl;
  CollectTypeRot();
  cout << "Collecting TDC Corruption Data " << endl;
  CollectTDCCor();
  cout << "Lets Find the Octets " << endl;
  CalcSimplSuper();
  cout << "Let's plot superratio of live times" << endl;
  Plot_Timing();
  cout << "Let's fill timing histos" << endl;
  Fill_Timing();
  cout << "Plot Energy Chis" << endl;
  Plot_E_Chis();
  cout << "Plot Type 2/3 Difference" << endl;
  Plot_MonRate();
  cout << "Getting Octets " << endl;  
  Collect_Octets();
  cout << "Finished Octets" << endl;
  Average_A();
  Collect_TvsE();
  Collect_23Anode();
  Collect_Stuff();
  Collect_Rad();
  Collect_Energy_Spectra();
  GammaBack();
  Collect_Gammas();
  Collect_TDCDiff();
  Plot_ChiDis();
  Plot_Multi();
  TrackAnodeMPV();
  TrackStats();
  average_type1();
  PlotRunTimes();
   */
}

/*void AnalysisTaskRefill(Int_t nrun, std::vector<Beta_Run*> btr, vector<Int_t> nrunl){
	TString refillFile;
	Double_t dummy9823;

  for(Int_t i=0; i< (int)nrunl.size();i++){
	refillFile="/home/jwwexler/robbyWork/histOut/hists_";
	refillFile+=nrunl[i];
	refillFile+=".root";

	cout << "Opening " << refillFile;

	TFile *fRefill = new TFile(refillFile,"READ");

	cout << endl;
	cout << Form("/%s_Analysis_%d_%d/hpe_%d",getenv("UCNA_ANA_AUTHOR"),50,3,nrunl[i]) << "= name of histo in " << refillFile << endl;
//	dummy9823 = (TH2F*)fRefill->Get(Form("/%s_Analysis_%d_%d/hpe_%d",getenv("UCNA_ANA_AUTHOR"),50,3,nrunl[i]))->GetEntries();

	btr[i]->hpe = (TH2F*)fRefill->Get(Form("/%s_Analysis_%d_%d/hpe_%d",getenv("UCNA_ANA_AUTHOR"),50,3,nrunl[i]));

	btr[i]->hpw = (TH2F*)fRefill->Get(Form("/%s_Analysis_%d_%d/hpw_%d",getenv("UCNA_ANA_AUTHOR"),50,3,nrunl[i]));  // Normally 50, 3 are RAD and rotated?
	cout << "i, hpe, hpw " << i << " " << btr[i]->hpe->GetEntries() << " " << btr[i]->hpw->GetEntries() << endl;
	delete fRefill;
  }

}*/

//---------------------------------------------------------------------------
Int_t analyze_background_runs(Int_t n, std::vector<Bck_Run*>bk,vector<Int_t> nrunl,Int_t remake)
{
// Analysis flow for background runs.
  for(Int_t i = 0 ; i < (int)nrunl.size() ; i++){
    if(bk[i]->OpenRun(nrunl[i])){
      bk[i]->Draw_2d(nrunl[i],i);
      bk[i]->Fill(i,remake,sep23);
      bk[i]->Scale2Time(1,1);
      bk[i]->Calculate_Backscatter(1,1);
      bk[i]->Diagnosis_Run();
    }
    bk[i]->DeleteHistos();
  }

 return 0;
}
//---------------------------------------------------------------------------
Int_t analyze_beta_runs(Int_t n, std::vector<Beta_Run*>bta,vector<Int_t> nrunl,Int_t remake)
{
 
  for(Int_t i = 0 ; i < (int)nrunl.size() ; i++){
    if(bta[i]->OpenRun(nrunl[i])){
      bta[i]->Draw_2d(nrunl[i],i);
      bta[i]->Fill(i,remake,sep23);
      bta[i]->Scale2Time(1,1);
      bta[i]->Diagnosis_Run();
    }
    bta[i]->DeleteHistos();
  }

 return 0;
}
//---------------------------------------------------------------------------
Int_t Subtract_Backgrounds(std::vector< Beta_Run *>bta,std::vector<Bck_Run*>bk, Int_t nb,Int_t nbk)
{
  
  // Loop through the beta runs and subtract off the background spectra
  Int_t i,j,k;
  
  for(Int_t n = 0 ; n < nb ; n++){
    i = bta[n]->GetBackgroundRun();
    j = 0;
    k = -1;
    if( i != 0 ) {
      do {
        k++;
	if(k == nbk) break;
        j = bk[k]->GetRunNumber();
      } while( i != j && k < nbk);
      if(i == j){
	bta[n]->Bkg_index = k;
//	bta[n]->SubBck(bk[k]);
//	bta[n]->GetEnergyChi();
//	bta[n]->Make_Kurie();
      }
    }
  }
  
  return 0;
  
}
//---------------------------------------------------------------------------
Int_t analyze_octets()
{
  fstream runlist;
  runlist.open("input_files/octet_list_D.txt",fstream::in);
  // Here you would build the Octets ......
  Int_t noctets;
  noct = 0;
  
  if(runlist.is_open()){
    do{
      runlist >> noctets >> noctstart >> noctstop;
      cout << noctets << "\t" << noctstart << "\t" << noctstop << endl;      
      if(noctstart >= runstart && noctets != 1111 && noctstart !=0){
        // Open the Octet object
        octet.push_back(new Octet(noct,noctstart,noctstop,radcut));
        // define histograms
        octet[noct]->Initialize_Histo(noct);
        // Use the MySQL database to determine the runs in each octet and fill the octet
        // defined runs.  If there are multiple DAQ runs for any such octet run they will
        // be summed such that the times are correctly counted and the rates are interms of
        // the total counts of all re
        octet[noct]->Find_Runs(btr,bckr,nbeta);
	// Need to adjust the octect calculation for the full octet analysis.
	octet[noct]->Calc_Super();
        octet[noct]->Calc_Super_Time();
	octet[noct]->Calc_A_sum();
	octet[noct]->Calc_A_multi();
	octet[noct]->Calc_A_sum_Bin();
	octet[noct]->Get_Rad_A();
        // Make and output file in the "output" directory with all relative octet rate info.
	octet[noct]->Debugger(noct);
	//octet[noct]->OutPutToDB();
	noct++;
      }
    }while(!(runlist.eof()));
  }
  noct--;
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

  nrun.push_back(irun);
  ntype.push_back(rty);

  if(rty == 2){
    nbeta++;
    nbck++;
    cout << "Enter Background :  " ;
    cin >> rb;
    nrun.push_back(rb);
    ntype.push_back(3);
    nbetas.push_back(irun);
    nrunb[0]  = rb;
    runtot    = 2;
    nbcks.push_back(rb);
  }else if(rty == 3){
    nbck++;
    runtot=1;
    nbcks.push_back(irun);
    nrunb[0] = irun;
  }

  return;
}
//---------------------------------------------------------------------------
Int_t ParseMPMOctetList(vector<Int_t> &RunListMPM)
{
  fstream mpmlist;
  char runline[500];
  char *seg;
  Int_t nn = 0;
  Int_t nocts = 0;
  vector<Int_t> current_oct (100);
  fstream foct,frlist;
  foct.open("input_files/octet_list_D.txt",fstream::out);
  frlist.open("input_files/parsed_run_list.txt",fstream::out);
  mpmlist.open("input_files/OctetList_Dummy.txt",fstream::in);

  if(mpmlist.is_open()){
    do{
      nn = 0;
      mpmlist.getline(runline,500);
      //cout << runline << endl;
      if(strncmp(runline,"#",1)!=0){
        if(strncmp(runline,"Octet",5)==0)nocts++;
	seg = strtok(runline,"= \t");
	while ( seg != NULL){
          //use this if statement to cut out an octet with nocts.
	  if(strlen(seg) == 5 && nocts != 9){
          //if(strlen(seg) == 5){
                  if(atoi(seg)!=0){
		  current_oct[nn] = atoi(seg);
		  nn++;
		  RunListMPM.push_back(atoi(seg));
		}
	  }
	  seg = strtok(NULL,"= \t");         
	}
	if(*min_element(current_oct.begin(),current_oct.begin()+nn) !=0){
	  foct << nocts << "\t" << *min_element(current_oct.begin(),current_oct.begin()+nn) << "\t";
	  foct << *max_element(current_oct.begin(),current_oct.begin()+nn) << endl;
	  cout << nocts << "\t" << *min_element(current_oct.begin(),current_oct.begin()+nn);
	  cout << "\t" << *max_element(current_oct.begin(),current_oct.begin()+nn) << endl;
	}
	for(int ii = 0 ; ii < nn ; ii++){
	  frlist << current_oct[ii] << endl;
	  current_oct[ii]=0.;
	}
      }
      cout << nocts << endl;
    } while(!(mpmlist.eof()));
  } else {
    cout << "Run list not found !!!" << endl;
    return -1;
  }
  cout << "Size of Run list " << (int)RunListMPM.size() << 
	  "\t" << *max_element(RunListMPM.begin(),RunListMPM.begin()+nn) << endl;
  mpmlist.close();
  foct.close();

  return (Int_t)RunListMPM.size();
};
//---------------------------------------------------------------------------
bool GetListofRuns(Int_t n1,Int_t n2)
{
   
  vector<Int_t> RunListMPM;
  fstream f;
  Int_t rn = 0;
  TSQLResult *res;
  TSQLRow *row;
  char buffer[500];

  Int_t nruns = ParseMPMOctetList(RunListMPM);
  if(nruns == -1)return kFALSE;
  // The user name and password is included in the 
  char *db = Form("mysql://%s/%s",getenv("UCNADBADDRESS"),getenv("UCNADB"));
  char *dbuser = Form("%s",getenv("UCNADBUSER"));
  char *dbpwd  = Form("%s",getenv("UCNADBPASS"));
  // Make a connection to the server at caltech
  sql = TSQLServer::Connect(db,dbuser,dbpwd);
  // Check if connected
  if(!(sql->IsConnected())){
    while(!(sql->IsConnected())){
      sql->Connect(db,dbuser,dbpwd);
    }
  }
  
  for(Int_t nn = 0 ; nn < nruns ; nn++){
  // ask for the run numbers and the gate valve state
  sprintf(buffer,"select t1.run_number,t1.gate_valve,t2.live_time_e from run as t1, analysis as t2 " 
	         "where t1.run_number between %d and %d and t1.run_number = %d and t1.run_type = 'Asymmetry' "
                 "and t1.run_number = t2.run_number ORDER by t1.run_number "
                 ,n1,n2,RunListMPM[nn]); 
  // send the query
  res = (TSQLResult*)sql->Query(buffer);
  //quick list
  if(res->GetRowCount() != 0){
    while((row = (TSQLRow*)res->Next())){
      rn = atoi(row->GetField(0));
      //make sure the run is not in the bad list............
      if(rn != 2536 && rn != 7837 && rn != 9812 && rn != 9808){
	  if( !strncmp(row->GetField(1),"Open",5) && atoi(row->GetField(2)) > 0){ 
	    //-----------------------------------------------------------------------
	    // Open a beta decay run
            if(checkruns(rn,nbetas,runtot)){
              nbetas.push_back(atoi(row->GetField(0)));
	      nrun.push_back(nbetas[nbeta]);
	      ntype.push_back(2);
	      btr.push_back( new Beta_Run(nbetas[nbeta],2,sql));
	      nbeta++;
	      runtot++;
	      //-----------------------------------------------
	      // Open a background run for the above beta run
              if(checkruns(btr[nbeta-1]->GetBackgroundRun(),nbcks,nbck)){
                nbcks.push_back(btr[nbeta-1]->GetBackgroundRun());
	        nrun.push_back(nbcks[nbck]);
	        nrunb[nbck]   = nbcks[nbck];
	        ntype.push_back(3);
	        bckr.push_back( new Bck_Run(nbcks[nbck],3,sql));
	        nbck++;
	        runtot++;
	      }
            }
	  }
      }
      delete row; 
     }
  } 
  delete res;
  }
  return kTRUE;
}
//---------------------------------------------------------------------------
void SetROOTOpt()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  TH1::SetDefaultSumw2(kTRUE);
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

  // Define arrays to hold Backscattering fractions and errorbars
  vector<Double_t> xrune,xrunw,xrunw1,xrunw23,xrune1,xrune23;
  vector<Double_t> ertt,erte,wrtt,wrte,ert1,er1e,wrt1,wr1e;
  vector<Double_t> ert23,er23e,wrt23,wr23e,erttf,ertef,wrttf,wrtef;
  vector<Double_t> ert1f,er1ef,wrt1f,wr1ef,ert23f,er23ef,wrt23f,wr23ef;
  vector<Double_t> xrunef,xrunwf,xrunw1f,xrunw23f,xrune1f,xrune23f;
   
  // Loop through background subtracted beta runs for 
  for(Int_t jrun = 0 ; jrun < nbeta ; jrun++){
    if(btr[jrun]->flipperOn == 1){
      SetRateV(xrune,ertt,erte,btr[jrun]->GetRunNumber(),btr[jrun]->bscat.e_all,TMath::Abs(btr[jrun]->bscat.e_alle),0);
      SetRateV(xrunw,wrtt,wrte,btr[jrun]->GetRunNumber(),btr[jrun]->bscat.w_all,TMath::Abs(btr[jrun]->bscat.w_alle),0);
           
      erttt.push_back(ertt[(int)ertt.size()-1] + wrtt[(int)wrtt.size()-1]);
      ertet.push_back(sqrt( wrte[(int)wrte.size()-1]*wrte[(int)wrte.size()-1] 
	                  + erte[(int)erte.size()-1]*erte[(int)erte.size()-1]));

      SetRateV(xrune1,ert1,er1e,btr[jrun]->GetRunNumber(),
             (btr[jrun]->bscat.etype_1 / (btr[jrun]->bscat.e_all))*100.,
                    Rate_Error(btr[jrun]->bscat.etype_1,
                                 btr[jrun]->bscat.e_all,
                                 btr[jrun]->rtime_e,
                                 btr[jrun]->bscat.etype_1_bck,
                                 btr[jrun]->bscat.e_all_bck),1);

      SetRateV(xrunw1,wrt1,wr1e,btr[jrun]->GetRunNumber(),
             (btr[jrun]->bscat.wtype_1 / (btr[jrun]->bscat.w_all)) * 100.,
                     Rate_Error(btr[jrun]->bscat.wtype_1,
			          btr[jrun]->bscat.w_all,
                                  btr[jrun]->rtime_w,
                                  btr[jrun]->bscat.wtype_1_bck,
                                  btr[jrun]->bscat.w_all_bck),1);

      SetRateV(xrune23,ert23,er23e,btr[jrun]->GetRunNumber(),
             (btr[jrun]->bscat.etype_23 / (btr[jrun]->bscat.e_all))*100.,
                       Rate_Error(btr[jrun]->bscat.etype_23,
                                     btr[jrun]->bscat.e_all,
                                     btr[jrun]->rtime_e,
                                     btr[jrun]->bscat.etype_23_bck,
                                     btr[jrun]->bscat.e_all_bck),1);

      SetRateV(xrunw23,wrt23,wr23e,btr[jrun]->GetRunNumber(),
             (btr[jrun]->bscat.wtype_23 / (btr[jrun]->bscat.w_all)) * 100.,
                    Rate_Error(btr[jrun]->bscat.wtype_23,
                                  btr[jrun]->bscat.w_all,
                                  btr[jrun]->rtime_w,
                                  btr[jrun]->bscat.wtype_23_bck,
                                  btr[jrun]->bscat.w_all_bck),1);

    } else if(btr[jrun]->flipperOn == 0){
      SetRateV(xrunef,erttf,ertef,btr[jrun]->GetRunNumber(),btr[jrun]->bscat.e_all,TMath::Abs(btr[jrun]->bscat.e_alle),0);
      SetRateV(xrunwf,wrttf,wrtef,btr[jrun]->GetRunNumber(),btr[jrun]->bscat.w_all,TMath::Abs(btr[jrun]->bscat.w_alle),0);
      
      ertttf.push_back(erttf[(int)erttf.size()-1] + wrttf[(int)wrttf.size()-1]);
      ertetf.push_back(sqrt( wrtef[(int)wrttf.size()-1]*wrtef[(int)wrttf.size()-1] 
                           + ertef[(int)ertef.size()-1]*ertef[(int)ertef.size()-1]));
     
      SetRateV(xrune1f,ert1f,er1ef,btr[jrun]->GetRunNumber(),
             (btr[jrun]->bscat.etype_1 / (btr[jrun]->bscat.e_all))*100.,
                     Rate_Error(btr[jrun]->bscat.etype_1,
                                  btr[jrun]->bscat.e_all,
                                  btr[jrun]->rtime_e,
                                  btr[jrun]->bscat.etype_1_bck,
                                  btr[jrun]->bscat.e_all_bck),1);

      SetRateV(xrunw1f,wrt1f,wr1ef,btr[jrun]->GetRunNumber(),
             (btr[jrun]->bscat.wtype_1 / (btr[jrun]->bscat.w_all)) * 100.,
                      Rate_Error(btr[jrun]->bscat.wtype_1,
                                  btr[jrun]->bscat.w_all,
                                  btr[jrun]->rtime_w,
                                  btr[jrun]->bscat.wtype_1_bck,
                                  btr[jrun]->bscat.w_all_bck),1);

      SetRateV(xrune23f,ert23f,er23ef,btr[jrun]->GetRunNumber(),
             (btr[jrun]->bscat.etype_23 / (btr[jrun]->bscat.e_all))*100.,
                       Rate_Error(btr[jrun]->bscat.etype_23,
                                     btr[jrun]->bscat.e_all,
                                     btr[jrun]->rtime_e,
                                     btr[jrun]->bscat.etype_23_bck,
                                     btr[jrun]->bscat.e_all_bck),1);
      
      SetRateV(xrunw23f,wrt23f,wr23ef,btr[jrun]->GetRunNumber(),
             (btr[jrun]->bscat.wtype_23 / (btr[jrun]->bscat.w_all)) * 100.,
                    Rate_Error(btr[jrun]->bscat.wtype_23,
                                  btr[jrun]->bscat.w_all,
                                  btr[jrun]->rtime_w,
                                  btr[jrun]->bscat.wtype_23_bck,
                                  btr[jrun]->bscat.w_all_bck),1);
    }
     
  }
  Double_t xre[(Int_t)xrune.size()];
  copy(xrune.begin(),xrune.end(),xre);

  // Raw Rate output 
  grtot = new TGraphErrors((Int_t)xrune.size(),xre,&erttt[0],0,&ertet[0]); 
  gre   = new TGraphErrors((Int_t)xrune.size(),xre,&ertt[0] ,0,&erte[0]);
  grw   = new TGraphErrors((Int_t)xrunw.size(),&xrunw[0],&wrtt[0] ,0,&wrte[0]);
  gre1  = new TGraphErrors((Int_t)xrune1.size(),&xrune1[0],&ert1[0] ,0,&er1e[0]);
  grw1  = new TGraphErrors((Int_t)xrunw1.size(),&xrunw1[0],&wrt1[0] ,0,&wr1e[0]);
  gre23 = new TGraphErrors((Int_t)xrune23.size(),&xrune23[0],&ert23[0],0,&er23e[0]);
  grw23 = new TGraphErrors((Int_t)xrunw23.size(),&xrunw23[0],&wrt23[0],0,&wr23e[0]);
  // Flipper on
  grtotf = new TGraphErrors((Int_t)xrunef.size(),&xrunef[0],&ertttf[0]  ,0,&ertetf[0]); 
  gref   = new TGraphErrors((Int_t)xrunef.size(),&xrunef[0],&erttf[0]   ,0,&ertef[0]);
  grwf   = new TGraphErrors((Int_t)xrunwf.size(),&xrunwf[0],&wrttf[0]   ,0,&wrtef[0]);
  gre1f  = new TGraphErrors((Int_t)xrune1f.size(),&xrune1f[0],&ert1f[0]   ,0,&er1ef[0]);
  grw1f  = new TGraphErrors((Int_t)xrunw1f.size(),&xrunw1f[0],&wrt1f[0]   ,0,&wr1ef[0]);
  gre23f = new TGraphErrors((Int_t)xrune23f.size(),&xrune23f[0],&ert23f[0]  ,0,&er23ef[0]);
  grw23f = new TGraphErrors((Int_t)xrunw23f.size(),&xrunw23f[0],&wrt23f[0]  ,0,&wr23ef[0]);

}
//------------------------------------------------------------------
void DrawRates()
{
  
  TF1 *fe1 = new TF1("fe1","[0]",runstart,runstop);
  TF1 *fe2 = new TF1("fe2","[0]",runstart,runstop);
  TF1 *fw1 = new TF1("fw1","[0]",runstart,runstop);
  TF1 *fw2 = new TF1("fw2","[0]",runstart,runstop);
  
  TF1 *ffe1 = new TF1("ffe1","[0]",runstart,runstop);
  TF1 *ffe2 = new TF1("ffe2","[0]",runstart,runstop);
  TF1 *ffw1 = new TF1("ffw1","[0]",runstart,runstop);
  TF1 *ffw2 = new TF1("ffw2","[0]",runstart,runstop);

  fe1->SetLineColor(3);
  fe2->SetLineColor(3);
  fw1->SetLineColor(4);
  fw2->SetLineColor(4);
  ffe1->SetLineColor(2);
  ffe2->SetLineColor(2);
  ffw1->SetLineColor(1);
  ffw2->SetLineColor(1);
  
  TCanvas *cgr = new TCanvas("cgr","Total Rates");
  cgr->cd();
  // flipper off
  ColorGraphic(gre,2,20,2,1,"Total Event Rates","Run Number","Total Rate s^{-1}");
  ColorGraphic(grw,1,20,2);
  ColorGraphic(gre1,3,20,2,1,"Type 1 Event Fractions (Flipper On)","Run Number","Fraction");
  ColorGraphic(grw1,4,20,2);
  ColorGraphic(gre23,3,20,2,1,"Type 2/3 Event Fractions (Flipper On)","Run Number","Fraction");
  ColorGraphic(grw23,4,20,2);
  // flipper on 
  ColorGraphic(gref,2,4,2);
  ColorGraphic(grwf,1,4,2);
  ColorGraphic(gre1f,2,4,2,1,"Type 1 Event Fractions (Flipper off)","Run Number","Fraction");
  ColorGraphic(grw1f,1,4,2);
  ColorGraphic(gre23f,2,4,2,1,"Type 2/3 Event Fractions (Flipper On)","Run Number","Fraction");
  ColorGraphic(grw23f,1,4,2);

  gre->Draw("AP");
  gre->GetYaxis()->SetRangeUser(0,50);
  grw->Draw("P");
  grwf->Draw("P");
  gref->Draw("P");
  
  TLegend *lRates = new TLegend(0.6,0.7,0.9,0.9);
  lRates->AddEntry(gre,"East Rate","lp");
  lRates->AddEntry(grw,"West Rate","lp");
  lRates->AddEntry(gref,"East Rate flipper off","lp");
  lRates->AddEntry(grwf,"West Rate flipper off","lp");
  lRates->Draw();
  cgr->Print("output_files/rates.pdf");

  TCanvas *cBckFractions = new TCanvas("cBckFractions","Backscatter Fractions");
  cBckFractions->Clear();
  cBckFractions->Divide(2,2);
  // Draw East Backscattering fractions............
  TLegend *lb1= new TLegend(0.5,0.6,0.9,0.9);
  DrawBackFractionPad(cBckFractions,1,gre1,grw1,fe1,fw1,lb1,8);
  TLegend *lb3 = new TLegend(0.5,0.6,0.9,0.9);
  DrawBackFractionPad(cBckFractions,2,gre1f,grw1f,ffe1,ffw1,lb3,8);
  TLegend *lb2 = new TLegend(0.5,0.6,0.9,0.9);
  DrawBackFractionPad(cBckFractions,3,gre23,grw23,fe2,fw2,lb2,6);
  TLegend *lb4 = new TLegend(0.5,0.6,0.9,0.9);
  DrawBackFractionPad(cBckFractions,4,gre23f,grw23f,ffe2,ffw2,lb4,6);
  cBckFractions->Print("output_files/backscatter_fractions.pdf");
 
  delete fe1;
  delete fe2;
  delete fw1;
  delete fw2;
  delete ffe1;
  delete ffe2;
  delete ffw1;
  delete ffw2;
  
}
//---------------------------------------------------------------------------
void CollectTypeRot()
{
  
  const int rbins = 150;
  Float_t rad     = 75.;
  Int_t last      = 0.  ;
  Double_t ecount[rbins][rbins],wcount[rbins][rbins];
  //Define some collection histograms------------------------------------
  DefineRotationCollectionHistos(rbins,rad);
  //--------------------------------------------
  // Clean the array,prepare to fill with Total events
  CleanArray(ecount,rbins,rbins);
  CleanArray(wcount,rbins,rbins);

	cout << "Pre-CTR Loop" << endl;

  // Loop through all events and fill the total array
  for(Int_t i = 0 ; i < nbeta ; i++){
    if(i == nbeta - 1)last = 1;

    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    CollectAllRot(ecount,rbins,rbins,btr[i]->hRote,hTotRote,last);
    CollectAllRot(wcount,rbins,rbins,btr[i]->hRotw,hTotRotw,last);
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  last = 0; //Reset last counter
  // Clean the array,prepare to fill with Type 1 events
  CleanArray(ecount,rbins,rbins);
  CleanArray(wcount,rbins,rbins);
  for(Int_t i = 0 ; i < nbeta ; i++){
    if(i == nbeta - 1)last = 1;
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    CollectAllRot(ecount,rbins,rbins,btr[i]->hRoteI,hTotRoteI,last);
    CollectAllRot(wcount,rbins,rbins,btr[i]->hRotwI,hTotRotwI,last);
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  last = 0;

	cout << "Halfway through loops" << endl;

  //-------------------------------------------------------------------
  CleanArray(ecount,rbins,rbins);
  CleanArray(wcount,rbins,rbins);
  for(Int_t i = 0 ; i <nbeta ; i++){
    if(i == nbeta - 1)last = 1;
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    CollectAllRot(ecount,rbins,rbins,btr[i]->hRote23,hTotRote23,last);
    CollectAllRot(wcount,rbins,rbins,btr[i]->hRotw23,hTotRotw23,last);
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  last = 0;
  //-------------------------------------------------------------------
  CleanArray(ecount,rbins,rbins);
  CleanArray(wcount,rbins,rbins);
  for(Int_t i = 0 ; i <nbeta ; i++){
    if(i == nbeta - 1)last = 1;
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    CollectAllRot(ecount,rbins,rbins,btr[i]->hRoteI23,hTotRoteI23,last);
    CollectAllRot(wcount,rbins,rbins,btr[i]->hRotwI23,hTotRotwI23,last);
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  last = 0;

  cout << "Post-Fill loops" << endl;

  //------------------------------------------------------------------
  // finished filling now draw
  
  hTotRoteI23->GetZaxis()->SetRangeUser(0,hTotRoteI23->GetBinContent(hTotRoteI23->GetMaximumBin())*1.1);
  hTotRotwI23->GetZaxis()->SetRangeUser(0,hTotRotwI23->GetBinContent(hTotRotwI23->GetMaximumBin())*1.1);
  hTotRote23->GetZaxis()->SetRangeUser(0,hTotRote23->GetBinContent(hTotRote23->GetMaximumBin())*1.1);
  hTotRotw23->GetZaxis()->SetRangeUser(0,hTotRotw23->GetBinContent(hTotRotw23->GetMaximumBin())*1.1);
 
  hTotRote->GetZaxis()->SetRangeUser(0,hTotRote->GetBinContent(hTotRote->GetMaximumBin())*1.1);
  hTotRotw->GetZaxis()->SetRangeUser(0,hTotRotw->GetBinContent(hTotRotw->GetMaximumBin())*1.1);
  hTotRoteI->GetZaxis()->SetRangeUser(0,hTotRoteI->GetBinContent(hTotRoteI->GetMaximumBin())*1.1);
  hTotRotwI->GetZaxis()->SetRangeUser(0,hTotRotwI->GetBinContent(hTotRotwI->GetMaximumBin())*1.1);

 
  TCanvas *crot = new TCanvas("crot","crot");
  crot->cd();

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
  crot->Print(Form("output_files/type_1_rotation_%d.pdf",btr[0]->GetGeo()));
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
  
  // ==================================================================
  // Now look at the Rotation directions............
  //Average_Displacement();
  TH2F *hAveEast1XRot = new TH2F("hAveEast1XRot","East Type 1 X Rot",75,-75,75,75,-75,75);
  TH2F *hAveEast1YRot = new TH2F("hAveEast1YRot","East Type 1 X Rot",75,-75,75,75,-75,75);
  TH2F *hAveWest1XRot = new TH2F("hAveWest1XRot","East Type 1 X Rot",75,-75,75,75,-75,75);
  TH2F *hAveWest1YRot = new TH2F("hAveWest1YRot","East Type 1 X Rot",75,-75,75,75,-75,75);
  //--------------------------------------------------------------------------------------

	cout << "Pre-Ultra loop" << endl;

/*  for(Int_t ii = 1 ; ii < 76 ; ii++){

	cout << ii << "th outer loop" << endl;

      for(Int_t jj = 1 ; jj < 76 ; jj++){
          for(Int_t irun = 0 ; irun < nbeta ; irun++){
    	      btr[irun]->Load_Histograms(bckr[btr[irun]->Bkg_index],0);
              hAveEast1XRot->Fill((ii*2.)-75.,(2.*jj)-75,btr[irun]->hEastType1XRot->GetBinContent(ii,jj));
              hAveEast1YRot->Fill((ii*2.)-75.,(2.*jj)-75,btr[irun]->hEastType1YRot->GetBinContent(ii,jj));
              hAveWest1XRot->Fill((ii*2.)-75.,(2.*jj)-75,btr[irun]->hWestType1XRot->GetBinContent(ii,jj));
              hAveWest1YRot->Fill((ii*2.)-75.,(2.*jj)-75,btr[irun]->hWestType1YRot->GetBinContent(ii,jj));
    	      btr[irun]->Remove_Histograms(bckr[btr[irun]->Bkg_index]);
          }
          hAveEast1XRot->SetBinContent(ii,jj,hAveEast1XRot->GetBinContent(ii,jj)/nbeta);
          hAveEast1YRot->SetBinContent(ii,jj,hAveEast1YRot->GetBinContent(ii,jj)/nbeta);
          hAveWest1XRot->SetBinContent(ii,jj,hAveWest1XRot->GetBinContent(ii,jj)/nbeta);
          hAveWest1YRot->SetBinContent(ii,jj,hAveWest1YRot->GetBinContent(ii,jj)/nbeta);
      }
  }*/

  for(Int_t irun = 0 ; irun < nbeta ; irun++){
    btr[irun]->Load_Histograms(bckr[btr[irun]->Bkg_index],0);
    for(Int_t ii = 1 ; ii < 76 ; ii++){
      for(Int_t jj = 1 ; jj < 76 ; jj++){

        hAveEast1XRot->Fill((ii*2.)-75.,(2.*jj)-75,btr[irun]->hEastType1XRot->GetBinContent(ii,jj));
        hAveEast1YRot->Fill((ii*2.)-75.,(2.*jj)-75,btr[irun]->hEastType1YRot->GetBinContent(ii,jj));
        hAveWest1XRot->Fill((ii*2.)-75.,(2.*jj)-75,btr[irun]->hWestType1XRot->GetBinContent(ii,jj));
        hAveWest1YRot->Fill((ii*2.)-75.,(2.*jj)-75,btr[irun]->hWestType1YRot->GetBinContent(ii,jj));
	}
      }
    btr[irun]->Remove_Histograms(bckr[btr[irun]->Bkg_index]);
  }
  cout << "Average Histos Filled" << endl;
  for(Int_t ii = 1 ; ii < 76 ; ii++){
      for(Int_t jj = 1 ; jj < 76 ; jj++){
          hAveEast1XRot->SetBinContent(ii,jj,hAveEast1XRot->GetBinContent(ii,jj)/nbeta);
          hAveEast1YRot->SetBinContent(ii,jj,hAveEast1YRot->GetBinContent(ii,jj)/nbeta);
          hAveWest1XRot->SetBinContent(ii,jj,hAveWest1XRot->GetBinContent(ii,jj)/nbeta);
          hAveWest1YRot->SetBinContent(ii,jj,hAveWest1YRot->GetBinContent(ii,jj)/nbeta);
      }
  }
  cout << "Average calculated" << endl;
  crot->Clear();

  TArrow *lEastDis[5625],*lWestDis[5625];
  Double_t Xinter,Yinter;
  TLine  *lEastNorm[5625],*lWestNorm[5625];
  hTotRote->Draw("colz");
  //------------------------------------------------------------------------------------------------------------------------
  // Define 2d histograms to hold the line intersection points.
  TH2F *hEastIntersection = new TH2F("hEastIntersection","East Point of Intersection;X(mm);Y(mm);"
                                     ,100,-100,100,100,-50,150);
  TH2F *hWestIntersection = new TH2F("hWestIntersection","West Point of Intersection;X(mm);Y(mm);"
                                     ,100,-100,100,100,-50,150);

  TH1D *hEastAngle = new TH1D("hEastAngle","East Rotation Angle; Angle; Counts",500,0,0.5);

  TH2F *hPosShiftedTot   = new TH2F("hPosShiftedTot"  ,"Shifted Position Difference "   ,150,-25,25,150,-25,25);
  TH2F *hPosUnShiftedTot = new TH2F("hPosUnShiftedTot"," UnShifted Position Difference ",150,-25,25,150,-25,25);

	cout << "Whatever the hell this is?" << endl;

  for(Int_t i = 0 ; i < nbeta ; i++){
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
     for(Int_t j = 1 ; j < hPosShiftedTot->GetNbinsX()*hPosShiftedTot->GetNbinsY() ; j++){
	
        hPosShiftedTot->SetBinContent(j  ,hPosShiftedTot->GetBinContent(j)   + btr[i]->hPosDiffShifted->GetBinContent(j));
        hPosUnShiftedTot->SetBinContent(j,hPosUnShiftedTot->GetBinContent(j) + btr[i]->hPosDiffUnShifted->GetBinContent(j));
     }
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  //-----------------------------------------------------------------------------------------------------------------------
  // Loop through all the normal lines and find the intersection points......
  for(Int_t iline = 0; iline < hAveEast1XRot->GetNbinsX()*hAveEast1XRot->GetNbinsY(); iline++)
          GetNorm(hAveEast1XRot,hAveEast1YRot,lEastNorm[iline],lEastDis[iline],crot,iline);
  // Get intersections and fill the histogram
  for(Int_t iline = 0; iline < 5625 ; iline++)
     for(Int_t jline = 0 ; jline < 5625 ; jline++)
        if(iline != jline)
	  if(GetInterSection(lEastNorm[iline],lEastNorm[jline],Xinter,Yinter) == 1)hEastIntersection->Fill(Xinter,Yinter,1);
  //------------------------------------------------------------------------------------------------------------------------

	cout << "Making plots" << endl;

  // Draw and print the canvas
  gPad->SetGrid(1);
  crot->Print("output_files/rotation_east_type1.pdf");
  crot->Clear();
  crot->Divide(1,2);
  crot->cd(1);
  hPosShiftedTot->Draw("colz");
  hPosShiftedTot->GetXaxis()->SetRangeUser(-25,25);
  hPosShiftedTot->GetYaxis()->SetRangeUser(-25,25);
  crot->cd(2);
  hPosUnShiftedTot->Draw("colz");
  hPosUnShiftedTot->GetXaxis()->SetRangeUser(-25,25);
  hPosUnShiftedTot->GetYaxis()->SetRangeUser(-25,25);
  crot->Print("output_files/position_shifts.pdf");
  crot->Clear();
  hEastIntersection->Draw("colz");
  crot->Print("output_files/intersection_points_east.pdf");
  Get2DGaussianFit(hEastIntersection,Xinter,Yinter);

  for(Int_t iline = 0; iline < 5625 ; iline++){
	Double_t angle = Return_Rotation_Angle(lEastDis[iline],Xinter,Yinter);
        if(angle > 0)hEastAngle->Fill(angle, 1.);
  }
 
  TF1 *fgausf = new TF1("fgausf","gaus",0.01,0.05);
  hEastAngle->Fit("fgausf","RME");

  cout << " East intersection : X = " << Xinter << "\t  Y = " << Yinter << endl;
  cout << " Average Rotation angle is " << fgausf->GetParameter(1) << " +/- " << fgausf->GetParError(1) << endl;//hEastAngle->GetBinCenter(hEastAngle->GetMaximumBin()) << endl;

  for(Int_t iline = 0 ; iline < 5625 ; iline++)
       if(hAveEast1XRot->GetBinContent(iline) > 0 || hAveEast1YRot->GetBinContent(iline) > 0)
          PositionRotate(hEastAngle->GetBinCenter(hEastAngle->GetMaximumBin()),Xinter,Yinter,lEastDis[iline]);

  crot->Clear();
  hEastAngle->Draw();
  crot->Print("output_files/rotation_angle_east.pdf");
  //------------------------------------------------------------------------------------------------------------------------
  // repeat process for the west side
  crot->Clear();
  hTotRotw->Draw("colz");

  for(Int_t iline = 0; iline < hAveWest1XRot->GetNbinsX()*hAveWest1XRot->GetNbinsY(); iline++)
          GetNorm(hAveWest1XRot,hAveWest1YRot,lWestNorm[iline],lWestDis[iline],crot,iline);
   //------------------------------------------------------------------------------------------------------------------------
   for(Int_t iline = 0; iline < 5625 ; iline++)
     for(Int_t jline = 0 ; jline < 5625 ; jline++)
        if(iline != jline)
          if(GetInterSection(lWestNorm[iline],lWestNorm[jline],Xinter,Yinter) == 1)hWestIntersection->Fill(Xinter,Yinter,1);

  gPad->SetGrid(1);
  crot->Print("output_files/rotation_west_type1.pdf"); 
  crot->Clear();
  hWestIntersection->Draw("colz");
  crot->Print("output_files/intersection_points_west.pdf");
  //-------------------------------------------------------------------------------------------------------------------------
}
//---------------------------------------------------------------------------
void CollectTDCCor()
{
  const int rbins = 200;
  Double_t ecount2[rbins][rbins];
  std::vector<Double_t> x1,srate,srater,brate,brater,sut,sute;
 
  //CleanArray(ecount2,rbins,rbins);
  // Set ecount to zero
  for(Int_t ibin = 0 ; ibin < rbins ; ibin++){
    for(Int_t jbin = 0 ; jbin < rbins ; jbin++){
      ecount2[ibin][jbin] = 0.;
    }
  }
  
  // Fill state array with the rate of TDC header footer failures
  for(Int_t i = 0 ; i < nbeta ; i++){
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    x1.push_back(btr[i]->GetRunNumber());
    
    srate.push_back(btr[i]->hETDC_cor->Integral());
    srater.push_back(sqrt(btr[i]->rtime_e*srate[i])/btr[i]->rtime_e);
    
    brate.push_back(btr[i]->hETDC_cor_passed->Integral());
    brater.push_back(sqrt(btr[i]->rtime_e*brate[i])/btr[i]->rtime_e);
    
    sut.push_back(srate[i] + brate[i] + erttt[i]);
    sute.push_back(sqrt( srater[i]*srater[i] + brater[i]*brater[i] + ertet[i]*ertet[i]));
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
 
  // Collect those failures into a Graph
  
  TGraphErrors *gcorr = new TGraphErrors((int)x1.size(),&x1[0],&srate[0],0,&srater[0]);
  TGraphErrors *gcorb = new TGraphErrors((int)x1.size(),&x1[0],&brate[0],0,&brater[0]);

  grstot = new TGraphErrors((int)x1.size(),&x1[0],&sut[0],0,&sute[0]);
  
  ColorGraphic(gcorr,2,20,2,1.,"Events with no Header/Footer","Run Number",
	       "Summed East + West Rate s^{-1}");
  ColorGraphic(gcorb,2,20,2);
  ColorGraphic(grstot,4,20,2);
  ColorGraphic(grtot,3,20,2);
  
  TCanvas *cHBrate = new TCanvas("cHBrate","cHBrate");
  cHBrate->cd();
  // Draw the TDC corrupted events
  gcorr->Draw("AP");
  grtot->Draw("P");
  gcorb->Draw("P");
  //  grstot->Draw("P");
  gcorr->GetYaxis()->SetRangeUser(0,20);
  cHBrate->Print("output_files/badTDC.pdf");

  TCanvas *cTDC = new TCanvas("cTDC","cTDC");
  hEvWtot  = new TH2F("hEvWtot","East vs West No Header/Footer;East(ns);West(ns)"
		      ,rbins,0,rbins,rbins,0,rbins);
  hEvWtotG = new TH2F("hEvWtotG","East vs West Good;East (ns);West(ns)"
		      ,rbins,0,rbins,rbins,0,rbins);
  // Add up the hEvWTDC histogram for all beta runs
  for(Int_t i = 0 ; i <nbeta ; i++){
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
	ecount2[ibin-1][jbin-1] += btr[i]->hEvWTDC_cor->GetBinContent(ibin,jbin);
      }
    }
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
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
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++)
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++)
          ecount2[ibin-1][jbin-1] += btr[i]->hEvWTDC->GetBinContent(ibin,jbin);
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++)
    for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++)
        hEvWtotG->SetBinContent(ibin,jbin,ecount2[ibin-1][jbin-1]);
  
  cTDC->Divide(2,1);
  cTDC->cd(1);
  hEvWtot->Draw("colz");
  cTDC->cd(2);
  hEvWtotG->Draw("colz");

  cTDC->Print("output_files/timingcurve.pdf");
  
  // delete cTDC;
  //delete cHBrate;
}
//------------------------------------------------------------------------------
void CollectAsym()
{
  using namespace TMath;
  
  std::vector<Double_t> Asym_Tot,Asym_Tot_Er,X,W,Temp,X1;

  for(Int_t j = 0 ; j < nbeta ; j++){
    // Reset the statistical weight and the value of the asymmetry
    W.push_back(0.);
    Temp.push_back(0.);
    // Loop through to calculate the error weighted average
    for(Int_t i = 0 ; i < nbeta ; i++){
        W[j]    += 1./Power(btr[i]->fAsym_Run_Err[j],2);
	Temp[j] += W[j]*btr[i]->fAsym_Run[j];
    }
    Asym_Tot.push_back(Temp[j] / W[j]); 
    Asym_Tot_Er.push_back(sqrt(1./W[j]));
    X.push_back((Double_t)(j+1)*10.);
  }
  
  for(Int_t i = 0 ; i < nbeta; i++){
    X1.push_back(btr[i]->GetRunNumber());
     W[i]    = btr[i]->fAsym_Run_Err[20];
     Temp[i] = btr[i]->fAsym_Run[20];
  }
  
  TGraphErrors *gAsym = new TGraphErrors((int)X.size(),&X[0],&Asym_Tot[0],0,&Asym_Tot_Er[0]);
  TGraphErrors *gChck = new TGraphErrors((int)X1.size(),&X1[0],&Temp[0],0,&W[0]);
  ColorGraphic(gAsym,2,20,2);
  ColorGraphic(gChck,4,20,4);
  
  TCanvas *cAsyms = new TCanvas("cAsyms","Asymmetry Calculations");
  cAsyms->Divide(1,2);
  
  cAsyms->cd(1);
  gAsym->Draw("AP");
  
  cAsyms->cd(2);
  gChck->Draw("Ap");

  cAsyms->Print("output_files/Asymmetry_1.pdf");

  delete gAsym;
  delete gChck;
  delete cAsyms;

}
//--------------------------------------------------------------------------------
void Plot_E_Chis()
{
  return; 
  // Plot a random runs background subtracted Beta energy spectrum against the 
  // appropiate PENELOPE prediction.
  
  TCanvas *cEchis = new TCanvas("cEchis","Energy and Backgrounds");
  TRandom3 *xrand = new TRandom3();
  xrand->SetSeed(0);

  Double_t x1[MAXRUNS],y1[MAXRUNS],y1er[MAXRUNS];
  Double_t x1w[MAXRUNS],y1w[MAXRUNS],y1erw[MAXRUNS];
  
  Int_t octn = xrand->Integer(nbeta);
  octn = 1;
  cEchis->Divide(2,1);
  cEchis->cd(1);
  
  cout << "Print Histogram for beta run " << btr[octn]->GetRunNumber() << endl;

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
 
  Double_t Redis[MAXRUNS],Rediser[MAXRUNS];
  // Find Regular SuperRatios
  fsuprat1.open("output_files/super_ratio.txt",fstream::out);

  for(Int_t k = 0 ; k < nbck ; k++){ // ORIGINAL
    bckr[k]->Load_Histograms();
//    Get_Base_Super(k,(char*)"A1",(char*)"A4");  // Original Get_Base_Super
//    Get_Base_Super(k,(char*)"A9",(char*)"A12"); //  uses btr, wrote
//    Get_Base_Super(k,(char*)"B1",(char*)"B4"); //   Get_Base_Super_Back
//    Get_Base_Super(k,(char*)"B9",(char*)"B12"); //  to read bckr

    Get_Base_Super_Back(k,(char*)"A1",(char*)"A4");  // Original Get_Base_Super
    Get_Base_Super_Back(k,(char*)"A9",(char*)"A12"); //  uses btr, wrote
    Get_Base_Super_Back(k,(char*)"B1",(char*)"B4"); //   Get_Base_Super_Back
    Get_Base_Super_Back(k,(char*)"B9",(char*)"B12"); //  to read bckr

    BkgdWestOne[k] = bckr[k]->BetasWest + 1.0;
    BkgdEastOne[k] = bckr[k]->BetasEast + 1.0;
    FirstWBkgTimeReal[k] = bckr[k]->CountTimeWFirstBeta;
    FirstEBkgTimeReal[k] = bckr[k]->CountTimeEFirstBeta;
    FirstWBkgTime[k] = bckr[k]->CountTimeWFirstBeta - 1.0;
    FirstEBkgTime[k] = bckr[k]->CountTimeEFirstBeta - 1.0;
    LastWBkgTime[k] = bckr[k]->CountTimeWBeta + 1.0;
    LastEBkgTime[k] = bckr[k]->CountTimeEBeta + 1.0;
    LastWBkgTimeReal[k] = bckr[k]->CountTimeWBeta;
    LastEBkgTimeReal[k] = bckr[k]->CountTimeEBeta;
    FirstWBkgTimeFive[k] = bckr[k]->CountTimeWFirstBeta + 5.0;
    FirstEBkgTimeFive[k] = bckr[k]->CountTimeEFirstBeta + 5.0;
    LastWBkgTimeFive[k] = bckr[k]->CountTimeWBeta - 5.0;
    LastEBkgTimeFive[k] = bckr[k]->CountTimeEBeta - 5.0;
    FirstWBkgTimeMore[k] = bckr[k]->CountTimeWFirstBeta + 5.1;
    FirstEBkgTimeMore[k] = bckr[k]->CountTimeEFirstBeta + 5.1;
    LastWBkgTimeMore[k] = bckr[k]->CountTimeWBeta - 5.1;
    LastEBkgTimeMore[k] = bckr[k]->CountTimeEBeta - 5.1;
    FirstWBkgTimePThree[k] = bckr[k]->CountTimeWFirstBeta + 0.3;
    FirstEBkgTimePThree[k] = bckr[k]->CountTimeEFirstBeta + 0.3;
    LastWBkgTimePThree[k] = bckr[k]->CountTimeWBeta - 0.3;
    LastEBkgTimePThree[k] = bckr[k]->CountTimeEBeta - 0.3;
    for(Int_t m = 2; m < BkgdWestOne[k]; m++){
    WestClockBkgArray[k][m] = bckr[k]->WestClock[m];
    }
    for(Int_t m = 2; m < BkgdEastOne[k]; m++){
    EastClockBkgArray[k][m] = bckr[k]->EastClock[m];
    }
    bckr[k]->Remove_Histograms();
  }
 
  for(Int_t i = 0 ; i < nbeta ; i++){
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    Get_Base_Super(i,(char*)"A2",(char*)"A5");
    Get_Base_Super(i,(char*)"A7",(char*)"A10");
    Get_Base_Super(i,(char*)"B2",(char*)"B5");
    Get_Base_Super(i,(char*)"B7",(char*)"B10");
    BetasWestOne[i] = btr[i]->BetasWest + 1.0;
    BetasEastOne[i] = btr[i]->BetasEast + 1.0;
    FirstWTimeReal[i] = btr[i]->CountTimeWFirstBeta;
    FirstETimeReal[i] = btr[i]->CountTimeEFirstBeta;
    FirstWTime[i] = btr[i]->CountTimeWFirstBeta - 1.0;
    FirstETime[i] = btr[i]->CountTimeEFirstBeta - 1.0;
    LastWTime[i] = btr[i]->CountTimeWBeta + 1.0;
    LastETime[i] = btr[i]->CountTimeEBeta + 1.0;
    LastWTimeReal[i] = btr[i]->CountTimeWBeta;
    LastETimeReal[i] = btr[i]->CountTimeEBeta;
    FirstWTimeFive[i] = btr[i]->CountTimeWFirstBeta + 5.0;
    FirstETimeFive[i] = btr[i]->CountTimeEFirstBeta + 5.0;
    LastWTimeFive[i] = btr[i]->CountTimeWBeta - 5.0;
    LastETimeFive[i] = btr[i]->CountTimeEBeta - 5.0;
    FirstWTimeMore[i] = btr[i]->CountTimeWFirstBeta + 5.1;
    FirstETimeMore[i] = btr[i]->CountTimeEFirstBeta + 5.1;
    LastWTimeMore[i] = btr[i]->CountTimeWBeta - 5.1;
    LastETimeMore[i] = btr[i]->CountTimeEBeta - 5.1;
    FirstWTimePThree[i] = btr[i]->CountTimeWFirstBeta + 0.3;
    FirstETimePThree[i] = btr[i]->CountTimeEFirstBeta + 0.3;
    LastWTimePThree[i] = btr[i]->CountTimeWBeta - 0.3;
    LastETimePThree[i] = btr[i]->CountTimeEBeta - 0.3;
    for(Int_t j = 2; j < BetasWestOne[i]; j++){
    WestClockArray[i][j] = btr[i]->WestClock[j];
    }
    for(Int_t j = 2; j < BetasEastOne[i]; j++){
    EastClockArray[i][j] = btr[i]->EastClock[j];
    } 
    WestDifference[i] = GetDifference(btr[i]->CountTimeWAll,btr[i]->CountTimeWBeta);
    EastDifference[i] = GetDifference(btr[i]->CountTimeEAll,btr[i]->CountTimeEBeta);
    WestDiffFirst[i] = GetDifference(btr[i]->CountTimeWFirstBeta,btr[i]->CountTimeWFirst);
    EastDiffFirst[i] = GetDifference(btr[i]->CountTimeEFirstBeta,btr[i]->CountTimeEFirst);
    /*LiveTimeBetaWest[i] = GetDifference(btr[i]->CountTimeWBeta,btr[i]->CountTimeWFirstBeta);
    LiveTimeBetaEast[i] = GetDifference(btr[i]->CountTimeEBeta,btr[i]->CountTimeEFirstBeta);
    LiveTimeAllWest[i] = GetDifference(btr[i]->CountTimeWAll,btr[i]->CountTimeWFirst);
    LiveTimeAllEast[i] = GetDifference(btr[i]->CountTimeEAll,btr[i]->CountTimeEFirst);*/
    DiffWestAllBeta[i] = GetDiff(btr[i]->CountTimeWAll,btr[i]->CountTimeWFirst,btr[i]->CountTimeWBeta,btr[i]->CountTimeWFirstBeta);
    DiffEastAllBeta[i] = GetDiff(btr[i]->CountTimeEAll,btr[i]->CountTimeEFirst,btr[i]->CountTimeEBeta,btr[i]->CountTimeEFirstBeta);
    if(btr[i]->flipperOn == 1){
    DiffBetaUp[i] = GetDiff(btr[i]->CountTimeEBeta,btr[i]->CountTimeEFirstBeta,btr[i]->CountTimeWBeta,btr[i]->CountTimeWFirstBeta);
    //cout << "Livetime for up-betas in east in run number " << btr[i]->GetRunNumber() << " is " << LiveTimeBetaEast[i] << " sec." << endl;
    //cout << "Livetime for up-betas in west in run number " << btr[i]->GetRunNumber() << " is " << LiveTimeBetaWest[i] << " sec." << endl;
    }
    if(btr[i]->flipperOn == 0) {
    DiffBetaDown[i] = GetDiff(btr[i]->CountTimeEBeta,btr[i]->CountTimeEFirstBeta,btr[i]->CountTimeWBeta,btr[i]->CountTimeWFirstBeta);
    //cout << "Livetime for down-betas in east in run number " << btr[i]->GetRunNumber() << " is " << LiveTimeBetaEast[i] << " sec." << endl;
    //cout << "Livetime for down-betas in west in run number " << btr[i]->GetRunNumber() << " is " << LiveTimeBetaWest[i] << " sec." << endl;
    }
    /*if(LiveTimeBetaWest[i] > 4200) {
    cout << "Livetime for betas in west in run number " << btr[i]->GetRunNumber() << " is " << LiveTimeBetaWest[i] << " sec." << endl;
    }
    if(LiveTimeAllWest[i] > 4200) {
    cout << "Livetime for all in west in run number " << btr[i]->GetRunNumber() << " is " << LiveTimeAllWest[i] << " sec." << endl;
    }
    if(LiveTimeBetaEast[i] > 4200) {
    cout << "Livetime for betas in east in run number " << btr[i]->GetRunNumber() << " is " << LiveTimeBetaEast[i] << " sec." << endl;
    }
    if(LiveTimeAllEast[i] > 4200) {
    cout << "Livetime for all in east in run number " << btr[i]->GetRunNumber() << " is " << LiveTimeAllEast[i] << " sec." << endl;
    }*/
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }

  fsuprat1.close();
 
  //========================================  

  TCanvas *cSuper = new TCanvas("cSuper","Super Ratios");
  cSuper->Divide(1,2);
  cSuper->cd(1);
  
  TGraphErrors *gSup = new TGraphErrors(Nsuper-1,Xrun,Super,0,SuperE);
  TF1 *fSup = new TF1("fSup","[0]",0.,Nsuper+1.);
  fSup->SetParameter(0,0.044);
  fSup->SetLineWidth(2);
  fSup->SetLineColor(2);
  
  ColorGraphic(gSup,4,20,2,1.,"Raw Super Ratio 0-600 keV (analysis choice A)",
		"Quartet Pair",  "(1 - #sqrt{S}) / (1 + #sqrt{S})");
  
  gSup->Draw("AP");
  gSup->Fit("fSup","RMEQ","");
  gPad->SetGrid();
 
  TGraphErrors *gSupC = new TGraphErrors(Nsuper-1,Xrun,SuperC,0,SuperCE);
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
  
  for(int i = 0 ; i < Nsuper-1 ; i++){
    Redis[i]   = (Super[i] - fSup->GetParameter(0))/SuperE[i];
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
  Double_t SupPos[4],SupPosE[4],Xpos[4];
 
  SupPos[0] = fLSup1->GetParameter(0);
  SupPos[1] = fLSup2->GetParameter(0);
  SupPos[2] = fLSup3->GetParameter(0);
  SupPos[3] = fLSup4->GetParameter(0);
  
  SupPosE[0] = fLSup1->GetParError(0);
  SupPosE[1] = fLSup2->GetParError(0);
  SupPosE[2] = fLSup3->GetParError(0);
  SupPosE[3] = fLSup4->GetParError(0);
  
  for(int i = 0; i < 4; i++)Xpos[i] = i+1;

  TGraphErrors *gPave = new TGraphErrors(4,Xpos,SupPos,0,SupPosE);
  ColorGraphic(gPave,4,20,2);
  gPave->Draw("AP");

  cPave->Print("output_files/positionave.pdf");

  // delete stuff
  delete cl;
  delete gRes;
  delete gSup1;
  delete gSup2;
  delete gSup3;
  delete gSup4;
  delete cSup2d;
  delete cSuper;
  delete fLSup1;
  delete fLSup2;
  delete fLSup3;
  delete fLSup4;  
  delete gPave;
  delete cPave;

  return;
}
//==========================================================================
void Get_Base_Super(Int_t i, char oct1[4],char oct2[4])
{
  
  // Eliminates the potential double counting of Quartets.

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
    }while(strcmp(btr[i+lookahead]->octtype,oct2) && lookahead < 6 && lookahead+1+i < nbeta);
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

    btr[i+lookahead]->Load_Histograms(bckr[btr[i+lookahead]->Bkg_index],0);
      fsuprat1 <<"Super ratios " <<  btr[i]->GetRunNumber() << "\t" << btr[i+lookahead]->GetRunNumber() << endl;
      // if the quartet is found calculate the Asymmetry from the
      // super ratio.
      // This needs to be adjusted to add in the cases where there are
      // multiple files for the same type of octet run.
      // Get Super Ratio

      SuperTiming[Nsuper-1] = GetSuperRatioEight(btr[i]->CountTimeWBeta,btr[i]->CountTimeWFirstBeta,btr[i]->CountTimeEBeta,btr[i]->CountTimeEFirstBeta,btr[i+lookahead]->CountTimeWBeta,btr[i+lookahead]->CountTimeWFirstBeta,btr[i+lookahead]->CountTimeEBeta,btr[i+lookahead]->CountTimeEFirstBeta);
      SuperTimingA[Nsuper-1] = GetSuperRatioEightAss(btr[i]->CountTimeWBeta,btr[i]->CountTimeWFirstBeta,btr[i]->CountTimeEBeta,btr[i]->CountTimeEFirstBeta,btr[i+lookahead]->CountTimeWBeta,btr[i+lookahead]->CountTimeWFirstBeta,btr[i+lookahead]->CountTimeEBeta,btr[i+lookahead]->CountTimeEFirstBeta);
      Super[Nsuper-1]  = GetSuperRatioAsymmetry(btr[i]->hwq,btr[i]->heq,
				     btr[i+lookahead]->hwq,btr[i+lookahead]->heq,nlow,nhigh);
      // Get Super Ratio error
      SuperE[Nsuper-1] = GetSuperRatioError(btr[i]->hwq,btr[i]->heq,
					  btr[i+lookahead]->hwq,btr[i+lookahead]->heq,
                                          btr[i]->rtime_e,btr[i+lookahead]->rtime_e,nlow,nhigh);
      // Get Analysis C Super Ratios
      SuperC[Nsuper-1]  = GetSuperRatioAsymmetry(btr[i]->hwqC,btr[i]->heqC,
				       btr[i+lookahead]->hwqC,btr[i+lookahead]->heqC,nlow,nhigh);
      // Get Super Ratio error
      SuperCE[Nsuper-1] = GetSuperRatioError(btr[i]->hwqC,btr[i]->heqC,
					    btr[i+lookahead]->hwqC,btr[i+lookahead]->heqC,
	                                    btr[i]->rtime_e,btr[i+lookahead]->rtime_e,nlow,nhigh);
      // Get Position Dependent Super Ratios
      
      SuperPos1[Nsuper-1] = GetSuperRatioAsymmetry(btr[i]->hAposw->GetBinContent(1,1),
					  btr[i]->hApose->GetBinContent(1,1),
					  btr[i+lookahead]->hAposw->GetBinContent(1,1),
					  btr[i+lookahead]->hApose->GetBinContent(1,1));
      
      SuperEPos1[Nsuper-1] = GetSuperRatioError(btr[i]->hAposw->GetBinContent(1,1),
					    btr[i]->hApose->GetBinContent(1,1),
					    btr[i+lookahead]->hAposw->GetBinContent(1,1),
					    btr[i+lookahead]->hApose->GetBinContent(1,1),
	                                    btr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      //----------------------------------------------------------------------------------------
      SuperPos2[Nsuper-1] = GetSuperRatioAsymmetry(btr[i]->hAposw->GetBinContent(2,1),
					  btr[i]->hApose->GetBinContent(2,1),
					  btr[i+lookahead]->hAposw->GetBinContent(2,1),
					  btr[i+lookahead]->hApose->GetBinContent(2,1));
      
      SuperEPos2[Nsuper-1] = GetSuperRatioError(btr[i]->hAposw->GetBinContent(2,1),
	                                        btr[i]->hApose->GetBinContent(2,1),
					        btr[i+lookahead]->hAposw->GetBinContent(2,1),
					        btr[i+lookahead]->hApose->GetBinContent(2,1),
						btr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      //----------------------------------------------------------------------------------------
      
      SuperPos3[Nsuper-1] = GetSuperRatioAsymmetry(btr[i]->hAposw->GetBinContent(1,2),
					  btr[i]->hApose->GetBinContent(1,2),
					  btr[i+lookahead]->hAposw->GetBinContent(1,2),
				          btr[i+lookahead]->hApose->GetBinContent(1,2));
      
      SuperEPos3[Nsuper-1] = GetSuperRatioError(btr[i]->hAposw->GetBinContent(1,2),
	                                        btr[i]->hApose->GetBinContent(1,2),
					        btr[i+lookahead]->hAposw->GetBinContent(1,2),
					        btr[i+lookahead]->hApose->GetBinContent(1,2),
						btr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      //--------------------------------------------------------------------------------------
      SuperPos4[Nsuper-1] = GetSuperRatioAsymmetry(btr[i]->hAposw->GetBinContent(2,2),
					  btr[i]->hApose->GetBinContent(2,2),
					  btr[i+lookahead]->hAposw->GetBinContent(2,2),
				          btr[i+lookahead]->hApose->GetBinContent(2,2));
      
      SuperEPos4[Nsuper-1] = GetSuperRatioError(btr[i]->hAposw->GetBinContent(2,2),
	                                        btr[i]->hApose->GetBinContent(2,2),
					        btr[i+lookahead]->hAposw->GetBinContent(2,2),
					        btr[i+lookahead]->hApose->GetBinContent(2,2),
						btr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      NOctet = Nsuper;
    btr[i+lookahead]->Remove_Histograms(bckr[btr[i+lookahead]->Bkg_index]);
      
      // Increment counters.
      Xrun[Nsuper-1] = Nsuper;
      Nsuper++;
    // reset lookahead variable to -4 so that we can look ahead and behind the current run
    // just incase.
    lookahead = -4;
    jstar = 0 ;
     
   }  
  return;
}
//--------------------------------------------------------------

void Get_Base_Super_Back(Int_t i, char oct1[4],char oct2[4])
{
  
  // Eliminates the potential double counting of Quartets.

  Int_t jstar = 0;
  
  if(i>0)
    
     if(!strcmp(bckr[i-1]->octtype,bckr[i]->octtype) && bckr[i-1]->rtime_e > bckr[i]->rtime_e) return;
  if(i < nbeta-1)
    if(!strcmp(bckr[i+1]->octtype,bckr[i]->octtype) && bckr[i+1]->rtime_e > bckr[i]->rtime_e) return;
  // Check for the existance of the quartet pair
  if(!strcmp(bckr[i]->octtype,oct1)){
    // prevent from looking a negative array indices
    if(i < 4) lookahead = 0;
    // Loop around to find the correct run
    do{      
      lookahead++;     
    }while(strcmp(bckr[i+lookahead]->octtype,oct2) && lookahead < 6 && lookahead+1+i < nbeta);
    // if 5 then its a bad data set
    if(lookahead == 6){
      cout << "Matching Run for Run " << bckr[i]->GetRunNumber() << " Not Found!! " << Nsuper << endl;
      if(bckr[i]->GetRunNumber() == 11090){
         do{
	   jstar++;
	 }while(bckr[jstar-1]->GetRunNumber() != 11093);
	 lookahead = jstar - i -1;
      }
      if(bckr[i]->GetRunNumber() == 11084){
	do{
	   jstar++;
	 }while(bckr[jstar-1]->GetRunNumber() != 11088);
	 lookahead = jstar-1 - i;
      }
    }// else {

	bckr[i+lookahead]->Load_Histograms();

      fsuprat1 <<"Super ratios " <<  bckr[i]->GetRunNumber() << "\t" << bckr[i+lookahead]->GetRunNumber() << endl;

      // if the quartet is found calculate the Asymmetry from the
      // super ratio.
      // This needs to be adjusted to add in the cases where there are
      // multiple files for the same type of octet run.
      // Get Super Ratio
      SuperTiming[Nsuper-1] = GetSuperRatioEight(bckr[i]->CountTimeWBeta,bckr[i]->CountTimeWFirstBeta,bckr[i]->CountTimeEBeta,bckr[i]->CountTimeEFirstBeta,bckr[i+lookahead]->CountTimeWBeta,bckr[i+lookahead]->CountTimeWFirstBeta,bckr[i+lookahead]->CountTimeEBeta,bckr[i+lookahead]->CountTimeEFirstBeta);
      SuperTimingA[Nsuper-1] = GetSuperRatioEightAss(bckr[i]->CountTimeWBeta,bckr[i]->CountTimeWFirstBeta,bckr[i]->CountTimeEBeta,bckr[i]->CountTimeEFirstBeta,bckr[i+lookahead]->CountTimeWBeta,bckr[i+lookahead]->CountTimeWFirstBeta,bckr[i+lookahead]->CountTimeEBeta,bckr[i+lookahead]->CountTimeEFirstBeta);
      Super[Nsuper-1]  = GetSuperRatioAsymmetry(bckr[i]->hwq,bckr[i]->heq,
				     bckr[i+lookahead]->hwq,bckr[i+lookahead]->heq,nlow,nhigh);
      // Get Super Ratio error
      SuperE[Nsuper-1] = GetSuperRatioError(bckr[i]->hwq,bckr[i]->heq,
					  bckr[i+lookahead]->hwq,bckr[i+lookahead]->heq,
                                          bckr[i]->rtime_e,bckr[i+lookahead]->rtime_e,nlow,nhigh);
      // Get Analysis C Super Ratios
      SuperC[Nsuper-1]  = GetSuperRatioAsymmetry(bckr[i]->hwqC,bckr[i]->heqC,
				       bckr[i+lookahead]->hwqC,bckr[i+lookahead]->heqC,nlow,nhigh);
      // Get Super Ratio error
      SuperCE[Nsuper-1] = GetSuperRatioError(bckr[i]->hwqC,bckr[i]->heqC,
					    bckr[i+lookahead]->hwqC,bckr[i+lookahead]->heqC,
	                                    bckr[i]->rtime_e,bckr[i+lookahead]->rtime_e,nlow,nhigh);
      // Get Position Dependent Super Ratios
      
      SuperPos1[Nsuper-1] = GetSuperRatioAsymmetry(bckr[i]->hAposw->GetBinContent(1,1),
					  bckr[i]->hApose->GetBinContent(1,1),
					  bckr[i+lookahead]->hAposw->GetBinContent(1,1),
					  bckr[i+lookahead]->hApose->GetBinContent(1,1));
      
      SuperEPos1[Nsuper-1] = GetSuperRatioError(bckr[i]->hAposw->GetBinContent(1,1),
					    bckr[i]->hApose->GetBinContent(1,1),
					    bckr[i+lookahead]->hAposw->GetBinContent(1,1),
					    bckr[i+lookahead]->hApose->GetBinContent(1,1),
	                                    bckr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      //----------------------------------------------------------------------------------------
      SuperPos2[Nsuper-1] = GetSuperRatioAsymmetry(bckr[i]->hAposw->GetBinContent(2,1),
					  bckr[i]->hApose->GetBinContent(2,1),
					  bckr[i+lookahead]->hAposw->GetBinContent(2,1),
					  bckr[i+lookahead]->hApose->GetBinContent(2,1));
      
      SuperEPos2[Nsuper-1] = GetSuperRatioError(bckr[i]->hAposw->GetBinContent(2,1),
	                                        bckr[i]->hApose->GetBinContent(2,1),
					        bckr[i+lookahead]->hAposw->GetBinContent(2,1),
					        bckr[i+lookahead]->hApose->GetBinContent(2,1),
						bckr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      //----------------------------------------------------------------------------------------
      
      SuperPos3[Nsuper-1] = GetSuperRatioAsymmetry(bckr[i]->hAposw->GetBinContent(1,2),
					  bckr[i]->hApose->GetBinContent(1,2),
					  bckr[i+lookahead]->hAposw->GetBinContent(1,2),
				          bckr[i+lookahead]->hApose->GetBinContent(1,2));
      
      SuperEPos3[Nsuper-1] = GetSuperRatioError(bckr[i]->hAposw->GetBinContent(1,2),
	                                        bckr[i]->hApose->GetBinContent(1,2),
					        bckr[i+lookahead]->hAposw->GetBinContent(1,2),
					        bckr[i+lookahead]->hApose->GetBinContent(1,2),
						bckr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      //--------------------------------------------------------------------------------------
      SuperPos4[Nsuper-1] = GetSuperRatioAsymmetry(bckr[i]->hAposw->GetBinContent(2,2),
					  bckr[i]->hApose->GetBinContent(2,2),
					  bckr[i+lookahead]->hAposw->GetBinContent(2,2),
				          bckr[i+lookahead]->hApose->GetBinContent(2,2));
      
      SuperEPos4[Nsuper-1] = GetSuperRatioError(bckr[i]->hAposw->GetBinContent(2,2),
	                                        bckr[i]->hApose->GetBinContent(2,2),
					        bckr[i+lookahead]->hAposw->GetBinContent(2,2),
					        bckr[i+lookahead]->hApose->GetBinContent(2,2),
						bckr[i]->rtime_e,btr[i+lookahead]->rtime_e);
      NOctet = Nsuper;
    
	bckr[i+lookahead]->Remove_Histograms();
  
      // Increment counters.
      Xrun[Nsuper-1] = Nsuper;
      Nsuper++;
    // reset lookahead variable to -4 so that we can look ahead and behind the current run
    // just incase.
    lookahead = -4;
    jstar = 0 ;
     
   }  
  return;
}
//--------------------------------------------------------------
//------------------------------------------------------------------------------------
void Plot_Timing()
{
  
  TH1F *hTimeDiffEndE = new TH1F(Form("hTimeDiffEndE"),"End of Run Difference, East",100,-0.05,0.95);
  TH1F *hTimeDiffEndW = new TH1F(Form("hTimeDiffEndW"),"End of Run Difference, West",100,-0.05,0.95);
  TH1F *hTimeDiffFirstE = new TH1F(Form("hTimeDiffFirstE"),"Begin of Run Difference, East",100,-0.05,6.95);
  TH1F *hTimeDiffFirstW = new TH1F(Form("hTimeDiffFirstW"),"Begin of Run Difference, West",100,-0.05,6.95);
  /*TH1F *hLiveTimeRecW = new TH1F(Form("hLiveTimeRecW"),"Recorded Live Time West",110,0,4400);
  TH1F *hLiveTimeRecE = new TH1F(Form("hLiveTimeRecE"),"Recorded Live Time West",110,0,4400); 
  TH1F *hLiveTimeBetaW = new TH1F(Form("hLiveTimeBetaW"),"True Live Time West (Betas)",110,0,4400);
  TH1F *hLiveTimeBetaE = new TH1F(Form("hLiveTimeBetaW"),"True Live Time East (Betas)",110,0,4400);
  TH1F *hLiveTimeAllW = new TH1F(Form("hLiveTimeAllW"),"True Live Time West (All)",110,0,4400);
  TH1F *hLiveTimeAllE = new TH1F(Form("hLiveTimeAllE"),"True Live Time East (All)",110,0,4400);*/
  TH1F *hLiveTimeDiffW = new TH1F(Form("hLiveTimeDiffW"),"Event/Beta Live Time Diff., W",100,-0.05,6.95);
  TH1F *hLiveTimeDiffE = new TH1F(Form("hLiveTimeDiffE"),"Event/Beta Live Time Diff., E",100,-0.05,6.95);
  TH1F *hDiffBetaUp = new TH1F(Form("hDiffBetaUp"),"Beta Live Time Difference (E-W) Spin Up",100,-10.0,10.0);
  TH1F *hDiffBetaDown = new TH1F(Form("hDiffBetaDown"),"Beta Live Time Difference (E-W) Spin Down",100,-10.0,10.0);

  for(Int_t i = 0; i < nbeta; i++) {

      //Fill histos.
      hTimeDiffEndE->Fill(EastDifference[i]);
      hTimeDiffEndW->Fill(WestDifference[i]);
      hTimeDiffFirstE->Fill(EastDiffFirst[i]);
      hTimeDiffFirstW->Fill(WestDiffFirst[i]);
  }
  TCanvas *timing = new TCanvas("superratiotiming", "superratiotiming",10,10,600,600);
  timing->Divide(2,2);
  timing->cd(1);
  hTimeDiffEndE->Draw("hist");
  hTimeDiffEndE->GetXaxis()->SetTitle("Time (sec)");
  hTimeDiffEndE->GetXaxis()->CenterTitle();
  timing->Update();
  timing->cd(2);
  hTimeDiffEndW->Draw("hist");
  hTimeDiffEndW->GetXaxis()->SetTitle("Time (sec)");
  hTimeDiffEndW->GetXaxis()->CenterTitle();
  timing->Update();
  timing->cd(3);
  hTimeDiffFirstE->Draw("hist");
  hTimeDiffFirstE->GetXaxis()->SetTitle("Time (sec)");
  hTimeDiffFirstE->GetXaxis()->CenterTitle();
  timing->Update();
  timing->cd(4);
  hTimeDiffFirstW->Draw("hist");
  hTimeDiffFirstW->GetXaxis()->SetTitle("Time (sec)");
  hTimeDiffFirstW->GetXaxis()->CenterTitle();
  timing->Print("output_files/timediff.pdf");
  
/*  for(Int_t i = 0; i < nbeta; i++) {

      //Fill histos.
      hLiveTimeRecE->Fill(btr[i]->rtime_e);
      hLiveTimeRecW->Fill(btr[i]->rtime_w);
      hLiveTimeAllE->Fill(LiveTimeAllEast[i]);
      hLiveTimeAllW->Fill(LiveTimeAllWest[i]);
      hLiveTimeBetaE->Fill(LiveTimeBetaEast[i]);
      hLiveTimeBetaW->Fill(LiveTimeBetaWest[i]);
  }

  TCanvas *livetime = new TCanvas("superratiolivetime", "superratiolivetime",10,10,600,600);
  livetime->Divide(2,3);
  livetime->cd(1);
  hLiveTimeAllE->Draw("hist");
  hLiveTimeAllE->GetXaxis()->SetTitle("Time (sec)");
  hLiveTimeAllE->GetXaxis()->CenterTitle();
  livetime->Update();
  livetime->cd(2);
  hLiveTimeAllW->Draw("hist");
  hLiveTimeAllW->GetXaxis()->SetTitle("Time (sec)");
  hLiveTimeAllW->GetXaxis()->CenterTitle();
  livetime->Update();
  livetime->cd(3);
  hLiveTimeBetaE->Draw("hist");
  hLiveTimeBetaE->GetXaxis()->SetTitle("Time (sec)");
  hLiveTimeBetaE->GetXaxis()->CenterTitle();
  livetime->Update();
  livetime->cd(4);
  hLiveTimeBetaW->Draw("hist");
  hLiveTimeBetaW->GetXaxis()->SetTitle("Time (sec)");
  hLiveTimeBetaW->GetXaxis()->CenterTitle();
  livetime->Update();
  livetime->cd(5);
  hLiveTimeRecE->Draw("hist");
  hLiveTimeRecE->GetXaxis()->SetTitle("Time (sec)");
  hLiveTimeRecE->GetXaxis()->CenterTitle();
  livetime->Update();
  livetime->cd(6);
  hLiveTimeRecW->Draw("hist");
  hLiveTimeRecW->GetXaxis()->SetTitle("Time (sec)");
  hLiveTimeRecW->GetXaxis()->CenterTitle();
  livetime->Print("output_files/livetimes.pdf");*/
 
  for(Int_t i = 0; i < nbeta; i++) {

      //Fill histos.
      hLiveTimeDiffW->Fill(DiffWestAllBeta[i]);
      hLiveTimeDiffE->Fill(DiffEastAllBeta[i]);
      if(btr[i]->flipperOn == 1){
      hDiffBetaUp->Fill(DiffBetaUp[i]);
      }
      if(btr[i]->flipperOn == 0){
      hDiffBetaDown->Fill(DiffBetaDown[i]);
      }
  }

  TCanvas *ltimediff = new TCanvas("superratioltimediff", "superratioltimediff",10,10,600,600);
  ltimediff->Divide(2,2);
  ltimediff->cd(1);
  hLiveTimeDiffW->Draw("hist");
  hLiveTimeDiffW->GetXaxis()->SetTitle("Time (sec)");
  hLiveTimeDiffW->GetXaxis()->CenterTitle();
  ltimediff->Update();
  ltimediff->cd(2);
  hLiveTimeDiffE->Draw("hist");
  hLiveTimeDiffE->GetXaxis()->SetTitle("Time (sec)");
  hLiveTimeDiffE->GetXaxis()->CenterTitle();
  ltimediff->Update();
  ltimediff->cd(3);
  hDiffBetaUp->Draw("hist");
  hDiffBetaUp->GetXaxis()->SetTitle("Time (sec)");
  hDiffBetaUp->GetXaxis()->CenterTitle();
  ltimediff->Update();
  ltimediff->cd(4);
  hDiffBetaDown->Draw("hist");
  hDiffBetaDown->GetXaxis()->SetTitle("Time (sec)");
  hDiffBetaDown->GetXaxis()->CenterTitle();
  ltimediff->Print("output_files/ltimediff.pdf");
}
//------------------------------------------------------------------------------------
void Fill_Timing()
{
//	cout << "Total number of events = " << numevents << endl;
  //      cout << "Total number of betas = " << numbetas << endl;
}
//------------------------------------------------------------------------------------
void Collect_Pos()
{
  // Declare things
  const int rbins = 150;
  Float_t rad = 75.;
  Double_t ecount[rbins][rbins],wcount[rbins][rbins];
  Int_t last = 0;
  
  hTotEPos = new TH2F("hTotEPos","East Type 1 ; X_{west} (mm) ; Y_{west} (mm)",
		      rbins,-rad,rad,rbins,-rad,rad);
  
  hTotWPos = new TH2F("hTotWPos","West Type  ; X_{east} (mm) ; Y_{east} (mm)",
		      rbins,-rad,rad,rbins,-rad,rad);
  //-------------------------------------------------------------------
  CleanArray(ecount,rbins,rbins);
  CleanArray(wcount,rbins,rbins);
  
// cout << "Analysis Refill ";
//  if(nbeta > 0)  AnalysisTaskRefill(nbeta, btr, nbetas);
// cout << "Succeeded!" << endl;

  for(Int_t i = 0 ; i < nbeta ; i++){
    if(i == nbeta - 1)last = 1;

    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);

    CollectAllRot(ecount,rbins,rbins,btr[i]->hpe,hTotEPos,last);
    CollectAllRot(wcount,rbins,rbins,btr[i]->hpw,hTotWPos,last);

    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);


  }
  last = 0;
  // Draw the Position of All Events---------------------------------------
  TCanvas *cPosAll = new TCanvas("cPosAll","Position All",800,400);
  cPosAll->Divide(2,1);
  cPosAll->cd(1);
  TEllipse *el1 = new TEllipse(0,0,radcut,radcut);
  TEllipse *el1r = new TEllipse(0,0,radcut-5,radcut-5);
  hTotEPos->GetZaxis()->SetRangeUser(0,1.1*hTotEPos->GetBinContent(hTotEPos->GetMaximumBin()));
  hTotWPos->GetZaxis()->SetRangeUser(0,1.1*hTotWPos->GetBinContent(hTotWPos->GetMaximumBin()));
  hTotEPos->Draw("colz");
  el1->SetLineColor(2);
  el1r->SetLineColor(2);
  el1->SetFillStyle(0);
  el1r->SetLineStyle(2);
  el1r->SetFillStyle(0);
  el1->Draw();
  el1r->Draw();
  cPosAll->cd(2);
  TEllipse *el2 = new TEllipse(0,0,radcut,radcut);
  TEllipse *el2r = new TEllipse(0,0,radcut-5,radcut-5);

  hTotWPos->Draw("colz");
  el2->SetLineColor(2);
  el2r->SetLineColor(2);
  el2r->SetLineStyle(2);
  el2->SetFillStyle(0); 
  el2r->SetFillStyle(0);
  el2->Draw();
  el2r->Draw();
  //------------------------------------------------------------------------
  // Do a height variation check
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
  
  cPosAll->Print(Form("output_files/position_all_%d.pdf",btr[0]->GetGeo()));
}
//--------------------------------------------------------------------
void Plot_23_Diff()
{
  TCanvas *cboring = new TCanvas("cboring","");
  cboring->Divide(2,1);
  cboring->cd(1);
  btr[1]->hPsDfET23->Draw("colz");
  
  std::vector<Double_t> nxr,nt2f;
  
  for(int i = 0 ; i < nbeta ; i++){
    nxr.push_back(btr[i]->heq->Integral()+btr[i]->hwq->Integral());
    nt2f.push_back(btr[i]->T2_frac);
  }

  cboring->cd(2);
  TGraph *gT2 = new TGraph((int)nxr.size(),&nxr[0],&nt2f[0]);
  gT2->SetMarkerColor(2);
  gT2->SetMarkerStyle(20);
  gT2->Draw("AP");

  cboring->Print("output_files/23_diff.pdf");
  
  delete gT2;
  delete cboring;
}
//-------------------------------------------------------------------------  
void Plot_MonRate()
{
  
  Double_t xt[MAXRUNS],y[MAXRUNS],yer[MAXRUNS];
  Double_t ek[MAXRUNS],wk[MAXRUNS];
  Double_t eker[MAXRUNS],wker[MAXRUNS];
  Double_t emean[MAXRUNS],wmean[MAXRUNS];
  
  TH1F *hEastEndPoint = new TH1F("hEastEndPoint","East End Point ; E_{vis} (keV) ; Counts",
				 100,700,900);
  TH1F *hWestEndPoint = new TH1F("hWestEndPoint","West End Point ; E_{vis} (keV) ; Counts",
				 100,700,900);				  
  
  for(int i = 0 ; i < nbeta ; i++){
    //--------------------------------------------------------------------------------------
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    xt[i]   = btr[i]->GetRunNumber();
    if(btr[i]->Mon_Rate > 0 ){
	    y[i]    = (btr[i]->heq->Integral() + btr[i]->hwq->Integral())/btr[i]->Mon_Rate;
            yer[i]  = y[i]*sqrt( 1./(btr[i]->rtime_e*btr[i]->heq->Integral()) 
                               + 1./(btr[i]->rtime_w*btr[i]->hwq->Integral())
                               + 1./(btr[i]->rtime_w*btr[i]->Mon_Rate));
    } else {
            y[i]    = 0.;
	    yer[i]  = 1.;
    }
    //-----------------------------------------------------------------------------------
    ek[i]   = btr[i]->E_Endpoint;
    wk[i]   = btr[i]->W_Endpoint;
    eker[i] = btr[i]->E_EndError;
    wker[i] = btr[i]->W_EndError;
    hEastEndPoint->Fill(btr[i]->E_Endpoint);
    hWestEndPoint->Fill(btr[i]->W_Endpoint);
    emean[i] = btr[i]->heq->GetMean();
    wmean[i] = btr[i]->hwq->GetMean();
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  
  TGraphErrors *gBdM = new TGraphErrors(nbeta,xt,y,0,yer);
  ColorGraphic(gBdM,2,20,1,1.,"Normalized #beta Rate","Run Number","#beta -Rate / Monitor 2 ");
 
  TCanvas  *cMon = new TCanvas("cMon","Beta-Rate Divided by Monitor 2");
  cMon->cd();
  gBdM->Draw("AP");
   
  cMon->Print("output_files/monitors.pdf"); 

  //delete gBdM;
  // delete cMon;

  TCanvas *cHom = new TCanvas("cHom","Stuff");
  cHom->cd();

  TGraphErrors *gKurieE = new TGraphErrors(nbeta,xt,ek,0,eker);
  TGraphErrors *gKurieW = new TGraphErrors(nbeta,xt,wk,0,wker);

  TGraph *gMeanE = new TGraph(nbeta,xt,emean);
  TGraph *gMeanW = new TGraph(nbeta,xt,wmean);
  
  ColorGraphic(gKurieE,2,20,1,1.,"East End Point","Run Number","Extracted End Point (keV)");
  ColorGraphic(gKurieW,4,20,1);
  ColorGraphic(gMeanE,2,20,1,1.,"Mean Energy","Run Number","Mean Energy (keV)");
  ColorGraphic(gMeanW,4,20,1);
  
  cHom->Divide(2,2);
  cHom->cd(1);
  gKurieE->Draw("AP");
  gKurieE->GetYaxis()->SetRangeUser(300,1000);
  gKurieW->Draw("P");
  
  cHom->cd(2);
  
  ColorGraphic(hEastEndPoint,2,20,2);
  ColorGraphic(hWestEndPoint,4,20,2);
  hEastEndPoint->GetYaxis()->SetTitle("N Counts");
  hEastEndPoint->GetXaxis()->SetTitle("End Point (keV)");
  hEastEndPoint->Draw();
  hWestEndPoint->Draw("same");
 
  cHom->cd(3);

  gMeanE->Draw("AP");
  gMeanW->Draw("P");

  cHom->Print("output_files/MeanEnergies.pdf");
  /*
  delete gKurieW;
  delete gKurieE;
  delete gMeanE;
  delete gMeanW;
  delete hEastEndPoint;
  delete hWestEndPoint;
  delete cHom;*/
  
}
//============================================================================================
void Collect_Octets()
{
   
  gStyle->SetOptStat(kTRUE); 

  using namespace TMath;
  
  TF1 *flA = new TF1("flA","[0]",0,noct);
  TF1 *flB = new TF1("flB","[0]",0,noct);

  flA->SetLineColor(2);
  flB->SetLineColor(4);

  // Analysis choice arrays..
  
  Asym_t A_sum_A,A_sum_B,A_multi_A,A_multi_B;
  vector<Asym_t> A_Aves (9), A_Aves_s (9);
  Asym_t A_quat,A_aveB,A_aveC,A_aveD,A_aveE,A_aveF,A_aveG,A_aveH,A_aveI,A_aves;
  Asym_t A_avesB,A_avesC,A_avesD,A_avesE,A_avesF,A_avesG,A_avesH,A_avesI;
  Asym_t A_ave,A_super;
  Double_t A_ave_tot[9],A_ave_toter[9],A_aves_tot[9],A_aves_toter[9];

  Int_t NoctC   = 0;
  Int_t FullOct = 0;
  Int_t NQuat   = 0;

  Double_t TotalCounts = 0.;
  vector<Double_t> CurrentCounts;
  vector<Double_t> OctNumber,OctDate;

  fstream foct,fsup,fquart;

  foct.open(Form("output_files/oct_list_%d.txt",btr[0]->GetGeo()),fstream::out);
  fsup.open(Form("output_files/super_list_%d.txt",btr[0]->GetGeo()),fstream::out);
  fquart.open(Form("output_files/quart_list_%d.txt",btr[0]->GetGeo()),fstream::out);

  foct <<   "Noct\t Runs \t A " << endl;
  fquart << "Noct\t Runs \t A quartet A \t\t A quartet B" << endl;
  fsup <<   "Noct\t Runs \t A sup A1 \t\t A super A2 \t\t Asuper B1 \t\t Asuper B2"<<endl;  
 
  TH1F *hTimeratio = new TH1F(Form("hTimeratio"),"Beta Superratio - 1",100,-0.002,0.002);
  TH1F *hTimeratioall = new TH1F(Form("hTimeratioall"),"Beta Superratio (All)  - 1",100,-0.002,0.002);
  TH1F *hSupsim = new TH1F(Form("hSupsim"),"Beta Superratio - 1",100,-0.002,0.002);  
  TH1F *hSupalsim = new TH1F(Form("hSupalsim"),"Beta Superratio - 1",100,-0.002,0.002);
  TH1F *hLWSimalEvents = new TH1F(Form("hLWSimEvt"),"NBetaEventsW last 5 s",100,0,200); 
  TH1F *hLESimalEvents = new TH1F(Form("hLESimEvt"),"NBetaEventsE last 5 s",100,0,200);
  TH1F *hFWSimalEvents = new TH1F(Form("hFWSimEvt"),"NBetaEventsW 1st 5 s",100,0,200);
  TH1F *hFESimalEvents = new TH1F(Form("hFESimEvt"),"NBetaEventsE 1st 5 s",100,0,200); 
  TH1F *hSupsims = new TH1F(Form("hSupsims"),"Beta Superratio - 1",100,-0.002,0.002);
  TH1F *hSupalsims = new TH1F(Form("hSupalsims"),"Beta Superratio - 1",100,-0.002,0.002);
  TH1F *hLWSimalEventss = new TH1F(Form("hLWSimEvts"),"NBetaEventsW last 5 s",100,0,200);
  TH1F *hLESimalEventss = new TH1F(Form("hLESimEvts"),"NBetaEventsE last 5 s",100,0,200);
  TH1F *hFWSimalEventss = new TH1F(Form("hFWSimEvts"),"NBetaEventsW 1st 5 s",100,0,200);
  TH1F *hFESimalEventss = new TH1F(Form("hFESimEvts"),"NBetaEventsE 1st 5 s",100,0,200);
  TH1F *hFirstPThreeE = new TH1F(Form("hFirstPThreeE"),"NBetaEvents East 1st 0.3 sec",16,-1,15);
  TH1F *hFirstPThreeW = new TH1F(Form("hFirstPThreeW"),"NBetaEvents West 1st 0.3 sec",16,-1,15);
  TH1F *hFirstPThreeEBkg = new TH1F(Form("hFirstPThreeEBkg"),"NBkgEvents East 1st 0.3 sec",10,-1,9);
  TH1F *hFirstPThreeWBkg = new TH1F(Form("hFirstPThreeWBkg"),"NBkgEvents West 1st 0.3 sec",10,-1,9);
  TH1F *hLastPThreeE = new TH1F(Form("hLastPThreeE"),"NBetaEvents East last 0.3 sec",16,-1,15);
  TH1F *hLastPThreeW = new TH1F(Form("hLastPThreeW"),"NBetaEvents West last 0.3 sec",16,-1,15);
  TH1F *hLastPThreeEBkg = new TH1F(Form("hLastPThreeEBkg"),"NBkgEvents East last 0.3 sec",10,-1,9);
  TH1F *hLastPThreeWBkg = new TH1F(Form("hLastPThreeWBkg"),"NBkgEvents West last 0.3 sec",10,-1,9);
  TH1F *hFirstFiveE = new TH1F(Form("hFirstFiveE"),"NBetaEventsE 1st 5 s",100,0,200);
  TH1F *hFirstFiveW = new TH1F(Form("hFirstFiveW"),"NBetaEventsW 1st 5 s",100,0,200);
  TH1F *hFirstFiveEBkg = new TH1F(Form("hFirstFiveEBkg"),"NBkgEvents East 1st 5.0 sec",16,-1,15);
  TH1F *hFirstFiveWBkg = new TH1F(Form("hFirstFiveWBkg"),"NBkgEvents West 1st 5.0 sec",16,-1,15);
  TH1F *hLastFiveE = new TH1F(Form("hLastFiveE"),"NBetaEventsE last 5 s",100,0,200);
  TH1F *hLastFiveW = new TH1F(Form("hLastFiveW"),"NBetaEventsW last 5 s",100,0,200);
  TH1F *hLastFiveEBkg = new TH1F(Form("hLastFiveEBkg"),"NBkgEvents East last 5.0 sec",16,-1,15);
  TH1F *hLastFiveWBkg = new TH1F(Form("hLastFiveWBkg"),"NBkgEvents West last 5.0 sec",16,-1,15);
  
  TH1F *hLastFiveEast[2000];
  TH1F *hLastFiveWest[2000];
  TH1F *hLastFiveBkgEast[2000];
  TH1F *hLastFiveBkgWest[2000];
  TH1F *hFirstFiveEast[2000];
  TH1F *hFirstFiveWest[2000];
  TH1F *hFirstFiveBkgEast[2000];
  TH1F *hFirstFiveBkgWest[2000];
// 
  char hnameeast[2000];
  char htitleeast[2000];
  char hnamewest[2000];
  char htitlewest[2000];
  /*char cname[2000];
  char ctitle[2000];
  char ftitle[2000];
  char eastevnt[2000];
  char westevnt[2000];
  char eastevntbkg[2000];
  char westevntbkg[2000];*/
  char hnameeastf[2000];
  char htitleeastf[2000];
  char hnamewestf[2000];
  char htitlewestf[2000];
  /*char eastevntf[2000];
  char westevntf[2000];
  char eastevntbkgf[2000];
  char westevntbkgf[2000];
  TCanvas *lastfive[2000];
  TPaveText *lwb[2000];
  TPaveText *leb[2000];
  TPaveText *le[2000];
  TPaveText *lw[2000];
  TPaveText *fwb[2000];
  TPaveText *feb[2000];
  TPaveText *fe[2000];
  TPaveText *fw[2000];*/
  char hnamebkgeast[2000];
  char htitlebkgeast[2000];
  char hnamebkgwest[2000];
  char htitlebkgwest[2000];
  char hnamebkgeastf[2000];
  char htitlebkgeastf[2000];
  char hnamebkgwestf[2000];
  char htitlebkgwestf[2000];  

  for(Int_t i = 0; i < nbeta; i++){
  sprintf(hnameeast,"hLastFiveEast%d",i);
  sprintf(htitleeast,"LastFiveSecEast%d",i);
  sprintf(hnamewest,"hLastFiveWest%d",i);
  sprintf(htitlewest,"LastFiveSecWest%d",i);
  sprintf(hnameeastf,"hFirstFiveEast%d",i);
  sprintf(htitleeastf,"FirstFiveSecEast%d",i);
  sprintf(hnamewestf,"hFirstFiveWest%d",i);
  sprintf(htitlewestf,"FirstFiveSecWest%d",i);
  hLastFiveEast[i] = new TH1F(hnameeast,htitleeast,305,LastETimeMore[i],LastETime[i]);
  hLastFiveWest[i] = new TH1F(hnamewest,htitlewest,305,LastWTimeMore[i],LastWTime[i]);
  hFirstFiveEast[i] = new TH1F(hnameeastf,htitleeastf,305,FirstETime[i],FirstETimeMore[i]);
  hFirstFiveWest[i] = new TH1F(hnamewestf,htitlewestf,305,FirstWTime[i],FirstWTimeMore[i]);
  for(Int_t j = 2; j < BetasEastOne[i]; j++) {
  hLastFiveEast[i]->Fill(EastClockArray[i][j]);
  hFirstFiveEast[i]->Fill(EastClockArray[i][j]);
  }
  hFirstFiveEast[i]->Fill(FirstETimeReal[i]);
  for(Int_t j = 2; j < BetasWestOne[i]; j++) {
  hLastFiveWest[i]->Fill(WestClockArray[i][j]);
  hFirstFiveWest[i]->Fill(WestClockArray[i][j]);
  }
  hFirstFiveWest[i]->Fill(FirstWTimeReal[i]);
  lowbine[i] = hLastFiveEast[i]->FindBin(LastETimeFive[i]);
  lowbinw[i] = hLastFiveWest[i]->FindBin(LastWTimeFive[i]);
  lowbineast[i] = hLastFiveEast[i]->FindBin(LastETimePThree[i]);
  lowbinwest[i] = hLastFiveWest[i]->FindBin(LastWTimePThree[i]);
  highbineast[i] = hLastFiveEast[i]->FindBin(LastETimeReal[i]);
  highbinwest[i] = hLastFiveWest[i]->FindBin(LastWTimeReal[i]);
  SumLastPThreeE[i] = hLastFiveEast[i]->Integral(lowbineast[i],highbineast[i]);
  SumLastPThreeW[i] = hLastFiveWest[i]->Integral(lowbinwest[i],highbinwest[i]);
  SumLastFiveE[i] = hLastFiveEast[i]->Integral(lowbine[i],highbineast[i]);
  SumLastFiveW[i] = hLastFiveWest[i]->Integral(lowbinw[i],highbineast[i]);
  RateLastFiveE[i] = SumLastFiveE[i]/(5.);
  RateLastFiveW[i] = SumLastFiveW[i]/(5.);
  ProbLastFiveE[i] = RateLastFiveE[i]*(0.000001);
  ProbLastFiveW[i] = RateLastFiveW[i]*(0.000001);
  lowbineastf[i] = hFirstFiveEast[i]->FindBin(FirstETimeReal[i]);
  lowbinwestf[i] = hFirstFiveWest[i]->FindBin(FirstWTimeReal[i]);
  highbineastf[i] = hFirstFiveEast[i]->FindBin(FirstETimePThree[i]);
  highbinwestf[i] = hFirstFiveWest[i]->FindBin(FirstWTimePThree[i]);
  highbinef[i] = hFirstFiveEast[i]->FindBin(FirstETimeFive[i]);
  highbinwf[i] = hFirstFiveWest[i]->FindBin(FirstWTimeFive[i]);
  SumFirstPThreeE[i] = hFirstFiveEast[i]->Integral(lowbineastf[i],highbineastf[i]);
  SumFirstPThreeW[i] = hFirstFiveWest[i]->Integral(lowbinwestf[i],highbinwestf[i]);
  SumFirstFiveE[i] = hFirstFiveEast[i]->Integral(lowbineastf[i],highbinef[i]);
  SumFirstFiveW[i] = hFirstFiveWest[i]->Integral(lowbinwestf[i],highbinwf[i]);
  RateFirstFiveE[i] = SumFirstFiveE[i]/(5.);
  RateFirstFiveW[i] = SumFirstFiveW[i]/(5.);
  ProbFirstFiveE[i] = RateFirstFiveE[i]*(0.000001);
  ProbFirstFiveW[i] = RateFirstFiveW[i]*(0.000001);
  SumLastFiveES[i] = SumLastFiveE[i];
  SumLastFiveWS[i] = SumLastFiveW[i];
  SumFirstFiveES[i] = SumFirstFiveE[i];
  SumFirstFiveWS[i] = SumFirstFiveW[i];
  /*sprintf(cname,"fivesecstudy%d",i);
  sprintf(ctitle,"fivesecstudy%d",i);
  lastfive[i] = new TCanvas(cname,ctitle,10,10,600,600);
  lastfive[i]->Divide(2,4);
  lastfive[i]->cd(1);
  hLastFiveEast[i]->Draw("hist");
  hLastFiveEast[i]->GetXaxis()->SetTitle("Time (sec)");
  hLastFiveEast[i]->GetXaxis()->CenterTitle();
  sprintf(eastevnt,"%.1f",SumLastPThreeE[i]);
  le[i] = new TPaveText(.73,.62,0.98,0.74,"NDC");
  le[i]->AddText("E Det Evnt (last 0.3 s):");
  le[i]->AddText(eastevnt);
  le[i]->SetTextSize(0.04);
  le[i]->Draw();
  lastfive[i]->Update();
  lastfive[i]->cd(2);
  hLastFiveWest[i]->Draw("hist");
  hLastFiveWest[i]->GetXaxis()->SetTitle("Time (sec)");
  hLastFiveWest[i]->GetXaxis()->CenterTitle();
  sprintf(westevnt,"%.1f",SumLastPThreeW[i]);
  lw[i] = new TPaveText(.73,.62,0.98,0.74,"NDC");
  lw[i]->AddText("W Det Evnt (last 0.3 s):");
  lw[i]->AddText(westevnt);
  lw[i]->SetTextSize(0.04);
  lw[i]->Draw();
  lastfive[i]->Update();
  lastfive[i]->cd(3);
  hFirstFiveEast[i]->Draw("hist");
  hFirstFiveEast[i]->GetXaxis()->SetTitle("Time (sec)");
  hFirstFiveEast[i]->GetXaxis()->CenterTitle();
  sprintf(eastevntf,"%.1f",SumFirstPThreeE[i]);
  fe[i] = new TPaveText(.73,.62,0.98,0.74,"NDC");
  fe[i]->AddText("E Det Evnt (1st 0.3 s):");
  fe[i]->AddText(eastevntf);
  fe[i]->SetTextSize(0.04);
  fe[i]->Draw();
  lastfive[i]->Update();
  lastfive[i]->cd(4);
  hFirstFiveWest[i]->Draw("hist");
  hFirstFiveWest[i]->GetXaxis()->SetTitle("Time (sec)");
  hFirstFiveWest[i]->GetXaxis()->CenterTitle();
  sprintf(westevntf,"%.1f",SumFirstPThreeW[i]);
  fw[i] = new TPaveText(.73,.62,0.98,0.74,"NDC");
  fw[i]->AddText("W Det Evnt (1st 0.3 s):");
  fw[i]->AddText(westevntf);
  fw[i]->SetTextSize(0.04);
  fw[i]->Draw();*/
  }

  for(Int_t k = 0; k < nbck; k++){
  sprintf(hnamebkgeast,"hLstFiveBkgEst%d",k);
  sprintf(htitlebkgeast,"LstFiveSecBkgEst%d",k);
  sprintf(hnamebkgwest,"hLstFiveBkgWst%d",k);
  sprintf(htitlebkgwest,"LstFiveSecBkgWst%d",k);
  sprintf(hnamebkgeastf,"hFrstFiveBkgEst%d",k);
  sprintf(htitlebkgeastf,"FrstFiveSecBkgEst%d",k);
  sprintf(hnamebkgwestf,"hFrstFiveBkgWst%d",k);
  sprintf(htitlebkgwestf,"FrstFiveSecBkgWst%d",k);
  hLastFiveBkgEast[k] = new TH1F(hnamebkgeast,htitlebkgeast,305,LastEBkgTimeMore[k],LastEBkgTime[k]);
  hLastFiveBkgWest[k] = new TH1F(hnamebkgwest,htitlebkgwest,305,LastWBkgTimeMore[k],LastWBkgTime[k]);
  hFirstFiveBkgEast[k] = new TH1F(hnamebkgeastf,htitlebkgeastf,305,FirstEBkgTime[k],FirstEBkgTimeMore[k]);
  hFirstFiveBkgWest[k] = new TH1F(hnamebkgwestf,htitlebkgwestf,305,FirstWBkgTime[k],FirstWBkgTimeMore[k]);
  for(Int_t m = 2; m < BkgdEastOne[k]; m++) {
  hLastFiveBkgEast[k]->Fill(EastClockBkgArray[k][m]);
  hFirstFiveBkgEast[k]->Fill(EastClockBkgArray[k][m]);
  }
  hFirstFiveBkgEast[k]->Fill(FirstEBkgTimeReal[k]);
  for(Int_t m = 2; m < BkgdWestOne[k]; m++) {
  hLastFiveBkgWest[k]->Fill(WestClockBkgArray[k][m]);
  hFirstFiveBkgWest[k]->Fill(WestClockBkgArray[k][m]);
  }
  hFirstFiveBkgWest[k]->Fill(FirstWBkgTimeReal[k]);
  lowbinebkg[k] = hLastFiveBkgEast[k]->FindBin(LastEBkgTimeFive[k]);
  lowbinwbkg[k] = hLastFiveBkgWest[k]->FindBin(LastWBkgTimeFive[k]);
  lowbineastbkg[k] = hLastFiveBkgEast[k]->FindBin(LastEBkgTimePThree[k]);
  lowbinwestbkg[k] = hLastFiveBkgWest[k]->FindBin(LastWBkgTimePThree[k]);
  highbineastbkg[k] = hLastFiveBkgEast[k]->FindBin(LastEBkgTimeReal[k]);
  highbinwestbkg[k] = hLastFiveBkgWest[k]->FindBin(LastWBkgTimeReal[k]);
  SumLastPThreeEBkg[k] = hLastFiveBkgEast[k]->Integral(lowbineastbkg[k],highbineastbkg[k]);
  SumLastPThreeWBkg[k] = hLastFiveBkgWest[k]->Integral(lowbinwestbkg[k],highbinwestbkg[k]);
  SumLastFiveEBkg[k] = hLastFiveBkgEast[k]->Integral(lowbinebkg[k],highbineastbkg[k]);
  SumLastFiveWBkg[k] = hLastFiveBkgWest[k]->Integral(lowbinwbkg[k],highbinwestbkg[k]);
  lowbineastbkgf[k] = hFirstFiveBkgEast[k]->FindBin(FirstEBkgTimeReal[k]);
  lowbinwestbkgf[k] = hFirstFiveBkgWest[k]->FindBin(FirstWBkgTimeReal[k]);
  highbineastbkgf[k] = hFirstFiveBkgEast[k]->FindBin(FirstEBkgTimePThree[k]);
  highbinwestbkgf[k] = hFirstFiveBkgWest[k]->FindBin(FirstWBkgTimePThree[k]);
  highbinebkgf[k] = hFirstFiveBkgEast[k]->FindBin(FirstEBkgTimeFive[k]);
  highbinwbkgf[k] = hFirstFiveBkgWest[k]->FindBin(FirstWBkgTimeFive[k]);
  SumFirstPThreeEBkg[k] = hFirstFiveBkgEast[k]->Integral(lowbineastbkgf[k],highbineastbkgf[k]);
  SumFirstPThreeWBkg[k] = hFirstFiveBkgWest[k]->Integral(lowbinwestbkgf[k],highbinwestbkgf[k]);
  SumFirstFiveEBkg[k] = hFirstFiveBkgEast[k]->Integral(lowbineastbkgf[k],highbinebkgf[k]);
  SumFirstFiveWBkg[k] = hFirstFiveBkgWest[k]->Integral(lowbinwestbkgf[k],highbinwbkgf[k]);
  /*sprintf(ftitle,"output_files/lastfive%d.pdf",k);
  lastfive[k]->Update();
  lastfive[k]->cd(5);
  hLastFiveBkgEast[k]->SetLineColor(6);
  hLastFiveBkgEast[k]->Draw("hist");
  hLastFiveBkgEast[k]->GetXaxis()->SetTitle("Time (sec)");
  hLastFiveBkgEast[k]->GetXaxis()->CenterTitle();
  sprintf(eastevntbkg,"%.1f",SumLastPThreeEBkg[k]);
  leb[k] = new TPaveText(.73,.62,0.98,0.74,"NDC");
  leb[k]->AddText("E Det Bkg Evnt (last 0.3 s):");
  leb[k]->AddText(eastevntbkg);
  leb[k]->Draw();
  lastfive[k]->Update();
  lastfive[k]->cd(6);
  hLastFiveBkgWest[k]->SetLineColor(6);
  hLastFiveBkgWest[k]->Draw("hist");
  hLastFiveBkgWest[k]->GetXaxis()->SetTitle("Time (sec)");
  hLastFiveBkgWest[k]->GetXaxis()->CenterTitle();
  sprintf(westevntbkg,"%.1f",SumLastPThreeWBkg[k]);
  lwb[k] = new TPaveText(.73,.62,0.98,0.74,"NDC");
  lwb[k]->AddText("W Det Bkg Evnt (last 0.3 s):");
  lwb[k]->AddText(westevntbkg);
  lwb[k]->SetTextSize(0.04);
  lwb[k]->Draw();
  lastfive[k]->Update();
  lastfive[k]->cd(7);
  hFirstFiveBkgEast[k]->SetLineColor(6);
  hFirstFiveBkgEast[k]->Draw("hist");
  hFirstFiveBkgEast[k]->GetXaxis()->SetTitle("Time (sec)");
  hFirstFiveBkgEast[k]->GetXaxis()->CenterTitle();
  sprintf(eastevntbkgf,"%.1f",SumFirstPThreeEBkg[k]);
  feb[k] = new TPaveText(.73,.62,0.98,0.74,"NDC");
  feb[k]->AddText("E Det Bkg Evnt (1st 0.3 s):");
  feb[k]->AddText(eastevntbkgf);
  feb[k]->SetTextSize(0.04);
  feb[k]->Draw();
  lastfive[k]->Update();
  lastfive[k]->cd(8);
  hFirstFiveBkgWest[k]->SetLineColor(6);
  hFirstFiveBkgWest[k]->Draw("hist");
  hFirstFiveBkgWest[k]->GetXaxis()->SetTitle("Time (sec)");
  hFirstFiveBkgWest[k]->GetXaxis()->CenterTitle();
  sprintf(westevntbkgf,"%.1f",SumFirstPThreeWBkg[k]);
  fwb[k] = new TPaveText(.73,.62,0.98,0.74,"NDC");
  fwb[k]->AddText("W Det Bkg Evnt (1st 0.3 s):");
  fwb[k]->AddText(westevntbkgf);
  fwb[k]->SetTextSize(0.04);
  fwb[k]->Draw();
  lastfive[k]->Print(ftitle);*/
  }

  for(Int_t i = 0; i < nbeta; i++) {
  hFirstPThreeE->Fill(SumFirstPThreeE[i]);
  hFirstPThreeW->Fill(SumFirstPThreeW[i]);
  hLastPThreeE->Fill(SumLastPThreeE[i]);
  hLastPThreeW->Fill(SumLastPThreeW[i]);
  hFirstFiveE->Fill(SumFirstFiveE[i]);
  hFirstFiveW->Fill(SumFirstFiveW[i]);
  hLastFiveE->Fill(SumLastFiveE[i]);
  hLastFiveW->Fill(SumLastFiveW[i]);
  }

  for(Int_t k = 0; k < nbck; k++) {
  hFirstPThreeEBkg->Fill(SumFirstPThreeEBkg[k]);
  hFirstPThreeWBkg->Fill(SumFirstPThreeWBkg[k]);
  hLastPThreeEBkg->Fill(SumLastPThreeEBkg[k]);
  hLastPThreeWBkg->Fill(SumLastPThreeWBkg[k]);
  hFirstFiveEBkg->Fill(SumFirstFiveEBkg[k]);
  hFirstFiveWBkg->Fill(SumFirstFiveWBkg[k]);
  hLastFiveEBkg->Fill(SumLastFiveEBkg[k]);
  hLastFiveWBkg->Fill(SumLastFiveWBkg[k]);
  }

  TCanvas *pthreef = new TCanvas("pthreef", "pthreef",10,10,600,600);
  pthreef->Divide(2,2);
  pthreef->cd(1);
  hFirstPThreeE->Draw("hist");
  hFirstPThreeE->GetXaxis()->SetTitle("NBetaEvents (1st 0.3 sec)");
  hFirstPThreeE->GetXaxis()->CenterTitle();
  hFirstPThreeE->GetYaxis()->SetTitle("NBetaRuns");
  hFirstPThreeE->GetYaxis()->CenterTitle();
  pthreef->Update();
  pthreef->cd(2);
  hFirstPThreeW->Draw("hist");
  hFirstPThreeW->GetXaxis()->SetTitle("NBetaEvents (1st 0.3 sec)");
  hFirstPThreeW->GetXaxis()->CenterTitle();
  hFirstPThreeW->GetYaxis()->SetTitle("NBetaRuns");
  hFirstPThreeW->GetYaxis()->CenterTitle();
  pthreef->Update();
  pthreef->cd(3);
  hFirstPThreeEBkg->SetLineColor(6);
  hFirstPThreeEBkg->Draw("hist");
  hFirstPThreeEBkg->GetXaxis()->SetTitle("NBkgEvents (1st 0.3 sec)");
  hFirstPThreeEBkg->GetXaxis()->CenterTitle();
  hFirstPThreeEBkg->GetYaxis()->SetTitle("NBkgRuns");
  hFirstPThreeEBkg->GetYaxis()->CenterTitle();
  pthreef->Update();
  pthreef->cd(4);
  hFirstPThreeWBkg->SetLineColor(6);
  hFirstPThreeWBkg->Draw("hist");
  hFirstPThreeWBkg->GetXaxis()->SetTitle("NBkgEvents (1st 0.3 sec)");
  hFirstPThreeWBkg->GetXaxis()->CenterTitle();
  hFirstPThreeWBkg->GetYaxis()->SetTitle("NBkgRuns");
  hFirstPThreeWBkg->GetYaxis()->CenterTitle();
  pthreef->Print("output_files/pthreef.pdf");

  //Make plots of first 5.0 sec per run of east/west detectors.

  TCanvas *pfivef = new TCanvas("pfivef", "pfivef",10,10,600,600);
  pfivef->Divide(2,2);
  pfivef->cd(1);
  hFirstFiveE->Draw("hist");
  gPad->Update();
  hFirstFiveE->GetXaxis()->SetTitle("NBetaEvents (1st 5.0 sec)");
  hFirstFiveE->GetXaxis()->CenterTitle();
  hFirstFiveE->GetYaxis()->SetTitle("NBetaRuns");
  hFirstFiveE->GetYaxis()->CenterTitle();
  pfivef->Update();
  pfivef->cd(2);
  hFirstFiveW->Draw("hist");
  gPad->Update();
  hFirstFiveW->GetXaxis()->SetTitle("NBetaEvents (1st 5.0 sec)");
  hFirstFiveW->GetXaxis()->CenterTitle();
  hFirstFiveW->GetYaxis()->SetTitle("NBetaRuns");
  hFirstFiveW->GetYaxis()->CenterTitle();
  pfivef->Update();
  pfivef->cd(3);
  hFirstFiveEBkg->SetLineColor(6);
  hFirstFiveEBkg->Draw("hist");
  hFirstFiveEBkg->GetXaxis()->SetTitle("NBkgEvents (1st 5.0 sec)");
  hFirstFiveEBkg->GetXaxis()->CenterTitle();
  hFirstFiveEBkg->GetYaxis()->SetTitle("NBkgRuns");
  hFirstFiveEBkg->GetYaxis()->CenterTitle();
  pfivef->Update();
  pfivef->cd(4);
  hFirstFiveWBkg->SetLineColor(6);
  hFirstFiveWBkg->Draw("hist");
  hFirstFiveWBkg->GetXaxis()->SetTitle("NBkgEvents (1st 5.0 sec)");
  hFirstFiveWBkg->GetXaxis()->CenterTitle();
  hFirstFiveWBkg->GetYaxis()->SetTitle("NBkgRuns");
  hFirstFiveWBkg->GetYaxis()->CenterTitle();
  pfivef->Print("output_files/pfivef.pdf");

//  TPaveStats *tFirstFiveE = (TPaveStats*) hFirstFiveE->FindObject("stats");
  TPaveStats *tFirstFiveE;
  tFirstFiveE = (TPaveStats*) hFirstFiveE->FindObject("stats");
  tFirstFiveE->SetName("FirstFiveE Stats");
  double XEF1 = tFirstFiveE->GetX1NDC();
  double YEF1 = tFirstFiveE->GetY1NDC();
  double XEF2 = tFirstFiveE->GetX2NDC();
  double YEF2 = tFirstFiveE->GetY2NDC();
  TPaveStats *tFirstFiveW = (TPaveStats*) hFirstFiveW->FindObject("stats");
  tFirstFiveW->SetName("FirstFiveW Stats");
  double XWF1 = tFirstFiveW->GetX1NDC();
  double YWF1 = tFirstFiveW->GetY1NDC();
  double XWF2 = tFirstFiveW->GetX2NDC();
  double YWF2 = tFirstFiveW->GetY2NDC();

  TCanvas *pthreel = new TCanvas("pthreel", "pthreel",10,10,600,600);
  pthreel->Divide(2,2);
  pthreel->cd(1);
  hLastPThreeE->Draw("hist");
  hLastPThreeE->GetXaxis()->SetTitle("NBetaEvents (last 0.3 sec)");
  hLastPThreeE->GetXaxis()->CenterTitle();
  hLastPThreeE->GetYaxis()->SetTitle("NBetaRuns");
  hLastPThreeE->GetYaxis()->CenterTitle();
  pthreel->Update();
  pthreel->cd(2);
  hLastPThreeW->Draw("hist");
  hLastPThreeW->GetXaxis()->SetTitle("NBetaEvents (last 0.3 sec)");
  hLastPThreeW->GetXaxis()->CenterTitle();
  hLastPThreeW->GetYaxis()->SetTitle("NBetaRuns");
  hLastPThreeW->GetYaxis()->CenterTitle();
  pthreel->Update();
  pthreel->cd(3);
  hLastPThreeEBkg->SetLineColor(6);
  hLastPThreeEBkg->Draw("hist");
  hLastPThreeEBkg->GetXaxis()->SetTitle("NBkgEvents (last 0.3 sec)");
  hLastPThreeEBkg->GetXaxis()->CenterTitle();
  hLastPThreeEBkg->GetYaxis()->SetTitle("NBkgRuns");
  hLastPThreeEBkg->GetYaxis()->CenterTitle();
  pthreel->Update();
  pthreel->cd(4);
  hLastPThreeWBkg->SetLineColor(6);
  hLastPThreeWBkg->Draw("hist");
  hLastPThreeWBkg->GetXaxis()->SetTitle("NBkgEvents (last 0.3 sec)");
  hLastPThreeWBkg->GetXaxis()->CenterTitle();
  hLastPThreeWBkg->GetYaxis()->SetTitle("NBkgRuns");
  hLastPThreeWBkg->GetYaxis()->CenterTitle();
  pthreel->Print("output_files/pthreel.pdf");

  //Make plots of last 5.0 sec per run in east/west detectors.

  TCanvas *pfivel = new TCanvas("pfivel", "pfivel",10,10,600,600);
  pfivel->Divide(2,2);
  pfivel->cd(1);
  hLastFiveE->Draw("hist");
  gPad->Update();
  hLastFiveE->GetXaxis()->SetTitle("NBetaEvents (last 5.0 sec)");
  hLastFiveE->GetXaxis()->CenterTitle();
  hLastFiveE->GetYaxis()->SetTitle("NBetaRuns");
  hLastFiveE->GetYaxis()->CenterTitle();
  pfivel->Update();
  pfivel->cd(2);
  hLastFiveW->Draw("hist");
  gPad->Update();
  hLastFiveW->GetXaxis()->SetTitle("NBetaEvents (last 5.0 sec)");
  hLastFiveW->GetXaxis()->CenterTitle();
  hLastFiveW->GetYaxis()->SetTitle("NBetaRuns");
  hLastFiveW->GetYaxis()->CenterTitle();
  pfivel->Update();
  pfivel->cd(3);
  hLastFiveEBkg->SetLineColor(6);
  hLastFiveEBkg->Draw("hist");
  hLastFiveEBkg->GetXaxis()->SetTitle("NBkgEvents (last 5.0 sec)");
  hLastFiveEBkg->GetXaxis()->CenterTitle();
  hLastFiveEBkg->GetYaxis()->SetTitle("NBkgRuns");
  hLastFiveEBkg->GetYaxis()->CenterTitle();
  pfivel->Update();
  pfivel->cd(4);
  hLastFiveWBkg->SetLineColor(6);
  hLastFiveWBkg->Draw("hist");
  hLastFiveWBkg->GetXaxis()->SetTitle("NBkgEvents (last 5.0 sec)");
  hLastFiveWBkg->GetXaxis()->CenterTitle();
  hLastFiveWBkg->GetYaxis()->SetTitle("NBkgRuns");
  hLastFiveWBkg->GetYaxis()->CenterTitle();
  pfivel->Print("output_files/pfivel.pdf");
  TPaveStats *tLastFiveE = (TPaveStats*) hLastFiveE->FindObject("stats");
  tLastFiveE->SetName("LastFiveE Stats");
  double XEL1 = tLastFiveE->GetX1NDC();
  double YEL1 = tLastFiveE->GetY1NDC();
  double XEL2 = tLastFiveE->GetX2NDC();
  double YEL2 = tLastFiveE->GetY2NDC();
  TPaveStats *tLastFiveW = (TPaveStats*) hLastFiveW->FindObject("stats");
  tLastFiveW->SetName("LastFiveW Stats");
  double XWL1 = tLastFiveW->GetX1NDC();
  double YWL1 = tLastFiveW->GetY1NDC();
  double XWL2 = tLastFiveW->GetX2NDC();
  double YWL2 = tLastFiveW->GetY2NDC();
  
  //Simulate event distribution in first/last 5.0 sec per run using 
  //a random number generation routine for a fixed population flat distribution.

  nbetatwo = nbeta/(2.);

  TRandom3 *rtest = new TRandom3();
  rtest->SetSeed(0);
  for(Int_t i = 0; i < nbetatwo; i++) {
  FFESIMPMin[2*i] = 5.1;
  FFESIMLMin[2*i+1] = 5.1;
  FFWSIMPMin[2*i] = 5.1;
  FFWSIMLMin[2*i+1] = 5.1;
  LFESIMPMax[2*i] = 0.;
  LFESIMLMax[2*i+1] = 0.;
  LFWSIMPMax[2*i] = 0.;
  LFWSIMLMax[2*i+1] = 0.;
  for(Int_t j = 0; j < 27; j++) {
  FirstFiveESIMP = rtest->Uniform(0,5);
  FirstFiveESIML = rtest->Uniform(0,5);
  if (FirstFiveESIMP < FFESIMPMin[2*i]) FFESIMPMin[2*i] = FirstFiveESIMP;
  if (FirstFiveESIML < FFESIMLMin[2*i+1]) FFESIMLMin[2*i+1] = FirstFiveESIML;
  }
  for(Int_t j = 0; j < 27; j++) {
  FirstFiveWSIMP = rtest->Uniform(0,5);
  FirstFiveWSIML = rtest->Uniform(0,5);
  if (FirstFiveWSIMP < FFWSIMPMin[2*i]) FFWSIMPMin[2*i] = FirstFiveWSIMP;
  if (FirstFiveWSIML < FFWSIMLMin[2*i+1]) FFWSIMLMin[2*i+1] = FirstFiveWSIML;
  }
  for(Int_t j = 0; j < 70; j++) {
  LastFiveESIMP = rtest->Uniform(3595,3600);
  LastFiveESIML = rtest->Uniform(3595,3600);
  if (LastFiveESIMP > LFESIMPMax[2*i]) LFESIMPMax[2*i] = LastFiveESIMP;
  if (LastFiveESIML > LFESIMLMax[2*i+1]) LFESIMLMax[2*i+1] = LastFiveESIML;
  }
  for(Int_t j = 0; j < 69; j++) {
  LastFiveWSIMP = rtest->Uniform(3595,3600);
  LastFiveWSIML = rtest->Uniform(3595,3600);
  if (LastFiveWSIMP > LFWSIMPMax[2*i]) LFWSIMPMax[2*i] = LastFiveWSIMP;
  if (LastFiveWSIML > LFWSIMLMax[2*i+1]) LFWSIMLMax[2*i+1] = LastFiveWSIML;
  }
  TSIMWP[2*i] = LFWSIMPMax[2*i] - FFWSIMPMin[2*i];
  TSIMEP[2*i] = LFESIMPMax[2*i] - FFESIMPMin[2*i];
  TSIMWL[2*i+1] = LFWSIMLMax[2*i+1] - FFWSIMLMin[2*i+1];
  TSIMEL[2*i+1] = LFESIMLMax[2*i+1] - FFESIMLMin[2*i+1];
  SupSIM[i] = (TSIMEP[2*i]*TSIMWL[2*i+1])/(TSIMWP[2*i]*TSIMEL[2*i+1]) - 1.0;
  hSupsim->Fill(SupSIM[i]);
  }
  
  //Simulate event distribution in first/last 5.0 sec per run using 
  //a random number generation routine for a variable population from 
  //a sample probability.

  for(Int_t i = 0; i < nbetatwo; i++) {
  FEAlP[2*i] = 5.1;
  FEAlL[2*i+1] = 5.1;
  FWAlP[2*i] = 5.1;
  FWAlL[2*i+1] = 5.1;
  LEAlP[2*i] = 0.;
  LEAlL[2*i+1] = 0.;
  LWAlP[2*i] = 0.;
  LWAlL[2*i+1] = 0.;
  mm[2*i] = 0;
  p[2*i] = 0;
  q[2*i] = 0;
  r[2*i] = 0;
  ml[2*i+1] = 0;
  pl[2*i+1] = 0;
  ql[2*i+1] = 0;
  rl[2*i+1] = 0;

  for(Int_t j = 1; j < 5000001; j++) {
  FirstFiveEAlP = rtest->Uniform(0,1);
  FirstFiveEAlL = rtest->Uniform(0,1);
  
  if (FirstFiveEAlP < 0.0000054) {
  FEAlPSel = j/(1000000.);
  mm[2*i]++;
  }

  if (FirstFiveEAlL < 0.0000054) {
  FEAlLSel = j/(1000000.);
  ml[2*i+1]++;
  }

  if (0. < FEAlPSel && FEAlPSel < FEAlP[2*i]) FEAlP[2*i] = FEAlPSel;

  if (0. < FEAlLSel && FEAlLSel < FEAlL[2*i+1]) FEAlL[2*i+1] = FEAlLSel;
  }

  for(Int_t j = 1; j < 5000001; j++) {
  FirstFiveWAlP = rtest->Uniform(0,1);
  FirstFiveWAlL = rtest->Uniform(0,1);
  if (FirstFiveWAlP < 0.0000054) {
  FWAlPSel = j/(1000000.);
  p[2*i]++;
  }

  if (FirstFiveWAlL < 0.0000054) {
  FWAlLSel = j/(1000000.);
  pl[2*i+1]++;
  }

  if (0. < FWAlPSel && FWAlPSel < FWAlP[2*i]) FWAlP[2*i] = FWAlPSel;

  if (0. < FWAlLSel && FWAlLSel < FWAlL[2*i+1]) FWAlL[2*i+1] = FWAlLSel;
  }

  for(Int_t j = 1; j < 5000001; j++) {
  LastFiveEAlP = rtest->Uniform(0,1);
  LastFiveEAlL = rtest->Uniform(0,1);
  if (LastFiveEAlP < 0.000014) {
  LEAlPSel = (j/(1000000.)) + 3595.;
  q[2*i]++;
  }

  if (LastFiveEAlL < 0.000014) {
  LEAlLSel = (j/(1000000.)) + 3595.;
  ql[2*i+1]++;
  }

  if (3600. > LEAlPSel && LEAlPSel > LEAlP[2*i]) LEAlP[2*i] = LEAlPSel;

  if (3600. > LEAlLSel && LEAlLSel > LEAlL[2*i+1]) LEAlL[2*i+1] = LEAlLSel;
  }

  for(Int_t j = 1; j < 5000001; j++) {
  LastFiveWAlP = rtest->Uniform(0,1);
  LastFiveWAlL = rtest->Uniform(0,1);
  if (LastFiveWAlP < 0.0000137) {
  LWAlPSel = (j/(1000000.)) + 3595.;
  r[2*i]++;
  }

  if (LastFiveWAlL < 0.0000137) {
  LWAlLSel = (j/(1000000.)) + 3595.;
  rl[2*i+1]++;
  }

  if (3600. > LWAlPSel && LWAlPSel > LWAlP[2*i]) LWAlP[2*i] = LWAlPSel;

  if (3600. > LWAlLSel && LWAlLSel > LWAlL[2*i+1]) LWAlL[2*i+1] = LWAlLSel;
  }
  TAlWP[2*i] = LWAlP[2*i] - FWAlP[2*i];
  TAlEP[2*i] = LEAlP[2*i] - FEAlP[2*i];
  TAlWL[2*i+1] = LWAlL[2*i+1] - FWAlL[2*i+1];
  TAlEL[2*i+1] = LEAlL[2*i+1] - FEAlL[2*i+1];
  SupAlSIM[i] = (TAlEP[2*i]*TAlWL[2*i+1])/(TAlWP[2*i]*TAlEL[2*i+1]) - 1.0;
  hSupalsim->Fill(SupAlSIM[i]);
  hLWSimalEvents->Fill(r[2*i]);
  hLESimalEvents->Fill(q[2*i]);
  hFWSimalEvents->Fill(p[2*i]);
  hFESimalEvents->Fill(mm[2*i]);
  hLWSimalEvents->Fill(rl[2*i+1]);
  hLESimalEvents->Fill(ql[2*i+1]);
  hFWSimalEvents->Fill(pl[2*i+1]);
  hFESimalEvents->Fill(ml[2*i+1]);

  }
 
  //Make plots of simulations of number of events in first/last 5.0 sec. of each detector.

  TCanvas *pfivesim = new TCanvas("pfivesim", "pfivesim",10,10,600,600);
  pfivesim->Divide(2,2);
  pfivesim->cd(1);
  hLWSimalEvents->SetLineColor(4);
  hLWSimalEvents->Draw("hist");
  gPad->Update();
  hLWSimalEvents->GetXaxis()->SetTitle("NSimEventsWest (last 5.0 sec)");
  hLWSimalEvents->GetXaxis()->CenterTitle();
  hLWSimalEvents->GetYaxis()->SetTitle("NSimRuns");
  hLWSimalEvents->GetYaxis()->CenterTitle();
  pfivesim->Update();
  pfivesim->cd(2);
  hLESimalEvents->SetLineColor(4);
  hLESimalEvents->Draw("hist");
  gPad->Update();
  hLESimalEvents->GetXaxis()->SetTitle("NSimEventsEast (last 5.0 sec)");
  hLESimalEvents->GetXaxis()->CenterTitle();
  hLESimalEvents->GetYaxis()->SetTitle("NSimRuns");
  hLESimalEvents->GetYaxis()->CenterTitle();
  pfivesim->Update();
  pfivesim->cd(3);
  hFWSimalEvents->SetLineColor(4);
  hFWSimalEvents->Draw("hist");
  gPad->Update();
  hFWSimalEvents->GetXaxis()->SetTitle("NSimEventsWest (1st 5.0 sec)");
  hFWSimalEvents->GetXaxis()->CenterTitle();
  hFWSimalEvents->GetYaxis()->SetTitle("NSimRuns");
  hFWSimalEvents->GetYaxis()->CenterTitle();
  pfivesim->Update();
  pfivesim->cd(4);
  hFESimalEvents->SetLineColor(4);
  hFESimalEvents->Draw("hist");
  gPad->Update();
  hFESimalEvents->GetXaxis()->SetTitle("NSimEventsEast (1st 5.0 sec)");
  hFESimalEvents->GetXaxis()->CenterTitle();
  hFESimalEvents->GetYaxis()->SetTitle("NSimRuns");
  hFESimalEvents->GetYaxis()->CenterTitle();
  pfivesim->Print("output_files/pfivesim.pdf");
  TPaveStats *tFESimalEvents = (TPaveStats*) hFESimalEvents->FindObject("stats");
  tFESimalEvents->SetX1NDC(XEF1);
  tFESimalEvents->SetX2NDC(XEF2);
  tFESimalEvents->SetY1NDC(YEF1-(YEF2-YEF1));
  tFESimalEvents->SetY2NDC(YEF1);
  TPaveStats *tFWSimalEvents = (TPaveStats*) hFWSimalEvents->FindObject("stats");
  tFWSimalEvents->SetX1NDC(XWF1);
  tFWSimalEvents->SetX2NDC(XWF2);
  tFWSimalEvents->SetY1NDC(YWF1-(YWF2-YWF1));
  tFWSimalEvents->SetY2NDC(YWF1);
  TPaveStats *tLESimalEvents = (TPaveStats*) hLESimalEvents->FindObject("stats");
  tLESimalEvents->SetX1NDC(XEL1);
  tLESimalEvents->SetX2NDC(XEL2);
  tLESimalEvents->SetY1NDC(YEL1-(YEL2-YEL1));
  tLESimalEvents->SetY2NDC(YEL1);
  TPaveStats *tLWSimalEvents = (TPaveStats*) hLWSimalEvents->FindObject("stats");
  tLWSimalEvents->SetX1NDC(XWL1);
  tLWSimalEvents->SetX2NDC(XWL2);
  tLWSimalEvents->SetY1NDC(YWL1-(YWL2-YWL1));
  tLWSimalEvents->SetY2NDC(YWL1);

  for(Int_t i = 0; i < nbetatwo; i++) {
  FFESIMPMins[2*i] = 5.1;
  FFESIMLMins[2*i+1] = 5.1;
  FFWSIMPMins[2*i] = 5.1;
  FFWSIMLMins[2*i+1] = 5.1;
  LFESIMPMaxs[2*i] = 0.;
  LFESIMLMaxs[2*i+1] = 0.;
  LFWSIMPMaxs[2*i] = 0.;
  LFWSIMLMaxs[2*i+1] = 0.;
  for(Int_t j = 0; j < SumFirstFiveES[2*i]; j++) {
  FirstFiveESIMPS = rtest->Uniform(0,5);
  if (FirstFiveESIMPS < FFESIMPMins[2*i]) FFESIMPMins[2*i] = FirstFiveESIMPS;
  }
  for(Int_t j = 0; j < SumFirstFiveES[2*i+1]; j++) {
  FirstFiveESIMLS = rtest->Uniform(0,5);
  if (FirstFiveESIMLS < FFESIMLMins[2*i+1]) FFESIMLMins[2*i+1] = FirstFiveESIMLS;
  }
  for(Int_t j = 0; j < SumFirstFiveWS[2*i]; j++) {
  FirstFiveWSIMPS = rtest->Uniform(0,5);
  if (FirstFiveWSIMPS < FFWSIMPMins[2*i]) FFWSIMPMins[2*i] = FirstFiveWSIMPS;
  }
  for(Int_t j = 0; j < SumFirstFiveWS[2*i+1]; j++) {
  FirstFiveWSIMLS = rtest->Uniform(0,5);
  if (FirstFiveWSIMLS < FFWSIMLMins[2*i+1]) FFWSIMLMins[2*i+1] = FirstFiveWSIMLS;
  }
  for(Int_t j = 0; j < SumLastFiveES[2*i]; j++) {
  LastFiveESIMPS = rtest->Uniform(3595,3600);
  if (LastFiveESIMPS > LFESIMPMaxs[2*i]) LFESIMPMaxs[2*i] = LastFiveESIMPS;
  }
  for(Int_t j = 0; j < SumLastFiveES[2*i+1]; j++) {
  LastFiveESIMLS = rtest->Uniform(3595,3600);
  if (LastFiveESIMLS > LFESIMLMaxs[2*i+1]) LFESIMLMaxs[2*i+1] = LastFiveESIMLS;
  }
  for(Int_t j = 0; j < SumLastFiveWS[2*i]; j++) {
  LastFiveWSIMPS = rtest->Uniform(3595,3600);
  if (LastFiveWSIMPS > LFWSIMPMaxs[2*i]) LFWSIMPMaxs[2*i] = LastFiveWSIMPS;
  }
  for(Int_t j = 0; j < SumLastFiveWS[2*i+1]; j++) {
  LastFiveWSIMLS = rtest->Uniform(3595,3600);
  if (LastFiveWSIMLS > LFWSIMLMaxs[2*i+1]) LFWSIMLMaxs[2*i+1] = LastFiveWSIMLS;
  }
  TSIMWPS[2*i] = LFWSIMPMaxs[2*i] - FFWSIMPMins[2*i];
  TSIMEPS[2*i] = LFESIMPMaxs[2*i] - FFESIMPMins[2*i];
  TSIMWLS[2*i+1] = LFWSIMLMaxs[2*i+1] - FFWSIMLMins[2*i+1];
  TSIMELS[2*i+1] = LFESIMLMaxs[2*i+1] - FFESIMLMins[2*i+1];
  SupSIMS[i] = (TSIMEPS[2*i]*TSIMWLS[2*i+1])/(TSIMWPS[2*i]*TSIMELS[2*i+1]) - 1.0;
  hSupsims->Fill(SupSIMS[i]);
  }

  for(Int_t i = 0; i < nbetatwo; i++) {
  FEAlPS[2*i] = 5.1;
  FEAlLS[2*i+1] = 5.1;
  FWAlPS[2*i] = 5.1;
  FWAlLS[2*i+1] = 5.1;
  LEAlPS[2*i] = 0.;
  LEAlLS[2*i+1] = 0.;
  LWAlPS[2*i] = 0.;
  LWAlLS[2*i+1] = 0.;
  ms[2*i] = 0;
  ps[2*i] = 0;
  qs[2*i] = 0;
  rs[2*i] = 0;
  msl[2*i+1] = 0;
  psl[2*i+1] = 0;
  qsl[2*i+1] = 0;
  rsl[2*i+1] = 0;

  for(Int_t j = 1; j < 5000001; j++) {
  FirstFiveEAlPS = rtest->Uniform(0,1);
  FirstFiveEAlLS = rtest->Uniform(0,1);
  
  if (FirstFiveEAlPS < ProbFirstFiveE[2*i]) {
  FEAlPSelS = j/(1000000.);
  ms[2*i]++;
  }

  if (FirstFiveEAlLS < ProbFirstFiveE[2*i+1]) {
  FEAlLSelS = j/(1000000.);
  msl[2*i+1]++;
  }

  if (0. < FEAlPSelS && FEAlPSelS < FEAlPS[2*i]) FEAlPS[2*i] = FEAlPSelS;

  if (0. < FEAlLSelS && FEAlLSelS < FEAlLS[2*i+1]) FEAlLS[2*i+1] = FEAlLSelS;
  }

  for(Int_t j = 1; j < 5000001; j++) {
  FirstFiveWAlPS = rtest->Uniform(0,1);
  FirstFiveWAlLS = rtest->Uniform(0,1);
  if (FirstFiveWAlPS < ProbFirstFiveW[2*i]) {
  FWAlPSelS = j/(1000000.);
  ps[2*i]++;
  }

  if (FirstFiveWAlLS < ProbFirstFiveW[2*i+1]) {
  FWAlLSelS = j/(1000000.);
  psl[2*i+1]++;
  }

  if (0. < FWAlPSelS && FWAlPSelS < FWAlPS[2*i]) FWAlPS[2*i] = FWAlPSelS;

  if (0. < FWAlLSelS && FWAlLSelS < FWAlLS[2*i+1]) FWAlLS[2*i+1] = FWAlLSelS;
  }

  for(Int_t j = 1; j < 5000001; j++) {
  LastFiveEAlPS = rtest->Uniform(0,1);
  LastFiveEAlLS = rtest->Uniform(0,1);
  if (LastFiveEAlPS < ProbLastFiveE[2*i]) {
  LEAlPSelS = (j/(1000000.)) + 3595.;
  qs[2*i]++;
  }

  if (LastFiveEAlLS < ProbLastFiveE[2*i+1]) {
  LEAlLSelS = (j/(1000000.)) + 3595.;
  qsl[2*i+1]++;
  }

  if (3600. > LEAlPSelS && LEAlPSelS > LEAlPS[2*i]) LEAlPS[2*i] = LEAlPSelS;

  if (3600. > LEAlLSelS && LEAlLSelS > LEAlLS[2*i+1]) LEAlLS[2*i+1] = LEAlLSelS;
  }

  for(Int_t j = 1; j < 5000001; j++) {
  LastFiveWAlPS = rtest->Uniform(0,1);
  LastFiveWAlLS = rtest->Uniform(0,1);
  if (LastFiveWAlPS < ProbLastFiveW[2*i]) {
  LWAlPSelS = (j/(1000000.)) + 3595.;
  rs[2*i]++;
  }

  if (LastFiveWAlLS < ProbLastFiveW[2*i+1]) {
  LWAlLSelS = (j/(1000000.)) + 3595.;
  rsl[2*i+1]++;
  }
 
  if (3600. > LWAlPSelS && LWAlPSelS > LWAlPS[2*i]) LWAlPS[2*i] = LWAlPSelS;

  if (3600. > LWAlLSelS && LWAlLSelS > LWAlLS[2*i+1]) LWAlLS[2*i+1] = LWAlLSelS;
  }
  TAlWPS[2*i] = LWAlPS[2*i] - FWAlPS[2*i];
  TAlEPS[2*i] = LEAlPS[2*i] - FEAlPS[2*i];
  TAlWLS[2*i+1] = LWAlLS[2*i+1] - FWAlLS[2*i+1];
  TAlELS[2*i+1] = LEAlLS[2*i+1] - FEAlLS[2*i+1];
  SupAlSIMS[i] = (TAlEPS[2*i]*TAlWLS[2*i+1])/(TAlWPS[2*i]*TAlELS[2*i+1]) - 1.0;
  hSupalsims->Fill(SupAlSIMS[i]);
  hLWSimalEventss->Fill(rs[2*i]);
  hLESimalEventss->Fill(qs[2*i]);
  hFWSimalEventss->Fill(ps[2*i]);
  hFESimalEventss->Fill(ms[2*i]);
  hLWSimalEventss->Fill(rsl[2*i+1]);
  hLESimalEventss->Fill(qsl[2*i+1]);
  hFWSimalEventss->Fill(psl[2*i+1]);
  hFESimalEventss->Fill(msl[2*i+1]);
  }

  TCanvas *pfivesims = new TCanvas("pfivesims", "pfivesims",10,10,600,600);
  pfivesims->Divide(2,2);
  pfivesims->cd(1);
  hLWSimalEventss->SetLineColor(3);
  hLWSimalEventss->Draw("hist");
  gPad->Update();
  hLWSimalEventss->GetXaxis()->SetTitle("NSimsEventsWest (last 5.0 sec)");
  hLWSimalEventss->GetXaxis()->CenterTitle();
  hLWSimalEventss->GetYaxis()->SetTitle("NSimsRuns");
  hLWSimalEventss->GetYaxis()->CenterTitle();
  pfivesims->Update();
  pfivesims->cd(2);
  hLESimalEventss->SetLineColor(3);
  hLESimalEventss->Draw("hist");
  gPad->Update();
  hLESimalEventss->GetXaxis()->SetTitle("NSimsEventsEast (last 5.0 sec)");
  hLESimalEventss->GetXaxis()->CenterTitle();
  hLESimalEventss->GetYaxis()->SetTitle("NSimsRuns");
  hLESimalEventss->GetYaxis()->CenterTitle();
  pfivesims->Update();
  pfivesims->cd(3);
  hFWSimalEventss->SetLineColor(3);
  hFWSimalEventss->Draw("hist");
  gPad->Update();
  hFWSimalEventss->GetXaxis()->SetTitle("NSimsEventsWest (1st 5.0 sec)");
  hFWSimalEventss->GetXaxis()->CenterTitle();
  hFWSimalEventss->GetYaxis()->SetTitle("NSimsRuns");
  hFWSimalEventss->GetYaxis()->CenterTitle();
  pfivesims->Update();
  pfivesims->cd(4);
  hFESimalEventss->SetLineColor(3);
  hFESimalEventss->Draw("hist");
  gPad->Update();
  hFESimalEventss->GetXaxis()->SetTitle("NSimsEventsEast (1st 5.0 sec)");
  hFESimalEventss->GetXaxis()->CenterTitle();
  hFESimalEventss->GetYaxis()->SetTitle("NSimsRuns");
  hFESimalEventss->GetYaxis()->CenterTitle();
  pfivesims->Print("output_files/pfivesims.pdf");
  TPaveStats *tFESimalEventss = (TPaveStats*) hFESimalEventss->FindObject("stats");
  tFESimalEventss->SetX1NDC(XEF1);
  tFESimalEventss->SetX2NDC(XEF2);
  tFESimalEventss->SetY1NDC(YEF1-(YEF2-YEF1)-(YEF2-YEF1));
  tFESimalEventss->SetY2NDC(YEF1-(YEF2-YEF1));
  TPaveStats *tFWSimalEventss = (TPaveStats*) hFWSimalEventss->FindObject("stats");
  tFWSimalEventss->SetX1NDC(XWF1);
  tFWSimalEventss->SetX2NDC(XWF2);
  tFWSimalEventss->SetY1NDC(YWF1-(YWF2-YWF1)-(YWF2-YWF1));
  tFWSimalEventss->SetY2NDC(YWF1-(YWF2-YWF1));
  TPaveStats *tLESimalEventss = (TPaveStats*) hLESimalEventss->FindObject("stats");
  tLESimalEventss->SetX1NDC(XEL1);
  tLESimalEventss->SetX2NDC(XEL2);
  tLESimalEventss->SetY1NDC(YEL1-(YEL2-YEL1)-(YEL2-YEL1));
  tLESimalEventss->SetY2NDC(YEL1-(YEL2-YEL1));
  TPaveStats *tLWSimalEventss = (TPaveStats*) hLWSimalEventss->FindObject("stats");
  tLWSimalEventss->SetX1NDC(XWL1);
  tLWSimalEventss->SetX2NDC(XWL2);
  tLWSimalEventss->SetY1NDC(YWL1-(YWL2-YWL1)-(YWL2-YWL1));
  tLWSimalEventss->SetY2NDC(YWL1-(YWL2-YWL1));

  //Make plots and overlays of data and simulations for first 5.0 sec per run of each detector.

  TCanvas *pfivefover = new TCanvas("pfivefover", "pfivefover",10,10,600,600);
  pfivefover->Divide(2,2);
  pfivefover->cd(1);
  hFWSimalEvents->SetLineColor(4);
  hFWSimalEvents->Draw("HIST");
  hFWSimalEvents->GetXaxis()->SetTitle("NEventsWest (1st 5.0 sec)");
  hFWSimalEvents->GetXaxis()->CenterTitle();
  hFWSimalEvents->GetYaxis()->SetTitle("NRuns");
  hFWSimalEvents->GetYaxis()->CenterTitle();
  hFWSimalEventss->SetLineColor(3);
  hFWSimalEventss->Draw("HIST SAME");
  hFWSimalEventss->GetXaxis()->SetTitle("NEventsWest (1st 5.0 sec)");
  hFWSimalEventss->GetXaxis()->CenterTitle();
  hFWSimalEventss->GetYaxis()->SetTitle("NRuns");
  hFWSimalEventss->GetYaxis()->CenterTitle();
  hFirstFiveW->Draw("HIST SAME");
  hFirstFiveW->GetXaxis()->SetTitle("NEventsWest (1st 5.0 sec)");
  hFirstFiveW->GetXaxis()->CenterTitle();
  hFirstFiveW->GetYaxis()->SetTitle("NRuns");
  hFirstFiveW->GetYaxis()->CenterTitle();
  tFWSimalEvents->Draw("SAME");
  tFWSimalEventss->Draw("SAME");
  tFirstFiveW->Draw("SAME");
  TLegend *legtwf = new TLegend(0.1,0.75,0.3,0.9);
   legtwf->AddEntry(hFirstFiveW,"Beta Run Data","l");
   legtwf->AddEntry(hFWSimalEvents,"Sim No. 1","l");
   legtwf->AddEntry(hFWSimalEventss,"Sim No. 2","l");
   legtwf->Draw();
  pfivefover->Update();
  pfivefover->cd(2);
  hFESimalEvents->SetLineColor(4);
  hFESimalEvents->Draw("HIST");
  hFESimalEvents->GetXaxis()->SetTitle("NEventsEast (1st 5.0 sec)");
  hFESimalEvents->GetXaxis()->CenterTitle();
  hFESimalEvents->GetYaxis()->SetTitle("NRuns");
  hFESimalEvents->GetYaxis()->CenterTitle();
  hFESimalEventss->SetLineColor(3);
  hFESimalEventss->Draw("HIST SAME");
  hFESimalEventss->GetXaxis()->SetTitle("NEventsEast (1st 5.0 sec)");
  hFESimalEventss->GetXaxis()->CenterTitle();
  hFESimalEventss->GetYaxis()->SetTitle("NRuns");
  hFESimalEventss->GetYaxis()->CenterTitle();
  hFirstFiveE->Draw("HIST SAME");
  hFirstFiveE->GetXaxis()->SetTitle("NEventsEast (1st 5.0 sec)");
  hFirstFiveE->GetXaxis()->CenterTitle();
  hFirstFiveE->GetYaxis()->SetTitle("NRuns");
  hFirstFiveE->GetYaxis()->CenterTitle();
  tFESimalEvents->Draw("SAME");
  tFESimalEventss->Draw("SAME");
  tFirstFiveE->Draw("SAME");
  TLegend *legtef = new TLegend(0.1,0.75,0.3,0.9);
   legtef->AddEntry(hFirstFiveE,"Beta Run Data","l");
   legtef->AddEntry(hFESimalEvents,"Sim No. 1","l");
   legtef->AddEntry(hFESimalEventss,"Sim No. 2","l");
   legtef->Draw();
  pfivefover->Update();
  pfivefover->cd(3);
  hFirstFiveEBkg->SetLineColor(6);
  hFirstFiveEBkg->Draw("hist");
  hFirstFiveEBkg->GetXaxis()->SetTitle("NEventsEast (1st 5.0 sec)");
  hFirstFiveEBkg->GetXaxis()->CenterTitle();
  hFirstFiveEBkg->GetYaxis()->SetTitle("NRuns");
  hFirstFiveEBkg->GetYaxis()->CenterTitle();
  TLegend *legtefb = new TLegend(0.1,0.8,0.3,0.9);
   legtefb->AddEntry(hFirstFiveEBkg,"Bkg Run Data","l");
   legtefb->Draw();
  pfivefover->Update();
  pfivefover->cd(4);
  hFirstFiveWBkg->SetLineColor(6);
  hFirstFiveWBkg->Draw("hist");
  hFirstFiveWBkg->GetXaxis()->SetTitle("NEventsWest (1st 5.0 sec)");
  hFirstFiveWBkg->GetXaxis()->CenterTitle();
  hFirstFiveWBkg->GetYaxis()->SetTitle("NRuns");
  hFirstFiveWBkg->GetYaxis()->CenterTitle();
  TLegend *legtwfb = new TLegend(0.1,0.8,0.3,0.9);
   legtwfb->AddEntry(hFirstFiveWBkg,"Bkg Run Data","l");
   legtwfb->Draw();
  pfivefover->Print("output_files/pfivefover.pdf");

  //Make plots and overlays of data and simulations for last 5.0 sec per run of each detector.

  TCanvas *pfivelover = new TCanvas("pfivelover", "pfivelover",10,10,600,600);
  pfivelover->Divide(2,2);
  pfivelover->cd(1);
  hLWSimalEvents->SetLineColor(4);
  hLWSimalEvents->Draw("HIST");
  hLWSimalEvents->GetXaxis()->SetTitle("NEventsWest (last 5.0 sec)");
  hLWSimalEvents->GetXaxis()->CenterTitle();
  hLWSimalEvents->GetYaxis()->SetTitle("NRuns");
  hLWSimalEvents->GetYaxis()->CenterTitle();
  hLWSimalEventss->SetLineColor(3);
  hLWSimalEventss->Draw("HIST SAME");
  hLWSimalEventss->GetXaxis()->SetTitle("NEventsWest (last 5.0 sec)");
  hLWSimalEventss->GetXaxis()->CenterTitle();
  hLWSimalEventss->GetYaxis()->SetTitle("NRuns");
  hLWSimalEventss->GetYaxis()->CenterTitle();
  hLastFiveW->Draw("HIST SAME");
  hLastFiveW->GetXaxis()->SetTitle("NEventsWest (last 5.0 sec)");
  hLastFiveW->GetXaxis()->CenterTitle();
  hLastFiveW->GetYaxis()->SetTitle("NRuns");
  hLastFiveW->GetYaxis()->CenterTitle();
  tLWSimalEvents->Draw("SAME");
  tLWSimalEventss->Draw("SAME");
  tLastFiveW->Draw("SAME");
  TLegend *legtwl = new TLegend(0.1,0.75,0.3,0.9);
   legtwl->AddEntry(hLastFiveW,"Beta Run Data","l");
   legtwl->AddEntry(hLWSimalEvents,"Sim No. 1","l");
   legtwl->AddEntry(hLWSimalEventss,"Sim No. 2","l");
   legtwl->Draw();
  pfivelover->Update();
  pfivelover->cd(2);
  hLESimalEvents->SetLineColor(4);
  hLESimalEvents->Draw("HIST");
  hLESimalEvents->GetXaxis()->SetTitle("NEventsEast (last 5.0 sec)");
  hLESimalEvents->GetXaxis()->CenterTitle();
  hLESimalEvents->GetYaxis()->SetTitle("NRuns");
  hLESimalEvents->GetYaxis()->CenterTitle();
  hLESimalEventss->SetLineColor(3);
  hLESimalEventss->Draw("HIST SAME");
  hLESimalEventss->GetXaxis()->SetTitle("NEventsEast (last 5.0 sec)");
  hLESimalEventss->GetXaxis()->CenterTitle();
  hLESimalEventss->GetYaxis()->SetTitle("NRuns");
  hLESimalEventss->GetYaxis()->CenterTitle();
  hLastFiveE->Draw("HIST SAME");
  hLastFiveE->GetXaxis()->SetTitle("NEventsEast (last 5.0 sec)");
  hLastFiveE->GetXaxis()->CenterTitle();
  hLastFiveE->GetYaxis()->SetTitle("NRuns");
  hLastFiveE->GetYaxis()->CenterTitle();
  tLESimalEvents->Draw("SAME");
  tLESimalEventss->Draw("SAME");
  tLastFiveE->Draw("SAME");
  TLegend *legtel = new TLegend(0.1,0.75,0.3,0.9);
   legtel->AddEntry(hLastFiveE,"Beta Run Data","l");
   legtel->AddEntry(hLESimalEvents,"Sim No. 1","l");
   legtel->AddEntry(hLESimalEventss,"Sim No. 2","l");
   legtel->Draw();
  pfivelover->Update();
  pfivelover->cd(3);
  hLastFiveEBkg->SetLineColor(6);
  hLastFiveEBkg->Draw("hist");
  hLastFiveEBkg->GetXaxis()->SetTitle("NEventsEast (last 5.0 sec)");
  hLastFiveEBkg->GetXaxis()->CenterTitle();
  hLastFiveEBkg->GetYaxis()->SetTitle("NRuns");
  hLastFiveEBkg->GetYaxis()->CenterTitle();
  TLegend *legtelb = new TLegend(0.1,0.8,0.3,0.9);
   legtelb->AddEntry(hLastFiveEBkg,"Bkg Run Data","l");
   legtelb->Draw();
  pfivelover->Update();
  pfivelover->cd(4);
  hLastFiveWBkg->SetLineColor(6);
  hLastFiveWBkg->Draw("hist");
  hLastFiveWBkg->GetXaxis()->SetTitle("NEventsWest (last 5.0 sec)");
  hLastFiveWBkg->GetXaxis()->CenterTitle();
  hLastFiveWBkg->GetYaxis()->SetTitle("NRuns");
  hLastFiveWBkg->GetYaxis()->CenterTitle();
  TLegend *legtwlb = new TLegend(0.1,0.8,0.3,0.9);
   legtwlb->AddEntry(hLastFiveWBkg,"Bkg Run Data","l");
   legtwlb->Draw();
  pfivelover->Print("output_files/pfivelover.pdf");

  //Fill histogram with superratio data.

  for(Int_t i = 0; i < noct ; i++){
    OctNumber.push_back((Double_t)i);

    cout << "Dave's superratio for A2+A5 is " << octet[i]->Asupertime1 << endl;
    cout << "Dave's superratio for A7+A10 is " << octet[i]->Asupertime2 << endl;
    cout << "Dave's superratio for B2+B5 is " << octet[i]->Bsupertime1 << endl;
    cout << "Dave's superratio for B7+B10 is " << octet[i]->Bsupertime2 << endl;
    hTimeratio->Fill(octet[i]->Asupertime1);
    hTimeratio->Fill(octet[i]->Asupertime2);
    hTimeratio->Fill(octet[i]->Bsupertime1);
    hTimeratio->Fill(octet[i]->Bsupertime2);
    hTimeratioall->Fill(octet[i]->Asupertime1a);
    hTimeratioall->Fill(octet[i]->Asupertime2a);
    hTimeratioall->Fill(octet[i]->Bsupertime1a);
    hTimeratioall->Fill(octet[i]->Bsupertime2a);
    foct << i << "\t" << octet[i]->GetFirst() << "-" << octet[i]->GetLast() << "\t" << octet[i]->A_sum[2];
    foct << "\t" << octet[i]->A_sumer[2] <<  endl;
    fquart << i<< "\t" << octet[i]->GetFirst() << "-" << octet[i]->GetLast() << "\t" << octet[i]->A_sum_A[2];
    fquart << "\t" <<octet[i]->A_sum_Aer[2] << "\t" << octet[i]->A_sum_B[2] << "\t"<< octet[i]->A_sum_Ber[2] << endl;
    fsup << i<< "\t" << octet[i]->GetFirst() << "-" << octet[i]->GetLast() << "\t";
    fsup << octet[i]->Asuper1[2] << "\t" << octet[i]->Asuper1er[2] << "\t" << octet[i]->Asuper2[2];
    fsup << "\t" << octet[i]->Asuper2er[2] << "\t";
    fsup << octet[i]->Bsuper1[2] << "\t" << octet[i]->Bsuper1er[2] << "\t" <<octet[i]->Bsuper2[2];
    fsup << "\t" <<octet[i]->Bsuper2er[2] << endl;
    
    // Output the oct list used in the analysis.
    if(IsNaN(octet[i]->Asuper2[2]) || IsNaN(octet[i]->Asuper1[2]) || IsNaN(octet[i]->Bsuper1[2]) || IsNaN(octet[i]->Bsuper2[2]) || 
     octet[i]->Asuper2[2] == 0 ||octet[i]->Asuper1[2] == 0 || octet[i]->Bsuper2[2] == 0 || octet[i]->Bsuper1[2] == 0   ){
      if(!(IsNaN(octet[i]->Asuper1[2])) && octet[i]->Asuper1[2] != 0  ){
          Fill_Asymmetry_Vector(A_super,octet[i]->Asuper1[2],octet[i]->Asuper1er[2],NoctC);
      }
      if(!(IsNaN(octet[i]->Asuper2[2])) && octet[i]->Asuper2[2] != 0  ){
          Fill_Asymmetry_Vector(A_super,octet[i]->Asuper2[2],octet[i]->Asuper2er[2],NoctC);
      }
      
      if(!(IsNaN(octet[i]->Bsuper1[2])) && octet[i]->Bsuper1[2] != 0 ){
          Fill_Asymmetry_Vector(A_super,octet[i]->Bsuper1[2],octet[i]->Bsuper1er[2],NoctC);
      }
      
      if(!(IsNaN(octet[i]->Bsuper2[2])) && octet[i]->Bsuper2[2] != 0 ){
         Fill_Asymmetry_Vector(A_super,octet[i]->Bsuper2[2],octet[i]->Bsuper2er[2],NoctC);
      }
    } else {
      Fill_Asymmetry_Vector(A_super,octet[i]->Asuper1[2],octet[i]->Asuper1er[2],NoctC);
      Fill_Asymmetry_Vector(A_super,octet[i]->Asuper2[2],octet[i]->Asuper2er[2],NoctC);
      Fill_Asymmetry_Vector(A_super,octet[i]->Bsuper1[2],octet[i]->Bsuper1er[2],NoctC);
      Fill_Asymmetry_Vector(A_super,octet[i]->Bsuper2[2],octet[i]->Bsuper2er[2],NoctC);
    }
    if(octet[i]->A_sum_A[2] != 0 ){
        Fill_Asymmetry_Vector(A_quat,octet[i]->A_sum_A[2],octet[i]->A_sum_Aer[2],NQuat);
    }
    if(octet[i]->A_sum_B[2] != 0 ){
        Fill_Asymmetry_Vector(A_quat,octet[i]->A_sum_B[2],octet[i]->A_sum_Ber[2],NQuat);
    }
      
      TotalCounts          += octet[i]->TotCounts;
      CurrentCounts.push_back(TotalCounts);
      OctDate.push_back(btr[octet[i]->GetIndex()]->RunDate->GetDate());

      for(Int_t ii = 0 ; ii < 9 ; ii++){
	      Fill_Asymmetry_Vector_Oct(A_Aves[ii],octet[i]->A_multi[ii],octet[i]->A_multier[ii],FullOct);
	      Fill_Asymmetry_Vector_Oct(A_Aves_s[ii],octet[i]->A_sum[ii],octet[i]->A_sumer[ii],FullOct);
      }
            
      Fill_Asymmetry_Vector_Oct(A_sum_A,octet[i]->A_sum_A[2],octet[i]->A_sum_Aer[2],FullOct);
      Fill_Asymmetry_Vector_Oct(A_sum_B,octet[i]->A_sum_B[2],octet[i]->A_sum_Ber[2],FullOct);
      Fill_Asymmetry_Vector_Oct(A_multi_A,octet[i]->A_multi_A[2],octet[i]->A_multi_Aer[2],FullOct);
      Fill_Asymmetry_Vector_Oct(A_multi_B,octet[i]->A_multi_B[2],octet[i]->A_multi_Ber[2],FullOct);

      FullOct++;
  } 
  
  //Make plots and overlays of superratio data and simulations.

//   TCanvas *supersimal = new TCanvas("superratiosimal", "superratiosimal",10,10,600,600);
//   hSupalsim->SetLineColor(4);
//   hSupalsim->GetXaxis()->SetLabelSize(0.03);
//   hSupalsim->GetYaxis()->SetLabelSize(0.03);
//   hSupalsim->Draw("HIST");
//   gPad->Update();
//   supersimal->Draw("output_files/timeratiosimal.pdf"); 
//   TPaveStats *tSupalsim = (TPaveStats*) hSupalsim->FindObject("stats");
//   tSupalsim->SetName("Sim2 Stats");
//   double X1 = tSupalsim->GetX1NDC();
//   double Y1 = tSupalsim->GetY1NDC();
//   double X2 = tSupalsim->GetX2NDC();
//   double Y2 = tSupalsim->GetY2NDC(); 
// 
//   TCanvas *supersim = new TCanvas("superratiosim", "superratiosim",10,10,600,600);
//   hSupsim->SetLineColor(2);
//   hSupsim->GetXaxis()->SetLabelSize(0.03);
//   hSupsim->GetYaxis()->SetLabelSize(0.03);
//   hSupsim->Draw("HIST");
//   gPad->Update();
//   supersim->Print("output_files/timeratio_sim.pdf");
//   TPaveStats *tSupsim = (TPaveStats*) hSupsim->FindObject("stats");
//   tSupsim->SetX1NDC(X1);
//   tSupsim->SetX2NDC(X2);
//   tSupsim->SetY1NDC(Y1-(Y2-Y1));
//   tSupsim->SetY2NDC(Y1);
// 
//   TCanvas *clocks = new TCanvas("superratioclocks", "superratioclocks",10,10,600,600);
//   hTimeratio->GetXaxis()->SetLabelSize(0.03);
//   hTimeratio->GetYaxis()->SetLabelSize(0.03);
//   hTimeratio->Draw("HIST");
//   gPad->Update();
//   clocks->Print("output_files/timeratio.pdf");
//   TPaveStats *tTimeratio = (TPaveStats*) hTimeratio->FindObject("stats");
//   tTimeratio->SetX1NDC(X1-(X2-X1));
//   tTimeratio->SetX2NDC(X1);
//   tTimeratio->SetY1NDC(Y1);
//   tTimeratio->SetY2NDC(Y2);
// 
//   TCanvas *timeoverlay = new TCanvas("timeoverlay", "timeoverlay",10,10,600,600);
//   hSupalsim->SetLineColor(4);
//   hSupalsim->GetXaxis()->SetLabelSize(0.03);
//   hSupalsim->GetYaxis()->SetLabelSize(0.03);
//   hSupalsim->Draw("HIST");
//   hSupsim->SetLineColor(2);
//   hSupsim->GetXaxis()->SetLabelSize(0.03);
//   hSupsim->GetYaxis()->SetLabelSize(0.03);
//   hSupsim->Draw("HIST SAME");
//   hTimeratio->GetXaxis()->SetLabelSize(0.03);
//   hTimeratio->GetYaxis()->SetLabelSize(0.03);
//   hTimeratio->Draw("HIST SAME");
//   tSupalsim->Draw("SAME");
//   tSupsim->Draw("SAME");
//   tTimeratio->Draw("SAME");
//   TLegend *leg = new TLegend(0.1,0.75,0.3,0.9);
//    leg->AddEntry(hTimeratio,"Beta Run Data","l");
//    leg->AddEntry(hSupsim,"Sim No. 1, flat","l");
//    leg->AddEntry(hSupalsim,"Sim No. 1","l");
//    leg->Draw();
//   timeoverlay->Print("output_files/timeratio_overlay.pdf");
// 
//   TCanvas *supersimals = new TCanvas("superratiosimals", "superratiosimals",10,10,600,600);
//   hSupalsims->SetLineColor(3);
//   hSupalsims->GetXaxis()->SetLabelSize(0.03);
//   hSupalsims->GetYaxis()->SetLabelSize(0.03);
//   hSupalsims->Draw("HIST");
//   gPad->Update();
//   supersimals->Draw("output_files/timeratiosimals.pdf");
//   TPaveStats *tSupalsims = (TPaveStats*) hSupalsims->FindObject("stats");
//   tSupalsims->SetName("Sim2 Stats");
//   double XS1 = tSupalsims->GetX1NDC();
//   double YS1 = tSupalsims->GetY1NDC();
//   double XS2 = tSupalsims->GetX2NDC();
//   double YS2 = tSupalsims->GetY2NDC();
// 
//   TCanvas *supersims = new TCanvas("superratiosims", "superratiosims",10,10,600,600);
//   hSupsims->SetLineColor(2);
//   hSupsims->GetXaxis()->SetLabelSize(0.03);
//   hSupsims->GetYaxis()->SetLabelSize(0.03);
//   hSupsims->Draw("HIST");
//   gPad->Update();
//   supersims->Print("output_files/timeratio_sims.pdf");
//   TPaveStats *tSupsims = (TPaveStats*) hSupsims->FindObject("stats");
//   tSupsims->SetX1NDC(XS1);
//   tSupsims->SetX2NDC(XS2);
//   tSupsims->SetY1NDC(YS1-(YS2-YS1));
//   tSupsims->SetY2NDC(YS1);
// 
//   TCanvas *clockss = new TCanvas("superratioclockss", "superratioclockss",10,10,600,600);
//   hTimeratio->GetXaxis()->SetLabelSize(0.03);
//   hTimeratio->GetYaxis()->SetLabelSize(0.03);
//   hTimeratio->Draw("HIST");
//   gPad->Update();
//   TPaveStats *tTimeratios = (TPaveStats*) hTimeratio->FindObject("stats");
//   tTimeratios->SetX1NDC(XS1-(XS2-XS1));
//   tTimeratios->SetX2NDC(XS1);
//   tTimeratios->SetY1NDC(YS1);
//   tTimeratios->SetY2NDC(YS2);
//   
//   TCanvas *clocksall = new TCanvas("clocksall", "clocksall",10,10,600,600);
//   hTimeratioall->GetXaxis()->SetLabelSize(0.03);
//   hTimeratioall->GetYaxis()->SetLabelSize(0.03);
//   hTimeratioall->Draw("HIST");
//   gPad->Update();
//   clocksall->Print("output_files/timeratioall.pdf");
//   TPaveStats *tTimeratioall = (TPaveStats*) hTimeratioall->FindObject("stats");
//   tTimeratioall->SetX1NDC(XS1-(XS2-XS1));
//   tTimeratioall->SetX2NDC(XS1);
//   tTimeratioall->SetY1NDC(YS1);
//   tTimeratioall->SetY2NDC(YS2);
// 
//   TCanvas *timeoverlays = new TCanvas("timeoverlays", "timeoverlays",10,10,600,600);
//   hSupalsims->SetLineColor(3);
//   hSupalsims->GetXaxis()->SetLabelSize(0.03);
//   hSupalsims->GetYaxis()->SetLabelSize(0.03);
//   hSupalsims->Draw("HIST");
//   hSupsims->SetLineColor(9);
//   hSupsims->GetXaxis()->SetLabelSize(0.03);
//   hSupsims->GetYaxis()->SetLabelSize(0.03);
//   hSupsims->Draw("HIST SAME");
//   hTimeratio->GetXaxis()->SetLabelSize(0.03);
//   hTimeratio->GetYaxis()->SetLabelSize(0.03);
//   hTimeratio->Draw("HIST SAME");
//   tSupalsims->Draw("SAME");
//   tSupsims->Draw("SAME");
//   tTimeratios->Draw("SAME");
//   TLegend *legs = new TLegend(0.1,0.75,0.3,0.9);
//    legs->AddEntry(hTimeratio,"Beta Run Data","l");
//    legs->AddEntry(hSupsims,"Sim No. 2, flat","l");
//    legs->AddEntry(hSupalsims,"Sim No. 2","l");
//    legs->Draw();
//   timeoverlays->Print("output_files/timeratio_overlays.pdf");
// 
//   TCanvas *timeoverlaysall = new TCanvas("timeoverlaysall", "timeoverlaysall",10,10,600,600);
//   hSupalsims->SetLineColor(3);
//   hSupalsims->GetXaxis()->SetLabelSize(0.03);
//   hSupalsims->GetYaxis()->SetLabelSize(0.03);
//   hSupalsims->Draw("HIST");
//   hSupsims->SetLineColor(9);
//   hSupsims->GetXaxis()->SetLabelSize(0.03);
//   hSupsims->GetYaxis()->SetLabelSize(0.03);
//   hSupsims->Draw("HIST SAME");
//   hTimeratioall->GetXaxis()->SetLabelSize(0.03);
//   hTimeratioall->GetYaxis()->SetLabelSize(0.03);
//   hTimeratioall->Draw("HIST SAME");
//   tSupalsims->Draw("SAME");
//   tSupsims->Draw("SAME");
//   tTimeratioall->Draw("SAME");
//   TLegend *legall = new TLegend(0.1,0.75,0.3,0.9);
//    legall->AddEntry(hTimeratioall,"Beta Run Data (All)","l");
//    legall->AddEntry(hSupsims,"Sim No. 2, flat","l");
//    legall->AddEntry(hSupalsims,"Sim No. 2","l");
//    legall->Draw();
//   timeoverlaysall->Print("output_files/timeratio_overlaysall.pdf");

  foct.close();
  fquart.close();
  fsup.close();

  for(Int_t ii = 0 ; ii < 9 ; ii++){
      Average_Array_Vector(A_Aves[ii],FullOct, A_ave_tot[ii], A_ave_toter[ii]);
      Average_Array_Vector(A_Aves_s[ii],FullOct, A_aves_tot[ii], A_aves_toter[ii]);
  }

  // Output asymmetry averages;
  fstream aveGout;
  aveGout.open("output_files/check_on_g.txt",fstream::out);
  for(int i = 0; i < FullOct ; i++){
     aveGout << A_Aves[6].A_ave[i] << " +/- " << A_Aves[6].A_error[i] << "\t\t"  << octet[i]->A_sum_A[6];
     aveGout << " +/- "   << octet[i]->A_sum_Aer[6]  << "\t\t"  << octet[i]->A_sum_B[6];
     aveGout << " +/- " << octet[i]->A_sum_Ber[6];
     aveGout << "\t\t"    << octet[i]->A_multi_A[6]  << " +/- " << octet[i]->A_multi_Aer[6];
     aveGout << "\t\t"    << octet[i]->A_multi_B[6]  << " +/- " << octet[i]->A_multi_Ber[6] << endl;
  }
  aveGout.close();
  //------------------------------------------------------------------------------------------------  
  Double_t xchoice[9];
  for(Int_t i = 0; i < 9 ; i++){
    xchoice[i] = i+1;
    cout << "Average A for Analysis " << i << "  " << A_ave_tot[i] << " +/- " << A_ave_toter[i];
    cout << "  " << A_aves_tot[i] << " +/- " << A_aves_toter[i]<< endl;
  }

  TF1 *fbeta = new TF1("fbeta",PoverE,0,2000,1);
  fbeta->SetParameter(0,1.);
  Double_t avebeta = fbeta->Integral(225,675)/450.;
  
  cout << "The average value for beta from 225,675 is " << avebeta << endl;
  cout << "Total Octet counts are : " << TotalCounts;

  for(Int_t i = 0; i < 9 ; i++){
    cout << "Average corrected A for Analysis " << i;
    cout << "  " << A_aves_tot[i]/(0.5*avebeta) << " +/- " << A_aves_toter[i]/(0.5*avebeta) <<  endl;
  }
  //----------------------------------------------------------------
  // Define the graphs to display the asymmetries ..................
  TGraphErrors *gsumA     = new TGraphErrors((int)A_sum_A.run_number.size()-1,&A_sum_A.run_number[0],&A_sum_A.A_ave[0],0,&A_sum_A.A_error[0]);
  TGraphErrors *gsumB     = new TGraphErrors((int)A_sum_B.run_number.size()-1,&A_sum_B.run_number[0],&A_sum_B.A_ave[0],0,&A_sum_B.A_error[0]);
  TGraphErrors *gmultiA   = new TGraphErrors((int)A_multi_A.run_number.size()-1,&A_multi_A.run_number[0],&A_multi_A.A_ave[0],0,&A_multi_A.A_error[0]);
  TGraphErrors *gmultiB   = new TGraphErrors((int)A_multi_B.run_number.size()-1,&A_multi_B.run_number[0],&A_multi_B.A_ave[0],0,&A_multi_B.A_error[0]);
  TGraphErrors *gOctAave  = new TGraphErrors((int)A_Aves[2].run_number.size()-1,&A_Aves[2].run_number[0],&A_Aves[2].A_ave[0],0,&A_Aves[2].A_error[0]);
  TGraphErrors *gSuperA   = new TGraphErrors((int)A_super.run_number.size()-1,&A_super.run_number[0],&A_super.A_ave[0],0,&A_super.A_error[0]);
  TGraphErrors *gQuartetA = new TGraphErrors((int)A_quat.run_number.size()-1,&A_quat.run_number[0],&A_quat.A_ave[0],0,&A_quat.A_error[0]);
  TGraph *gCnts           = new TGraph(noct-1,&OctDate[0],&CurrentCounts[0]);
  TGraphErrors *gAcho     = new TGraphErrors(9,xchoice,A_aves_tot,0,A_aves_toter);
  //-------------------------------------------------------------------------------
  //Push results to mpm db base for full comparison...
  //OutPutToMPMDB(noct);
  //---------------------------------------------------------------------------------
  ColorGraphic(gAcho,4,20,2,0.8,"Asymmetry by Analysis Choice","Analysis Choice","Raw Asymmetry");
  //-------------------------------------------------------------------------------------------------------------
  ColorGraphic(gsumA,2,20,2);
  ColorGraphic(gsumB,3,20,2);
  ColorGraphic(gmultiA,4,20,2);
  ColorGraphic(gmultiB,1,20,2);
  //--------------------------------------------------------------------------------------------------------------  
  ColorGraphic(gOctAave ,2,20,2,0.6,"Octet Integral Asymmetry 175 - 775 keV","Octet Number","Asym_{raw}");
  ColorGraphic(gSuperA  ,2,20,2,0.6,"Super-Ratio Integral Asymmetry 175 - 775 keV","Super-Ratio Number","Asym_{raw}");
  ColorGraphic(gQuartetA,2,20,2,0.6,"Quartet Integral Asymmetry 175 - 775 keV","Quartet Number","Asym_{raw}");
  // --------------------------------------------------------------------------------
  // Plot some shit..........
  Plot_RawAsymmetries(gOctAave,gQuartetA,gSuperA,FullOct,NoctC,NQuat);
  // 
  TCanvas *cdiff = new TCanvas("cdiff","Sum Product Difference");
  cdiff->Divide(1,3);
  cdiff->cd(1);
  gsumA->Draw("AP");
  gmultiA->Draw("P");
  cdiff->cd(2);
  gsumB->Draw("AP");
  gmultiB->Draw("P");
  cdiff->cd(3);
  ColorGraphic(gCnts,2,20,2);
  gCnts->GetXaxis()->SetTimeDisplay(1);
  gCnts->GetXaxis()->SetTimeFormat("%d / %m");
  gCnts->Draw("AP");
  cdiff->Print("output_files/sumproductdiffernece.pdf");

  // Plot the Asymmetry my analysis choice
  TCanvas *cAnalysisChoice = new TCanvas("cAnalysisChoice","Analysis Choice");
  cAnalysisChoice->cd();
  fstream fchmonte;
  fstream fchanalysis;
  fchmonte.open("input_files/analysischoice.txt",fstream::in);
  fchanalysis.open("output_files/analysis_choice_data.txt",fstream::out);
  Double_t MonteA[9];
  for(Int_t i = 0 ; i < 9 ; i++){
	fchmonte >> MonteA[i];
	MonteA[i] = TMath::Abs(MonteA[i]);
        fchanalysis << xchoice[i] << "\t" << A_aves_tot[i] << "\t" << A_aves_toter[i] << endl;
  }
  
  fchmonte.close();
  fchanalysis.close();

  TGraph *gMonte = new TGraph(9,xchoice,MonteA);
  ColorGraphic(gMonte,4,20,2);
  gAcho->Draw("AP");
  gMonte->Draw("P");
  gAcho->GetXaxis()->Set(1000,0,10);
  gAcho->GetXaxis()->SetBinLabel(100,"A");
  gAcho->GetXaxis()->SetBinLabel(200,"B");
  gAcho->GetXaxis()->SetBinLabel(300,"C");
  gAcho->GetXaxis()->SetBinLabel(400,"D");
  gAcho->GetXaxis()->SetBinLabel(500,"E");
  gAcho->GetXaxis()->SetBinLabel(600,"F");
  gAcho->GetXaxis()->SetBinLabel(700,"G");
  gAcho->GetXaxis()->SetBinLabel(800,"H");
  gAcho->GetXaxis()->SetBinLabel(900,"I");
  gAcho->GetXaxis()->LabelsOption("uh");
  gPad->SetGrid();
  cAnalysisChoice->Print("output_files/AnalysisChoice.pdf");
  /*
  delete gsumA;
  delete gsumB;
  delete fbeta;
  delete gmultiA;
  delete gmultiB;
  delete gOctAave;
  delete gSuperA;
  delete gQuartetA;
  delete gCnts;
  delete cdiff;
  */
  gStyle->SetOptStat(kFALSE); 
};

void OutPutToMPMDB(Int_t noct)
{

/*
  AnaCutSpec cuts;
  cuts.emin   = octet[0]->hAsyA[2]->GetBinCenter(nlow) - octet[0]->hAsyA[2]->GetBinWidth(nlow)/2.;
  cuts.emax   = octet[0]->hAsyA[2]->GetBinCenter(nhigh)+ octet[0]->hAsyA[2]->GetBinWidth(nhigh)/2.;
  cuts.radius = radcut;

  AnalysisDB *adb = AnalysisDB::getADB();
  AnaResult Asym;
  
  Asym.anatp     = AnaResult::ANA_ASYM;
  Asym.datp      = AnaResult::REAL_DATA;
  Asym.etypes.insert(TYPE_0_EVENT);
  Asym.etypes.insert(TYPE_I_EVENT);
  Asym.etypes.insert(TYPE_II_EVENT);
  Asym.etypes.insert(TYPE_III_EVENT);
 //------------------------------------------
  Asym.anach     = ANCHOICE_C;
  Asym.s         = BOTH;
  Asym.afp       = AFP_OTHER;
 //-----------------------------------------
 // Find and Delete Previous entries........
  vector<AnaResult> adbm = adb->findMatching(Asym);
  for(Int_t i = 0 ; i < (int)adbm.size(); i++)
       adb->deleteAnaResult(adbm[i].arid);
  // output the asymmetry results...........
//  for(AnalysisChoice j = ANCHOICE_C ; j <= ANCHOICE_C ; ++j){
  	for(Int_t ii = 0 ; ii < noct ; ii++){
		//------------------------------------------
		// Octets ----------------------------------
  		Asym.startRun  = octet[ii]->GetFirst();
  		Asym.endRun    = octet[ii]->GetLast();
		Asym.grouping  = AnaResult::GROUP_OCTET;
		//----------------------------------------
		// set and upload values, errors, and cuts
  		Asym.value     = octet[ii]->A_sum[2];
  		Asym.err       = octet[ii]->A_sumer[2];
  		Asym.csid      = adb->uploadCutSpec(cuts);
  		adb->uploadAnaResult(Asym);
		//=========================================
		Asym.grouping  = AnaResult::GROUP_QUARTET;
		Asym.startRun  = octet[ii]->GetQuartetAStart();
		Asym.endRun    = octet[ii]->GetQuartetAEnd();
		Asym.value     = octet[ii]->A_sum_A[2];
		Asym.err       = octet[ii]->A_sum_Aer[2];
		Asym.csid      = adb->uploadCutSpec(cuts);
		adb->uploadAnaResult(Asym);
		//=================================================
		Asym.startRun = octet[ii]->GetQuartetBStart();
		Asym.endRun   = octet[ii]->GetQuartetBEnd();
		Asym.value    = octet[ii]->A_sum_B[2];
		Asym.err      = octet[ii]->A_sum_Ber[2];
		Asym.csid      = adb->uploadCutSpec(cuts);
		adb->uploadAnaResult(Asym);
		//=========================================
		Asym.grouping = AnaResult::GROUP_PAIR;
		Asym.startRun = octet[ii]->GetSuperA1Start();
		Asym.endRun   = octet[ii]->GetSuperA1End();
		Asym.value    = octet[ii]->Asuper1[2];
		Asym.err      = octet[ii]->Asuper1er[2];
		Asym.csid      = adb->uploadCutSpec(cuts);
		adb->uploadAnaResult(Asym);
		//=================================================
		Asym.startRun = octet[ii]->GetSuperA2Start();
		Asym.endRun   = octet[ii]->GetSuperA2End();
		Asym.value    = octet[ii]->Asuper2[2];
		Asym.err      = octet[ii]->Asuper2er[2];
		Asym.csid      = adb->uploadCutSpec(cuts);
		adb->uploadAnaResult(Asym);
		//=================================================
		Asym.startRun = octet[ii]->GetSuperB1Start();
		Asym.endRun   = octet[ii]->GetSuperB1End();
		Asym.value    = octet[ii]->Bsuper1[2];
		Asym.err      = octet[ii]->Bsuper1er[2];
		Asym.csid      = adb->uploadCutSpec(cuts);
		adb->uploadAnaResult(Asym);
		//=================================================
		Asym.startRun = octet[ii]->GetSuperB2Start();
		Asym.endRun   = octet[ii]->GetSuperB2End();
		Asym.value    = octet[ii]->Bsuper2[2];
		Asym.err      = octet[ii]->Bsuper2er[2];
		Asym.csid      = adb->uploadCutSpec(cuts);
		adb->uploadAnaResult(Asym);

  	}
  //}
  */  
}
//---------------------------------------------------------------------------------------------------------------
void Fill_Asymmetry_Vector(Asym_t &A,Double_t oct,Double_t octer,Int_t &index)
{
    A.A_ave.push_back(TMath::Abs(oct));
    A.A_error.push_back(octer);
    A.run_number.push_back(index + 1.);
    index++;
};
//---------------------------------------------------------------------------------------------------------------
void Fill_Asymmetry_Vector_Oct(Asym_t &A,Double_t oct,Double_t octer,Int_t index)
{
    A.A_ave.push_back(TMath::Abs(oct));
    A.A_error.push_back(octer);
    A.run_number.push_back(index + 1.);
};
//----------------------------------------------------------------------------------------------------------------
void Plot_RawAsymmetries(TGraphErrors *g1,TGraphErrors *g2,TGraphErrors *g3,Int_t noct,Int_t NoctC,Int_t NQuat)
{
 
  TH1F *hOctResid = new TH1F("hOctResid","Octet Residuals; #Delta A / A ; Cnts",20,-10,10);
  TH1F *hQutResid = new TH1F("hQutResid","Quartet Residuals; #Delta A / A ; Cnts",20,-10,10);
  TH1F *hSupResid = new TH1F("hSupResid","Super Ratio Residuals; #Delta A / A ; Cnts",20,-10,10);
 
  Double_t Xasy,Asy,Asyer;

  //----------------------------------------------------------------------------------------------
  TCanvas *cass = new TCanvas("cass","Octet");
  cass->cd();
  // Define Pads
  TPad *pOctV = new TPad("pOctV","Octet Values",0,0.66,0.75,1.00);
  TPad *pQutV = new TPad("pQutV","Octet Values",0,0.33,0.75,0.66);
  TPad *pSupV = new TPad("pSupV","Octet Values",0,0.00,0.75,0.33);
  TPad *pOctR = new TPad("pOctR","Octet Values",0.75,0.66,1.,1.00);
  TPad *pQutR = new TPad("pQutR","Octet Values",0.75,0.33,1.,0.66);
  TPad *pSupR = new TPad("pSupR","Octet Values",0.75,0.00,1.,0.33);
  //-----------------------------------------------------------------------------------------------
  pOctV->Draw();
  pQutV->Draw();
  pSupV->Draw();
  pOctR->Draw();
  pQutR->Draw();
  pSupR->Draw();
  //-----------------------------------------------------------------------------------------------
  TF1 *flO = new TF1("flO","[0]",0,noct  + 1.);
  TF1 *flS = new TF1("flS","[0]",0,NoctC + 1.);
  TF1 *flQ = new TF1("flQ","[0]",0,NQuat + 1.);
  //-----------------------------------------------------------------------------------------------
  g1->Fit("flO","RMEQ","goff");
  g2->Fit("flQ","RMEQ","goff");
  g3->Fit("flS","RMEQ","goff");
  //-----------------------------------------------------------------------------------------------
  // Fill Residual histograms
  for(Int_t i = 0 ; i < noct ; i++){
        g1->GetPoint(i,Xasy,Asy);
        Asyer = g1->GetErrorY(i);
        hOctResid->Fill( (Asy - flO->GetParameter(0))/Asyer,1);
  }
  for(Int_t i = 0 ; i < NQuat ; i++){
        g2->GetPoint(i,Xasy,Asy);
        Asyer = g2->GetErrorY(i);
        hQutResid->Fill( (Asy - flQ->GetParameter(0))/Asyer,1);
  }
  for(Int_t i = 0 ; i < NoctC ; i++){
        g3->GetPoint(i,Xasy,Asy);
        Asyer = g3->GetErrorY(i);
        hSupResid->Fill( (Asy - flS->GetParameter(0))/Asyer,1);
  }
  // Draw the OctetAsymmetries----------------------------------------------------------------------
  pOctV->cd();
  pOctV->SetMargin(0.1,0.05,0.1,0.1);
  g1->GetYaxis()->SetTitleSize(0.08);
  g1->GetYaxis()->SetTitleOffset(0.5);
  g1->GetYaxis()->SetRangeUser(0.03,0.07);
  g1->Draw("AP");
  g1->Fit("flO","RMEQ","");
  //------------------------------------------------------------------------------------------
  TLegend *legO = new TLegend(0.5,0.7,0.95,0.9);
  legO->AddEntry(g1,Form("A = %6.5f #pm %6.5f",flO->GetParameter(0),flO->GetParError(0)),"lp");
  legO->AddEntry(g1,Form("#chi^{2}/#nu = %6.5f",flO->GetChisquare()/flO->GetNDF()),"");
  legO->AddEntry(g1,Form("Probability = %6.5f",flO->GetProb()),"");
  legO->SetFillColor(0);
  legO->Draw();
  pOctR->cd();
  hOctResid->Draw();
  // Draw the Quartet Asymmetries
  pQutV->cd();
  pQutV->SetMargin(0.1,0.05,0.1,0.1);
  g2->GetYaxis()->SetRangeUser(0.03,0.07);
  g2->GetYaxis()->SetTitleSize(0.08);
  g2->GetYaxis()->SetTitleOffset(0.5);
  g2->Draw("AP");
  g2->Fit("flQ","RMEQ","");
  TLegend *legQ = new TLegend(0.5,0.7,0.95,0.9);
  legQ->AddEntry(g2,Form("A = %6.5f #pm %6.5f",flQ->GetParameter(0),flQ->GetParError(0)),"lp");
  legQ->AddEntry(g2,Form("#chi^{2}/#nu = %6.5f",flQ->GetChisquare()/flQ->GetNDF()),"");
  legQ->AddEntry(g2,Form("Probability = %6.5f",flQ->GetProb()),"");
  legQ->SetFillColor(0);
  legQ->Draw();
  pQutR->cd();
  hQutResid->Draw();
  // Draw the Super Ratio Asymmetries
  pSupV->cd();
  pSupV->SetMargin(0.1,0.05,0.1,0.1);
  g3->GetYaxis()->SetRangeUser(0.03,0.07);
  g3->GetYaxis()->SetTitleSize(0.08);
  g3->GetYaxis()->SetTitleOffset(0.5);
  g3->Draw("AP");
  g3->Fit("flS","RMEQ","");
  TLegend *legS = new TLegend(0.5,0.7,0.95,0.9);
  legS->AddEntry(g3,Form("A = %6.5f #pm %6.5f",flS->GetParameter(0),flS->GetParError(0)),"lp");
  legS->AddEntry(g3,Form("#chi^{2}/#nu = %6.5f",flS->GetChisquare()/flS->GetNDF()),"");
  legS->AddEntry(g3,Form("Probability = %6.5f",flS->GetProb()),"");
  legS->SetFillColor(0);
  legS->Draw();
  pSupR->cd();
  hSupResid->Draw();
  // output the pdf

  cass->Print("output_files/asymmetries_out.pdf");
  cass->Print("output_files/asymmetries_out.C");
  /*
  delete flO;
  delete flQ;
  delete flS;
  delete legS;
  delete cass;
  */
}
//---------------------------------------------------------------------------------------------------------------------
void Average_A()
{
  using namespace TMath;
   
  vector<Double_t> x (200),y (200),yer (200);
  vector<Double_t> yb (200),yber (200);
  vector<Double_t> ytot (200),ytoter (200);
  vector<Double_t> y1 (200),y1er (200);
  vector<Double_t> y1b (200),y1ber (200);
  vector<Double_t> y1tot (200),y1toter (200);
  
  TF1 *fbeta   = new TF1("fbeta",PoverE,0,2000,1);
  fbeta->SetParameter(0,1.);
  TF1 *fbetaf  = new TF1("fbetaf",PoverE,200,700,1);
  TF1 *fbetaf2 = new TF1("fbetaf2",PoverE,0,2000,1);
  // Calculate the weighted average of the Asymmetry
  for(Int_t j = 0 ; j < noct ; j++){
//    if(octet[j]->nA2[2] != 0 || octet[j]->nA5[2] != 0 )
      cout << j << "\t" << octet[j]->GetFirst() << endl;
      if(octet[j]->GetFirst() < 14888){
	Average_All_Hists(octet[j]->hAsyA[2],y,yer);
        Average_All_Hists(octet[j]->hAsyB[2],yb,yber);
        Average_All_Hists(octet[j]->hAsyTot[2],ytot,ytoter);
      } else if(octet[j]->GetFirst() > 14887){
	Average_All_Hists(octet[j]->hAsyA[2],y1,y1er);
        Average_All_Hists(octet[j]->hAsyB[2],y1b,y1ber);
        Average_All_Hists(octet[j]->hAsyTot[2],y1tot,y1toter);
      }
  }
  // Set the bin center for all graphs
  for(Int_t i = 0 ; i <= octet[0]->hAsyA[2]->GetNbinsX()-1 ;i++)
       x[i]      = octet[0]->hAsyA[2]->GetBinCenter(i+1);
  //------------------------------------------------------------------------------
  // Invert the error weighted average and apply physics to the average raw 
  // asymmetries
  Return_Asymmetry(y,yer,octet[0]->hAsyA[2]->GetNbinsX());
  Return_Asymmetry(yb,yber,octet[0]->hAsyB[2]->GetNbinsX());
  Return_Asymmetry(ytot,ytoter,octet[0]->hAsyTot[2]->GetNbinsX());
  //============================================================================
  Return_Asymmetry(y1,y1er,octet[0]->hAsyA[2]->GetNbinsX());
  Return_Asymmetry(y1b,y1ber,octet[0]->hAsyB[2]->GetNbinsX());
  Return_Asymmetry(y1tot,y1toter,octet[0]->hAsyTot[2]->GetNbinsX());
  //============================================================================
  // Output to a text file the results 
  fstream fasye,fasye2;
  if(btr[0]->GetRunNumber() < 12000){
    fasye.open(Form("output_files/asymmetry_e_%d.txt",btr[0]->GetGeo()),fstream::out);
  } else {

    fasye.open("output_files/asymmetry_e_a0_16.txt",fstream::out);
    fasye2.open("output_files/asymmetry_e_a17_50.txt",fstream::out);
  }
  for(Int_t i = 0 ; i< octet[0]->hAsyTot[2]->GetNbinsX() ; i++){
    fasye << x[i] << "\t" << ytot[i] << "\t" << ytoter[i] << "\t" << y[i];
    fasye << "\t" << yer[i] << "\t" << yb[i] << "\t" << yber[i] << endl;
    fasye2 << x[i] << "\t" << y1tot[i] << "\t" << y1toter[i] << "\t" << y1[i];
    fasye2 << "\t" << y1er[i] << "\t" << y1b[i] << "\t" << y1ber[i] << endl;
  }
  fasye.close();
  fasye2.close();
  //------------------------------------------------------------------
  TGraphErrors *gAvsE    = new TGraphErrors((int)x.size(),&x[0],&y1[0],0,&y1er[0]);
  TGraphErrors *gAvsEb   = new TGraphErrors((int)x.size(),&x[0],&y1b[0],0,&y1ber[0]);
  TGraphErrors *gAvsEtot = new TGraphErrors((int)x.size(),&x[0],&y1tot[0],0,&y1toter[0]);

  ColorGraphic(gAvsE,2,20,2,1);
  ColorGraphic(gAvsEb,4,20,2,1);
  ColorGraphic(gAvsEtot,2,20,2,1,"Octet Asymmetry vs. Energy in 10 keV bins","Reconstructed Energy (keV)","Measured Asymmetry");

  TCanvas *cGreat = new TCanvas("cGreat","Energy vs. A",500,800);
  cGreat->Divide(1,2);
  cGreat->cd(1);
  
  TF1 *fAbeta = new TF1("fAbeta","[0]",225,675);
  TF1 *fAbeta2 = new TF1("fAbeta2","[0]",225,675);
  
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
  cGreat->Print("output_files/asymmetry_v_e.pdf");

  TH1F *hAsymErrs = new TH1F("hAsymErrs",
			     "Asymmetry Absolute Statistical Errors;Reconstructed Energy (keV);Absolute Error",
			     octet[0]->hAsyTot[2]->GetNbinsX(),0.,2000.);

  TCanvas *cAveErrs = new TCanvas("cAveErrs","Asymmetry Errors");
  vector<Double_t> A_per_errors,X_per;

  for(Int_t i = 0 ; i < 80 ; i++){
    X_per.push_back(0.);
    if(ytot[i] != 0)
      A_per_errors.push_back(100.*Abs(ytoter[i]/ytot[i]));
    else 
      A_per_errors.push_back(1.);
    
    hAsymErrs->SetBinContent(i,A_per_errors[i]);
    
  }  
  
  ColorGraphic(hAsymErrs,2,20,2);
  hAsymErrs->Draw("");
  
  cAveErrs->Print("output_files/errors_macro.C");
  cAveErrs->Print("output_files/errors_macro.pdf");
  /*
  delete hAsymErrs;
  delete cAveErrs; 
  delete cGreat;*/
}

void Collect_TvsE()
{

  // Calculate the average Energy vs. dt for the type 1 backscatters 
  // then compare it with the results from Monte Claro.

  const Int_t rbins = 100;
  Double_t ecount[rbins][rbins],wcount[rbins][rbins];
  
  hTotETimeVsE = new TH2F("hTotETimeVsE",
			  "East Type 1 E_{vis} vs. Timing ; E_{vis}(keV) ;"
			  "TDC Timing (ns)",rbins,0,1000,rbins,0,200);
  
  hTotWTimeVsE = new TH2F("hTotWTimeVsE",
			  "West Type 1 E_{vis} vs. Timing ; E_{vis} (keV) ;"
			  "TDC Timing (ns)",rbins,0,1000,rbins,0,200);
  
  // Fill the Collection 2d histogram

  for(Int_t ibin = 0 ; ibin < rbins ; ibin++){
    for(Int_t jbin = 0 ; jbin < rbins ; jbin++){
      ecount[ibin][jbin] = 0.;
      wcount[ibin][jbin] = 0.;
    }
  }
  
  for(Int_t i = 0 ; i <nbeta ; i++){
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    for(Int_t ibin = 1 ; ibin < rbins ; ibin++){
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
        if(btr[i]->hETimeVsE->GetBinContent(ibin,jbin) > 0)
	  ecount[ibin-1][jbin-1] += btr[i]->hETimeVsE->GetBinContent(ibin,jbin)*btr[i]->rtime_e;
	if(btr[i]->hWTimeVsE->GetBinContent(ibin,jbin) > 0)
	  wcount[ibin-1][jbin-1] += btr[i]->hWTimeVsE->GetBinContent(ibin,jbin)*btr[i]->rtime_w;
      }
    }
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  
  for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
    for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
      hTotETimeVsE->SetBinContent(ibin,jbin,ecount[ibin-1][jbin-1]);
      hTotWTimeVsE->SetBinContent(ibin,jbin,wcount[ibin-1][jbin-1]);
    }
  }
  // Finished Filling now compare...

  TCanvas *check = new TCanvas("check","check");
  check->Divide(3,2);
  check->cd(1);
  hTotETimeVsE->Draw("colz");
  check->cd(2);
  hTotWTimeVsE->Draw("colz");
  
  fstream fback;
  fback.open("output_files/backscatter_e_vs_t.txt",fstream::out);
  
  Double_t xe[100],ye[100],xw[100],yw[100],yet,ywt;
  Double_t yer[100],ywr[100];

  for(Int_t i = 0 ; i < rbins; i++){
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
    fbck2.open("input_files/evt_ave.txt",fstream::in);
    Toffset  = 0.;
    Toffsete = 5.; 
  } else if(btr[0]->GetGeo() == 1){
    fbck2.open("input_files/evt_ave.txt",fstream::in);
    Toffset = 10.;
  }
  else if(btr[0]->GetGeo() == 2 || btr[0]->GetGeo() == 3 ){
    fbck2.open("input_files/evt_ave.txt",fstream::in);
    Toffset = -4;
    Toffsete = 17.;
  }
 
  Int_t ntbins = 160;
  Double_t xmc[ntbins],emc[ntbins],emcer[ntbins],xmcw[ntbins];
  Double_t Ediff[ntbins];

  for(Int_t i = 0 ; i < ntbins ; i++){
    fbck2 >> xmc[i] >> emc[i] >> emcer[i];
 //   cout << "Monte time " << xmc[i] << " data time " << xe[i] << endl;
    xmcw[i] = xmc[i] - Toffset;
    xmc[i]  = xmc[i] + Toffsete;
  } 
  fbck2.close();

  for(Int_t i = 0 ; i < 100 ; i++){
       	if(ye[i] > 0){
		Ediff[i] = (ye[i] - emc[i+(Int_t)Toffset])/(sqrt(yer[i]*yer[i]+TMath::Power(emcer[i+(Int_t)Toffset],2)));
	}else 
  		Ediff[i] = 0.;
  }
  
  TGraphErrors *gMcTe  = new TGraphErrors(ntbins,xmc,emc,0,emcer);
  TGraphErrors *gMcTw  = new TGraphErrors(ntbins,xmcw,emc,0,emcer);
  TGraphErrors *gAveTe = new TGraphErrors(100,xe,ye,0,yer);
  TGraphErrors *gAveTw = new TGraphErrors(100,xw,yw,0,ywr);
  TGraph *gEdiffs = new TGraph(100,xe,Ediff);
  
  ColorGraphic(gAveTe,2,20,2,0.6,"East Average Detected Evis","TDC Timing (ns)","#LT E_{vis} #GT");
  ColorGraphic(gAveTw,2,20,2,0.6,"West Average Detected Evis","TDC Timing (ns)","#LT E_{vis} #GT");
  ColorGraphic(gMcTe,4,7,2,0.6);
  ColorGraphic(gMcTw,4,7,2,0.6); 
  ColorGraphic(gEdiffs,4,7,2);

  check->cd(3);
  gAveTe->Draw("AP");
  gAveTe->GetXaxis()->SetRangeUser(20,170);
  gMcTe->Draw("l");
  check->cd(4);
  gAveTw->Draw("AP");
  gAveTw->GetXaxis()->SetRangeUser(20,170);
  gMcTw->Draw("l");
  check->cd(5);
  gEdiffs->Draw("AP");

  check->Print("output_files/energy_vs_timing.pdf");
  /*
  delete gMcTe;
  delete gAveTe;
  delete gMcTw;
  delete gAveTw;
  delete check;*/
}
//----------------------------------------------------------------------------------------------------------
void Collect_23Anode()
{
   hTE23Anode2d = new TH2F("hTE23Anode2d",
		     "East Type 23 E vs. W Anode; East Anode ; West Anode"
		     ,100,0,3000,100,0,3000);
		     
   hTW23Anode2d = new TH2F("hTW23Anode2d",
		     "West Type 23 E vs. W Anode; East Anode ; West Anode"
		     ,100,0,3000,100,0,3000);
  
  const Int_t rbins = 150;     
  Double_t ecount[rbins][rbins],wcount[rbins][rbins];
  CleanArray(ecount,rbins,rbins);
  CleanArray(wcount,rbins,rbins);
 
  for(Int_t i = 0 ; i <nbeta ; i++){
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
      for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
        if(btr[i]->hE23Anode2d->GetBinContent(ibin,jbin) > 0)
	ecount[ibin-1][jbin-1] += btr[i]->hE23Anode2d->GetBinContent(ibin,jbin)*btr[i]->rtime_e;
	if(btr[i]->hW23Anode2d->GetBinContent(ibin,jbin) > 0)
	wcount[ibin-1][jbin-1] += btr[i]->hW23Anode2d->GetBinContent(ibin,jbin)*btr[i]->rtime_w;
      }
    }
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  
  for(Int_t ibin = 1 ; ibin < rbins+1 ; ibin++){
    for(Int_t jbin = 1 ; jbin < rbins+1 ; jbin++){
      hTE23Anode2d->SetBinContent(ibin,jbin,ecount[ibin-1][jbin-1]);
      hTW23Anode2d->SetBinContent(ibin,jbin,wcount[ibin-1][jbin-1]);
    }
  }
  
  TCanvas *cAnode = new TCanvas("cAnode","cAnode");
  cAnode->Divide(2,2);
  cAnode->cd(1);
  hTE23Anode2d->Draw("colz");
  cAnode->cd(2);
  hTW23Anode2d->Draw("colz");

  TH1F *hEAnode23Tot = new TH1F("hEAnode23Tot","East ",100,0,20); 
  TH1F *hWAnode23Tot = new TH1F("hWAnode23Tot","East ",100,0,20);
  
  for(Int_t i = 0 ; i < nbeta ; i++){
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
	for(Int_t j = 0 ; j < hEAnode23Tot->GetNbinsX() ; j++){
	hEAnode23Tot->SetBinContent(j,hEAnode23Tot->GetBinContent(j) + btr[i]->hEAnode23->GetBinContent(j)*btr[i]->rtime_e);
	hWAnode23Tot->SetBinContent(j,hWAnode23Tot->GetBinContent(j) + btr[i]->hWAnode23->GetBinContent(j)*btr[i]->rtime_w);
  	}
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }

  cAnode->cd(3);
  hEAnode23Tot->Draw();
  cAnode->cd(4);
  hWAnode23Tot->Draw();

  fstream fanode;
  fanode.open("output_files/anode_23_spec.txt",fstream::out);
 
  for(Int_t i = 1 ; i < hEAnode23Tot->GetNbinsX() ; i++){
	fanode << hEAnode23Tot->GetBinCenter(i) << "\t" << hEAnode23Tot->GetBinContent(i) << "\t" << hWAnode23Tot->GetBinContent(i) << endl;
  }
  fanode.close();
  cAnode->Print("output_files/anode.pdf");
  /*
  delete hEAnode23Tot;
  delete hWAnode23Tot;
  delete cAnode;
  */
  
}
//---------------------------------------------------------------------------
void Collect_Stuff()
{

  TCanvas *cStN = new TCanvas("cStN","Signal To Noise");
  cStN->cd(0);
  Double_t xS2N[MAXRUNS],yS2N[MAXRUNS],yS2Nw[MAXRUNS];
  
  
  for(int i = 0 ; i < nbeta ; i++){
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
       xS2N[i] = btr[i]->GetRunNumber();
       yS2N[i] = btr[i]->E_Sig_Nos;
       yS2Nw[i] = btr[i]->W_Sig_Nos;
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  
  TGraph *gS2Ne = new TGraph(nbeta,xS2N,yS2N);
  ColorGraphic(gS2Ne,2,20,2);
  TGraph *gS2Nw = new TGraph(nbeta,xS2N,yS2Nw);
  ColorGraphic(gS2Nw,4,20,2);
  
  gS2Ne->Draw("AP");
  gS2Nw->Draw("P");

  cStN->Print("output_files/SigToNoise.pdf");
  /*
  delete cStN;
  delete gS2Ne;
  delete gS2Nw;
  */ 
};
//-------------------------------------------------------------------------------
void Collect_Rad()
{

  using namespace TMath;
  
  Double_t x[20],y[20],yer[20];
  
  // average asymmetry as a function of radius from each octet....

  for(Int_t i = 0; i < 12 ;i++){
    x[i]      = 0.;
    y[i]      = 0.;
    yer[i]    = 0.;
  } 
  
  for(Int_t j = 0 ; j < noct ; j++){
    if(!(IsNaN(octet[j]->Asuper2[2]))){
      for(Int_t i = 0; i < 12 ; i++ ){
	// Insure that the asymmetry is calculated for that bin
	if(octet[j]->A_rader[i] != 1 && !(IsNaN(octet[j]->A_rader[i])) && octet[j]->A_rader[i]>0 ){
          // sum the values divided by the error squared...
	  //  cout << j << " " << i << "  "  << octet[j]->A_rad[i] << " +/- " << octet[j]->A_rader[i] << endl; 
	  y[i] += TMath::Abs(octet[j]->A_rad[i] / Power(octet[j]->A_rader[i],2));
	  // sum the weights......
	  yer[i] += 1./Power(octet[j]->A_rader[i],2);
	}
      }
    }
  }
  
  for(Int_t i = 0 ; i < 12 ;i++){
    x[i]   = (i+1)*5;

    if(yer[i] != 0.){
      y[i]   = y[i] / yer[i];
      yer[i] = sqrt( 1./yer[i]);
    }
    if(IsNaN(y[i])){
	y[i] = 0.0;
	yer[i] = 0.1;
    }
  }
  
  
  TGraphErrors *gAvsRad = new TGraphErrors(12,x,y,0,yer);
  ColorGraphic(gAvsRad,2,20,2,1.,"Raw Octet Asymmetry for R_{i} < R < R_{i+1}","Radius (mm)","A_{raw}");
  
  TF1 *frad1 = new TF1("frad1","[0]",0,60);
  frad1->SetLineColor(2);

  TCanvas *cRadA = new TCanvas("cRadA","Radi Dependent");
  cRadA->cd();

  gAvsRad->GetYaxis()->SetRangeUser(0.02,0.08);
  gAvsRad->Draw("AP");
  gAvsRad->GetYaxis()->SetRangeUser(0.03,0.06);
  gAvsRad->Fit("frad1","REMQ");

  fstream fradout;
  fradout.open("output_files/asym_rad_out.txt",fstream::out);
  for(Int_t i = 0 ; i < 12;i++){
	fradout << x[i] << "\t" << y[i] << "\t" << yer[i] << endl;
  }
  fradout.close();

  TLegend *lrado = new TLegend(0.55,0.75,0.9,0.9);
  lrado->AddEntry(gAvsRad,Form("A = %6.4f #pm %6.4f",frad1->GetParameter(0),
						    frad1->GetParError(0)),"lp");
  lrado->AddEntry(gAvsRad,Form("#chi^{2}/#nu = %6.4f",frad1->GetChisquare()/frad1->GetNDF()),"");
  lrado->SetFillColor(0);
  lrado->Draw();
  cRadA->Print("output_files/Radial_Dependence.pdf");

  DrawRadialCounts(yer,frad1->GetParameter(0));
  /*
  delete cRadA;
  delete lrado;
  delete gAvsRad;
  delete frad1;
  */
};
//-----------------------------------------------------------------------------------------------------
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
  Double_t missedcounts   = 0.;
  Double_t tloss          = 0.;
  //
  // Loop through and fill histograms for the signal and back ground spectra
  //
  Fill_Foreground(&west_time_on,&east_time_on,&west_time_off,&east_time_off,&tloss,&missedcounts);
  Fill_Background(&west_time_onb,&east_time_onb,&west_time_offb,&east_time_offb);
  // Sum total time
  Double_t totaltime = (west_time_on + east_time_on)/2. + (west_time_off + east_time_off)/2.;
  
  // Set error bars in the collected histograms
  fstream parbck;
  parbck.open(Form("output_files/back_ground_%d_on.txt",btr[0]->GetGeo()),fstream::out);
  fstream parbckoff;
  parbckoff.open(Form("output_files/back_ground_%d_off.txt",btr[0]->GetGeo()),fstream::out);
  parbck << "E bin\tRe\tRe_er\tRw\tRw_er\tRe_1\tRe_1er\tRw_1\tRw_1er\tRe_2\tRe_2er\tRw_2\tRw_2er\t" << endl;
  parbckoff << "E bin\tRe\tRe_er\tRw\tRw_er\tRe_1\tRe_1er\tRw_1\tRw_1er\tRe_2\tRe_2er\tRw_2\tRw_2er\t" << endl;

  fstream parfor;
  parfor.open(Form("output_files/fore_ground_%d_on.txt",btr[0]->GetGeo()),fstream::out);
  fstream parforoff;
  parforoff.open(Form("output_files/for_ground_%d_off.txt",btr[0]->GetGeo()),fstream::out);
  parfor << "E bin\tRe\tRe_er\tRw\tRw_er\tRe_1\tRe_1er\tRw_1\tRw_1er\tRe_2\tRe_2er\tRw_2\tRw_2er\t" << endl;
  parforoff << "E bin\tRe\tRe_er\tRw\tRw_er\tRe_1\tRe_1er\tRw_1\tRw_1er\tRe_2\tRe_2er\tRw_2\tRw_2er\t" << endl;
  
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
  
  //-----------------------------------------------------------------------------  
  cout << "East On Type 1 "   << Etype1on  << " +/- " << erEt1on <<   endl;
  cout << "East On Type 2/3 " << Etype23on << " +/- " << erEt23on << endl;
			      
  cout << "East Off Type 1 "   << Etype1off << " +/- " << erEt1off << endl;
  cout << "East Off Type 2/3 " << Etype23off<< " +/- " << erEt23off << endl;
			      
  cout << "West On Type 1 "   << Wtype1on << " +/- " << erWt1on<< endl;
  cout << "West On Type 2/3 " << Wtype23on << " +/- " << erWt23on<< endl;
			      
  cout << "West Off Type 1 "   << Wtype1off  << " +/- " << erWt1off<< endl;
  cout << "West Off Type 2/3 " << Wtype23off << " +/- " << erWt23off<< endl;

  cout << "-------------------------------------------------------------" << endl; 
  cout << " energy range is " << hEFlipperOn->GetBinCenter(nlow) << " - " << hEFlipperOn->GetBinCenter(nhigh) << endl;
  cout << "Total East Counts On  (0-800keV) " << hEFlipperOn->Integral(nlow,nhigh);
  cout << " Bck (0-800keV) " << hEFlipperOn_B->Integral(nlow,nhigh);
  cout << " Sig / Noise " << hEFlipperOn->Integral(nlow,nhigh)/ hEFlipperOn_B->Integral(nlow,nhigh) << endl;

  cout << "Total East Counts Off (0-800keV) " << hEFlipperOff->Integral(nlow,nhigh);
  cout << " Bck (0-800keV) " << hEFlipperOff_B->Integral(nlow,nhigh);
  cout << " Sig / Noise " << hEFlipperOff->Integral(nlow,nhigh)/ hEFlipperOff_B->Integral(nlow,nhigh) << endl;

  cout << "Total West Counts On  (0-800keV) " << hWFlipperOn->Integral(nlow,nhigh);
  cout << " Bck (0-800keV) " << hWFlipperOn_B->Integral(nlow,nhigh);
  cout << " Sig / Noise " << hWFlipperOn->Integral(nlow,nhigh)/ hWFlipperOn_B->Integral(nlow,nhigh)<< endl;

  cout << "Total West Counts Off (0-800keV) " << hWFlipperOff->Integral(nlow,nhigh);
  cout << " Bck (0-800keV) " << hWFlipperOff_B->Integral(nlow,nhigh);
  cout << " Sig / Noise " << hWFlipperOff->Integral(nlow,nhigh)/ hWFlipperOff_B->Integral(nlow,nhigh) << endl;

  cout << "Total Counts (0-800keV) " << (hEFlipperOn->Integral(nlow,nhigh)*east_time_on + 
    hEFlipperOff->Integral(nlow,nhigh)*east_time_off + hWFlipperOn->Integral(nlow,nhigh)*west_time_on + hWFlipperOff->Integral(nlow,nhigh)*west_time_off) << endl;

  cout << " Missed Time " << tloss << endl;
  cout << "Possible Missed Counts " << missedcounts << endl;
  cout << " Run time " << totaltime << endl;
  cout << " East time " << east_time_on << " off  " << east_time_off << endl;
  cout << " West Time " << west_time_on << " off  " << west_time_off << endl;
  
  
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
    
    parfor << i << "\t" << setprecision(4) << hEFlipperOn->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOn->GetBinError(i) << "\t";
    parfor << setprecision(4) << hWFlipperOn->GetBinContent(i)   << "\t"<< setprecision(4) << hWFlipperOn->GetBinError(i)   << "\t";
    parfor << setprecision(4) << hEFlipperOn_I->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOn_I->GetBinError(i) << "\t";
    parfor << setprecision(4) << hWFlipperOn_I->GetBinContent(i) << "\t"<< setprecision(4) << hWFlipperOn_I->GetBinError(i) << "\t";
    parfor << setprecision(4) << hEFlipperOn_2->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOn_2->GetBinError(i) << "\t";
    parfor << setprecision(4) << hWFlipperOn_2->GetBinContent(i) << "\t"<< setprecision(4) << hWFlipperOn_2->GetBinError(i) << endl;
    
    parforoff << i << "\t" << setprecision(4) << hEFlipperOff->GetBinContent(i) << "\t"<< hEFlipperOff->GetBinError(i) << "\t";
    parforoff << setprecision(4) << hWFlipperOff->GetBinContent(i)   << "\t"<< setprecision(4) << hWFlipperOff->GetBinError(i) << "\t";
    parforoff << setprecision(4) << hEFlipperOff_I->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOff_I->GetBinError(i) << "\t";
    parforoff << setprecision(4) << hWFlipperOff_I->GetBinContent(i) << "\t"<< setprecision(4) << hWFlipperOff_I->GetBinError(i) << "\t";
    parforoff << setprecision(4) << hEFlipperOff_2->GetBinContent(i) << "\t"<< setprecision(4) << hEFlipperOff_2->GetBinError(i) << "\t";
    parforoff << setprecision(4) << hWFlipperOff_2->GetBinContent(i) << "\t"<< setprecision(4) << hWFlipperOff_2->GetBinError(i) << endl;
  }
  
  parbckoff.close();
  parbck.close();
  
  parforoff.close();
  parfor.close();
   
  TCanvas *c1bk = new TCanvas("c1bk","One Clock back?",700,800);
  c1bk->Divide(3,4);
  TLegend *legSpc1,*legSpc2,*legSpc3,*legSpc4,*legSpc5,*legSpc6,*legSpc7,*legSpc8;
  TLegend *legSpc9,*legSpc10,*legSpc11,*legSpc12;


    btr[1]->Load_Histograms(bckr[btr[1]->Bkg_index],0);
    btr[3]->Load_Histograms(bckr[btr[3]->Bkg_index],0);
    btr[4]->Load_Histograms(bckr[btr[4]->Bkg_index],0);
    btr[5]->Load_Histograms(bckr[btr[5]->Bkg_index],0);
    btr[6]->Load_Histograms(bckr[btr[6]->Bkg_index],0);

  DrawEnerPanel(hWFlipperOff  ,btr[3]->hEERef , c1bk,1,legSpc1,west_time_off);
  DrawEnerPanel(hWFlipperOff_I,btr[3]->hEERef1, c1bk,2,legSpc2,west_time_off);
  DrawEnerPanel(hWFlipperOff_2,btr[3]->hEERef2, c1bk,3,legSpc3,west_time_off);
  DrawEnerPanel(hWFlipperOn   ,btr[4]->hEERef , c1bk,4,legSpc4,west_time_on);
  DrawEnerPanel(hWFlipperOn_I ,btr[4]->hEERef1, c1bk,5,legSpc5,west_time_on);
  DrawEnerPanel(hWFlipperOn_2 ,btr[4]->hEERef2, c1bk,6,legSpc6,west_time_on);
  DrawEnerPanel(hEFlipperOff  ,btr[5]->hEERef , c1bk,7,legSpc7,east_time_off);
  DrawEnerPanel(hEFlipperOff_I,btr[5]->hEERef1, c1bk,8,legSpc8,east_time_off);
  DrawEnerPanel(hEFlipperOff_2,btr[5]->hEERef2, c1bk,9,legSpc9,east_time_off);
  DrawEnerPanel(hEFlipperOn   ,btr[6]->hEERef , c1bk,10,legSpc10,east_time_on);
  DrawEnerPanel(hEFlipperOn_I ,btr[6]->hEERef1, c1bk,11,legSpc11,east_time_on);
  DrawEnerPanel(hEFlipperOn_2 ,btr[6]->hEERef2, c1bk,12,legSpc12,east_time_on);
  
  btr[1]->hEERef->Scale(hEFlipperOn->Integral(1,80)/btr[1]->hEERef->Integral(1,80));
  btr[1]->hEERef1->Scale(1./btr[1]->hEERef1->Integral(1,80));
  btr[1]->hEERef2->Scale(1./btr[1]->hEERef2->Integral(1,80));

    btr[1]->Remove_Histograms(bckr[btr[1]->Bkg_index]);
    btr[3]->Remove_Histograms(bckr[btr[3]->Bkg_index]);
    btr[4]->Remove_Histograms(bckr[btr[4]->Bkg_index]);
    btr[5]->Remove_Histograms(bckr[btr[5]->Bkg_index]);
    btr[6]->Remove_Histograms(bckr[btr[6]->Bkg_index]);

  c1bk->Print("output_files/e_spec_compares.pdf");

//   TCanvas *cTot = new TCanvas("cTot","Total Rates");
//   cTot->Divide(3,4);
//   // Draw the averaged energy spectra from monte carlo and data......
//   Double_t scaling = DrawEnergyPanel(hEFlipperOff  ,btr[1]->hEERef ,hEFlipperOff,btr[1]->hEERef,1,cTot);
//   DrawEnergyPanel(hEFlipperOff_I,btr[1]->hEERef1,hEFlipperOff_I,btr[1]->hEERef1,2,cTot,scaling);
//   DrawEnergyPanel(hEFlipperOff_2,btr[1]->hEERef2,hEFlipperOff_2,btr[1]->hEERef2,3,cTot,scaling);
// 
//   scaling = DrawEnergyPanel(hWFlipperOff  ,btr[1]->hEERef ,hWFlipperOff,btr[1]->hEWRef,4,cTot);
//   DrawEnergyPanel(hWFlipperOff_I,btr[1]->hEERef1,hWFlipperOff_I,btr[1]->hEERef1,5,cTot,scaling);
//   DrawEnergyPanel(hWFlipperOff_2,btr[1]->hEERef2,hWFlipperOff_2,btr[1]->hEERef2,6,cTot,scaling);
// 
//   scaling = DrawEnergyPanel(hEFlipperOn   ,btr[2]->hEERef ,hEFlipperOn ,btr[2]->hEERef,7,cTot);
//   DrawEnergyPanel(hEFlipperOn_I ,btr[2]->hEERef1,hEFlipperOn_I,btr[2]->hEERef1,8,cTot,scaling);
//   DrawEnergyPanel(hEFlipperOn_2 ,btr[2]->hEERef2,hEFlipperOn_2,btr[2]->hEERef2,9,cTot,scaling);
// 
//   scaling = DrawEnergyPanel(hWFlipperOn   ,btr[2]->hEERef ,hWFlipperOn ,btr[2]->hEERef,10,cTot); 
//   DrawEnergyPanel(hWFlipperOn_I ,btr[2]->hEERef1,hWFlipperOn_I ,btr[2]->hEERef1,11,cTot,scaling);
//   DrawEnergyPanel(hWFlipperOn_2 ,btr[2]->hEERef2,hWFlipperOn_2 ,btr[2]->hEERef2,12,cTot,scaling);
// 
//   // output the spectra
//   cTot->Print("output_files/TotalRates.pdf");

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
  Double_t x1[80],x1er[80],sig[80];
  
  for(Int_t i = 0 ; i < 80 ; i++){
    x1[i]   = hEFlipperOn->GetBinCenter(i+1);
    x1er[i] = 0.001;
    sig[i]  = 1.;
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

    btr[1]->Load_Histograms(bckr[btr[1]->Bkg_index],0);
   for(Int_t i = 0 ; i < 80 ; i++){
   	if(hEFlipperOn->GetBinError(i+1) > 0)
    		chiE0on[i] = Power((hEFlipperOn->GetBinContent(i+1) -
				  btr[1]->hEERef->GetBinContent(i+1))/
				  hEFlipperOn->GetBinError(i+1),1);			 
	if(hEFlipperOn_I->GetBinError(i+1) > 0)
   		chiE1on[i]  = Power((hEFlipperOn_I->GetBinContent(i+1) -
				btr[1]->hEERef1->GetBinContent(i+1))/
				hEFlipperOn_I->GetBinError(i+1),1);
   	if(hEFlipperOn_2->GetBinError(i+1) > 0)
		chiE2on[i] += Power((hEFlipperOn_2->GetBinContent(i+1) -
			btr[1]->hEERef2->GetBinContent(i+1))/
			hEFlipperOn_2->GetBinError(i+1),1);
   	if(hEFlipperOff->GetBinError(i+1) > 0)	
   		chiE0off[i] += Power((hEFlipperOff->GetBinContent(i+1) -
			btr[1]->hEERef->GetBinContent(i+1))/
			hEFlipperOff->GetBinError(i+1),1);
  	 if(hEFlipperOff_I->GetBinError(i+1) > 0)	
   		chiE1off[i] += Power((hEFlipperOff_I->GetBinContent(i+1) -
			btr[1]->hEERef1->GetBinContent(i+1))/
			hEFlipperOff_I->GetBinError(i+1),1);
  	if(hEFlipperOff_2->GetBinError(i+1) > 0)
	       chiE2off[i] += Power((hEFlipperOff_2->GetBinContent(i+1) -
			btr[1]->hEERef2->GetBinContent(i+1))/
			hEFlipperOff_2->GetBinError(i+1),1);
   	if(hWFlipperOn->GetBinError(i+1) > 0)			
   		chiW0on[i] += Power((hWFlipperOn->GetBinContent(i+1) -
			btr[1]->hEWRef->GetBinContent(i+1))/
			hWFlipperOn->GetBinError(i+1),1);
   	if(hWFlipperOn_I->GetBinError(i+1) > 0)			
   		chiW1on[i] += Power((hWFlipperOn_I->GetBinContent(i+1) -
			btr[1]->hEWRef1->GetBinContent(i+1))/
			hWFlipperOn_I->GetBinError(i+1),1); 
   	if(hWFlipperOn_2->GetBinError(i+1) > 0)
  		 chiW2on[i] += Power((hWFlipperOn_2->GetBinContent(i+1) -
			btr[1]->hEWRef2->GetBinContent(i+1))/
			hWFlipperOn_2->GetBinError(i+1),1);
   	if(hWFlipperOff->GetBinError(i+1) > 0)		
   		chiW0off[i] += Power((hWFlipperOff->GetBinContent(i+1) -
			btr[1]->hEWRef->GetBinContent(i+1))/
			hWFlipperOff->GetBinError(i+1),1);
  	if(hWFlipperOff_I->GetBinError(i+1) > 0)
		   chiW1off[i] += Power((hWFlipperOff_I->GetBinContent(i+1) -
			btr[1]->hEWRef1->GetBinContent(i+1))/
			hWFlipperOff_I->GetBinError(i+1),1);
	if(hWFlipperOff_2->GetBinError(i+1) > 0)
   		chiW2off[i] += Power((hWFlipperOff_2->GetBinContent(i+1) -
					btr[1]->hEWRef2->GetBinContent(i+1))/
				        hWFlipperOff_2->GetBinError(i+1),1);
   
  }
  
    btr[1]->Remove_Histograms(bckr[btr[1]->Bkg_index]);

  TGraphErrors *gEO0 = new TGraphErrors(80,x1,chiE0on,x1er,sig);
  ColorGraphic(gEO0,2,20,2,0.8,"East AFP On","Energy (keV)","#sigma");
  TGraphErrors *gEO1 = new TGraphErrors(80,x1,chiE1on,x1er,sig);
  ColorGraphic(gEO1,3,20,2,0.8);
  TGraphErrors *gEO2 = new TGraphErrors(80,x1,chiE2on,x1er,sig);
  ColorGraphic(gEO2,4,20,2,0.8);
  
  TGraphErrors *gEF0 = new TGraphErrors(80,x1,chiE0off,x1er,sig);
  ColorGraphic(gEF0,2,20,2,0.8,"East AFP Off","Energy (keV)","#sigma");
  TGraphErrors *gEF1 = new TGraphErrors(80,x1,chiE1off,x1er,sig);
  ColorGraphic(gEF1,3,20,2,0.8);
  TGraphErrors *gEF2 = new TGraphErrors(80,x1,chiE2off,x1er,sig);
  ColorGraphic(gEF2,4,20,2,0.8);
  
  TGraphErrors *gWO0 = new TGraphErrors(80,x1,chiW0on,x1er,sig);
  ColorGraphic(gWO0,2,20,2,0.8,"West AFP On","Energy (keV)","#sigma");
  TGraphErrors *gWO1 = new TGraphErrors(80,x1,chiW1on,x1er,sig);
  ColorGraphic(gWO1,3,20,2,0.8);
  TGraphErrors *gWO2 = new TGraphErrors(80,x1,chiW2on,x1er,sig);
  ColorGraphic(gWO2,4,20,2,0.8);
  
  TGraphErrors *gWF0 = new TGraphErrors(80,x1,chiW0off,x1er,sig);
  ColorGraphic(gWF0,2,20,2,0.8,"West AFP Off","Energy (keV)","#sigma");
  TGraphErrors *gWF1 = new TGraphErrors(80,x1,chiW1off,x1er,sig);
  ColorGraphic(gWF1,3,20,2,0.8);
  TGraphErrors *gWF2 = new TGraphErrors(80,x1,chiW2off,x1er,sig);
  ColorGraphic(gWF2,4,20,2,0.8);
   
  TCanvas *cETChi = new TCanvas("cETChi","Total Energy #chi^{2}");
  cETChi->Divide(2,2);
  
  cETChi->cd(1);
  gEO0->Draw("AP");
  gEO1->Draw("P");
  gEO2->Draw("P");
  
  cETChi->cd(2);
  gEF0->Draw("AP");
  gEF1->Draw("P");
  gEF2->Draw("P");
  
  cETChi->cd(3);

  gWO0->Draw("AP");
  gWO1->Draw("P");
  gWO2->Draw("P");
  
  cETChi->cd(4);

  gWF0->Draw("AP");
  gWF1->Draw("P");
  gWF2->Draw("P");
  
  cETChi->Print("output_files/EnergyChi.pdf");

  TCanvas *cFinal = new TCanvas("cFinal","Final");
  cFinal->Divide(2,1);
  cFinal->cd(1);
  hEFlipperOn->Draw();
  cFinal->cd(2);
  TH1F *hEFlipperClone = (TH1F*)hEFlipperOn->Clone("hEFlipperClone");
  hEFlipperClone->Draw();
  cFinal->Print("output_files/Final.pdf");

  //------------------------------------------------------
  // output the rough super ratio 
  
  Double_t R1p = hEFlipperOn_2->Integral(nlow,nhigh);
  Double_t R1m = hWFlipperOn_2->Integral(nlow,nhigh);
  Double_t R2p = hEFlipperOff_2->Integral(nlow,nhigh);
  Double_t R2m = hWFlipperOff_2->Integral(nlow,nhigh);
  
  Double_t SS = (R1p*R2m)/(R2p*R1m);  
  Double_t As_quick = (1.-sqrt(SS))/(1.+sqrt(SS));
  cout<<"Simple Super ratio for Type 2/3's is " << As_quick << endl;
  /*
  delete cFinal;
  delete cETChi;
  */
};
//------------------------------------------------------------------------------------------
void Define_E_Spec()
{
  Int_t nxbins = 200;

  hEFlipperOn = new TH1F("hEFlipperOn","East On Type 0; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOn = new TH1F("hWFlipperOn","West On Type 0; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hEFlipperOff = new TH1F("hEFlipperOff","East Off Type 0; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOff = new TH1F("hWFlipperOff","West Off Type 0; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
			 
  hEFlipperOn_I = new TH1F("hEFlipperOn_I","East On Type 1; Energy(keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hWFlipperOn_I = new TH1F("hWFlipperOn_I","West On Type 1; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hEFlipperOff_I = new TH1F("hEFlipperOff_I","East Off Type 1; Energy (keV) ;"
			    "Rate (Hz)",nxbins,0,2000);
  hWFlipperOff_I = new TH1F("hWFlipperOff_I","West Off Type 1; Energy (keV) ; Rate"
			 "(Hz)",nxbins,0,2000);
			 
  hEFlipperOn_2 = new TH1F("hEFlipperOn_2","East On Type 2/3; Energy(keV) ; Rate (Hz)"
		         ,nxbins,0,2000);
  hWFlipperOn_2 = new TH1F("hWFlipperOn_2","West On Type 2/3; Energy (keV) ; Rate (Hz)"
			 ,nxbins,0,2000);
  hEFlipperOff_2 = new TH1F("hEFlipperOff_2","East Off Type 2/3 ; Energy (keV) ;"
			    "Rate (Hz)",nxbins,0,2000);
  hWFlipperOff_2 = new TH1F("hWFlipperOff_2","West Off Type 2/3 ; Energy (keV) ; Rate"
			 "(Hz)",nxbins,0,2000);

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

}
void AddRates(Double_t *Cter,TH1F *hOut,TH1F *hIn,Int_t ibin)
{
  if(hIn->GetBinError(ibin)>0){
		Cter[ibin-1] += 1./TMath::Power(hIn->GetBinError(ibin),2);
		hOut->SetBinContent(ibin,hOut->GetBinContent(ibin) +
				   hIn->GetBinContent(ibin)/TMath::Power(hIn->GetBinError(ibin),2));  
  }else{
    Cter[ibin-1] += 0.;
  }
}
//---------------------------------------------------------------------------
void Fill_Foreground(Double_t *t1,Double_t *t2, Double_t *t3, Double_t *t4,Double_t *tl,Double_t *nmiss)
{
  Int_t nbinsX = 200;
  using namespace TMath;
  
  Double_t t1val = *t1;
  Double_t t2val = *t2;
  Double_t t3val = *t3;
  Double_t t4val = *t4;
  Double_t tloss = *tl;
  Double_t missc = *nmiss;

  Double_t nBinsX0=0;

  Double_t EastOn[nbinsX],EastOn_I[nbinsX],EastOn_2[nbinsX];
  Double_t EastOff[nbinsX],EastOff_I[nbinsX],EastOff_2[nbinsX];
  Double_t WestOn[nbinsX],WestOn_I[nbinsX],WestOn_2[nbinsX];
  Double_t WestOff[nbinsX],WestOff_I[nbinsX],WestOff_2[nbinsX];
  
  for(Int_t i = 0 ; i < nbinsX ; i++){
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
  missc = 0.;
  tloss = 0.;
  
  // Loop over all beta runs and sum up the background subtracted spectra
  fstream betalist;
  betalist.open("output_files/list_of_beta_run.txt",fstream::out);

  for(Int_t i = 0 ; i < nbeta ; i++){
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    betalist << i << "\t" << btr[i]->GetRunNumber();
    tloss += btr[i]->time_difference;
    missc += (btr[i]->heq->Integral(nlow,nhigh) + btr[i]->hwq->Integral(nlow,nhigh))*btr[i]->time_difference;
    betalist << "\t" << btr[i]->time_difference << "\t" << (btr[i]->rtime_w+btr[i]->rtime_e)/2. << "\t";
    betalist << (btr[i]->heq->Integral(nlow,nhigh) + btr[i]->hwq->Integral(nlow,nhigh)) << "\t";
    betalist << (btr[i]->heq->Integral(nlow,nhigh) + btr[i]->hwq->Integral(nlow,nhigh))* btr[i]->time_difference << "\t";
    betalist << btr[i]->heq->Integral(nlow,nhigh)*btr[i]->rtime_e + btr[i]->hwq->Integral(nlow,nhigh)*btr[i]->rtime_w << endl;
    
    for(Int_t bin = 1 ; bin <= hEFlipperOn->GetNbinsX(); bin++){
	if(btr[i]->flipperOn == 1){
	      AddRates(EastOn,hEFlipperOn,btr[i]->heq,bin);
              AddRates(WestOn,hWFlipperOn,btr[i]->hwq,bin);
	      AddRates(EastOn_I,hEFlipperOn_I,btr[i]->heqF,bin);
              AddRates(WestOn_I,hWFlipperOn_I,btr[i]->hwqF,bin);
	      AddRates(EastOn_2,hEFlipperOn_2,btr[i]->heqG,bin);
              AddRates(WestOn_2,hWFlipperOn_2,btr[i]->hwqG,bin);
	      // Incriment the run time
	      if(bin == 1){
		t2val += btr[i]->rtime_e;
		t1val += btr[i]->rtime_w;
	      }
	} else if(btr[i]->flipperOn == 0){ 
	      AddRates(EastOff,hEFlipperOff,btr[i]->heq,bin);
              AddRates(WestOff,hWFlipperOff,btr[i]->hwq,bin);
	      AddRates(EastOff_I,hEFlipperOff_I,btr[i]->heqF,bin);
              AddRates(WestOff_I,hWFlipperOff_I,btr[i]->hwqF,bin);
	      AddRates(EastOff_2,hEFlipperOff_2,btr[i]->heqG,bin);
              AddRates(WestOff_2,hWFlipperOff_2,btr[i]->hwqG,bin); 
	// Incriment the run time
	      if(bin == 1){
		  t4val += btr[i]->rtime_e;
		  t3val += btr[i]->rtime_w;
	      }
	}
    }

	if(i==0) nBinsX0=btr[0]->heq->GetNbinsX();

    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }

  betalist.close();

  for(Int_t i = 1 ; i <= nBinsX0 ; i++){
    
    hEFlipperOn->SetBinContent(i,   hEFlipperOn->GetBinContent(i)   / EastOn[i-1]);
    hEFlipperOn_I->SetBinContent(i, hEFlipperOn_I->GetBinContent(i) / EastOn_I[i-1]);
    hEFlipperOn_2->SetBinContent(i, hEFlipperOn_2->GetBinContent(i) / EastOn_2[i-1]);
    
    hEFlipperOff->SetBinContent(i,   hEFlipperOff->GetBinContent(i)   / EastOff[i-1]);
    hEFlipperOff_I->SetBinContent(i, hEFlipperOff_I->GetBinContent(i) / EastOff_I[i-1]);
    hEFlipperOff_2->SetBinContent(i, hEFlipperOff_2->GetBinContent(i) / EastOff_2[i-1]);
    
    hWFlipperOn->SetBinContent(i,   hWFlipperOn->GetBinContent(i)   / WestOn[i-1]);
    hWFlipperOn_I->SetBinContent(i, hWFlipperOn_I->GetBinContent(i) / WestOn_I[i-1]);
    hWFlipperOn_2->SetBinContent(i, hWFlipperOn_2->GetBinContent(i) / WestOn_2[i-1]);
    
    hWFlipperOff->SetBinContent(i,   hWFlipperOff->GetBinContent(i)   / WestOff[i-1]);
    hWFlipperOff_I->SetBinContent(i, hWFlipperOff_I->GetBinContent(i) / WestOff_I[i-1]);
    hWFlipperOff_2->SetBinContent(i, hWFlipperOff_2->GetBinContent(i) / WestOff_2[i-1]);
    
    hEFlipperOn->SetBinError(i, 1./Sqrt(EastOn[i-1]));
    hEFlipperOn_I->SetBinError(i, 1./Sqrt(EastOn_I[i-1]));
    hEFlipperOn_2->SetBinError(i, 1./Sqrt(EastOn_2[i-1]));
    hEFlipperOff->SetBinError(i, 1./Sqrt(EastOff[i-1]));
    hEFlipperOff_I->SetBinError(i, 1./Sqrt(EastOff_I[i-1]));
    hEFlipperOff_2->SetBinError(i, 1./Sqrt(EastOff_2[i-1]));
    
    hWFlipperOn->SetBinError(i, 1./Sqrt(WestOn[i-1]));
    hWFlipperOn_I->SetBinError(i, 1./Sqrt(WestOn_I[i-1]));
    hWFlipperOn_2->SetBinError(i, 1./Sqrt(WestOn_2[i-1]));
    hWFlipperOff->SetBinError(i, 1./Sqrt(WestOff[i-1]));
    hWFlipperOff_I->SetBinError(i, 1./Sqrt(WestOff_I[i-1]));
    hWFlipperOff_2->SetBinError(i, 1./Sqrt(WestOff_2[i-1])); 
    
  }
  
  *t1 = t1val;
  *t2 = t2val;
  *t3 = t3val;
  *t4 = t4val;
  *tl = tloss;
  *nmiss = missc;
 
}

//---------------------------------------------------------------------------
void Fill_Background(Double_t* t1,Double_t* t2, Double_t* t3, Double_t* t4)
{

  using namespace TMath;
  
  Int_t nbinsX = 200;
  
  Double_t t1val = *t1;
  Double_t t2val = *t2;
  Double_t t3val = *t3;
  Double_t t4val = *t4;

    btr[0]->Load_Histograms(bckr[btr[0]->Bkg_index],0);
  Double_t nBinsX0=btr[0]->heq->GetNbinsX();
    btr[0]->Remove_Histograms(bckr[btr[0]->Bkg_index]);

  Double_t EastOn[nbinsX],EastOn_I[nbinsX],EastOn_2[nbinsX];
  Double_t EastOff[nbinsX],EastOff_I[nbinsX],EastOff_2[nbinsX];
  Double_t WestOn[nbinsX],WestOn_I[nbinsX],WestOn_2[nbinsX];
  Double_t WestOff[nbinsX],WestOff_I[nbinsX],WestOff_2[nbinsX];
  
  for(Int_t i = 0; i < nbinsX ; i++){
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
    bckr[i]->Load_Histograms();
    if(bckr[i]->GetRunNumber() != last){
    for(Int_t bin = 1 ; bin <= hEFlipperOn_B->GetNbinsX(); bin++){
      if(bckr[i]->flipperOn == 1){
        AddRates(EastOn,hEFlipperOn_B,bckr[i]->heq,bin);
        AddRates(WestOn,hWFlipperOn_B,bckr[i]->hwq,bin);
	AddRates(EastOn_I,hEFlipperOn_B_I,bckr[i]->heqF,bin);
        AddRates(WestOn_I,hWFlipperOn_B_I,bckr[i]->hwqF,bin);
	AddRates(EastOn_2,hEFlipperOn_B_2,bckr[i]->heqG,bin);
        AddRates(WestOn_2,hWFlipperOn_B_2,bckr[i]->hwqG,bin);
	// Incriment the run time
	if(bin == 1){
	  t2val += bckr[i]->rtime_e;
	  t1val += bckr[i]->rtime_w;
	}
      } else {
	AddRates(EastOff,hEFlipperOff_B,bckr[i]->heq,bin);
        AddRates(WestOff,hWFlipperOff_B,bckr[i]->hwq,bin);
	AddRates(EastOff_I,hEFlipperOff_B_I,bckr[i]->heqF,bin);
        AddRates(WestOff_I,hWFlipperOff_B_I,bckr[i]->hwqF,bin);
	AddRates(EastOff_2,hEFlipperOff_B_2,bckr[i]->heqG,bin);
        AddRates(WestOff_2,hWFlipperOff_B_2,bckr[i]->hwqG,bin);
	// Incriment the run time
	if(bin == 1) {
	  t4val += bckr[i]->rtime_e;
	  t3val += bckr[i]->rtime_w;
 	}
       }
       last = bckr[i]->GetRunNumber();
      }
    }
    bckr[i]->Remove_Histograms();
  }
 
   for(Int_t i = 1 ; i <= nBinsX0 ; i++){
   // cout << i << "\t"<< EastOn[i-1] << "\t" << EastOff[i-1] << endl;
    hEFlipperOn_B->SetBinContent(i  , hEFlipperOn_B->GetBinContent(i)   / EastOn[i-1]);
    hEFlipperOn_B_I->SetBinContent(i, hEFlipperOn_B_I->GetBinContent(i) / EastOn_I[i-1]);
    hEFlipperOn_B_2->SetBinContent(i, hEFlipperOn_B_2->GetBinContent(i) / EastOn_2[i-1]);
   
    hEFlipperOff_B->SetBinContent(i  , hEFlipperOff_B->GetBinContent(i)   / EastOff[i-1]);
    hEFlipperOff_B_I->SetBinContent(i, hEFlipperOff_B_I->GetBinContent(i) / EastOff_I[i-1]);
    hEFlipperOff_B_2->SetBinContent(i, hEFlipperOff_B_2->GetBinContent(i) / EastOff_2[i-1]);
   
    hWFlipperOn_B->SetBinContent(i  , hWFlipperOn_B->GetBinContent(i)   / WestOn[i-1]);
    hWFlipperOn_B_I->SetBinContent(i, hWFlipperOn_B_I->GetBinContent(i) / WestOn_I[i-1]);
    hWFlipperOn_B_2->SetBinContent(i, hWFlipperOn_B_2->GetBinContent(i) / WestOn_2[i-1]);
   
    hWFlipperOff_B->SetBinContent(i,   hWFlipperOff_B->GetBinContent(i)   / WestOff[i-1]);
    hWFlipperOff_B_I->SetBinContent(i, hWFlipperOff_B_I->GetBinContent(i) / WestOff_I[i-1]);
    hWFlipperOff_B_2->SetBinContent(i, hWFlipperOff_B_2->GetBinContent(i) / WestOff_2[i-1]);
   
    hEFlipperOn_B->SetBinError(i, 1./Sqrt(EastOn[i-1]));
    hEFlipperOn_B_I->SetBinError(i, 1./Sqrt(EastOn_I[i-1]));
    hEFlipperOn_B_2->SetBinError(i, 1./Sqrt(EastOn_2[i-1]));
    hEFlipperOff_B->SetBinError(i, 1./Sqrt(EastOff[i-1]));
    hEFlipperOff_B_I->SetBinError(i, 1./Sqrt(EastOff_I[i-1]));
    hEFlipperOff_B_2->SetBinError(i, 1./Sqrt(EastOff_2[i-1]));
   
    hWFlipperOn_B->SetBinError(i, 1./Sqrt(WestOn[i-1]));
    hWFlipperOn_B_I->SetBinError(i, 1./Sqrt(WestOn_I[i-1]));
    hWFlipperOn_B_2->SetBinError(i, 1./Sqrt(WestOn_2[i-1]));
    hWFlipperOff_B->SetBinError(i, 1./Sqrt(WestOff[i-1]));
    hWFlipperOff_B_I->SetBinError(i, 1./Sqrt(WestOff_I[i-1]));
    hWFlipperOff_B_2->SetBinError(i, 1./Sqrt(WestOff_2[i-1])); 
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
//----------------------------------------------------------------------------------
void GammaBack()
{ 
  using namespace TMath;

  const Int_t nxbins = 200;
  //===================================================================================
  vector<Double_t> aveEO (nxbins)  ,aveWO (nxbins)  ,aveE1 (nxbins)  ,aveW1 (nxbins);
  vector<Double_t> aveEOe (nxbins) ,aveWOe (nxbins) ,aveE1e (nxbins) ,aveW1e (nxbins);
  vector<Double_t> avesEO (nxbins) ,avesWO (nxbins) ,avesE1 (nxbins) ,avesW1(nxbins);
  vector<Double_t> avesEOe (nxbins),avesWOe (nxbins),avesE1e (nxbins),avesW1e (nxbins);
  //===================================================================================
  vector<Double_t> x,REON,REOFF,RWON,RWOFF,REONe,REOFFe,RWONe,RWOFFe,xoff;
  //===================================================================================
  Double_t reonr=0.,reonre=0.,reofr=0.,reofre=0.,rwonr=0.,rwonre=0.,rwofr=0.,rwofre=0.;
  
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
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
      if(btr[i]->flipperOn == 1){
	  Average_All_Hists(btr[i]->hENoMWPC,aveEO,aveEOe);
          Average_All_Hists(btr[i]->hEMWPC,avesEO,avesEOe);
          Average_All_Hists(btr[i]->hWNoMWPC,aveWO,aveWOe);
          Average_All_Hists(btr[i]->hWMWPC,avesWO,avesWOe);
      } else if(btr[i]->flipperOn == 0){
          Average_All_Hists(btr[i]->hENoMWPC,aveE1,aveE1e);
          Average_All_Hists(btr[i]->hEMWPC,avesE1,avesE1e);
          Average_All_Hists(btr[i]->hWNoMWPC,aveW1,aveW1e);
          Average_All_Hists(btr[i]->hWMWPC,avesW1,avesW1e);
      }
  
    //-----------------------------------------------------------------
    // Look at individual runs
    //----------------------------------------------------------------
   
    if(btr[i]->flipperOn == 1 && btr[i]->MWPC_RatioE_e !=0 && btr[i]->MWPC_RatioW_e !=0){
      x.push_back(btr[i]->GetRunNumber());
      
      reonr +=  btr[i]->MWPC_RatioE / Power(btr[i]->MWPC_RatioE_e,2);
      rwonr +=  btr[i]->MWPC_RatioW / Power(btr[i]->MWPC_RatioW_e,2);

      reonre +=  1. / Power(btr[i]->MWPC_RatioE_e,2);
      rwonre +=  1. / Power(btr[i]->MWPC_RatioW_e,2);

      REON.push_back(btr[i]->MWPC_RatioE_fit);
      RWON.push_back(btr[i]->MWPC_RatioW_fit);

      REONe.push_back(2.2*btr[i]->MWPC_RatioE_fit_e);
      RWONe.push_back(2.2*btr[i]->MWPC_RatioW_fit_e);

    } else if(btr[i]->flipperOn == 0 && btr[i]->MWPC_RatioE_e !=0 && btr[i]->MWPC_RatioW_e !=0){

      xoff.push_back(btr[i]->GetRunNumber());

      reofr +=  btr[i]->MWPC_RatioE / Power(btr[i]->MWPC_RatioE_e,2);
      rwofr +=  btr[i]->MWPC_RatioW / Power(btr[i]->MWPC_RatioW_e,2);

      reofre +=  1. / Power(btr[i]->MWPC_RatioE_e,2);
      rwofre +=  1. / Power(btr[i]->MWPC_RatioW_e,2);
      
      REOFF.push_back(btr[i]->MWPC_RatioE_fit);
      RWOFF.push_back(btr[i]->MWPC_RatioW_fit);
      
      REOFFe.push_back(2.2*btr[i]->MWPC_RatioE_fit_e);
      RWOFFe.push_back(2.2*btr[i]->MWPC_RatioW_fit_e);

    }
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  
  if((int)REON.size() == 0) return;
  
  TGraphErrors *gRen = new TGraphErrors((int)x.size(),&x[0],&REON[0],0,&REONe[0]);
  ColorGraphic(gRen,2,20,1);
  
  TGraphErrors *gRef = new TGraphErrors((int)xoff.size(),&xoff[0],&REOFF[0],0,&REOFFe[0]);
  ColorGraphic(gRef,4,20,1);
   
  TGraphErrors *gRwn = new TGraphErrors((int)x.size(),&x[0],&RWON[0],0,&RWONe[0]);
  ColorGraphic(gRwn,2,20,1);
  
  TGraphErrors *gRwf = new TGraphErrors((int)xoff.size(),&xoff[0],&RWOFF[0],0,&RWOFFe[0]);
  ColorGraphic(gRwf,4,20,1);
  
  Double_t interEon = 0.;
  Double_t interEof = 0.;
  Double_t interWon = 0.;
  Double_t interWof = 0.;
  
  reonre = 1./sqrt(reonre);
  rwonre = 1./sqrt(rwonre);
  reofre = 1./sqrt(reofre);
  rwofre = 1./sqrt(rwofre);

  reonr  *= reonre*reonre;
  reofr  *= reofre*reofre;
  rwonr  *= rwonre*rwonre;
  rwofr  *= rwofre*rwofre;

  cout << "rwofre " << rwofre << "  rwonr "  << rwonr << endl;
  cout << "rwonre " << rwonre << "  rwofr "  << rwofr << endl;

  Return_Asymmetry(aveEO,aveEOe,hEFlipperOn_G->GetNbinsX());
  Return_Asymmetry(aveE1,aveE1e,hEFlipperOn_G->GetNbinsX());
  Return_Asymmetry(aveWO,aveWOe,hWFlipperOn_G->GetNbinsX());
  Return_Asymmetry(aveW1,aveW1e,hWFlipperOn_G->GetNbinsX());

  for(Int_t j = 1 ; j <=  hEFlipperOn_G->GetNbinsX(); j++){

    if( j >=nlow && j <=nhigh){
      interEon += aveEOe[j-1]*aveEOe[j-1];
      interEof += aveE1e[j-1]*aveE1e[j-1];
      interWon += aveWOe[j-1]*aveWOe[j-1];
      interWof += aveW1e[j-1]*aveW1e[j-1];
    }

    //cout << j<< "\t" << aveEO[j-1] << "\t" << aveEOe[j-1] << endl;  

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
    /*
    cout << "East On " << j << " " << hEFlipperOn_G->GetBinContent(j-1) << "  " << hEFlipperOn_G->GetBinError(j-1) << endl;
    cout << "East Off" << j << " " << hEFlipperOff_G->GetBinContent(j-1) << "  " << hEFlipperOff_G->GetBinError(j-1) << endl;
    cout << "West On " << j << " " << hWFlipperOn_G->GetBinContent(j-1) << "  " << hWFlipperOn_G->GetBinError(j-1) << endl;  
    cout << "West Off "<< j << " " << hWFlipperOff_G->GetBinContent(j-1) << "  " << hWFlipperOff_G->GetBinError(j-1) << endl;
    */
  }

  interEon = sqrt(Abs(interEon));
  interEof = sqrt(Abs(interEof));
  interWon = sqrt(Abs(interWon));
  interWof = sqrt(Abs(interWof));

  TCanvas *cNoMWPC = new TCanvas("cNoMWPC","No wire chamber");
  cNoMWPC->Divide(2,2);
  cNoMWPC->cd(1);

  TLegend *lEst;
  DrawMWPCPanel1(hEFlipperOn_G,hEFlipperOn_S,hEFlipperOff_G,hEFlipperOff_S,lEst,interEof,reofre,interEon,reonre);

  cNoMWPC->cd(2);

  TF1 *fEnl = new TF1("fEnl","[0]",runstart,runstop);
  TF1 *fEfl = new TF1("fEfl","[0]",runstart,runstop);
  
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
			  fEnl->GetParameter(0),fEnl->GetParError(0)),"lp");
  lEn->AddEntry(gRen,Form("#chi^{2} / #nu = %5.3f",fEnl->GetChisquare()/fEnl->GetNDF()),"lp");
  lEn->AddEntry(gRef,Form("(Off) Average Ratio = %5.3f #pm %5.3f",
			  fEfl->GetParameter(0),fEfl->GetParError(0)),"lp");
  lEn->AddEntry(gRef,Form("#chi^{2} / #nu = %5.3f",fEfl->GetChisquare()/fEfl->GetNDF()),"lp");
  lEn->SetFillColor(0);
  lEn->Draw();

  cNoMWPC->cd(3);

  TLegend *lWst;
  DrawMWPCPanel1(hWFlipperOn_G,hWFlipperOn_S,hWFlipperOff_G,hWFlipperOff_S,lWst,interWof,rwofre,interWon,rwonre);

  cNoMWPC->cd(4);

  TF1 *fWnl = new TF1("fWnl","[0]",runstart,runstop);
  TF1 *fWfl = new TF1("fWfl","[0]",runstart,runstop);

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
			  fWnl->GetParameter(0),fWnl->GetParError(0)),"lp");
  lWn->AddEntry(gRwn,Form("#chi^{2} / #nu = %5.3f",fWnl->GetChisquare()/fWnl->GetNDF()),"lp");
  lWn->AddEntry(gRwf,Form("(Off) Average Ratio = %5.3f #pm %5.3f",
			  fWfl->GetParameter(0),fWfl->GetParError(0)),"lp");
  lWn->AddEntry(gRwf,Form("#chi^{2} / #nu = %5.3f",fWfl->GetChisquare()/fWfl->GetNDF()),"lp");
  lWn->SetFillColor(0);
  lWn->Draw();
  
  cNoMWPC->Print("output_files/MWPC_gamma_check.pdf");
  /*
  TCanvas *cMWPC_Gamma = new TCanvas("cMWPC_Gamma","Gammas");
  cMWPC_Gamma->Divide(2,1);
  cMWPC_Gamma->cd(1);
  btr[0]->hWMWPC->SetTitle(Form("Run : %d",btr[7]->GetRunNumber()));
  btr[0]->hWNoMWPC->SetLineColor(2);
  btr[0]->hWMWPC->SetLineColor(4);
  btr[0]->hWMWPC->Draw();
  bckr[btr[0]->Bkg_index]->hWMWPC->Draw("same");
  bckr[btr[0]->Bkg_index]->hWNoMWPC->SetLineColor(3);
  bckr[btr[0]->Bkg_index]->hWNoMWPC->Draw("same");
  btr[0]->hWNoMWPC->Draw("same");

  cMWPC_Gamma->cd(2);
  btr[0]->hWMWPC->SetTitle(Form("Run : %d",btr[6]->GetRunNumber()));
  btr[0]->hWNoMWPC->SetLineColor(2);
  btr[0]->hWMWPC->SetLineColor(4);
  btr[0]->hWMWPC->Draw();
  bckr[btr[0]->Bkg_index]->hWMWPC->Draw("same");
  bckr[btr[0]->Bkg_index]->hWNoMWPC->SetLineColor(3);
  bckr[btr[0]->Bkg_index]->hWNoMWPC->Draw("same");
  btr[0]->hWNoMWPC->Draw("same");

  cMWPC_Gamma->Print("output_files/randomgammaback.pdf");
  
  delete mgRe;
  delete lEn;
  delete lWst;
  delete fWnl;
  delete fWfl;
  delete mgRw;
  delete lWn;
  delete cNoMWPC;
  delete cMWPC_Gamma;
  */
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
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    for(Int_t j = 1 ; j <= hTDCDiffTot_Off->GetNbinsX() ; j++){
      if(btr[i]->flipperOn == 0){
	hTDCDiffTot_Off->SetBinContent(j,hTDCDiffTot_Off->GetBinContent(j)+btr[i]->hTDCDiff->GetBinContent(j));
      } else {
	hTDCDiffTot_On->SetBinContent(j,hTDCDiffTot_On->GetBinContent(j)+ btr[i]->hTDCDiff->GetBinContent(j));
      }
    }
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  
  TCanvas *cTD = new TCanvas("cTD","TDC Difference");
  cTD->Divide(1,2);
  cTD->cd(1);
  hTDCDiffTot_On->Draw();
  cTD->cd(2);
  hTDCDiffTot_Off->Draw();
  
  cTD->Print("output_files/tdcdifference.pdf");
  /*
  delete hTDCDiffTot_On;
  delete hTDCDiffTot_Off;
  delete cTD;
  */
  return;
}
//------------------------------------------------------------------------------------------------
void Collect_Gammas()
{

  Double_t x[MAXRUNS],fge[MAXRUNS],fgw[MAXRUNS],fr[MAXRUNS];
  
  for(Int_t i = 0 ; i < nbeta ; i++){
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    x[i]   = btr[i]->GetRunNumber();
    fge[i] = btr[i]->hGammaCounts->GetBinContent(1) / btr[i]->hGammaCountsg->GetBinContent(1);
    fgw[i] = btr[i]->hGammaCounts->GetBinContent(2) / btr[i]->hGammaCountsg->GetBinContent(2);
    if(fgw[i] !=0)fr[i]  = fge[i] / fgw[i];
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
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

  cGamRat->Print("output_files/GammaRatio.pdf");
  /*
  delete gGaRat;
  delete gGaRatw;
  delete gRat;
  delete cGamRat;*/
}
//--------------------------------------------------------------------------------------------------
void Plot_ChiDis()
{
 
  TF1 *fbeta2 = new TF1("fbeta2",PoverE,0,2000,1);  
  fbeta2->SetParameter(0,1.);
  TF1 *fline  = new TF1("fline","[0]",0,1000);
  TCanvas *cTest = new TCanvas("cTest","Chi Square Test");
  cTest->cd();
  //Double_t chiA = 0.;

  fstream fch;
  
  if(btr[0]->GetGeo()!=2){
    fch.open(Form("output_files/chi_squares_%d.txt",btr[0]->GetGeo()),fstream::out);
  } else {
    if(btr[0]->GetRunNumber() <12000){
      fch.open("output_files/chi_squares_2.txt",fstream::out);
    } else {
      fch.open("output_files/chi_squares_3.txt",fstream::out);
    }
  }
  
  TH1F *hChiOct = new TH1F("hChiOct","#chi^{2} Distribution",40,0.,40.);
  TH1F *hTest   = new TH1F("hTest","testing",octet[0]->hAsyTot[2]->GetNbinsX()
			      		    ,octet[0]->hAsyTot[2]->GetBinCenter(1) - octet[0]->hAsyTot[2]->GetBinWidth(1)/2.
					    ,octet[0]->hAsyTot[2]->GetBinCenter(octet[0]->hAsyTot[2]->GetNbinsX()) + 
 			                     octet[0]->hAsyTot[2]->GetBinWidth(1)/2.);
  
  for(Int_t i = 0 ; i < noct ; i++){
    for(Int_t j = 1 ; j <= octet[i]->hAsyTot[2]->GetNbinsX() ; j++){
      hTest->SetBinContent(j,octet[i]->hAsyTot[2]->GetBinContent(j) / fbeta2->Eval(octet[i]->hAsyTot[2]->GetBinCenter(j)));
      hTest->SetBinError(j,octet[i]->hAsyTot[2]->GetBinError(j) / fbeta2->Eval(octet[i]->hAsyTot[2]->GetBinCenter(j)));
    }   hTest->Fit(fline,"REMQ","",octet[i]->hAsyTot[2]->GetBinCenter(nlow),octet[i]->hAsyTot[2]->GetBinCenter(nhigh));
    hChiOct->Fill(fline->GetChisquare());
    fch << i << "\t" << fline->GetChisquare() << endl;
  }
  fch.close();

  hChiOct->Draw();
  TF1 *fchi = new TF1("fchi",ChiSquareDis,0,40.,2);
  fchi->SetParameter(0,hChiOct->Integral());
  fchi->SetParameter(1,(nhigh-nlow-1));
  fchi->SetLineColor(2);
  fchi->Draw("same");

  cTest->Print("output_files/ChiTesting.pdf");
  /*
  delete fchi;
  delete hTest;
  delete hChiOct;
  delete fline;
  delete fbeta2;
  delete cTest; */   
}

void PrintPDF()
{

}

void DefineRotationCollectionHistos(Int_t rbins,Float_t rad)
{

  hTotRote = new TH2F("hTotRote",
		   "East Type 1 Backscattering Rotation ; X_{west} (mm) ; Y_{west} (mm)"
		      ,rbins,-rad,rad,rbins,-rad,rad);
  hTotRoteI = new TH2F("hTotRoteI",
		      "East Type Trigger Side Position ; X_{east} (mm) ; Y_{east} (mm)",
		       rbins,-rad,rad,rbins,-rad,rad);
  
  hTotRotwI = new TH2F("hTotRotwI",
		      "West Type Trigger Side Position ; X_{west} (mm) ; Y_{west} (mm)",
		       rbins,-rad,rad,rbins,-rad,rad);
  hTotRotw = new TH2F("hTotRotw",
		      "West Type 1 Backscattering Rotation ; X_{east} (mm) ; Y_{east} (mm)"
		      ,rbins,-rad,rad,rbins,-rad,rad);
 
  hTotRote23 = new TH2F("hTotRote23",
			"East Type 23 Backscattering Rotation ; X_{west} (mm) ; Y_{west} (mm)"
			,rbins,-rad,rad,rbins,-rad,rad);
  
  hTotRotw23 = new TH2F("hTotRotw23",
			"West Type 23 Backscattering Rotation ; X_{east} (mm) ; Y_{east} (mm)"
			,rbins,-rad,rad,rbins,-rad,rad);

  hTotRoteI23 = new TH2F("hTotRoteI23",
			 "East Type 23 trigger side Position ; X_{east} (mm) ; Y_{east} (mm)"
			 ,rbins,-rad,rad,rbins,-rad,rad);
  
  hTotRotwI23 = new TH2F("hTotRotwI23",
			 "West Type 23 trigger side position ; X_{west} (mm) ; Y_{west} (mm)"
			 ,rbins,-rad,rad,rbins,-rad,rad); 

  return;
}
//---------------------------------------------------------------------------------------------------------------
bool SetAnalysisType(Int_t *remake){
  // Determine the type of analysis to perform
  cout << "Single Run or from list (1 of single, !=1 for list)" << endl;
 //cin >> runopt;
  runopt = 2;

  // Get the list of Runs to analyze
  if(runopt == 1)
   GetSingleRun();
  else {
    do{
      cout << "Geometry A fuck   : 7659  - 9189  " << endl;
      cout << "Geometry B   : 9333  - 10333 " << endl;
      cout << "Geometry C   : 10404 - 11095 "<< endl;
      cout << "Geometry C09 : 12430 - 12700 " << endl;
      cout << "2010 Data    : 13642 - 15954 " << endl;
      cout << "Adding List of Runs"  << endl;
      cout << "Enter first run : ";
     // cin >> runstart;
      cout << "Enter last run : ";
      //cin >> runstop;
      runstart = NRUNFIRST;
      runstop  = NRUNLAST;
    }while(runstart > runstop);
  
    cout << "Remake Histograms ? ( 1 = Yes , 0 = No)" << endl;
    cin >> *remake;

    if(!(GetListofRuns(runstart,runstop)))return kFALSE;
  }
  return kTRUE;
}
//--------------------------------------------------------------------------
void Load23SepArray(){
  
  fstream ref;
  Int_t j1,j2,j3,j4;
  
  // this get the file for the type23 separation fits from Monte Carlo....
  
  if(runstart >= 7659 && runstop <= 9189){
    nlist = 1;
    ref.open("input_files/type_23_a.txt",fstream::in);
  }
  else if(runstart >= 9333 && runstop <= 10333){
    nlist = 2;
    ref.open("type_23_b.txt",fstream::in);
  } 
  else if(runstart >= 10400 && runstop <= 11095){
    nlist = 3;
    ref.open("type_23_c.txt",fstream::in);
  }
  else if(runstart >= 12420 && runstop <= 12701){
    ref.open("input_files/type_23_c.txt",fstream::in);
    nlist = 4;
  } else if(runstart > 13000){
    ref.open("input_files/type_23_c.txt",fstream::in);
    nlist = 5;
  }

  for(Int_t i = 0 ; i < 100 ; i++){
    ref >> j1 >> sep23[i] >> j2 >> j3 >> j4;
  }
  ref.close();

}
//-----------------------------------------------------------------------------------------------
void Plot_Multi()
{
  return;
  /*--------------------------------------------------------------------------------------------
    Plot the MWPC Multiplicities for the East and West detectors for each event type.
  */

  TH1F *hEType0_MultiT  = new TH1F("hEType0_MultiT","Type 0; Multiplicity;Cnts",20,0,20);
  TH1F *hEType1_MultiT  = new TH1F("hEType1_MultiT","Type 1",20,0,20);
  TH1F *hEType23_MultiT = new TH1F("hEType23_MultiT","Type 23",20,0,20);

  TH1F *hWType0_MultiT  = new TH1F("hWType0_MultiT","Type 0; Multiplicity;Cnts",20,0,20);
  TH1F *hWType1_MultiT  = new TH1F("hWType1_MultiT","Type 1",20,0,20);
  TH1F *hWType23_MultiT = new TH1F("hWType23_MultiT","Type 23",20,0,20);

  for(Int_t i = 0 ; i < nbeta ; i++){
    for(Int_t j = 1 ; j < hEType0_MultiT->GetNbinsX() ; j++){

        hEType0_MultiT->SetBinContent(j,  hEType0_MultiT->GetBinContent(j) + btr[i]->hEType0_Multi->GetBinContent(j));
        hEType1_MultiT->SetBinContent(j,  hEType1_MultiT->GetBinContent(j) + btr[i]->hEType1_Multi->GetBinContent(j));
        hEType23_MultiT->SetBinContent(j, hEType23_MultiT->GetBinContent(j)+ btr[i]->hEType23_Multi->GetBinContent(j));

        hWType0_MultiT->SetBinContent(j,  hWType0_MultiT->GetBinContent(j) + btr[i]->hWType0_Multi->GetBinContent(j));
        hWType1_MultiT->SetBinContent(j,  hWType1_MultiT->GetBinContent(j) + btr[i]->hWType1_Multi->GetBinContent(j));
        hWType23_MultiT->SetBinContent(j, hWType23_MultiT->GetBinContent(j)+ btr[i]->hWType23_Multi->GetBinContent(j));
    }
  }
  
  TCanvas *cPlotMulti = new TCanvas("cPlotMulti","Multiplicity");
  cPlotMulti->Divide(2,1);
  cPlotMulti->cd(1);

  hEType0_MultiT->DrawNormalized();
  ColorGraphic(hEType1_MultiT,2,20,1);
  hEType1_MultiT->DrawNormalized("same");
  ColorGraphic(hEType23_MultiT,4,20,1);
  hEType23_MultiT->DrawNormalized("same");
  
 
  cPlotMulti->cd(2);

  hWType0_MultiT->DrawNormalized();
  ColorGraphic(hWType1_MultiT,2,20,1);
  hWType1_MultiT->DrawNormalized("same");
  ColorGraphic(hWType23_MultiT,4,20,1);
  hWType23_MultiT->DrawNormalized("same");
 
  cPlotMulti->Print("output_files/multiplots.pdf");
  // delete the histograms......................................
  /*
  delete hEType0_MultiT;
  delete hEType1_MultiT;
  delete hEType23_MultiT;
  delete hWType0_MultiT;
  delete hWType1_MultiT;
  delete hWType23_MultiT;
  delete cPlotMulti;*/

}
//--------------------------------------------------------------------------------------
Bool_t checkruns(Int_t currentrun, vector<Int_t> openruns,Int_t nruns)
{
  //simple check to make sure a file isn't already open
  if(openruns.size() == 0)return kTRUE;

  for(Int_t ii = 0 ; ii < (Int_t)openruns.size() ; ii++)
	if(currentrun == openruns[ii])return kFALSE;
      
  return kTRUE;
}
//--------------------------------------------------------------------------------------------
void DrawBackFractionPad(TCanvas *cCan,Int_t nPad,TGraphErrors *gEast,TGraphErrors *gWest,TF1 *fEast,
                         TF1 *fWest,TLegend *legBack,Double_t Ymax)
{

  cCan->cd(nPad);

  gEast->Draw("AP");
  gEast->Fit(fEast->GetName(),"REMQ+");
  gWest->Draw("P");
  gWest->Fit(fWest->GetName(),"REMQ+");

  gEast->GetYaxis()->SetRangeUser(0,Ymax);
  
  legBack->AddEntry(gEast,Form("East %5.3f #pm %5.3f",
                         fEast->GetParameter(0),fEast->GetParError(0)),"lp");
  legBack->AddEntry(gWest,Form("West %5.3f #pm %5.3f",
                         fWest->GetParameter(0),fWest->GetParError(0)),"lp");
  legBack->SetFillColor(kWhite);
  legBack->Draw();

};
//----------------------------------------------------------------------------------------------
void DrawRadialCounts(Double_t *ARader,Double_t Aave)
{
  // Sum total Radial Counts
  vector<Double_t> AllRadCnts (Rbins);
  vector<Double_t> RadErrCnt (Rbins);
  vector<Double_t> RadErrUnSc (Rbins);
  vector<Double_t> Xrad;
  
  for(Int_t ii = 0 ; ii < noct ; ii++){
     for(Int_t nrad = 0 ; nrad < Rbins ; nrad++){
        Xrad.push_back(nrad*5.+2.5);
        AllRadCnts[nrad] += octet[ii]->TotRadCounts[nrad];
     }
  }

  for(Int_t nrad = 0 ; nrad < Rbins ; nrad++){
        for(Int_t jrad = 0 ; jrad <= nrad ; jrad++){
          RadErrCnt[nrad] += 1./(ARader[jrad]*ARader[jrad]);
        }
        RadErrUnSc[nrad] = 100.0/sqrt(RadErrCnt[nrad])/Aave;
        RadErrCnt[nrad]  = 100.e6/sqrt(RadErrCnt[nrad])/Aave;
  }
  

  TGraph *gRadCnt     = new TGraph(Rbins,&Xrad[0],&AllRadCnts[0]);
  TGraph *gRadErr     = new TGraph(Rbins,&Xrad[0],&RadErrCnt[0]);
  TGraph *gRadErrUnSc = new TGraph(Rbins,&Xrad[0],&RadErrUnSc[0]);
  ColorGraphic(gRadCnt,4,20,1);
  ColorGraphic(gRadErr,2,20,1);
  ColorGraphic(gRadErrUnSc,2,7,1);

  gRadCnt->SetTitle("Counts vs. Radius");
  gRadCnt->GetXaxis()->SetTitle("Radial Bin (5mm)");
  gRadCnt->GetYaxis()->SetTitle("Counts");
  gRadCnt->GetXaxis()->CenterTitle();
  gRadCnt->GetYaxis()->CenterTitle();

  TCanvas *cRadCnt = new TCanvas("cRadCnt","Radial Counter");
  cRadCnt->cd();
  gRadCnt->Draw("ap");
  gPad->SetGrid();
  gRadErr->Draw("p");  

  TGaxis *axis = new TGaxis(gRadCnt->GetXaxis()->GetXmax(),gRadCnt->GetYaxis()->GetXmin()
			   ,gRadCnt->GetXaxis()->GetXmax(),gRadCnt->GetYaxis()->GetXmax()
                           ,0,gRadCnt->GetYaxis()->GetXmax()/1e6,510,"+L");
  axis->SetLineColor(2);
  axis->SetLabelColor(2);
  axis->SetTitle("Cumulative Statistical Error [%]");
  axis->CenterTitle();
  axis->SetTitleColor(2);
  axis->Draw();

  TPad *gZoom = new TPad("gZoom","Zoomed in Pad",0.55,0.4,0.85,0.6,0,1,1);
  gZoom->SetGrid();
  gZoom->Draw();
  gZoom->cd();
  gRadErrUnSc->Draw("ap");
  gRadErrUnSc->GetXaxis()->SetRangeUser(35,60);
  gRadErrUnSc->GetYaxis()->SetRangeUser(0.4,0.70);
  gRadErrUnSc->GetYaxis()->SetNdivisions(503);
  gRadErrUnSc->GetYaxis()->SetLabelSize(0.08);
  gRadErrUnSc->GetXaxis()->SetLabelSize(0.07);
  gRadErrUnSc->GetXaxis()->SetTitle("Radial Bin (5mm)");
  gRadErrUnSc->GetYaxis()->SetTitle("Statistical Error [%]");
  gRadErrUnSc->SetTitle("Cumulative Stat. Err. 35-60mm");

  cRadCnt->Print("output_files/Radial_Counter.pdf");
  /*
  delete gZoom;
  delete axis;
  delete gRadErr;
  delete gRadCnt;
  delete cRadCnt;*/
};
//--------------------------------------------------------------------------------
void TrackAnodeMPV()
{
  //----------------------------------------------------------------------------------------
  // This compiles the average anode MPV from a landau fit verses run.
 
  vector<Double_t> xNrun,EastAnodeMPV,WestAnodeMPV,EastMPVer,WestMPVer;
    btr[0]->Load_Histograms(bckr[btr[0]->Bkg_index],0);

  TF1 *fLand = new TF1("fLand","landau",0.1,10);

  TH1F *hEastAnodeTot = new TH1F("hEastAnodeTot","East Anode Spectrum;Energy (keV);Counts",
				  btr[0]->hEAnode->GetNbinsX(),btr[0]->hWAnode->GetXaxis()->GetXmin(),
                                  btr[0]->hEAnode->GetXaxis()->GetXmax());

  TH1F *hWestAnodeTot = new TH1F("hWestAnodeTot","West Anode Spectrum;Energy (keV);Counts",
                                  btr[0]->hWAnode->GetNbinsX(),btr[0]->hWAnode->GetXaxis()->GetXmin(),
                                  btr[0]->hWAnode->GetXaxis()->GetXmax());

  TH1F *hEastAnodeTotO = new TH1F("hEastAnodeTotO","East Anode Spectrum;Energy (keV);Counts",
                                  btr[0]->hEAnode->GetNbinsX(),btr[0]->hWAnode->GetXaxis()->GetXmin(),
                                  btr[0]->hEAnode->GetXaxis()->GetXmax());

  TH1F *hWestAnodeTotO = new TH1F("hWestAnodeTotO","West Anode Spectrum;Energy (keV);Counts",
                                  btr[0]->hWAnode->GetNbinsX(),btr[0]->hWAnode->GetXaxis()->GetXmin(),
                                  btr[0]->hWAnode->GetXaxis()->GetXmax());
    btr[0]->Remove_Histograms(bckr[btr[0]->Bkg_index]);
  
  for(Int_t ii = 0 ; ii < nbeta ; ii++){
    btr[ii]->Load_Histograms(bckr[btr[ii]->Bkg_index],0);
     xNrun.push_back(btr[ii]->GetRunNumber());
  
     btr[ii]->hEAnode->Fit("fLand","REMQ");
  
     EastAnodeMPV.push_back(fLand->GetParameter(1));
     EastMPVer.push_back(fLand->GetParError(1));
  
     btr[ii]->hWAnode->Fit(fLand,"REMQ");
     WestAnodeMPV.push_back(fLand->GetParameter(1));
     WestMPVer.push_back(fLand->GetParError(1));

     if(btr[ii]->flipperOn == 0){
	     for(Int_t jj = 1 ; jj <= hEastAnodeTot->GetNbinsX(); jj++){ 
        	hEastAnodeTot->SetBinContent(jj,hEastAnodeTot->GetBinContent(jj) + btr[ii]->hEAnode->GetBinContent(jj));
        	hWestAnodeTot->SetBinContent(jj,hWestAnodeTot->GetBinContent(jj) + btr[ii]->hWAnode->GetBinContent(jj));
     	     }
     } else {
       	    for(Int_t jj = 1 ; jj <= hEastAnodeTot->GetNbinsX(); jj++){
                hEastAnodeTotO->SetBinContent(jj,hEastAnodeTotO->GetBinContent(jj) + btr[ii]->hEAnode->GetBinContent(jj));
                hWestAnodeTotO->SetBinContent(jj,hWestAnodeTotO->GetBinContent(jj) + btr[ii]->hWAnode->GetBinContent(jj));
             }
     }
    btr[ii]->Remove_Histograms(bckr[btr[ii]->Bkg_index]);
  }
  TGraphErrors *gEastAnodeMPV = new TGraphErrors(nbeta,&xNrun[0],&EastAnodeMPV[0],0,&EastMPVer[0]);
  TGraphErrors *gWestAnodeMPV = new TGraphErrors(nbeta,&xNrun[0],&WestAnodeMPV[0],0,&WestMPVer[0]);
  ColorGraphic(gEastAnodeMPV,2,20,1);
  ColorGraphic(gWestAnodeMPV,4,20,1);

  TCanvas *cAnodeTrack = new TCanvas("cAnodeTrack","Track the Anode MPV");
  cAnodeTrack->cd();
  gEastAnodeMPV->Draw("AP");
  gWestAnodeMPV->Draw("P");
  gEastAnodeMPV->GetYaxis()->SetRangeUser(0,4);
  cAnodeTrack->Print("output_files/AnodeMPV.pdf");

  fstream fanode_out;
  fanode_out.open("output_files/anode_spectrum_total.txt",fstream::out);
  for(Int_t ii = 1 ; ii <= hEastAnodeTot->GetNbinsX() ; ii++){
     fanode_out << hEastAnodeTot->GetBinCenter(ii)  << "\t";
     fanode_out << hEastAnodeTot->GetBinContent(ii) << "\t" << hWestAnodeTot->GetBinContent(ii);
     fanode_out << "\t" <<  hEastAnodeTotO->GetBinContent(ii) << "\t" << hWestAnodeTotO->GetBinContent(ii) << endl;
  }
  fanode_out.close();
  /*
  delete fLand; 
  delete gEastAnodeMPV;
  delete gWestAnodeMPV;
  delete hEastAnodeTot;
  delete hWestAnodeTot;
  delete hEastAnodeTotO;
  delete hWestAnodeTotO;
  delete cAnodeTrack;  */
};
//==========================================================================================================
void TrackStats()
{

  vector<Double_t> Cnts,Date;
  Int_t nValidRuns=0;

  fstream p;
  p.open("output_files/stattrack.dat",fstream::out);

  for(Int_t i = 0 ; i < nbeta; i++){
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
	if(btr[i]->heq->Integral() > 0){
		nValidRuns++;
		Cnts.push_back((btr[i]->heq->Integral()*btr[i]->rtime_e 
				 +   btr[i]->hwq->Integral()*btr[i]->rtime_w));
                if(nValidRuns > 1)Cnts[nValidRuns-1] += Cnts[nValidRuns-2];
		Date.push_back(btr[i]->RunDate->Convert());
                
      //          cout << Date[nValidRuns-1] << "\t" << Cnts[nValidRuns-1] <<"\t"<<btr[i]->RunDate->GetDate() <<"\t"<<btr[i]->RunDate->GetMonth()<< "\t" << btr[i]->RunDate->GetYear()<<endl;

		p << Date[nValidRuns-1] << "\t" << (btr[i]->heq->Integral()*btr[i]->rtime_e
                                 +   btr[i]->hwq->Integral()*btr[i]->rtime_w) << endl;
	}
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }
  p.close();

  TGraph *gStats = new TGraph((Int_t)Date.size(),&Date[0],&Cnts[0]);
  
  TCanvas *cStats = new TCanvas("cStats","Statistics",600,400);
  cStats->cd();
 
  ColorGraphic(gStats,2,20,1,1.,"Total Counts","Date","Counts");
  gStats->Draw("APL");
  gStats->GetXaxis()->SetTimeDisplay(1);
  gStats->GetXaxis()->SetTimeFormat("%d-%b");
  gStats->GetXaxis()->SetNdivisions(501);
  cStats->Print("output_files/TrackStats.pdf");

  
};

void average_type1()
{
  TH1F *hE1Type_Pr = new TH1F("hE1Type_Pr",";Energy;Cts",100,0,1);
  TH1F *hE1Type_Sc = new TH1F("hE1Type_Sc",";Energy;Cts",100,0,1);
  TH1F *hW1Type_Pr = new TH1F("hW1Type_Pr",";Energy;Cts",100,0,1);
  TH1F *hW1Type_Sc = new TH1F("hW1Type_Sc",";Eenrgy;Cts",100,0,1);

  for(Int_t i = 0 ; i < nbeta ; i++){
    btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
    for(Int_t j = 1 ; j <= 100 ; j++){
      double epk1 = btr[i]->hEType1_Primary->GetBinCenter(j);
      hE1Type_Pr->SetBinContent(j,hE1Type_Pr->GetBinContent(j) + btr[i]->hEType1_Primary->GetBinContent(j)*btr[i]->rtime_e);
      hE1Type_Sc->SetBinContent(j,hE1Type_Sc->GetBinContent(j) + btr[i]->hEType1_Secondary->GetBinContent(j)*btr[i]->rtime_e);
      hW1Type_Pr->SetBinContent(j,hW1Type_Pr->GetBinContent(j) + btr[i]->hWType1_Primary->GetBinContent(j)*btr[i]->rtime_w);
      hW1Type_Sc->SetBinContent(j,hW1Type_Sc->GetBinContent(j) + btr[i]->hWType1_Secondary->GetBinContent(j)*btr[i]->rtime_w);
    }
    btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
  }

  TCanvas *cType1 = new TCanvas("cType1","Type 1 Primary vs. Secondary",600,600);
  cType1->cd();
  fstream penin;
  penin.open("input_files/type_1_e_fraction_3.txt",fstream::in);
  TH1F *hSimT1p = new TH1F("hSimT1p","East",100,0,1);
  TH1F *hSimT2p = new TH1F("hSimT2P","West",100,0,1);
  for(Int_t i = 1 ; i <= 100 ; i++){
    double h1,h2,h3;
    penin >> h1 >> h2 >> h3;
    hSimT1p->SetBinContent(i,h2);
    hSimT2p->SetBinContent(i,h3);
  }
  penin.close();
  //hSimT2p->Sumw2();
  //hSimT1p->Sumw2();
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
  cType1->Print("output_files/type_1_comp_pr.pdf");

};
//-------------------------------------------------------------------------------------------------------------
void  PlotRunTimes()
{
	vector<Double_t> nRun,livetime,nRunb,livetimeb;
	vector<Double_t> livetimeW,livetimeWb;
	for(Int_t i = 0 ; i < nbeta ; i++){
    		btr[i]->Load_Histograms(bckr[btr[i]->Bkg_index],0);
		nRun.push_back(btr[i]->GetRunNumber());
                livetime.push_back(btr[i]->rtime_e);
		livetimeW.push_back(btr[i]->rtime_w);
		livetimeb.push_back(bckr[btr[i]->Bkg_index]->rtime_e);
		livetimeWb.push_back(bckr[btr[i]->Bkg_index]->rtime_w);
    		btr[i]->Remove_Histograms(bckr[btr[i]->Bkg_index]);
	}
	
	TGraph *gTimesE  = new TGraph((int)nRun.size(),&nRun[0],&livetime[0]);
	TGraph *gTimesW  = new TGraph((int)nRun.size(),&nRun[0],&livetimeW[0]);
	TGraph *gTimesEb = new TGraph((int)nRun.size(),&nRun[0],&livetimeb[0]);
	TGraph *gTimesWb = new TGraph((int)nRun.size(),&nRun[0],&livetimeWb[0]);

        TCanvas *cLiveTime = new TCanvas("cLiveTime","Live Time",800,600);
        cLiveTime->cd();

        ColorGraphic(gTimesE,1,24,1,1,"Live Time","Run Number","Live Time (s)");
        ColorGraphic(gTimesW,1,20,1,1,"Live Time","Run Number");
        ColorGraphic(gTimesEb,2,24,1,1,"Live Time","Run Number");
        ColorGraphic(gTimesWb,2,20,1,1,"Live Time","Run Number");

        gTimesE->Draw("ap");
	gTimesW->Draw("p");
	gTimesEb->Draw("p");
	gTimesWb->Draw("p");
        cLiveTime->Print("output_files/livetimes.pdf");
}

#endif
