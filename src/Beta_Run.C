#ifndef __beta_run_c__
#define __beta_run_c__

#include "Beta_Run.h"

using namespace std;
 
Double_t Int_Asym(Double_t x1, Double_t x2);
Double_t Asym_Err(Double_t x1,Double_t x2);
Double_t fermi(Double_t x);
Double_t Rate_Error(TH1F *h1, TH1F *h2, Float_t t1, Float_t t2);
Double_t Rate_Error(TH1F *h1, TH1F *h2, TH1F *h3,TH1F *h4,Float_t t1, Float_t t2);
Double_t GetIntError(TH1F *h1, TH1F *h2, TH1F *h3,TH1F *h4,Float_t t1, Float_t t2);

//ClassImp(Beta_Run);

Beta_Run::Beta_Run(Int_t n,Int_t m,TSQLServer *sql):Run(n,m,sql)
{
  if(GetBackGround(sql) == 0) cout << "failure" << endl;
  Find_TDC_Cut(n,sql);
  Load_Rotation_Matrix(GetGeo()); // Rotate the detector corrdinates
  Load_Background(GetRunNumber(),GetRunType(),GetGeo());
}

Beta_Run::~Beta_Run()
{

}

void Beta_Run::Load_Histograms(Bck_Run *br,Bool_t SUBBCK,Int_t run=0)
{
     run=this->GetRunNumber();
     //cout << "Load Histograms Run " << run << " " << br->GetRunNumber() << endl;
     this->GetHistograms(run);
     this->ScaleList(this->HEastAn,this->rtime_e);
     this->ScaleList(this->HWestAn,this->rtime_w);
     br->GetHistograms(br->GetRunNumber());
     br->ScaleList(br->HEastAn,br->rtime_e);
     br->ScaleList(br->HWestAn,br->rtime_w);
//     SUBBCK=kFALSE;
     if(SUBBCK)this->SubBck(br);

}

void Beta_Run::Remove_Histograms(Bck_Run *br)
{
   this->DeleteHistos();
   br->DeleteHistos();
}

Int_t Beta_Run::Draw_2d(Int_t nr,Int_t n)
{
  
  hpw = new TH2F(Form("btr[%i]->hpw",n),"West Position;X;Y",150,-75,75,150,-75,75);
  hpe = new TH2F(Form("btr[%i]->hpe",n),"East Position;X;Y",150,-75,75,150,-75,75);
  
  return 0;
}

Int_t Beta_Run::Fill(Int_t n,Int_t remake,Double_t *sep,Int_t nrun)
{
  //-----------------------------------------------------------------------------------+
  // Loop through events in the tree to fill histograms for analysis.  Ideally         |
  // this happens only on the first pass and the histograms are saved to the Analysis  |
  // directory created in the ROOT file.  Later analysis reads in the histograms       |
  // for analysis which greatly speeds things up.                                      |
  //-----------------------------------------------------------------------------------+

  Float_t East_Time_Cnt=0.;
  Float_t West_Time_Cnt=0.;

  Float_t OldTE = 0.;
  Float_t OldTW = 0.;  

  for(Int_t i = 0 ; i < 100 ; i ++){
    // fill the 2/3 separation array for use in the monte carlo type identification routine
    seppar[i] = sep[i];
  }

  if(remake == 1 || !AnalysisDirExist){
    Initialize_hist(0,1,1);
    TFile *f2 = new TFile(Form("%s/hists/spec_%d.root",getenv("UCNAOUTPUTDIR"),nrun),"READ");
    //cout << "Reading " << Form("%s/hists/spec_%d.root",getenv("UCNAOUTPUTDIR"),nrun) << endl;
    hmrIn = (TH1F*)f2->Get("UCN_Mon_4_Rate");
    for(Int_t MRbin=0; MRbin<hmrIn->GetNbinsX(); MRbin++){
     hmr1->Fill(hmrIn->GetBinCenter(MRbin),hmrIn->GetBinContent(MRbin));
    }
    delete hmrIn; 
    f2->Close();

    for(Int_t i = 0 ; i < t1->GetEntries() ; i++){ 
      t1->GetEntry(i);  // Get Enetry from the Tree;
      ScaleTDC();
      // Stupid live time counter.  Live times are up to date in the mysql database so 
      // these are being used to up date the clocks.

	TofE=Tof;
	TofW=Tof;

      if(TofE < OldTE)East_Time_Cnt++;
      OldTE = TofE;
      
      if(TofW < OldTW)West_Time_Cnt++;
      OldTW = TofW;
    
      //-------------------------------
      if(Check_Vetos()){
      // Passed Muon Veto
	Book_Raw(i);
      } else { 
	Book_Muons();
      }
    }

    hGammaCounts->SetBinContent(1,GammasEast);
    hGammaCounts->SetBinContent(2,GammasWest);
    hGammaCountsg->SetBinContent(1,GammasEastg);
    hGammaCountsg->SetBinContent(2,GammasWestg);

    if(GetGeo() == 1){
      // if geometry b correct the live times based on gamma counts
      rtime_e = rtime_e*(hGammaCountsg->GetBinContent(1))/(hGammaCounts->GetBinContent(1)); 
      rtime_w = rtime_w*(hGammaCountsg->GetBinContent(2))/(hGammaCounts->GetBinContent(2)); 
    }
       
    SaveHistograms(kTRUE);
  } else if( remake == 0){
    //GetHistograms(0);
    //cout << "Initializing Open Run num " << GetRunNumber() << endl;
    GetHistograms(GetRunNumber());
    if(GetRunNumber() == 9983){
      cout << "East live time uncorrected " << rtime_e << endl;
      cout << "West live time uncorrected " << rtime_w << endl;
    }
    if(GetGeo() == 1){
      // if geometry b correct the live times based on gamma counts
      rtime_e = rtime_e*(hGammaCountsg->GetBinContent(1))/(hGammaCounts->GetBinContent(1)); 
      rtime_w = rtime_w*(hGammaCountsg->GetBinContent(2))/(hGammaCounts->GetBinContent(2)); 
    }
    if(GetRunNumber() == 9983){
      cout << "East live time corrected " << rtime_e << endl;
      cout << "West live time corrected " << rtime_w << endl;
    }
    SaveHistograms(kFALSE);
  }

  return 0;
}

Int_t Beta_Run::Draw_Hists(Int_t n)
{
  c1->cd(1);
  hpe->Draw("colz");
  
  c1->cd(2);
  hpw->Draw("colz");

  c1->cd(3);
  heq->Draw();
    
  c2->cd(3);
  gPad->SetLogy();
  hte->Draw();
    
  c2->cd(1);
  hEtype_1->Draw();

  c1->cd(4);
  hwq->Draw();

  c2->cd(4);
  gPad->SetLogy();
  htw->Draw();
   
  c2->cd(2);
  hWtype_1->Draw();
  
  c2->cd(5);
  hRote->Draw("colz");
  
  c2->cd(6);
  hRotw->Draw("colz");
  
  delete c1;
  delete c2;
  
  return 0;
}

Int_t Beta_Run::SubBck(Bck_Run *br)
{
  
  //Subtract backgrounds from beta data, determine signal to noise ratios and 
  // do some other shit.
/*
   cout << "Subtracting background for Run " <<  GetRunNumber() << " using run " << br->GetRunNumber() << endl;
   cout << "Background run time " << br->rtime_e << "  " << br->rtime_w << endl;
   cout << "West rates : " <<  hwq->Integral() <<"\t"<<  br->hwq->Integral() << endl;
   cout << "East rates : " <<  heq->Integral() <<"\t"<<  br->heq->Integral() << endl;
*/
   btime_e = br->rtime_e;
   btime_w = br->rtime_w;
   TObjArrayIter *EastIter  = new TObjArrayIter(HEastAn);
   TObjArrayIter *EastIterb = new TObjArrayIter(br->HEastAn);

   TObjArrayIter *WestIter  = new TObjArrayIter(HWestAn);
   TObjArrayIter *WestIterb = new TObjArrayIter(br->HWestAn);
   TH1 *hS,*hB;

   for(Int_t ii = 0 ; ii < HEastAn->GetEntries() ; ii++){
     hS = (TH1*)EastIter->Next();
     hB = (TH1*)EastIterb->Next();
    /* cout << "Signal " << hS->GetName() << "\t bin : " << hS->GetNbinsX() << endl;
     cout << "Background " << hB->GetName() << "\t bin : " << hB->GetNbinsX() << endl;*/
     hS->Add(hB,-1);
     if(ii < 9){
	BkgRtE[ii] = hB->Integral(nlow,nhigh);
	Rate_Error((TH1F*)hS,(TH1F*)hB,rtime_e,br->rtime_e);
     }
   }
  for(Int_t ii = 0 ; ii < HWestAn->GetEntries() ; ii++){
     hS = (TH1*)WestIter->Next();
     hB = (TH1*)WestIterb->Next();
     hS->Add(hB,-1);
     if(ii < 9){
	BkgRtW[ii] = hB->Integral(nlow,nhigh);
	Rate_Error((TH1F*)hS,(TH1F*)hB,rtime_w,br->rtime_w);
     }
   }

  BkgTe = br->rtime_e;
  BkgTw = br->rtime_w;
  Rate_Error(hEAnode23,br->hEAnode23,rtime_e,br->rtime_e);
  Rate_Error(hWAnode23,br->hWAnode23,rtime_w,br->rtime_e);
  
  Rate_Error(hEMWPC,br->hEMWPC,rtime_e,br->rtime_e);
  
  Rate_Error(hENoMWPC,br->hENoMWPC,hEMWPC,br->hEMWPC,rtime_e,br->rtime_e);
  
  Rate_Error(hWMWPC  ,br->hWMWPC ,rtime_w,br->rtime_w);
  
  Rate_Error(hWNoMWPC,br->hWNoMWPC,hWMWPC,br->hWMWPC,rtime_w,br->rtime_w);

  Rate_Error(hEtype_1 ,br->hEtype_1 ,rtime_e,br->rtime_e);
  Rate_Error(hEtype_23,br->hEtype_23,rtime_e,br->rtime_e);
  Rate_Error(hWtype_1 ,br->hWtype_1 ,rtime_w,br->rtime_w);
  Rate_Error(hWtype_23,br->hWtype_23,rtime_w,br->rtime_w);
 
  TF1 *fMWPCL = new TF1("fMWPCL","[0]",100,1000);
  if(hEMWPC->Integral(nlow,nhigh) !=0){
    MWPC_RatioE   = hENoMWPC->Integral(nlow,nhigh) / hEMWPC->Integral(nlow,nhigh);
    MWPC_RatioE_e = GetIntError(hENoMWPC,br->hENoMWPC,hEMWPC,br->hEMWPC,rtime_e,br->rtime_e);
  
    hENoMWPC->Fit("fMWPCL","QRME","goff");

    MWPC_RatioE_fit   = (nhigh-nlow)*fMWPCL->GetParameter(0) / hEMWPC->Integral(nlow,nhigh);
    MWPC_RatioE_fit_e = (nhigh-nlow)*fMWPCL->GetParError(0) / hEMWPC->Integral(nlow,nhigh);
  }
  if(hWMWPC->Integral(nlow,nhigh) != 0 ){
    MWPC_RatioW   = hWNoMWPC->Integral(nlow,nhigh)/ hWMWPC->Integral(nlow,nhigh);
    MWPC_RatioW_e = GetIntError(hWNoMWPC,br->hWNoMWPC,hWMWPC,br->hWMWPC,rtime_w,br->rtime_w);
			      
    hWNoMWPC->Fit("fMWPCL","QRME","goff");		      
			      
    MWPC_RatioW_fit   = (nhigh-nlow)*fMWPCL->GetParameter(0) / hWMWPC->Integral(nlow,nhigh);
    MWPC_RatioW_fit_e = (nhigh-nlow)*fMWPCL->GetParError(0) / hWMWPC->Integral(nlow,nhigh);
  }
  for(int i = 1 ; i <= heq->GetNbinsX() ; i++){
      hESigNos->SetBinContent(i,heq->GetBinContent(i));
      hESigNos->SetBinError(i,heq->GetBinError(i));
      hWSigNos->SetBinContent(i,hwq->GetBinContent(i));
      hWSigNos->SetBinError(i,hwq->GetBinError(i));
  }
  
   hETimeVsE->Add(br->hETimeVsE,-1);
   hWTimeVsE->Add(br->hWTimeVsE,-1);
   
  E_Sig_Nos = heq->Integral(nlow,nhigh)/br->heq->Integral(nlow,nhigh);
  W_Sig_Nos = hwq->Integral(nlow,nhigh)/br->hwq->Integral(nlow,nhigh);
 
  hESigNos->Divide(br->heq);
  hWSigNos->Divide(br->hwq);
  
  hAposw->Add(br->hAposw,-1);
  hApose->Add(br->hApose,-1);
  
//   hETDC_cor->Add(br->hETDC_cor,-1);
//   hWTDC_cor->Add(br->hWTDC_cor,-1);
// 
//   hETimeVsE->Add(br->hETimeVsE,-1);
//   hWTimeVsE->Add(br->hWTimeVsE,-1);
//   
//   hEmuon->Add(br->hEmuon,-1);
//   hWmuon->Add(br->hWmuon,-1);
  
  for(int i = 0 ; i < 12 ; i++){
//     hERad[i]->Add(br->hERad[i],-1);
//     hWRad[i]->Add(br->hWRad[i],-1);
//     
//     Rate_Error(hERad[i],br->hERad[i],rtime_e,br->rtime_e);
//     Rate_Error(hWRad[i],br->hWRad[i],rtime_w,br->rtime_w);

    erad[i]    = hERad[i]->Integral(nlow,nhigh);
    wrad[i]    = hWRad[i]->Integral(nlow,nhigh);
	  
    erader[i]  = sqrt( erad[i]*rtime_e + 
		       br->hERad[i]->Integral(nlow,nhigh)*br->rtime_e)/rtime_e;
    wrader[i]  = sqrt( wrad[i]*rtime_w + 
		       br->hWRad[i]->Integral(nlow,nhigh)*br->rtime_w)/rtime_w;
  }
  
  emuon_rate = hEmuon->Integral();
  wmuon_rate = hWmuon->Integral();
 
  bscat.etype_1_bck  = br->heqF->Integral(nlow,nhigh)*br->rtime_e;
  bscat.etype_23_bck = br->heqG->Integral(nlow,nhigh)*br->rtime_e;
  
  bscat.wtype_1_bck  = br->hwqF->Integral(nlow,nhigh)*br->rtime_w;
  bscat.wtype_23_bck = br->hwqG->Integral(nlow,nhigh)*br->rtime_w;
  
  bscat.w_all_bck    = br->hwq->Integral(nlow,nhigh)*br->rtime_w;
  bscat.e_all_bck    = br->heq->Integral(nlow,nhigh)*br->rtime_e;
  
  delete WestIter;
  delete WestIterb;
  delete EastIter;
  delete EastIterb;
  
  return 0;
}

Bool_t Beta_Run::Check_Vetos()
{
 
  return kTRUE;
  
  if(TaggedTopE + TaggedBackE + TaggedDriftE 
     + TaggedBackW + TaggedDriftW != 0 ) return kFALSE;
  
  return kTRUE; 
}

Int_t Beta_Run::GetBackGround(TSQLServer *sql)
{
  // Open Connection to that analysis database
  if(!(sql->IsConnected()))ReConnect(sql);
  TSQLResult *res;
  TSQLRow *row;
  char buffer[300];
 
  // Fetch the octet type of run from DB
  sprintf(buffer,"select run_number,asym_oct from run where run_number = %d ",GetRunNumber());
  res = (TSQLResult*)sql->Query(buffer);
  if(res->GetRowCount() != 0){
    while((row = (TSQLRow*)res->Next())){
      sprintf(octtype,"%s",row->GetField(1));
      delete row;
    }
  }
  
  delete res;
  // Find the matching background run in the database
  // cout << "at run " << GetRunNumber() << "  " << octtype << endl;
  if(!strncmp(octtype,"A2",2)) {
    sprintf(buffer," select t1.run_number,t1.geometry,t1.asym_oct,t1.gate_valve " 
	    " from run as t1, analysis as t2 where t1.run_number between %d and %d and "
	    " t1.run_type = 'Asymmetry' and t1.asym_oct = 'A1'"/* and t2.live_time_e > 0"
	    " and t2.run_number=t1.run_number"*/,GetRunNumber()-5,GetRunNumber()+5);
  } else if(!strncmp(octtype,"A5",2)){
    sprintf(buffer," select t1.run_number,t1.geometry,t1.asym_oct,t1.gate_valve " 
	    " from run as t1, analysis as t2 where t1.run_number between %d and %d and "
	    " t1.run_type = 'Asymmetry' and t1.asym_oct = 'A4'"/*
	    " and t2.run_number = t1.run_number and t2.live_time_e >0 "*/,GetRunNumber()-5,GetRunNumber()+5);
  } else if(!strncmp(octtype,"A7",2)){
    sprintf(buffer," select t1.run_number,t1.geometry,t1.asym_oct,t1.gate_valve " 
	    " from run as t1, analysis as t2 where t1.run_number between %d and %d and "
	    " t1.run_type = 'Asymmetry' and t1.asym_oct = 'A9'"/*
	    " and t2.run_number = t1.run_number and t2.live_time_e >0 "*/,GetRunNumber()-5,GetRunNumber()+5);
  } else if(!strncmp(octtype,"A10",2)){
    sprintf(buffer," select t1.run_number,t1.geometry,t1.asym_oct,t1.gate_valve "
	    " from run as t1, analysis as t2 where t1.run_number between %d and %d and "
	    " t1.run_type = 'Asymmetry' and t1.asym_oct = 'A12'"/*
	    "and t2.run_number = t1.run_number and t2.live_time_e >0 "*/,GetRunNumber()-5,GetRunNumber()+5);
  } else if(!strncmp(octtype,"B2",2)){
    if((GetRunNumber() < 8277 && GetRunNumber() > 9189) || GetRunNumber() > 12400){
    sprintf(buffer," select t1.run_number,t1.geometry,t1.asym_oct,t1.gate_valve "
	    " from run as t1, analysis as t2 where t1.run_number between %d and %d and "
	    " t1.run_type = 'Asymmetry' and t1.asym_oct = 'B1'"/*
	    " and t1.run_number = t2.run_number and t2.live_time_e >0"*/,GetRunNumber()-7,GetRunNumber()+7);
    } else {
      sprintf(buffer," select t1.run_number,t1.geometry,t1.asym_oct,t1.gate_valve "
	    " from run as t1, analysis as t2 where t1.run_number between %d and %d and "
	    " t1.run_type = 'Asymmetry' and t1.asym_oct = 'A12'"/*
	    " and t1.run_number = t2.run_number and t2.live_time_e >0"*/,GetRunNumber()-5,GetRunNumber()+5);
    }
  } else if(!strncmp(octtype,"B5",2)){
    sprintf(buffer," select t1.run_number,t1.geometry,t1.asym_oct,t1.gate_valve "
	    " from run as t1, analysis as t2 where t1.run_number between %d and %d and "
	    " t1.run_type = 'Asymmetry' and t1.asym_oct = 'B4'"
	    /*" and t1.run_number = t2.run_number and t2.live_time_e >0"*/,GetRunNumber()-5,GetRunNumber()+5);
  } else if(!strncmp(octtype,"B7",2)){
    sprintf(buffer," select t1.run_number,t1.geometry,t1.asym_oct,t1.gate_valve "
	    " from run as t1, analysis as t2 where t1.run_number between %d and %d and "
	    " t1.run_type = 'Asymmetry' and t1.asym_oct = 'B9'"
	   /* " and t2.run_number = t1.run_number and t2.live_time_e >0 "*/,GetRunNumber()-5,GetRunNumber()+5);
  } else if(!strncmp(octtype,"B10",2)){
    sprintf(buffer," select t1.run_number,t1.geometry,t1.asym_oct,t1.gate_valve "
	    " from run as t1, analysis as t2 where t1.run_number between %d and %d and "
	    " t1.run_type = 'Asymmetry' and t1.asym_oct = 'B12'"
	    /*" and t2.run_number = t1.run_number and t2.live_time_e >0 "*/,GetRunNumber()-5,GetRunNumber()+5);
  }
  
  res = (TSQLResult*)sql->Query(buffer);
  runb = 0;
  if(res->GetRowCount() != 0){
    while((row = (TSQLRow*)res->Next())){
      runb = atoi(row->GetField(0));
      delete row;
    }
  }
  
  //cout << GetRunNumber() << " bk is " << runb << endl;

  // the following list is a fix for bad background runs...  
  if(GetRunNumber() == 15321) runb = 15320;
  if(GetRunNumber() == 16190) runb = 16189;
  if(GetRunNumber() == 14419) runb = 14418;
  if(GetRunNumber() > 15111 && GetRunNumber() < 15115) runb = 15102;
  if(GetRunNumber() == 13950) runb = 13943;
  if(GetRunNumber() == 12531) runb = 12524;
  if(GetRunNumber() == 12647) runb = 12646;
  if(GetRunNumber() == 12443) runb = 12442;
  if(GetRunNumber() == 11048) runb = 11058;
  if(GetRunNumber() == 11022) runb = 11032;
  if(GetRunNumber() == 11007) runb = 11020;
  if(GetRunNumber() == 10891) runb = 10905;
  if(GetRunNumber() == 10545) runb = 10539;
  if(GetRunNumber() == 10535) runb = 10552;
  if(GetRunNumber() == 10308) runb = 10321;
  if(GetRunNumber() == 10252) runb = 10263;
  if(GetRunNumber() == 10281) runb = 10291;
  if(GetRunNumber() == 10148) runb = 10160;
  if(GetRunNumber() == 10011) runb = 10022;
  if(GetRunNumber() == 9754)  runb = 9771;
  if(GetRunNumber() == 9163)  runb = 9172;
  if(GetRunNumber() == 9024)  runb = 9029;
  if(GetRunNumber() == 9788)  runb = 9798;
  if(GetRunNumber() == 9718)  runb = 9728;
  if(GetRunNumber() == 9577)  runb = 9568;
  if(GetRunNumber() == 9430)  runb = 9440;
  if(GetRunNumber() == 9085)  runb = 9087;
  if(GetRunNumber() == 8658)  runb = 8662;
  if(GetRunNumber() == 9861)  runb = 9858;
  if(GetRunNumber() == 9759)  runb = 9756;
  if(GetRunNumber() == 9610)  runb = 9603;
  if(GetRunNumber() == 9387)  runb = 9398;
  if(GetRunNumber() == 9596)  runb = 9593;
  if(GetRunNumber() == 9594)  runb = 9593;
  if(GetRunNumber() == 9052)  runb = 9046;
  if(GetRunNumber() == 8993)  runb = 8984;
  if(GetRunNumber() == 8864)  runb = 8853;
  if(GetRunNumber() == 8854)  runb = 8853;
  if(GetRunNumber() == 8740)  runb = 8744;
  if(GetRunNumber() == 8985)  runb = 8984;
  if(GetRunNumber() == 8814)  runb = 8826;
  if(GetRunNumber() == 8706)  runb = 8716;
  if(GetRunNumber() == 8595)  runb = 8605;
  if(GetRunNumber() == 8570)  runb = 8581;
  if(GetRunNumber() == 8518)  runb = 8529;
  if(GetRunNumber() == 8312)  runb = 8377;
  if(GetRunNumber() == 8278)  runb = 8292;
  if(GetRunNumber() == 8008)  runb = 8007;
  if(GetRunNumber() == 7942)  runb = 7941;
  if(GetRunNumber() == 7921)  runb = 7920;
  if(GetRunNumber() == 7879)  runb = 7878;
  if(GetRunNumber() == 7876)  runb = 7875;
  if(GetRunNumber() == 7873)  runb = 7872;
  if(GetRunNumber() == 7831)  runb = 7830;
  if(GetRunNumber() == 7729)  runb = 7728;
  if(GetRunNumber() == 7709)  runb = 7707;
  if(GetRunNumber() == 7702)  runb = 7701;
  if(GetRunNumber() == 7663)  runb = 7662;
  if(GetRunNumber() == 7687)  runb = 7686;
  
  if(runb == 0) cout << "failed to find background for run " << GetRunNumber() << endl;

  delete res;
    
  return runb;
}
//------------------------------------------------------------------------------------
Double_t Int_Asym(Double_t x1,Double_t x2)
{
  // Calculate the Simple Asymmetry 
  return (x1-x2)/(x1+x2);
}
//------------------------------------------------------------------------------------
Double_t Asym_Err(Double_t x1,Double_t x2)
{
  using namespace TMath;

  return 2.*sqrt(x1*x2 / Power(x1+x2,3));
}
//------------------------------------------------------------------------------------
Int_t Beta_Run::GetSimpleAsym()
{
  // Calculate the Simple Rate Difference Integral Asymmetry 
  using namespace TMath;
  
  for(Int_t i = 0 ; i < 40 ; i++){
     fAsym_Run[i]     = Int_Asym(heq->Integral(i+1,60),hwq->Integral(i+1,60));
     fAsym_Run_Err[i] = Asym_Err(heq->Integral(i+1,60),hwq->Integral(i+1,60));
  }
  
  return -1;
}
//-------------------------------------------------------------------------------------
Double_t Beta_Run::GetEnergyChi()
{
  using namespace TMath;
  
  TF1 *fline = new TF1("fline","[0]",900.,1500.);
  
  fstream fref;
  if(GetRunNumber() >=  7659 && GetRunNumber() <= 9189){
    fref.open("input_files/histogram_geoA.txt",fstream::in);
  } else if(GetRunNumber() >= 9356 && GetRunNumber() <= 10333) {
    fref.open("input_files/histogram_geoB.txt",fstream::in);
  } else if(GetRunNumber() >= 10404 && GetRunNumber() <= 11095) {
    fref.open("input_files/histogram_geoC.txt",fstream::in);
  } else if(GetRunNumber() > 13905){
    fref.open("input_files/histogram_e.txt",fstream::in);
  } else{ 
    fref.open("input_files/histogram_e.txt",fstream::in);
  }
  Double_t temp1,temp2,temp3,temp4;
  
  Int_t ibin = 1;
  do{
     fref >> temp1 >> temp2 >> temp3 >> temp4;
     hEERef->SetBinContent(ibin,temp2);
     hEERef1->SetBinContent(ibin,temp3);
     hEERef2->SetBinContent(ibin,temp4);

//	cout << ibin << "th Bin, contents " << temp2 << " " << temp3 << " " << temp4 << endl;

     ibin++;
  }while(ibin < 101);
  
  fref.close();
  
  heq->Fit(fline,"RMEQ","goff");
  Residual_Bkg[0] = fline->GetParameter(0);
  Resid_Bkger[0]  = fline->GetParError(0);
  
  hwq->Fit(fline,"RMEQ","goff");
  Residual_Bkg[1] = fline->GetParameter(0);
  Resid_Bkger[1]  = fline->GetParError(0);
  
  Chi = 0.;
  
  for(Int_t i = 1 ; i <= heq->GetNbinsX() ; i++){
    Chi += Diff[i-1]/80;
  }
   
  delete fline;
  
  return 1.;
}

Int_t Beta_Run::Make_Kurie()
{
  // Simple function to calculate the Kurie plots and find the
  // end point energy east and west sides.
  
  TF1 *fline = new TF1("fline","[0]+[1]*x",200,650);
  FillKurie(heq,hEKurie);  
  hEKurie->Fit(fline,"RMEQ","");
  E_Endpoint = TMath::Abs(fline->GetParameter(0) / fline->GetParameter(1));

  E_EndError = E_Endpoint*sqrt(
			       TMath::Power(fline->GetParError(0)
					    /fline->GetParameter(0),2)
			       +TMath::Power(fline->GetParError(1)/
					     fline->GetParameter(1),2)
			       );
  FillKurie(hwq,hWKurie);
  
  hWKurie->Fit(fline,"RMQE","");
  
  W_Endpoint = TMath::Abs(fline->GetParameter(0) / fline->GetParameter(1));
  
  W_EndError = W_Endpoint*sqrt(
			       TMath::Power(fline->GetParError(0)
					    /fline->GetParameter(0),2)
			       +TMath::Power(fline->GetParError(1)/
					     fline->GetParameter(1),2)
			       );
  
  
  if(W_Endpoint < 500 || W_Endpoint > 1000)
    cout << "Odd West End Point " << W_Endpoint << "\t for Run " << GetRunNumber() << endl;
  
  if(E_Endpoint < 500 || E_Endpoint > 1000)
    cout << "Odd East End Point " << E_Endpoint << "\t for Run " << GetRunNumber() << endl;
  
 delete fline;
 
 return -1;
}

void Beta_Run::FillKurie(TH1F *hEspec,TH1F *hKurie)
{
  Float_t mass = 510.9991; // Electron Rest mass
  Float_t xt   = 0.;
  Float_t pt   = 0.;
  Float_t et   = 0.;
  Float_t kt   = 0.;
  Float_t kter = 0.;

  for(Int_t i = 1 ; i <= hEspec->GetNbinsX() ; i++){
    xt = hEspec->GetBinContent(i);
    et = hEspec->GetBinCenter(i);
    pt = et*(et+2.*mass);
    if(xt > 0){
      kt   = sqrt(xt / (pt*fermi(hEspec->GetBinCenter(i))));
      kter = kt*(hEspec->GetBinError(i)/hEspec->GetBinContent(i));
    }else{
      kt   = 0;
      kter = 0.;
    }
    hKurie->SetBinContent(i,kt);
    hKurie->SetBinError(i,kter);
  }

}
//-------------------------------------------------------------------------
Double_t fermi(Double_t x)
{
  //  Expansion of the Fermi function from Wilkinson 1981 Nucl Physics A377.
  using namespace TMath;
  
  if(x == 0.) x+=50.;
  
  Double_t m = 510.9991; // electron rest mass;
  Double_t perm = 8.85418e-12; // 
  Double_t alpha = Qe()*Qe()/(Hbar()*C()*4.*Pi()*perm); // fine structure constant
  Double_t R = 1.e-15; // Proton Raduis
  
  Double_t W = x/m + 1.; // Total energy in terms of electron rest mass;
  Double_t p = sqrt(W*W-1.); // momentum 
  Double_t pi = Pi();  
  Double_t gammae = 0.577215; // Euler's constant
    
  //  Double_t RC = 0.0390;

  Double_t a0 = 1.;
  Double_t a1 = pi*W/p;
  Double_t a2 = 11./4. - gammae - log(2.*p*R) + pi*pi*(W/p)*(W/p);
  Double_t a3 = pi*(W/p)*(11./4. - gammae - log(2.*p*R));
  
  Double_t f = a0 + a1*alpha + a2*alpha*alpha + a3*alpha*alpha*alpha;

  if( IsNaN(f) == 1){
    f = 1.;
  } 
  
  return f;
}

//-----------------------------------------------------------------------
Double_t Rate_Error(TH1F *h1, TH1F *h2, Float_t t1,Float_t t2)
{
  
  // Rate Error in a background subtracted spectra.
  // Bin error should be correctly set in the Run::Scale2Time function
  // Thus for background subtraction, h1 and h2 are in rates and 
  // the errors per bin are in rates so dR_beta = sqrt(dR_s^2 + dR_b^2)
  //

  using namespace TMath;

  Double_t error = 0.;
  
  for(Int_t i = 1 ; i <= h1->GetNbinsX() ; i++){
  
    error = Sqrt(Abs(h1->GetBinError(i)*h1->GetBinError(i) + h2->GetBinError(i)*h2->GetBinError(i)));
    h1->SetBinError(i,error);
				          
  }
  
  return -1.e3;
}
//----------------------------------------------------------------------------------------
Double_t Rate_Error(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, Float_t t1,Float_t t2)
{
  
  //----------------------------------------------------------------------------------------------
  // Calculate the error in the MWPC on off test.  With St = (Snomwpc - Bnomwpc) - (Smwpc -Bmwpc)
  // These look complicated becasuse the h1 = St, so its all ready by background and Smwpc subtracted.
  //----------------------------------------------------------------------------------------------
  
  using namespace TMath;
  Double_t error = 0.;
  for(Int_t i = 1 ; i <= h1->GetNbinsX() ; i++){

    error =  Sqrt(Abs((h1->GetBinContent(i)+2.*h3->GetBinContent(i) +h2->GetBinContent(i) +
		       h4->GetBinContent(i))*t1 +
		      (h2->GetBinContent(i) + h4->GetBinContent(i))*t2))/t1;
    
    if(error <= 0) error = 1./t1;

    h1->SetBinError(i,error);
				          
  }
  
  return -1.e3;
}
//----------------------------------------------------------------------------------------

Double_t GetIntError(TH1F *h1, TH1F *h2, TH1F *h3,TH1F *h4,Float_t t1, Float_t t2)
{
  //---------------------------------------------------
  // Calculate the error in the mwpc on / off ratio
  //---------------------------------------------------
  
  using namespace TMath;

  Double_t S  = (h1->Integral(nlow,nhigh) + h2->Integral(nlow,nhigh) + h3->Integral(nlow,nhigh))*t1;
  Double_t Bb = h2->Integral(nlow,nhigh)*t2;
  Double_t Sn = (h3->Integral(nlow,nhigh)+h4->Integral(nlow,nhigh))*t1;
  Double_t Bn = h4->Integral(nlow,nhigh)*t2;

  Double_t Tm1 = Bb / Power(Sn-Bn,2);
  Double_t Tm2 = Sn / Power(Sn-Bn,2);
  Double_t Tm3 = Sn*Power( 1./(Sn-Bn) + (S-Sn+Bn-Bb)/Power(Sn-Bn,2),2);
  Double_t Tm4 = Bn*Power( 1./(Sn-Bn) + (S-Sn+Bn-Bb)/Power(Sn-Bn,2),2);

  Double_t er = Sqrt(Tm1+Tm2+Tm3+Tm4);
  
  return er;
		
};

#endif
