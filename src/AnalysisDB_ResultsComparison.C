#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSQLServer.h>
#include <TGraph.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <TF1.h>
 
struct OctResult_t{
  vector<Double_t> nRun;
  vector<Double_t> nArun;
  vector<Double_t> EType0;
  vector<Double_t> EType1;
  vector<Double_t> EType23;
  vector<Double_t> WType0;
  vector<Double_t> WType1;
  vector<Double_t> WType23;
  
  vector<Double_t> A_A;
  vector<Double_t> A_B;
  vector<Double_t> A_C;
  vector<Double_t> A_D;
  vector<Double_t> A_E;
  vector<Double_t> A_F;
  vector<Double_t> A_G;
  vector<Double_t> A_H;
  vector<Double_t> A_I;
  vector<Double_t> A_Aer;
  vector<Double_t> A_Ber;
  vector<Double_t> A_Cer;
  vector<Double_t> A_Der;
  vector<Double_t> A_Eer;
  vector<Double_t> A_Fer;
  vector<Double_t> A_Ger;
  vector<Double_t> A_Her;
  vector<Double_t> A_Ier;
};
//
Bool_t checkruns(Int_t currentrun, vector<Double_t> openruns,Int_t nruns)
{
  //simple check to make sure a file isn't already open
  //if(openruns.size() > 0)return kTRUE;

  for(Int_t ii = 0 ; ii < (Int_t)openruns.size() ; ii++){
	if(currentrun == openruns[ii]){
	  cout << currentrun << "\t" << openruns[ii] << endl;
	  return kFALSE;
	}
  }
      
  return kTRUE;
};
//======================================================================================
void SetDifference(OctResult_t Rp,OctResult_t Mp,OctResult_t &Df)
{
  Int_t rpIndex,mpIndex;
  
  for(Int_t ii = 0 ; ii < (int)Rp.nRun.size() ; ii++){
    mpIndex = -1;
    for(Int_t jj = 0 ; jj < (int)Mp.nRun.size() ; jj++){
      if(Rp.nRun[ii] == Mp.nRun[jj])mpIndex = jj;
    }
    if(mpIndex != -1){
      Df.nRun.push_back(Rp.nRun[ii]);
      Df.EType0.push_back(Rp.EType0[ii]   - Mp.EType0[mpIndex]);
      Df.EType1.push_back(Rp.EType1[ii]   - Mp.EType1[mpIndex]);
      Df.EType23.push_back(Rp.EType23[ii] - Mp.EType23[mpIndex]);
      Df.WType0.push_back(Rp.WType0[ii]   - Mp.WType0[mpIndex]);
      Df.WType1.push_back(Rp.WType1[ii]   - Mp.WType1[mpIndex]);
      Df.WType23.push_back(Rp.WType23[ii] - Mp.WType23[mpIndex]); 
    }
  }
  
  
};
void PrintSizes(OctResult_t Oct)
{
  cout << "Sizes " << endl;
  cout << "Number of Unique Runs found " << (int)Oct.nRun.size() << endl;
  cout << "Type 0 East :  " << (int)Oct.EType0.size() << endl;
  cout << "Type 1 East :  " << (int)Oct.EType1.size() << endl;
  cout << "Type 23 East :  " << (int)Oct.EType23.size() << endl;
  cout << "Type 0 West :  " << (int)Oct.WType0.size() << endl;
  cout << "Type 1 West :  " << (int)Oct.WType1.size() << endl;
  cout << "Type 23 West :  " << (int)Oct.WType23.size() << endl;
}
//-------------------------------------------------------------------------------------
void PlotRates(OctResult_t MPMres, OctResult_t RWPres, OctResult_t Diff)
{
   
  TGraph *gRWPT0 = new TGraph((int)RWPres.EType0.size(),&RWPres.nRun[0],
			      &RWPres.EType0[0]);
  TGraph *gRWPT1 = new TGraph((int)RWPres.EType1.size(),&RWPres.nRun[0],
			      &RWPres.EType1[0]);
  TGraph *gRWPT23 = new TGraph((int)RWPres.EType23.size(),&RWPres.nRun[0],
			      &RWPres.EType23[0]);
  
  TGraph *gMPMT0 = new TGraph((int)MPMres.EType0.size(),&MPMres.nRun[0],
			      &MPMres.EType0[0]);
  TGraph *gMPMT1 = new TGraph((int)MPMres.EType1.size(),&MPMres.nRun[0],
			      &MPMres.EType1[0]);
  TGraph *gMPMT23 = new TGraph((int)MPMres.EType23.size(),&MPMres.nRun[0],
			      &MPMres.EType23[0]);
  
  TGraph *gRWPT0w = new TGraph((int)RWPres.WType0.size(),&RWPres.nRun[0],
			      &RWPres.WType0[0]);
  TGraph *gRWPT1w = new TGraph((int)RWPres.WType1.size(),&RWPres.nRun[0],
			      &RWPres.WType1[0]);
  TGraph *gRWPT23w = new TGraph((int)RWPres.WType23.size(),&RWPres.nRun[0],
			      &RWPres.WType23[0]);
  
  TGraph *gMPMT0w = new TGraph((int)MPMres.WType0.size(),&MPMres.nRun[0],
			      &MPMres.WType0[0]);
  TGraph *gMPMT1w = new TGraph((int)MPMres.WType1.size(),&MPMres.nRun[0],
			      &MPMres.WType1[0]);
  TGraph *gMPMT23w = new TGraph((int)MPMres.WType23.size(),&MPMres.nRun[0],
			      &MPMres.WType23[0]);
  
  TGraph *gDiffT0 = new TGraph((int)Diff.EType0.size(),&Diff.nRun[0],
			      &Diff.EType0[0]);
  TGraph *gDiffT1 = new TGraph((int)Diff.EType1.size(),&Diff.nRun[0],
			      &Diff.EType1[0]);
  TGraph *gDiffT23 = new TGraph((int)Diff.EType23.size(),&Diff.nRun[0],
			      &Diff.EType23[0]);
  
  TGraph *gDiffT0w = new TGraph((int)Diff.WType0.size(),&Diff.nRun[0],
			      &Diff.WType0[0]);
  TGraph *gDiffT1w = new TGraph((int)Diff.WType1.size(),&Diff.nRun[0],
			      &Diff.WType1[0]);
  TGraph *gDiffT23w = new TGraph((int)Diff.WType23.size(),&Diff.nRun[0],
			      &Diff.WType23[0]);
  
  gRWPT0->SetMarkerStyle(20);
  gRWPT1->SetMarkerStyle(20);
  gRWPT23->SetMarkerStyle(20);
  
  gRWPT0w->SetMarkerStyle(20);
  gRWPT1w->SetMarkerStyle(20);
  gRWPT23w->SetMarkerStyle(20);
  
  gMPMT0->SetMarkerStyle(24);
  gMPMT0->SetMarkerColor(2);
  gMPMT1->SetMarkerStyle(24);
  gMPMT1->SetMarkerColor(2);
  gMPMT23->SetMarkerStyle(24);
  gMPMT23->SetMarkerColor(2);
  
  gMPMT0w->SetMarkerStyle(24);
  gMPMT0w->SetMarkerColor(2);
  gMPMT1w->SetMarkerStyle(24);
  gMPMT1w->SetMarkerColor(2);
  gMPMT23w->SetMarkerStyle(24);
  gMPMT23w->SetMarkerColor(2);
  
  
  TCanvas *cEvent = new TCanvas("cEvents","Events");
  cEvent->Divide(2,3);
  cEvent->cd(1);
  gDiffT0->Draw("ap");
  cEvent->cd(3);
  gDiffT1->Draw("ap");
  cEvent->cd(5);
  gDiffT23->Draw("ap");
  cEvent->cd(2);
  gDiffT0w->Draw("ap");
  cEvent->cd(4);
  gDiffT1w->Draw("ap");
  cEvent->cd(6);
  gDiffT23w->Draw("ap");
 
 
};

void FillRates(OctResult_t &R,char query[500],TSQLServer *sql)
{
  
  TSQLResult *res;
  TSQLRow *row;
  
  res = (TSQLResult*)sql->Query(query);
  cout << "Number of Rows found : " << res->GetRowCount() << endl;
  
  if(res->GetRowCount() != 0){
    Int_t e0=0,e1=0,e2=0;
    Int_t w0=0,w1=0,w2=0; 
    while((row = (TSQLRow*)res->Next())){
      if(strncmp(row->GetField(1),"East",4)==0){
	
	if((int)R.nRun.size() == 0)
		  R.nRun.push_back(atoi(row->GetField(2)));
	else if(checkruns(atoi(row->GetField(2)),R.nRun,(Int_t)R.nRun.size())){
	          R.nRun.push_back(atoi(row->GetField(2)));
		  e0=0;
		  e1=0;
		  e2=0;
		  w0=0;
		  w1=0;
		  w2=0;
	}
	
	if(strncmp(row->GetField(3),"0",2) == 0 && e0 == 0){
		  R.EType0.push_back(atof(row->GetField(4)));
		  e0=1;
	}else if(strncmp(row->GetField(3),"I",2) == 0 && e1 == 0){
		  R.EType1.push_back(atof(row->GetField(4)));
		  e1=1;
	}else if(strncmp(row->GetField(3),"II,III",6) == 0 && e2 == 0){
		  R.EType23.push_back(atof(row->GetField(4)));
		  e2=1;
	}
	
      }	else if(strncmp(row->GetField(1),"West",4)==0){
	if(strncmp(row->GetField(3),"0",2) == 0 && w0==0){
		  R.WType0.push_back(atof(row->GetField(4)));
		  w0=1;
	}else if(strncmp(row->GetField(3),"I",2) == 0 && w1 == 0){
		  R.WType1.push_back(atof(row->GetField(4)));
		  w1=1;
	}else if(strncmp(row->GetField(3),"II,III",6) == 0 && w2 == 0){
		  R.WType23.push_back(atof(row->GetField(4)));
		  w2=1;
	}
      }
      delete row;
    };
  };
  delete res;
  
}

void FillMPMRates(OctResult_t &R,char query[500],TSQLServer *sql)
{
  
  TSQLResult *res;
  TSQLRow *row;
  
  res = (TSQLResult*)sql->Query(query);
  cout << "Number of Rows found : " << res->GetRowCount() << endl;
  if(res->GetRowCount() != 0){
    Int_t e0=0,e1=0,e2=0;
    Int_t w0=0,w1=0,w2=0;
    while((row = (TSQLRow*)res->Next())){
      if(strncmp(row->GetField(1),"East",4)==0){
	
	if((int)R.nRun.size() == 0)
		  R.nRun.push_back(atoi(row->GetField(2)));
	else if(checkruns(atoi(row->GetField(2)),R.nRun,(Int_t)R.nRun.size())){
	          R.nRun.push_back(atoi(row->GetField(2)));
		  e0=0;
		  e1=0;
		  e2=0;
		  w0=0;
		  w1=0;
		  w2=0;
	}
	
	if(strncmp(row->GetField(3),"0",1) == 0 && e0==0){
		  R.EType0.push_back(atof(row->GetField(4)));
		  e0=1;
	}else if(strncmp(row->GetField(3),"I",2) == 0 && e1==0){
		  R.EType1.push_back(atof(row->GetField(4)));
		  e1=1;
	}else if(strncmp(row->GetField(3),"II",2) == 0 && e2==0){
		  R.EType23.push_back(atof(row->GetField(4)));
		  e2=1;
	}else if(strncmp(row->GetField(3),"III",3) == 0 && e2 == 1){
		  R.EType23.back() += atof(row->GetField(4));
		  e2=2;
	}
      }	else if(strncmp(row->GetField(1),"West",4)==0){
	if(strncmp(row->GetField(3),"0",1) == 0 && w0 == 0){
		  R.WType0.push_back(atof(row->GetField(4)));
		  w0=1;
	}else if(strncmp(row->GetField(3),"I",2) == 0 && w1 ==0){
		  R.WType1.push_back(atof(row->GetField(4)));
		  w1 =1;
	}else if(strncmp(row->GetField(3),"II",2) == 0 && w2 == 0){
		  R.WType23.push_back(atof(row->GetField(4)));
		  w2 =1;
	}else if(strncmp(row->GetField(3),"III",3) == 0 && w2==1){
		  R.WType23.back() += atof(row->GetField(4));
		  w2 =2;
	}
      }	
      delete row;
    };
  };
}

void FillMPMAsymmetries(OctResult_t &M,char query[500],TSQLServer *sql)
{
  TSQLResult *res;
  TSQLRow *row;
  
  res = (TSQLResult*)sql->Query(query);
  cout << "Number of Rows found : " << res->GetRowCount() << endl;
  if(res->GetRowCount() != 0){
    while((row = (TSQLRow*)res->Next())){
//       if((int)M.nRun.size() == 0){
// 	M.nArun.push_back(atof(row->GetField(1)));
// 	M.A_C.push_back(atof(row->GetField(3)));
// 	M.A_Cer.push_back(atof(row->GetField(4)));
//       } else if(checkruns(atoi(row->GetField(1)),M.nRun,(Int_t)M.nRun.size())){
	M.nArun.push_back(atof(row->GetField(1)));
	M.A_C.push_back(atof(row->GetField(3)));
	M.A_Cer.push_back(atof(row->GetField(4)));
   //   }
    }
    delete row;
  }
  delete res;
}

void FillRWPAsymmetries(OctResult_t &M,char query[500],TSQLServer *sql)
{
  TSQLResult *res;
  TSQLRow *row;
  
  res = (TSQLResult*)sql->Query(query);
  cout << "Number of Rows found : " << res->GetRowCount() << endl;
  if(res->GetRowCount() != 0){
    while((row = (TSQLRow*)res->Next())){
      if((int)M.nRun.size() == 0){
	M.nArun.push_back(atof(row->GetField(1)));
	M.A_C.push_back(atof(row->GetField(3)));
	M.A_Cer.push_back(atof(row->GetField(4)));
      } else if(checkruns(atoi(row->GetField(1)),M.nRun,(Int_t)M.nRun.size())){
	M.nArun.push_back(atof(row->GetField(1)));
	M.A_C.push_back(atof(row->GetField(3)));
	M.A_Cer.push_back(atof(row->GetField(4)));
      }
    }
    delete row;
  }
  delete res;
}

void FillRWPQAsymmetries(OctResult_t &M,char query[500],TSQLServer *sql)
{
  TSQLResult *res;
  TSQLRow *row;
  
  res = (TSQLResult*)sql->Query(query);
  cout << "Number of Rows found : " << res->GetRowCount() << endl;
  if(res->GetRowCount() != 0){
    while((row = (TSQLRow*)res->Next())){
      if(atof(row->GetField(1))>0){
	  if((int)M.nRun.size() == 0){
	    M.nArun.push_back(atof(row->GetField(1)));
	    M.A_C.push_back(TMath::Abs(atof(row->GetField(3))));
	    M.A_Cer.push_back(atof(row->GetField(4)));
	  } else if(checkruns(atoi(row->GetField(1)),M.nRun,(Int_t)M.nRun.size())){
	    M.nArun.push_back(atof(row->GetField(1)));
	    M.A_C.push_back(TMath::Abs(atof(row->GetField(3))));
	    M.A_Cer.push_back(atof(row->GetField(4)));
	  }
      }
    }
    delete row;
  }
  delete res;
}

void ReadAnalsyisDB()
{
  char query[500];
  char *db     = Form("mysql://%s/%s",getenv("UCNAANARESDBADDRESS"),getenv("UCNAANADB"));
  char *dbuser = Form("%s",getenv("UCNADBUSER"));
  char *dbpwd  = Form("%s",getenv("UCNADBPASS"));
  // Make a connection to the server at caltech
  TSQLServer *sql = TSQLServer::Connect(db,dbuser,dbpwd);
  // Check if connected
  if(!(sql->IsConnected())){
    while(!(sql->IsConnected())){
      sql->Connect(db,dbuser,dbpwd);
    }
  }

 
  OctResult_t MPMres,RWPres,Diff,RWPQres,RWPSres,MPMQres,MPMSres;
  
  sprintf(query,"select t1.author,t1.side,t1.start_run,t1.event_type,t1.value,t1.err,t2.cut_spec_id from "
	       " analysis_results as t1, cut_spec as t2 where t1.type='Counts' and t1.source='Data' and t1.author='RWP' "
	       " and t1.ana_choice = '' and t1.cut_spec_id = t2.cut_spec_id and t2.radius=50 order by t1.start_run");
  FillRates(RWPres,query,sql);

  
  sprintf(query,"select t1.author,t1.side,t1.start_run,t1.event_type,t1.value,t1.err,t2.cut_spec_id from "
	        " analysis_results as t1, cut_spec as t2 where t1.type='Counts' and t1.source='Data' and t1.author='MPM' "
	        "and t1.ana_choice = 'C' and t1.grouping='pair' and t2.cut_spec_id = t1.cut_spec_id and t2.radius = 50 and t2.energy_min = 220 order by t1.start_run");
  
  FillMPMRates(MPMres,query,sql);
  
  sprintf(query,"select author,start_run,ana_choice,value,err,cut_spec_id from analysis_results where type='Asymmetry' " 
                 "and source='Data' and author='MPM' and gate_valve='Open' and value>0 and end_run-start_run < 200 and "
		 " grouping = 'octet' order by start_run");
  FillMPMAsymmetries(MPMres,query,sql);
  
  sprintf(query,"select author,start_run,ana_choice,value,err,cut_spec_id from analysis_results where type='Asymmetry' " 
                 "and source='Data' and event_type='0,I,II,III' and author='RWP' and value>0 and grouping ='octet' and end_run-start_run < 200 "
		 "and ana_choice='C' order by start_run");
  FillRWPAsymmetries(RWPres,query,sql);
  
  sprintf(query,"select author,start_run,ana_choice,value,err,cut_spec_id from analysis_results where type='Asymmetry' " 
                 "and source='Data' and author='RWP'and event_type='0,I,II,III' and grouping ='quartet' and end_run-start_run < 200 "
		 "and ana_choice='C' order by start_run");
  
  FillRWPQAsymmetries(RWPQres,query,sql);
  sprintf(query,"select author,start_run,ana_choice,value,err,cut_spec_id from analysis_results where type='Asymmetry' " 
                 "and source='Data' and author='RWP'and event_type='0,I,II,III' and grouping ='pair' and end_run-start_run < 200 "
		 "and ana_choice='C' order by start_run");
  FillRWPQAsymmetries(RWPSres,query,sql);
  
  sprintf(query,"select author,start_run,ana_choice,value,err,cut_spec_id from analysis_results where type='Asymmetry' " 
                 "and source='Data' and author='MPM' and gate_valve='Open' and value>0 and end_run-start_run < 200 and "
		 " grouping = 'quartet' order by start_run");
  FillMPMAsymmetries(MPMQres,query,sql);
  
  sprintf(query,"select author,start_run,ana_choice,value,err,cut_spec_id from analysis_results where type='Asymmetry' " 
                 "and source='Data' and author='MPM' and gate_valve='Open' and value>0 and end_run-start_run < 200 and "
		 " grouping = 'pair' order by start_run");
  FillMPMAsymmetries(MPMSres,query,sql);
  
  
  SetDifference(RWPres,MPMres,Diff);
  PrintSizes(RWPres);
  PrintSizes(MPMres);
  PrintSizes(Diff);
  
  PlotRates(MPMres,RWPres,Diff);
  
  TCanvas *cAsym = new TCanvas("cAsym","Asymmetries",500,700);
  cAsym->cd();
  TGraphErrors *gMPMA = new TGraphErrors((int)MPMres.nArun.size(),&MPMres.nArun[0],&MPMres.A_C[0],0,&MPMres.A_Cer[0]);
  TGraphErrors *gMPMAQ = new TGraphErrors((int)MPMQres.nArun.size(),&MPMQres.nArun[0],&MPMQres.A_C[0],0,&MPMQres.A_Cer[0]);
  TGraphErrors *gMPMAS = new TGraphErrors((int)MPMSres.nArun.size(),&MPMSres.nArun[0],&MPMSres.A_C[0],0,&MPMSres.A_Cer[0]);

  TGraphErrors *gRWPA = new TGraphErrors((int)RWPres.nArun.size()-1,&RWPres.nArun[0],&RWPres.A_C[0],0,&RWPres.A_Cer[0]);
  TGraphErrors *gRWPAQ = new TGraphErrors((int)RWPQres.nArun.size()-1,&RWPQres.nArun[0],&RWPQres.A_C[0],0,&RWPQres.A_Cer[0]);
  TGraphErrors *gRWPAS = new TGraphErrors((int)RWPSres.nArun.size()-1,&RWPSres.nArun[0],&RWPSres.A_C[0],0,&RWPSres.A_Cer[0]);
  
  cout << "Sizes " << endl;
  cout << "Octets   : " << (int)MPMres.nArun.size() << "\t" << (int)RWPres.nArun.size() << endl;
  cout << "Quartets : " << (int)MPMQres.nArun.size() << "\t" << (int)RWPQres.nArun.size() << endl;
  cout << "SR-Pair  : " << (int)MPMSres.nArun.size() << "\t" << (int)RWPSres.nArun.size() << endl;
  
  fstream f1;
  f1.open("octet_res.txt",fstream::out);
  
  for(Int_t i = 0 ; i < (int)RWPres.nArun.size() ; i++){
    f1 << RWPres.nArun[i] << "\t" << RWPres.A_C[i];
    if(i<(int)MPMres.nArun.size()) 
      f1 << "\t" << MPMres.nArun[i] << "\t" << MPMres.A_C[i];
    f1<< endl;
  }
  
  f1 << "===========================" << endl;
  
  for(Int_t i = 0 ; i < (int)RWPQres.nArun.size() ; i++){
    f1 << RWPQres.nArun[i] << "\t" << RWPQres.A_C[i];
    if(i<(int)MPMQres.nArun.size()) 
      f1 << "\t" << MPMQres.nArun[i] << "\t" << MPMQres.A_C[i];
    f1<< endl;
  }
  
  f1 << "===========================" << endl;
  
   for(Int_t i = 0 ; i < (int)RWPSres.nArun.size() ; i++){
    f1 << RWPSres.nArun[i] << "\t" << RWPSres.A_C[i];
    if(i<(int)MPMSres.nArun.size()) 
      f1 << "\t" << MPMSres.nArun[i] << "\t" << MPMSres.A_C[i];
    f1<< endl;
  }
  
  f1 << "===========================" << endl;
  f1.close();


  TF1 *fNO = new TF1("fNO","pol0",14000,16250);
  TF1 *fNQ = new TF1("fNQ","pol0",14000,16250);
  TF1 *fNS = new TF1("fNS","pol0",14000,16250);

  TF1 *fCO = new TF1("fC0","pol0",14000,16250);
  TF1 *fCQ = new TF1("fCQ","pol0",14000,16250);
  TF1 *fCS = new TF1("fCS","pol0",14000,16250);
  
   fNO->SetLineColor(2);
   fNO->SetLineStyle(2);
   fNO->SetLineWidth(1);
   fNQ->SetLineColor(2);
   fNQ->SetLineStyle(2);
   fNQ->SetLineWidth(1);
   fNS->SetLineColor(2);
   fNS->SetLineStyle(2);
   fNS->SetLineWidth(1);
   
   fCO->SetLineStyle(2);
   fCO->SetLineWidth(1);
   fCQ->SetLineStyle(2);
   fCQ->SetLineWidth(1);
   fCS->SetLineStyle(2);
   fCS->SetLineWidth(1);
   
  cAsym->Divide(1,3);
  cAsym->cd(1);
  gMPMA->SetMarkerStyle(24);
  gRWPA->SetMarkerStyle(20);
  gRWPA->SetMarkerColor(2);
  gRWPA->SetLineColor(2);
  gRWPA->Draw("ap");
  gRWPA->Fit(fNO,"RMEQ","same");
  gMPMA->Draw("p");
  gMPMA->Fit(fCO,"RMEQ","same");
  
  TLegend *lg1 = new TLegend(0.2,0.60,0.6,0.88);
  lg1->SetLineColor(0);
  lg1->SetFillColor(0);
  lg1->AddEntry(gRWPA,Form("NCSU : A_{raw} = %6.5f #pm %6.5f",
			     fNO->GetParameter(0),fNO->GetParError(0)),"lp");
  lg1->AddEntry(gMPMA,Form("CIT  : A_{raw} = %6.5f #pm %6.5f",
			     fCO->GetParameter(0),fCO->GetParError(0)),"lp");
  lg1->Draw();
  cAsym->cd(2);
  gMPMAQ->SetMarkerStyle(24);
  gRWPAQ->SetMarkerStyle(20);
  gRWPAQ->SetMarkerColor(2);
  gRWPAQ->SetLineColor(2);
  gRWPAQ->Draw("ap");
  gRWPAQ->Fit(fNQ,"RMEQ","same");
  gMPMAQ->Draw("p");
  gMPMAQ->Fit(fCQ,"REMQ","same");
  
  TLegend *lg2 = new TLegend(0.2,0.60,0.6,0.88);
  lg2->SetLineColor(0);
  lg2->SetFillColor(0);
  lg2->AddEntry(gRWPAQ,Form("NCSU : A_{raw} = %6.5f #pm %6.5f",
			     fNQ->GetParameter(0),fNQ->GetParError(0)),"lp");
  lg2->AddEntry(gMPMAQ,Form("CIT  : A_{raw} = %6.5f #pm %6.5f",
			     fCQ->GetParameter(0),fCQ->GetParError(0)),"lp");
  lg2->Draw();
  
  cAsym->cd(3);
  gRWPAS->SetTitle("");
  gMPMAS->SetMarkerStyle(24);
  gRWPAS->SetMarkerStyle(20);
  gRWPAS->SetMarkerColor(2);
  gRWPAS->SetLineColor(2);
  gRWPAS->Draw("ap");
  gMPMAS->Draw("p");
  gRWPAS->Fit(fNS,"RMEQ","same");
  gMPMAS->Fit(fCS,"RMEQ","same");
  
  TLegend *lg3 = new TLegend(0.2,0.60,0.6,0.88);
  lg3->SetLineColor(0);
  lg3->SetFillColor(0);
  lg3->AddEntry(gRWPAS,Form("NCSU : A_{raw} = %6.5f #pm %6.5f",
			     fNS->GetParameter(0),fNS->GetParError(0)),"lp");
  lg3->AddEntry(gMPMAS,Form("CIT  : A_{raw} = %6.5f #pm %6.5f",
			     fCS->GetParameter(0),fCS->GetParError(0)),"lp");
  lg3->Draw();

};