#include "Bck_Run.h"

using namespace std;

Bck_Run::Bck_Run(Int_t n,Int_t m,TSQLServer *sql):Run(n,m,sql)
{
  Find_TDC_Cut(n,sql);
  Load_Rotation_Matrix(GetGeo()); // Rotate the detector corrdinates
  Load_Background(GetRunNumber(),GetRunType(),GetGeo()); // Load Parent Background directory 
}

Bck_Run::~Bck_Run()
{

}

Int_t Bck_Run::Draw_2d(Int_t nr,Int_t n)
{
  
  hpw = new TH2F(Form("bckr[%i]->hpw",n),"West Position;X;Y",150,-75,75,150,-75,75);
  hpe = new TH2F(Form("bckr[%i]->hpe",n),"East Position;X;Y",150,-75,75,150,-75,75);
  
  // Draw 2-d Position Plots of the sources
  
  return 0;

}

void Bck_Run::Load_Histograms(Int_t run=0)
{
     run=this->GetRunNumber();
     //cout << "Load_Histograms Run " << run << endl;
     this->GetHistograms(run);
     this->ScaleList(this->HEastAn,this->rtime_e);
     this->ScaleList(this->HWestAn,this->rtime_w);
}

void Bck_Run::Remove_Histograms()
{
    this->DeleteHistos();
}

//--------------------------------------------------------------------------------------------------------------------
Int_t Bck_Run::Fill(Int_t n,Int_t remake,Double_t *sep,Int_t nrun)
{

  Float_t East_Time_Cnt=0.;
  Float_t West_Time_Cnt=0.;

  Float_t OldTE = 0.;
  Float_t OldTW = 0.;  

  for(Int_t i = 0 ; i < 100 ; i ++){
    seppar[i] = sep[i];
  }

  //cout << "run num = " << GetRunNumber() << endl;
  if(remake == 1 || !AnalysisDirExist){
    Initialize_hist(0,1,1);
    cout << "Reading " << Form("%s/hists/spec_%d.root",getenv("UCNAOUTPUTDIR"),nrun) << endl;
    TFile *f2 = new TFile(Form("%s/hists/spec_%d.root",getenv("UCNAOUTPUTDIR"),nrun),"READ"); 
    hmrIn = (TH1F*) f2->Get("UCN_Mon_4_Rate");
    for(Int_t MRbin=0; MRbin<hmrIn->GetNbinsX(); MRbin++){
	hmr1->Fill(hmrIn->GetBinCenter(MRbin),hmrIn->GetBinContent(MRbin));
    }
    delete hmrIn; 
    f2->Close();
    //cout << "run num = " << GetRunNumber() << endl;
    for(Int_t i = 0 ; i < t1->GetEntries() ; i++){

      t1->GetEntry(i);  // Get Enetry from the Tree;
      //cout << i << "\t" << "run num = " << GetRunNumber() << endl;
      ScaleTDC();
      // Stupid live time counter.  Live times are up to date in the mysql database so 
      // these are being used to up date the clocks.
      if(TofE < OldTE)East_Time_Cnt++;
      OldTE = TofE;
      
      if(TofW < OldTW)West_Time_Cnt++;
      OldTW = TofW; 

      if(Check_Vetos()){
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
      rtime_e = rtime_e*(hGammaCountsg->GetBinContent(1))/(hGammaCounts->GetBinContent(1)); 
      rtime_w = rtime_w*(hGammaCountsg->GetBinContent(2))/(hGammaCounts->GetBinContent(2)); 
    }
    SaveHistograms(kTRUE);
  } else if(remake == 0){
    //GetHistograms(0);
    //cout << "Initializing Open Run num " << GetRunNumber() << endl;
    GetHistograms(GetRunNumber());
    if(GetGeo() ==1){    
      rtime_e = rtime_e*(hGammaCountsg->GetBinContent(1))/(hGammaCounts->GetBinContent(1)); 
      rtime_w = rtime_w*(hGammaCountsg->GetBinContent(2))/(hGammaCounts->GetBinContent(2)); 
    }
    SaveHistograms(kFALSE);

  }

  return 0;

}

Int_t Bck_Run::Draw_Hists(Int_t n)
{
    c1->cd(1);
    hpe->Draw("colz");
    
    c1->cd(2);
    hpw->Draw("colz");
    
    
    c1->cd(3);
    heq->Draw();
    
    c2->cd(1);
    gPad->SetLogy();
    hte->Draw();
    
    c2->cd(2);
    hEtype_1->Draw();

    c1->cd(4);
    hwq->Draw();

    c2->cd(4);
    gPad->SetLogy();
    htw->Draw();
   
    c2->cd(3);
    hWtype_1->Draw();
  
    delete c1;
    delete c2;
    
   return 0;
}

Bool_t Bck_Run::Check_Vetos()
{
  return kTRUE;
  
  if(TaggedTopE || TaggedBackE || TaggedDriftE ||  TaggedBackW || TaggedDriftW ) return kFALSE;
  
  return kTRUE; 
}

