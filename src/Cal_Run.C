#include "Cal_Run.h"

using namespace std;

//ClassImp(Cal_Run);

//----------------------------------------------------------------------------------------//
Cal_Run::Cal_Run(Int_t n,Int_t m,Int_t b):Run(n,m,b)
{

  runnum  = n;
  runtype = m;
  const char* filepath = getenv("UCNAOUTPUTDIR"); 
  
  f1 = new TFile(Form("%s/hists/spec_%d.root",filepath,n),"UPDATE");
 
  if(!(f1->IsOpen()))
        cout << " Failed to open file "<<endl;
  t1 = (TTree*)f1->Get("phys");

  if(!(f1->cd("Analysis")))
        f1->mkdir("Analysis");

  gROOT->cd();
  dHists = new TDirectory("dHists",Form("Histograms_%d",n));
  dHists->cd();
  dHists->Append(t1);

  SetBranches();
  if(GetBackGround(sql) == 0) cout << "failure" << endl;

  Load_Rotation_Matrix(GetGeo()); // Rotate the detector corrdinates
  Load_Background(runnum,runtype,GetGeo());

/*
  f1 = new TFile(Form("../reduced_trees/cal_%4i.root",n),"READ");
  if(f1->IsOpen()){
    cout << "File cal_"<<n<<".root has been successfully openned" << endl;
  } else {
    cout << "Failed to open file"<<endl;
  }
   // set run number
  t1 = (TTree*)f1->Get("tr");
  t1->SetBranchAddress("prout",&chb);*/

}
//----------------------------------------------------------------------------------------//
Cal_Run::~Cal_Run()
{
 
}
//----------------------------------------------------------------------------------------//
Int_t Cal_Run::get_peak_pos(peak_t *fred,TH1D* h1,Int_t &npeaks)
{
  // -- Find X coordinate of peaks ----------
  TSpectrum *spos = new TSpectrum();
  TF1 *fg1 = new TF1("fg1","gaus",-80,80);
  Int_t peak_ctr  = 0;
  
  // Use TSpectrum to search the for the peaks.
  npeaks = spos->Search(h1,2,"",0.15);

  //extract peaks
  Float_t *xpeaks = spos->GetPositionX();

  cout << "TSpectrum found : " << npeaks << " peaks." << endl;  
    
  // Fit all peaks
  while(peak_ctr < npeaks){
    fg1->FixParameter(1,xpeaks[peak_ctr]);
    h1->Fit("fg1","RQME","",xpeaks[peak_ctr]-20.,xpeaks[peak_ctr]+20.);
    fred[peak_ctr].sig = fg1->GetParameter(2);
    cout << "Peak " << peak_ctr+1 << " is at ";
    cout << xpeaks[peak_ctr] << " and width "<< fg1->GetParameter(2) << endl;
    fred[peak_ctr].x   = xpeaks[peak_ctr];
    fred[peak_ctr].sig = (fg1->GetParameter(2) < 2.) ? fg1->GetParameter(2) : 2.;
    peak_ctr++;
  } 
  
  delete spos;
  delete fg1;
  
  return 0;
}
//----------------------------------------------------------------------------------------//
Int_t Cal_Run::Draw_2d(Int_t nr,Int_t n)
{
  
  hpw = new TH2F(Form("cr1[%i]->hpw",n),"West Position;X;Y",160,-80,80,160,-80,80);
  hpe = new TH2F(Form("cr1[%i]->hpe",n),"East Position;X;Y",160,-80,80,160,-80,80);
  
  // Draw 2-d Position Plots of the sources
  
  t1->Draw(Form("wypos:wxpos >> cr1[%i]->hpw",n),
	   Form("wbvetotdc<200 && etdc < 200 && wymulti>0 && westa>%f && wxmulti>0 && run == %i",
		w_acut,nr)
	   ,"goff");
   
  t1->Draw(Form("eypos:expos >> cr1[%i]->hpe",n),
	   Form("eymulti > 0 && ebvetotdc < 200 && wtdc < 200 &&  easta>%f && exmulti>0 && run == %i",
		e_acut,nr)
	   ,"goff");
  return 0;
}
//----------------------------------------------------------------------------------------//
Int_t Cal_Run::Fill(Int_t n)
{
  Int_t peakl=0;

  e_acut = 30.;
  w_acut = 30.;

  for(Int_t i = 0 ; i < t1->GetEntries() ; i++){
    if(i%100000 == 0)
        cout << "At event " << i << " of " << t1->GetEntries() << endl;
    
    t1->GetEntry(i);  // Get Enetry from the Tree;
    
    if(chb.ebvetotdc < 200 && chb.etoptdc < 200 && chb.wbvetotdc < 200  ){
      if(chb.easta > e_acut && chb.westa < w_acut){
        while(!Cut_on_peak(chb.expos,chb.eypos,expeaks[peakl].sig,eypeaks[peakl].sig,
			   expeaks[peakl].x,eypeaks[peakl].x)){
	 peakl++;
	 if(peakl >= nex){
	   peakl = -1;
	   break;
    	}
       }
        if(peakl >= 0){
	 heq[peakl]->Fill(chb.easte);
	 he1[peakl]->Fill(chb.e1);
	 he2[peakl]->Fill(chb.e2);
	 he3[peakl]->Fill(chb.e3);
	 he4[peakl]->Fill(chb.e4);
        }
      } else if(chb.easta < e_acut && chb.westa > w_acut){
        while(!Cut_on_peak(chb.wxpos,chb.wypos,wxpeaks[peakl].sig,wypeaks[peakl].sig,
			 wxpeaks[peakl].x,wypeaks[peakl].x)){
	 peakl++; 
	 if(peakl >= nwx){
	    peakl = -1;
	   break;
	 }
        }
        if(peakl >= 0){
	 hwq[peakl]->Fill(chb.weste);
	 hw1[peakl]->Fill(chb.w1);
	 hw2[peakl]->Fill(chb.w2);
	 hw3[peakl]->Fill(chb.w3);
	 hw4[peakl]->Fill(chb.w4);
        }
      } else if(chb.easta > e_acut && chb.westa > w_acut){
        while(!(Cut_on_peak(chb.expos,chb.eypos,expeaks[peakl].sig,eypeaks[peakl].sig,
			 expeaks[peakl].x,eypeaks[peakl].x)) &&
	    !(Cut_on_peak(chb.wxpos,chb.wypos,wxpeaks[peakl].sig,wypeaks[peakl].sig,
			 wxpeaks[peakl].x,wypeaks[peakl].x))) {
	peakl++;
	if(peakl >= nex){
	  peakl = -1;
	  break;
	}
      }
      if(peakl >= 0){
	htw[peakl]->Fill(chb.wtdc);
	hte[peakl]->Fill(chb.etdc);
	if(chb.easte>0 && chb.weste>0){
	 if(chb.etdc > 200 && chb.etdc < 2900){
	    // East Type 1 backscatters
	    hEtype_1[peakl]->Fill(chb.easte + chb.weste);
	    heq[peakl]->Fill(chb.easte      + chb.weste);
       	 }else if(chb.wtdc > 200 && chb.wtdc < 2900){
	   // West Type 1 backscatters
	    hWtype_1[peakl]->Fill(chb.easte + chb.weste);
	    hwq[peakl]->Fill(chb.easte      + chb.weste);
	 }
	} else if(chb.easte > 0 && chb.weste <= 0){
	  hEtype_23[peakl]->Fill(chb.easte);
	  heq[peakl]->Fill(chb.easte);
	} else if(chb.weste > 0 && chb.easte <= 0){
	  hWtype_23[peakl]->Fill(chb.weste);
	  hwq[peakl]->Fill(chb.weste);
	}
      }
    }
    peakl=0;
   }
  }
return 0;
}
//----------------------------------------------------------------------------------------//
Int_t Cal_Run::Fit_Response(Int_t n)
{
  // Currently this function performs a simple fit to a Gaussain around the mean of the
  // the response peak.  In the future this should be expanded to include the semi-analytic
  // model or a Gamma Distribution Ask Leah.

  // Define Gaussain Fitting Function
  TF1 *fg = new TF1("fg","gaus",1000,10000);
  fg->SetLineWidth(2);
  
  //Open a holder canvas
  TCanvas *c9 = new TCanvas();
  
  //Define Spectrum object and holder variables
  
  TSpectrum *spos = new TSpectrum();
  Int_t npeaks,prvp;
  Float_t meanh[10],sigh[10],err[10];
  Int_t totpeaks=0;
  Float_t *peaks;

  prvp = 0;
  
  for(Int_t i = 0 ; i < nex; i++){
    npeaks    = spos->Search(heq[i],6,"",0.2);
    cout << npeaks << " Found for position peak " << i << endl;
    totpeaks += npeaks;

    peaks = (Float_t*)malloc(npeaks*sizeof(Float_t));
    peaks = spos->GetPositionX();

    for(int j = prvp ; j < totpeaks ; j++){
      heq[i]->Fit("fg","QRME","",
		  (Double_t)peaks[j-prvp] - 500.,
		  (Double_t)peaks[j-prvp] + 500.);

      meanh[j] = fg->GetParameter(1);
      sigh[j]  = fg->GetParameter(2);
      err[j]   = fg->GetParError(1);
    }
    prvp += npeaks;
  }

  ecal = (calp_t*)malloc(totpeaks*sizeof(calp_t));
  
  for(int i = 0 ; i < totpeaks ; i++){
    ecal[i].adc = meanh[i];
    ecal[i].sig = sigh[i];
    ecal[i].err = err[i];
  }
  nep = totpeaks;

  totpeaks = 0;
  prvp = 0;
  
  for(Int_t i = 0 ; i < nwx; i++){
    npeaks    = spos->Search(hwq[i],6,"",0.2);
    cout << npeaks << " Found for position peak " << i << endl;
    totpeaks += npeaks;
    
    peaks = (Float_t*)malloc(npeaks*sizeof(Float_t));
    peaks = spos->GetPositionX();

    for(int j = prvp ; j < totpeaks ; j++){
      hwq[i]->Fit("fg","QRME","",
		  (Double_t)peaks[j-prvp]-500.,
		  (Double_t)peaks[j-prvp]+500.);

      meanh[j] = fg->GetParameter(1);
      sigh[j]  = fg->GetParameter(2);
      err[j]   = fg->GetParError(1);
    }
    prvp += npeaks;
  }

  wcal = (calp_t*)malloc(totpeaks*sizeof(calp_t));

  for(int i = 0 ; i < totpeaks ; i++){
    wcal[i].adc = meanh[i];
    wcal[i].sig = sigh[i];
    wcal[i].err = err[i];
  }
  nwp = totpeaks;

  // Delete fitting function, Tspectrum object and the holder canvas
  delete fg;
  delete spos;
  delete c9;

  return 0;
}
//----------------------------------------------------------------------------------------//
Int_t Cal_Run::FillGraph(Int_t n)
{
  
  Double_t yepeaks[nep],ywpeaks[nwp],xe[nep],xw[nwp],xer[nep],xwr[nwp];
  Double_t sek[nep],swk[nwp],sekr[nep],swkr[nep];
  Double_t yepeaksr[nep],ywpeaksr[nwp];
  
  Double_t evis_tin,evis_bi1,evis_bi2,evis_cd;
  
// Set the Evis scale for the calibration sources in the various geometries...
  
  if(GetRunNumber() < 9200){
    evis_tin = 335.;
    evis_bi1 = 462.;
    evis_bi2 = 965.;
    evis_cd  = 35.;
  } else if(GetRunNumber() > 9200 && GetRunNumber()< 10400) { 
    evis_tin = 335.;
    evis_bi1 = 462.;
    evis_bi2 = 965.;
    evis_cd  = 35.;
  } else if(GetRunNumber() > 10400) { 
    evis_tin = 335.;
    evis_bi1 = 462.;
    evis_bi2 = 965.;
    evis_cd  = 35.;
  }
  
  if(nep == 1){ 
      xe[0] = evis_tin;
  } else if( nep == 3){
      xe[0] = evis_tin;   
      xe[1] = evis_bi1;
      xe[2] = evis_bi2;
  } else if( nep == 4){
      xe[0] = evis_cd;
      xe[1] = evis_tin;
      xe[2] = evis_bi1;
      xe[3] = evis_bi2;
  }

   if(nwp == 1){
    xw[0] = evis_tin;
  } else if( nwp == 3){
    xw[0] = evis_tin;
    xw[1] = evis_bi1;
    xw[2] = evis_bi2;
  } else if( nwp == 4){
    xw[0] = evis_cd;
    xw[1] = evis_tin;
    xw[2] = evis_bi1;
    xw[3] = evis_bi2;
  }
  
  for(Int_t i = 0; i < nep ; i++){
    xer[i] = 0.1;
    yepeaks[i] = (Double_t)ecal[i].adc;
    yepeaksr[i] = (Double_t)ecal[i].err; 
  }

  Int_t ct = 0;

  for(Int_t i = 0;i<nep;i++){
    ct = 0;
    for(Int_t j = 0; j < nep ; j++){
      if(yepeaks[i] > yepeaks[j])ct++;
    }
    sek[ct] = yepeaks[i];
    sekr[ct] = yepeaksr[i];
  }
    
  for(Int_t i = 0; i < nwp ; i++){
    xwr[i] = 0.1;
    ywpeaks[i] = (Double_t)wcal[i].adc;
    ywpeaksr[i] = (Double_t)wcal[i].err;
  }

  for(Int_t i = 0;i<nwp;i++){
    ct = 0;
    for(Int_t j = 0; j < nwp ; j++){
      if(ywpeaks[i] > ywpeaks[j])ct++;
    }
    swk[ct] = ywpeaks[i];
    swkr[ct] = ywpeaksr[i];
  }

  gcale = new TGraphErrors(nep,sek,xe,sekr,xer);
  gcalw = new TGraphErrors(nwp,swk,xw,swkr,xwr);
  
  return 0;
}
//----------------------------------------------------------------------------------------//
Int_t Cal_Run::Locate_peaks(Int_t n,Int_t nrun)
{
 // Project 2d plot in X to find the x-pos of the sources
  TF1 *fg = new TF1("fg","gaus",-80,80);

  TH1D* hecx = hpe->ProjectionX();
  TH1D* hwcx = hpw->ProjectionX();
 
  // Initialize peak holding structs ----------------------------------
  
  expeaks = (peak_t*)malloc(5.*sizeof(peak_t));
  wxpeaks = (peak_t*)malloc(5.*sizeof(peak_t));
  
  get_peak_pos(expeaks,hecx,nex);
  get_peak_pos(wxpeaks,hwcx,nwx);
  
  eypeaks = (peak_t*)malloc(nex*sizeof(peak_t));
  wypeaks = (peak_t*)malloc(nwx*sizeof(peak_t));

  //------------------------------------------------------------------
  // Find y-pos of peaks

  for(Int_t i = 0 ; i < nex ; i++){
    cout << "Analyzing Peak " << i << endl;
    hecy[i] = new TH1F(Form("cr1[%i]->hecy[%i]",n,i),
			    Form("cr1[%i]->hecy[%i]",n,i),160,-80,80);
    
    t1->Draw(Form("eypos >> cr1[%i]->hecy[%i]",n,i),
	     Form("ebvetotdc <200 && etoptdc < 200 && expos > %f && expos < %f && easta > 40 && run == %i",
		  expeaks[i].x - expeaks[i].sig,
		  expeaks[i].x + expeaks[i].sig
		  ,nrun),"goff");
    
    hecy[i]->Fit("fg","QRME","goff");
    
    eypeaks[i].x   = fg->GetParameter(1);
    eypeaks[i].sig = ( fg->GetParameter(2) < 5. ) ? fg->GetParameter(2) : 5.;
    
  } 
  
  for(int i= 0 ;i < nwx ; i++){
    
    cout << "Analyzing Peak " << i << endl;
    
    hwcy[i] = new TH1F(Form("cr1[%i]->hwcy[%i]",n,i),
			    Form("cr1[%i]->hwcy[%i]",n,i),160,-80,80);
    
    t1->Draw(Form("wypos >> cr1[%i]->hwcy[%i]",n,i),
	     Form("wbvetotdc<200 && wxpos>%f && wxpos<%f && westa>30 && run == %i",
		  wxpeaks[i].x - wxpeaks[i].sig,
		  wxpeaks[i].x + wxpeaks[i].sig
		  ,nrun),"");
    
    hwcy[i]->Fit("fg","RMEQ","goff");
    
    wypeaks[i].x   = fg->GetParameter(1);
    wypeaks[i].sig = ( fg->GetParameter(2) > 5. ) ? 2. : fg->GetParameter(2);
    
  }  
  
  //----------------------------------------------------------------------+
  // Now apply those peaks that have been found to the position cuts for  |
  // creation of the Qadc Histograms                                      |
  //----------------------------------------------------------------------+
  nsig = 3.;
  
  for(int i = 0 ; i< nex ;i++){
    cout << expeaks[i].x << "\t " << nsig*expeaks[i].sig << endl;
    epw[i] = new TEllipse(expeaks[i].x,
			  eypeaks[i].x,
			  nsig*expeaks[i].sig,
			  nsig*eypeaks[i].sig);
    epw[i]->SetFillStyle(0);
  }
  for(int i = 0 ; i< nwx ;i++){
    wpw[i] = new TEllipse(wxpeaks[i].x,
			  wypeaks[i].x,
			  nsig*wxpeaks[i].sig,
			  nsig*wypeaks[i].sig);
    wpw[i]->SetFillStyle(0);
  }


  delete fg;
  
  return 0;
}
//----------------------------------------------------------------------------------------//
Int_t Cal_Run::Draw_Hists(Int_t n)
{
  for(Int_t i = 0 ; i < nex ; i++){
    c1->cd(3);
    
    if(i==0){
      heq[i]->Draw();
      gPad->SetLogy();
    } else {
      heq[i]->Draw("same");
    }
  /*  
    c1->cd(8);
    if(i==0){
      gPad->SetLogy();
      hte[i]->Draw();
    } else {
      hte[i]->Draw("same");
    }
    c1->cd(6);
    if(i==0){
      hEtype_1[i]->Draw();
    } else {
      hEtype_1[i]->Draw("same");
    }*/
  }

  for(int i = 0;i< nwx ; i++){
    c1->cd(4);
    if(i==0){
      // gPad->SetLogy();
      hwq[i]->Draw();
    } else {
      hwq[i]->Draw("same");
    }
/*
    c1->cd(7);
    if(i==0){
      gPad->SetLogy();
      htw[i]->Draw();
    } else {
      htw[i]->Draw("same");
    }
    c1->cd(5);
    if(i==0){
      hWtype_1[i]->Draw();
    } else {
      hWtype_1[i]->Draw("same");
    }*/
  }
  return 0;
}
//----------------------------------------------------------------------------------------//
Bool_t Cal_Run::Cut_on_peak(Float_t x,Float_t y,Float_t sig1,Float_t sig2,Float_t xp,Float_t yp)
{
  using namespace TMath;

  if(Power(x-xp,2)/Power(sig1,2) + Power(y-yp,2)/Power(sig2,2) < 3)
    return kTRUE;
  else 
    return kFALSE;
}

Int_t Cal_Run::Draw_Pos(Int_t nex, Int_t nwx)
{
  c1->cd(2);
  hpw->Draw("colz");

  for(Int_t j = 0; j < nwx ; j++){
    wpw[j]->SetLineColor(2+j);
    wpw[j]->SetLineWidth(2);
    wpw[j]->Draw();
  }

  c1->cd(1);
  hpe->Draw("colz");
  for(Int_t j = 0 ; j < nex ; j++){
    epw[j]->SetLineColor(2+j);
    epw[j]->SetLineWidth(2);
    epw[j]->Draw();
  }
  return 0;
  
}
//-----------------------------------------------------------------------------------------------//
Int_t Cal_Run::SubBck(Bck_Run *br)
{

  for(Int_t i = 0 ; i < nex ;i++){
    heq[i]->Add(br->heq[0],-1);
  }

  for(Int_t i = 0 ; i< nwx ; i++){
    hwq[i]->Add(br->hwq[0],-1);
  }
   return 0;
}
