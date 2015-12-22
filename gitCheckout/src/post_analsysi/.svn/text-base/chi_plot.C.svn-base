#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <fstream>
#include <iostream>

#define nlow 10
#define nhigh 27

Double_t ChiSquareDis(Double_t *x,Double_t *par);

void dostuff()
{
   Double_t col,chi;
   
   fstream f1;
   
   TH1F *hChi = new TH1F("hChi","#chi^{2} Distribution ; #chi^{2} ; Counts " , 10,0,40);
   
   f1.open("chi_squares_0.txt",fstream::in);
   
   do{
     f1 >> col >> chi;
     hChi->Fill(chi,1);
   }while(!(f1.eof()));
   
   f1.close();
   f1.open("chi_squares_1.txt",fstream::in);
   
   do{
     f1 >> col >> chi;
     hChi->Fill(chi,1);
   }while(!(f1.eof()));
   
   f1.close();
   f1.open("chi_squares_2.txt",fstream::in);
   
   do{
     f1 >> col >> chi;
     hChi->Fill(chi,1);
   }while(!(f1.eof()));
   
   f1.close();
   f1.open("chi_squares_3.txt",fstream::in);
   
   do{
     f1 >> col >> chi;
     hChi->Fill(chi,1);
   }while(!(f1.eof()));
   
   f1.close();
   
   
   TCanvas *c1 =  new TCanvas("c1","chi");
   c1->cd();
   hChi->GetXaxis()->CenterTitle();
   hChi->GetYaxis()->CenterTitle();
   hChi->SetMarkerStyle(20);
   hChi->Draw("E1X0");
   TF1 *fchi = new TF1("fchi",ChiSquareDis,0,40.,2);
   fchi->SetParameter(0,4.*hChi->Integral());
   fchi->SetParameter(1,17);
   fchi->SetLineColor(2);
   fchi->Draw("same");
   
   c1->Print("chi_data.eps");
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