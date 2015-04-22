#ifndef __MathFunc_C__
#define __MathFunc_C__

#include "analysisMathFunc.h"
 
Double_t ChiSquareDis(Double_t *x,Double_t *par)
{
  
  using namespace TMath;
  
  Double_t xx  = x[0];   // chi^2
  Double_t nd  = par[1]; // degrees of freedom
  Double_t scl = par[0]; // scale
  
  Double_t chip = scl*(Power(xx,nd/2.-1)*exp(-xx/2.)/(Gamma(nd/2.)*Power(2,nd/2.)));
  return chip;
  
}

//-----------------------------------------------------------------------------------------
void CollectAllRot(Double_t x[][150],Int_t xbin,Int_t ybin,TH2F *h,TH2F *h2,Int_t last)
{
  
  for(Int_t i = 1 ; i < xbin + 1 ; i++){
    for(Int_t j = 1 ; j < ybin + 1 ; j++){
      x[i-1][j-1] += h->GetBinContent(i,j);
    }
  }

  if(last == 1){
    for(Int_t ibin = 1 ; ibin < xbin+1 ; ibin++){
      for(Int_t jbin = 1 ; jbin < ybin+1 ; jbin++){
	h2->SetBinContent(ibin,jbin,x[ibin-1][jbin-1]);
      }
    }
  }
  
  return;
}
//--------------------------------------------------------------
void CleanArray(Double_t x[][150],Int_t xbin,Int_t ybin)
{
  

  for(Int_t i = 0 ; i < xbin ; i++){
    for(Int_t j = 0 ; j < ybin ; j++){
      x[i][j] = 0.;
    }
  }

  return;
}
//----------------------------------------------------------------
void Average_Array_Vector(Asym_t A,Int_t n, Double_t &y_Average,Double_t &y_Average_er)
{
  // 
  //   y_Average    = Sum_0^n( x[i]i/ xer[i]^2) / Sum_0^n(1./xer[i]^2)
  //   y_Average_er = Sqrt(Sum(1./xer[i]^2))
  //

  Double_t ytemp   = 0.;
  Double_t ytemper = 0.;

  for(Int_t j = 0 ; j < n ; j++){
    Double_t x   = A.A_ave[j];
    Double_t xer = A.A_error[j];
    if(!(TMath::IsNaN(x)) && xer > 0 && xer != 1){
          ytemp    += x / TMath::Power(xer,2);
          ytemper  += 1./TMath::Power(xer,2);
    }
  }

  y_Average    = ytemp/ytemper;  // Error weighted Average of the Array passed
  y_Average_er = sqrt(1./ytemper);

  return;
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
//----------------------------------------------------------------------------
Double_t Bin_Errors(Double_t r1,Double_t t1,Double_t r2,Double_t t2)
{
  using namespace TMath;
  using namespace std;

  Double_t er = 0.;
  if(IsNaN(r1) != 0 ) r1 = 0.;
  
  er = Sqrt(Abs( (r1 + r2) + r2));
  
  if(IsNaN(er) != 0){
    cout << "Error calculation is bad \t " << r1 << "\t" << r2 <<endl;
    er = 0.;
  }
  return er;
}

//---------------------------------------------------------------------------
Double_t Rate_Error(Float_t r1,Float_t r2,Float_t t1,Float_t b1,Float_t b2)
{
  using namespace TMath;
  
  Double_t er;
//  std::cout << r1 << " " << b1 << "  " << r2 << "  " << b2 << std::endl;
  er  = (r1*t1 + b1)/Power((r1*t1 - b1),2);
  er += (r2*t1 + b2)/Power((r2*t1 - b2),2);
     
  er = sqrt(Abs(er));
  
 return er;
}
//----------------------------------------------------------------------------------------
void Average_All_Hists(TH1F *h1,std::vector<Double_t> &ave,std::vector<Double_t> &avee)
{
   using namespace TMath;
  
   for(Int_t nbin = 1; nbin <= h1->GetNbinsX() ; nbin++){
        if(h1->GetBinError(nbin) > 0){
		ave[nbin-1]  += h1->GetBinContent(nbin) / Power(h1->GetBinError(nbin),2);
        	avee[nbin-1] += 1./Power(h1->GetBinError(nbin),2);
	}
   }	

}
//---------------------------------------------------------------------------
void Return_Asymmetry(std::vector<Double_t> &ave,std::vector<Double_t> &avee,Int_t bins)
{

   TF1 *fb = new TF1("fb",PoverE,0,2000,1);
   fb->SetParameter(0,0.5);
  
   for(Int_t i = 0 ; i < bins ; i++){
       if(avee[i] > 0){
           ave[i] = ave[i] / (avee[i] * /*fb->Eval((double)(i+1)*10));*/fb->Integral((double)i*10,(double)(i+1)*10)/10.);
           avee[i] = 1./(sqrt(avee[i]) * /*fb->Eval((double)(i+1)*10));*/fb->Integral((double)i*10,(double)(i+1)*10)/10.);
       } else {
           ave[i]  = 0.;
           avee[i] = 0.1;
       }
   }
   delete fb;
};

Double_t PoverE(Double_t *x,Double_t *par)
{
  // simple p/e = \beta function

  Double_t MassE = 510.998910;
  return par[0]*sqrt(x[0]*(x[0]+2.*MassE))/(MassE+x[0]);

};

Double_t BetaSpec(Double_t *x,Double_t *par)
{
  using namespace TMath;
  // the input x is the kinetic energy of the electron
  //-----------------------------------------------------
  //Values from 2010 PDG---------------------------------
  Double_t MassE = 510.998910;
  Double_t Q     = 1293.33214;
  Double_t pe    = Sqrt(x[0]*(x[0]+2.*MassE));  
  Double_t perm = 8.85418e-12; // 
  Double_t alpha = Qe()*Qe()/(Hbar()*C()*4.*Pi()*perm); // fine structure constant
  Double_t R = 0.8768e-15; // Proton Raduis
  Double_t W = x[0]/MassE + 1.; // Total energy in terms of electron rest mass;
  Double_t p = sqrt(W*W-1.); // momentum 
  Double_t pi = Pi();
  Double_t gammae = 0.577215; // Euler's constant
  Double_t a0 = 1.;
  Double_t a1 = pi*W/p;
  Double_t a2 = 11./4. - gammae - log(2.*p*R) + pi*pi*(W/p)*(W/p);
  Double_t a3 = pi*(W/p)*(11./4. - gammae - log(2.*p*R));
  Double_t f = a0 + a1*alpha + a2*alpha*alpha + a3*alpha*alpha*alpha;
  if( IsNaN(f) == 1) f = 1.;
  Double_t b_fierz = 0.;//par[0];
  Double_t Fierz = (1. + b_fierz*(MassE/(MassE+x[0])));
 
  return Power(Q-(x[0]+MassE),2)*(x[0]+MassE)*pe*f*Fierz;
};
//------------------------------------------------------------------------------
Double_t GetSuperRatioAsymmetry(TH1F *R1M, TH1F *R1P,TH1F *R2M, TH1F *R2P,Int_t n1,Int_t n2)
{

  // Get the Super Ratio Asymmetry straight forward..

  using namespace TMath;

  Double_t R = R1M->Integral(n1,n2)*R2P->Integral(n1,n2) /
    (R1P->Integral(n1,n2)*R2M->Integral(n1,n2));
  Double_t S = Abs((1. - sqrt(R)) / ( 1. + sqrt(R)));

  return S;
};
//------------------------------------------------------------------------------
Double_t GetSuperRatio(TH1F *R1M, TH1F *R1P,TH1F *R2M, TH1F *R2P,Int_t n1,Int_t n2)
{

  // Get the Super Ratio Asymmetry straight forward..

  using namespace TMath;

  Double_t R = R1M->Integral(n1,n2)*R2P->Integral(n1,n2) /
    (R1P->Integral(n1,n2)*R2M->Integral(n1,n2));
  
  return R;
};
//-----------------------------------------------------------------------------
Float_t GetSuperRatioEightAss(Float_t R1M, Float_t R1MI, Float_t R1P, Float_t R1PI, Float_t R2M, Float_t R2MI, Float_t R2P, Float_t R2PI)
{

  // Get the Super Ratio Assymmetry straight forward.. with Double values instead of from 
  // histograms

  using namespace TMath;

  Float_t R = ((R1M - R1MI)*(R2P - R2PI) / ((R1P - R1PI)*(R2M - R2MI)));
  Float_t S = (1. - sqrt(R)) / ( 1. + sqrt(R));

  return S;
};
//------------------------------------------------------------------------------
Double_t GetSuperRatio(Double_t R1M, Double_t R1P,Double_t R2M, Double_t R2P)
{

  // Get the Super Ratio - 1 straight forward.. with Double values instead of from 
  // histograms

  using namespace TMath;

  Double_t R = (R1M*R2P / (R1P*R2M)) - 1.;

  return R;
};
//------------------------------------------------------------------------------
Float_t GetSuperRatioEight(Float_t R1M, Float_t R1MI, Float_t R1P, Float_t R1PI, Float_t R2M, Float_t R2MI, Float_t R2P, Float_t R2PI)
{

  // Get the Super Ratio - 1 straight forward.. with Double values instead of from
  // histograms

  using namespace TMath;

  Float_t R = ((R1M - R1MI)*(R2P - R2PI) / ((R1P - R1PI)*(R2M - R2MI))) - 1.;
 
  return R;
};
//-----------------------------------------------------------------------------
Double_t GetSuperRatioAsymmetry(Double_t R1M, Double_t R1P, Double_t R2M, Double_t R2P)
{

  // Get the Super Ratio Assymmetry straight forward.. with Double values instead of from
  // histograms

  using namespace TMath;

  Double_t R = (R1M*R2P / (R1P*R2M));
  Double_t S = (1. - sqrt(R)) / ( 1. + sqrt(R));

  return S;
};
//-----------------------------------------------------------------------------
Double_t GetSuperRatioError(TH1F *R1M, TH1F *R1P,TH1F *R2M, TH1F *R2P,Double_t time1,Double_t time2,Int_t n1,Int_t n2)
{

  // ----------------------------------------------------------------------------+
  // Calculate the Error in the Super Ratio Asymmetry which will be inflated     |
  // once the background numbers are added in.  For now I'm just trying to get   |
  // a number.                                                                   |
  //-----------------------------------------------------------------------------+
  using namespace TMath;

  Double_t R = R1M->Integral(n1,n2)*R2P->Integral(n1,n2) /
    (R1P->Integral(n1,n2)*R2M->Integral(n1,n2));

  Double_t K1 = sqrt(R) / Power(1.+sqrt(R),2);

  Double_t T1 = sqrt(R1M->Integral(n1,n2)/time1) / R1M->Integral(n1,n2);
  Double_t T2 = sqrt(R1P->Integral(n1,n2)/time1) / R1P->Integral(n1,n2);

  Double_t T3 = sqrt(R2M->Integral(n1,n2)/time2) / R2M->Integral(n1,n2);
  Double_t T4 = sqrt(R2P->Integral(n1,n2)/time2) / R2P->Integral(n1,n2);

  Double_t Error = K1*sqrt(T1*T1+T2*T2+T3*T3+T4*T4);

  return Error;
};
//-------------------------------------------------------------------------------------------------------------
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
};
//-------------------------------------------------------------------------------------------------------------
Float_t GetDiff(Float_t R1M, Float_t R1P, Float_t R2M, Float_t R2P)
{

  // Get the Difference straight forward..

  using namespace TMath;

  Float_t D = (R1M - R1P) - (R2M - R2P);

  return D;
};
//-------------------------------------------------------------------------------------------------------------
Float_t GetDifference(Float_t R1M, Float_t R1P)
{

  // Get the Difference straight forward..

  using namespace TMath;

  Float_t D = R1M - R1P;

  return D;
};
//-----------------------------------------------------------------------------------------------------
void ColorGraphic(TH1 *h, Int_t color, Int_t Style, Int_t Width,Int_t nmarker,Double_t msize)
{
  h->SetMarkerColor(color);
  h->SetLineColor(color);
  h->SetLineWidth(Width);
  h->SetLineStyle(Style);
  h->SetMarkerStyle(nmarker);
  h->SetMarkerSize(msize);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}
//---------------------------------------------------------------------------------------------------------------
void ColorGraphic(TH1 *h, Int_t color, Int_t Style, Int_t Width,Double_t Size)
{
  h->SetLineColor(color);
  h->SetLineWidth(Width);
  h->SetLineStyle(Style);
  h->SetMarkerSize(Size);
  h->SetMarkerColor(color);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
};
//------------------------------------------------------------------------------------
void ColorGraphic(TGraphErrors *h, Int_t color, Int_t Style, Int_t Width,Float_t Size,
                  const char *Title,const char *XTitle,const char *YTitle)
{
  h->SetLineColor(color);
  h->SetMarkerStyle(Style);
  h->SetMarkerSize(Size);
  h->SetMarkerColor(color);
  h->SetTitle(Title);
  h->GetXaxis()->SetTitle(XTitle);
  h->GetYaxis()->SetTitle(YTitle);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
};
//------------------------------------------------------------------------------------
void ColorGraphic(TGraph *h, Int_t color, Int_t Style, Int_t Width,Float_t Size,
                  const char *Title,const char *XTitle,const char *YTitle)
{
  h->SetLineColor(color);
  h->SetMarkerStyle(Style);
  h->SetMarkerColor(color);
  h->SetTitle(Title);
  h->SetMarkerSize(Size);
  h->GetXaxis()->SetTitle(XTitle);
  h->GetYaxis()->SetTitle(YTitle);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
};
//=====================================================================================
void SetRate(Double_t *xR,Double_t *Rate,Double_t *RateEr,Double_t Run,
             Double_t mRate,Double_t mRateEr,Int_t &index)
{ 
   xR[index]     = Run;
   Rate[index]   = mRate;
   RateEr[index] = mRateEr;
  
   if(!(TMath::IsNaN(mRateEr)))index++;

};
//------------------------------------------------------------------------------------
void SetRateV(std::vector<Double_t> &xR,std::vector<Double_t> &Rate,std::vector<Double_t> &RateEr,
              Double_t Run,Double_t mRate,Double_t mRateEr,Int_t nBck)
{


   if(TMath::IsNaN(mRateEr))return;

   xR.push_back(Run);
   Rate.push_back(mRate);
   if(nBck == 0)
  	 RateEr.push_back(mRateEr);
   else if(nBck == 1 )
         RateEr.push_back(mRateEr*mRate);
  
//   if(!(TMath::IsNaN(mRateEr)))index++;
}; 

void GetNorm(TH2F *hAveXRot,TH2F *hAveYRot,TLine *&lNorm,TArrow *&aDis,TCanvas *cCan,Int_t nline)
{
     //-----------------------------------------------------------------------------
     // Select the passed canvas.
     Double_t radius = 75.;
     //-------------------------------------------------------------
     TVector3 DisVec,PosVec,NormVec; 
     TVector3 zhat(0,0,1); // z unit vector
     TVector3 xhat(1,0,0); // x unit vector
     Int_t binx,biny,binz;
     // choice the canvas
     cCan->cd();
     // Get the bin we are looking at
     hAveXRot->GetBinXYZ(nline+1,binx,biny,binz);
     // Find position and displacement..................
     Double_t xpos   = 2.*(Double_t)binx - radius;
     Double_t ypos   = 2.*(Double_t)biny - radius;
     Double_t xdis   = hAveXRot->GetBinContent(binx,biny);
     Double_t ydis   = hAveYRot->GetBinContent(binx,biny);
     // define the diplacement arrow...................
     aDis = new TArrow(xpos,ypos,xpos-xdis,ypos-ydis,0.004,"|>");
     //-----------------------------------------------------------------------------
     DisVec.SetXYZ(xdis,ydis,0);   // Displacement vector
     PosVec.SetXYZ(xpos,ypos,0);   // Position vector
     NormVec = DisVec.Cross(zhat); // Normal Vector
     NormVec *= 50./NormVec.Mag();
     //-----------------------------------------------------------------------------
     Double_t slope = (NormVec(0) > 0) ? NormVec(1)/NormVec(0) : 0;
     // intercept
     Double_t intercept = xpos*tan(NormVec.Angle(xhat)) + ypos;
     // ----------------------------------------------------------------------------
     aDis->SetLineWidth(1);
     aDis->SetFillColor(1);
     // 
     lNorm = new TLine(xpos,ypos,xpos+NormVec(0),ypos+NormVec(1));
     lNorm->SetLineColor(2);

     if(sqrt(xdis*xdis + ydis*ydis) > 2 && sqrt(xdis*xdis + ydis*ydis) < 10){
         aDis->Draw("");
         lNorm->Draw("");
     }

};

Int_t GetInterSection(TLine *Line1,TLine *Line2,Double_t &Xin,Double_t &Yin)
{
     using namespace std;
     using namespace TMath;

     if(IsNaN(Line1->GetX2()) || IsNaN(Line2->GetX2())){
        // if these lines are bad (don't have a real displacement vector 
        // ignore them
	Xin =  0.;
        Yin =  0.;
        return 0;
     } 

     // Define a few vectors...................     

     TVector3 Vec1(Line1->GetX2() - Line1->GetX1() , Line1->GetY2() - Line1->GetY1(),0);
     TVector3 Vec2(Line2->GetX2() - Line2->GetX1() , Line2->GetY2() - Line2->GetY1(),0);
     TVector3 xhat(1,0,0);

     // Get the angle between x-unit vector.... 
     Double_t angle1 = (Line1->GetX1() < 0) ? Vec1.Angle(xhat) : Vec1.Angle(-xhat);
     Double_t angle2 = (Line2->GetX1() < 0) ? Vec2.Angle(xhat) : Vec2.Angle(-xhat);
    
     if(Vec1(0) != 0 && Vec2(0) != 0 ){

	 Double_t slope1 = Vec1(1)/Vec1(0);
    	 Double_t slope2 = Vec2(1)/Vec2(0);

     	 Double_t inter1 = Abs(Line1->GetX1())*tan(angle1) + Line1->GetY1();
    	 Double_t inter2 = Abs(Line2->GetX1())*tan(angle2) + Line2->GetY1();
     
    	 Xin = (inter2 - inter1 ) / (slope1 - slope2);
     
     	 Double_t Y1 = Xin*slope1 + inter1;
     	 Double_t Y2 = Xin*slope2 + inter2;
         Yin = Y1;

         if(Abs(Yin) > 200 || Abs(Xin) >200){
	   Xin = 0.;
           Yin = 0.;
           return 0;
        }
     } else {
       Xin = 0.;
       Yin = 0.;
       return 0;
     };
   
    return 1;
};
//-------------------------------------------------------------------------------------------
void Get2DGaussianFit(TH2F *h2,Double_t &X,Double_t &Y)
{
  
   TCanvas *c = new TCanvas("c");
   Int_t binx,biny,binz;
  
   std::cout << "Intergral counts of the 2d histo is " << h2->Integral() << std::endl;
   h2->GetBinXYZ(h2->GetMaximumBin(),binx,biny,binz);
  
   std::cout << "Max bin "  << binx << " , " << biny << std::endl;
   Double_t Xmax = binx*(h2->GetXaxis()->GetXmax()-h2->GetXaxis()->GetXmin())/h2->GetNbinsX() + h2->GetXaxis()->GetXmin();
   Double_t Ymax = biny*(h2->GetYaxis()->GetXmax()-h2->GetYaxis()->GetXmin())/h2->GetNbinsY() + h2->GetYaxis()->GetXmin();

   std::cout << "Max position " << Xmax << " , " << Ymax << std::endl;

   TF2 *Gaus2d = new TF2("Gaus2d","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",-100,100,-50,150);

   Gaus2d->SetParLimits(0,h2->GetBinContent(h2->GetMaximumBin()),1e6);
   Gaus2d->SetParLimits(1,Xmax-10,Xmax+10);
   Gaus2d->SetParLimits(2,2,15);
   Gaus2d->SetParLimits(3,Ymax-10,Ymax+10);
   Gaus2d->SetParLimits(4,2,30);
   Gaus2d->SetParameter(1,Xmax);
   Gaus2d->SetParameter(3,Ymax);

   if(h2->Integral() > 0)h2->Fit("Gaus2d","EM");
   // Set the fit parameters
   X = Gaus2d->GetParameter(1);
   Y = Gaus2d->GetParameter(3);

   h2->Draw("colz");
   Gaus2d->Draw("l same");

   c->Print(Form("output/%s.pdf",h2->GetTitle()));

   delete Gaus2d;   
   delete c;

};
//=========================================================================================
Double_t Return_Rotation_Angle(TArrow *aDis,Double_t X,Double_t Y)
{
  //====================================================================================+
  // Define vectors from the average intersection point to the beginning and end of     |
  // the displacement vector.                                                           |
  //====================================================================================+
  TVector3 v1(aDis->GetX1()-X,aDis->GetY1()-Y);
  TVector3 v2(aDis->GetX2()-X,aDis->GetY2()-Y);
  // Use the vector library operations to get the angle between v1 and v2.
  Double_t angle = v1.Angle(v2);

  return angle;
};

Double_t PositionRotate(Double_t Angle,Double_t X,Double_t Y,TArrow *aDis)
{

  using namespace TMath;
  //---------------------------------------------------
  TMatrixD Rot;
  Rot.Clear();
  Rot.ResizeTo(2,2);
  Rot(0,0) = Cos(Angle);
  Rot(0,1) = Sin(Angle);
  Rot(1,0) = -Sin(Angle);
  Rot(1,1) = Cos(Angle);
  //---------------------------------------------------
  TVectorD Pos;
  Pos.ResizeTo(2);
  Pos(0) = aDis->GetX1();
  Pos(1) = aDis->GetX2();

  TVectorD Oprime;
  Oprime.ResizeTo(2);
  Oprime(0) = X;
  Oprime(1) = Y;
  //---------------------------------------------------
  // R' = |M|(R-T) + T
  Pos -= Oprime;
  Pos *= Rot; 
  Pos += Oprime;
  //---------------------------------------------------
  Double_t err = sqrt(Power(Pos(0) - aDis->GetX2(),2) +
                      Power(Pos(1) - aDis->GetY2(),2));

  return err;
};
//=========================================================================================

void DrawEnerPanel(TH1F *hData,TH1F *hMC,TCanvas *c,Int_t npad,TLegend *lleg,Double_t time)
{
  using namespace std;

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
    //  double MM = hMC->Integral(1,32);
    //  double NN = hData->Integral(1,32);
      double s1 = hData->GetBinError(i);
     // double s2 = hMC->GetBinError(i);
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


#endif
