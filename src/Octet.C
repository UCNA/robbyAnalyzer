#include "Octet.h"

using namespace std;

Octet::Octet(Int_t n,Int_t j, Int_t m,Double_t radius)
{
  Num   = n;
  First = j;
  Last  = m;
  radial_cut = radius;
}

Octet::~Octet()
{
  
}

void Octet::CalculateSuperTime(Float_t R1M, Float_t R1MI, Float_t R1P, Float_t R1PI, Float_t R2M, Float_t R2MI, Float_t R2P, Float_t R2PI, Float_t &R)
{

  // Get the Super Ratio - 1 straight forward.. with Floats.

  R = ((R1M - R1MI)*(R2P - R2PI) / ((R1P - R1PI)*(R2M - R2MI))) - 1.;

};

void Octet::CalculateSuper(Double_t R1,Double_t R2,Double_t R3,Double_t R4,
			   Double_t B1,Double_t B2,Double_t B3,Double_t B4,  
			   Double_t t1,Double_t t1b,Double_t t2,Double_t t2b,
			   Double_t t3,Double_t t3b,Double_t t4,Double_t t4b,
			   Double_t &A,Double_t &Aer)
{
    Double_t RR;
    Double_t er1 = ErrorMultiRate(R1  ,B1 ,t1 ,t1b ); 
    Double_t er2 = ErrorMultiRate(R2  ,B2 ,t2 ,t2b );
    Double_t er3 = ErrorMultiRate(R3  ,B3 ,t3 ,t3b );
    Double_t er4 = ErrorMultiRate(R4  ,B4 ,t4 ,t4b ); 
  
    Double_t ersum = TMath::Power(er1/R1,2) + TMath::Power(er2/R2,2) 
		    +TMath::Power(er3/R3,2) + TMath::Power(er4/R4,2);  

    RR = (R1*R2)/(R3*R4);
  
    if(TMath::IsNaN(RR)){
      A = 0.;
      Aer = 0.;
    } else {
      A   = (1-TMath::Sqrt(RR))/(1+TMath::Sqrt(RR));
      Aer = TMath::Sqrt(RR)/TMath::Power((1. + TMath::Sqrt(RR)),2) * TMath::Sqrt(ersum);
    }
};
//================================================================================
Double_t Octet::Calc_Super()
{
  //-------------------------------------------------------------------------
  // Calculate the Super Ratio 
  // inputs..........
  //
  // nAi(w)[j] where i=2,5,7,10 and j is the index of the analysis choice
  // 0 = A, 1 = B, 2 = C, ...,etc. We've chosen to look at method C.
  //
  // nBi(w)[j]... the same but for the B quartet.
  //
  // tAi,tBi  ... the run time for the signal
  // tAib,tBib... the run time for the background runs
  using namespace TMath;
  //Loop through analysis choices
  for(Int_t ich = 0 ; ich < 9 ; ich++){
      //==========================================================
      CalculateSuper(nA2[ich],nA5w[ich],nA5[ich],nA2w[ich]
		    ,nA2b[ich],nA5bw[ich],nA5b[ich],nA2bw[ich]
		    ,tA2,tA2b,tA5,tA5b,tA2w,tA2bw,tA5w,tA5bw
		    ,Asuper1[ich],Asuper1er[ich]);
      //==========================================================
      CalculateSuper(nA10[ich],nA7w[ich],nA7[ich],nA10w[ich]
		    ,nA10b[ich],nA7bw[ich],nA7b[ich],nA10bw[ich]
		    ,tA10,tA10b,tA7,tA7b,tA10w,tA10bw,tA7w,tA7bw
		    ,Asuper2[ich],Asuper2er[ich]);
      //==========================================================
      CalculateSuper(nB2[ich],nB5w[ich],nB5[ich],nB2w[ich]
		    ,nB2b[ich],nB5bw[ich],nB5b[ich],nB2bw[ich]
		    ,tB2,tB2b,tB5,tB5b,tB2w,tB2bw,tB5w,tB5bw
		    ,Bsuper1[ich],Bsuper1er[ich]);
      //==========================================================
      CalculateSuper(nB10[ich],nB7w[ich],nB7[ich],nB10w[ich]
		    ,nB10b[ich],nB7bw[ich],nB7b[ich],nB10bw[ich]
		    ,tB10,tB10b,tB7,tA7b,tB10w,tB10bw,tB7w,tB7bw
		    ,Bsuper2[ich],Bsuper2er[ich]);
      //==========================================================
  }
  return 0.;
}
//=====================================================================================
Float_t Octet::Calc_Super_Time()
{
  //-------------------------------------------------------------------------
  // Calculate the Super Ratio
  // inputs..........
  //
  //
  // tA(e,w)i,tB(e,w)i  ... the time stamp of first (beta) event in run
  // tA(e,w)f,tB(e,w)f  ... the time stamp of last (beta) event in run
  using namespace TMath;
  //Loop through analysis choices
      //==========================================================
      CalculateSuperTime(tA2ef,tA2ei,tA5ef,tA5ei,tA2wf,tA2wi,tA5wf,tA5wi
                    ,Asupertime1);
      //==========================================================
      CalculateSuperTime(tA10ef,tA10ei,tA7ef,tA7ei,tA10wf,tA10wi,tA7wf,tA7wi
                    ,Asupertime2);
      //==========================================================
      CalculateSuperTime(tB5ef,tB5ei,tB2ef,tB2ei,tB5wf,tB5wi,tB2wf,tB2wi
                    ,Bsupertime1);
      //==========================================================
      CalculateSuperTime(tB7ef,tB7ei,tB10ef,tB10ei,tB7wf,tB7wi,tB10wf,tB10wi
                    ,Bsupertime2);
      //==========================================================
  return 0.;
}
//=====================================================================================
Float_t Octet::Calc_Super_Time_All()
{
  //-------------------------------------------------------------------------
  // Calculate the Super Ratio
  // inputs..........
  //
  //
  // tA(e,w)i,tB(e,w)i  ... the time stamp of first event in run
  // tA(e,w)f,tB(e,w)f  ... the time stamp of last event in run
  using namespace TMath;
  //Loop through analysis choices
      //==========================================================
      CalculateSuperTime(tA2efa,tA2eia,tA5efa,tA5eia,tA2wfa,tA2wia,tA5wfa,tA5wia
                    ,Asupertime1a);
      //==========================================================
      CalculateSuperTime(tA10efa,tA10eia,tA7efa,tA7eia,tA10wfa,tA10wia,tA7wfa,tA7wia
                    ,Asupertime2a);
      //==========================================================
      CalculateSuperTime(tB5efa,tB5eia,tB2efa,tB2eia,tB5wfa,tB5wia,tB2wfa,tB2wia
                    ,Bsupertime1a);
      //==========================================================
      CalculateSuperTime(tB7efa,tB7eia,tB10efa,tB10eia,tB7wfa,tB7wia,tB10wfa,tB10wia
                    ,Bsupertime2a);
      //==========================================================
  return 0.;
}
//=====================================================================================
Double_t Octet::Calc_A_multi()
{
  
  using namespace TMath;
  
  Double_t RA1,RA2,RA3,RA4,RB1,RB2,RB3,RB4;
  Double_t erA1,erA2,erA3,erA4,erA5,erA6,erA7,erA8,ersumA;
  Double_t erB1,erB2,erB3,erB4,erB5,erB6,erB7,erB8,ersumB;
  Double_t RR;
  Double_t pow = 0.;
  
  for(Int_t j = 0 ; j < 9 ; j++){
  
    if(nA2[j]   == 0 )   nA2[j]  = 1.;
    if(nA2w[j]  == 0 )  nA2w[j]  = 1.;
    if(nA5[j]   == 0 )   nA5[j]  = 1.;
    if(nA5w[j]  == 0 )  nA5w[j]  = 1.;
    if(nA7[j]   == 0 )   nA7[j]  = 1.;
    if(nA7w[j]  == 0 )  nA7w[j]  = 1.;
    if(nA10[j]  == 0 )  nA10[j]  = 1.;
    if(nA10w[j] == 0 ) nA10w[j]  = 1.;
  
    // Calculate the rates 
    RA1 = nA2[j]*nA10[j];
    RA2 = nA5w[j]*nA7w[j];
    RA3 = nA5[j]*nA7[j];
    RA4 = nA2w[j]*nA10w[j];
  
    //Check for zero's from missing runs...
    // this allow for quartet, octet creation from sets 
    // where runs are missing
    RA1 = ( RA1 != 0 ) ? RA1 : 1.;
    RA2 = ( RA2 != 0 ) ? RA2 : 1.;
    RA3 = ( RA3 != 0 ) ? RA3 : 1.;
    RA4 = ( RA4 != 0 ) ? RA4 : 1.;
    
    erA1 = (nA2[j]   != 1.) ? ErrorMultiRate(nA2[j],nA2b[j],tA2,tA2b) : 0;
    erA2 = (nA5[j]   != 1.) ? ErrorMultiRate(nA5[j],nA5b[j],tA5,tA5b) : 0;
    erA3 = (nA7[j]   != 1.) ? ErrorMultiRate(nA7[j],nA7b[j],tA7,tA7b) : 0;
    erA4 = (nA10[j]  != 1.) ? ErrorMultiRate(nA10[j],nA10b[j],tA10,tA10b): 0;
    
    erA5 = (nA2w[j]  != 1.) ? ErrorMultiRate(nA2w[j],nA2bw[j],tA2w,tA2bw) : 0;
    erA6 = (nA5w[j]  != 1.) ? ErrorMultiRate(nA5w[j],nA5bw[j],tA5w,tA5bw) : 0;
    erA7 = (nA7w[j]  != 1.) ? ErrorMultiRate(nA7w[j],nA7bw[j],tA7w,tA7bw) : 0;
    erA8 = (nA10w[j] != 1.) ? ErrorMultiRate(nA10w[j],nA10bw[j],tA10w,tA10bw) : 0;
      
    ersumA = Power(erA1/nA2[j],2) + Power(erA2/nA5[j],2) + Power(erA3/nA7[j],2) + Power(erA4/nA10[j],2) 
           + Power(erA5/nA2w[j],2) + Power(erA6/nA5w[j],2) + Power(erA7/nA7w[j],2) + Power(erA8/nA10w[j],2); 
    
    if(nA2[j]  != 1)pow++;
    if(nA5[j]  != 1)pow++;
    if(nA7[j]  != 1)pow++;
    if(nA10[j] != 1)pow++;
    
    RR = Power((RA1*RA2)/(RA3*RA4),1./pow);
    
    A_multi_A[j]   = ( 1. - RR )/( 1.+ RR );
    A_multi_Aer[j] = (2./pow)*RR/Power((1. + RR),2) * sqrt(ersumA);
  //--------------------------------------------------------------------------
  // check to see if octet is incomplete
  //
    if(nB2[j]  == 0 )   nB2[j] = 1.;
    if(nB2w[j] == 0 )  nB2w[j] = 1.;
    if(nB5[j]  == 0 )   nB5[j] = 1.;
    if(nB5w[j] == 0 )  nB5w[j] = 1.;
    if(nB7[j]  == 0 )   nB7[j] = 1.;
    if(nB7w[j] == 0 )  nB7w[j] = 1.;
    if(nB10[j] == 0 )  nB10[j] = 1.;
    if(nB10w[j]== 0 ) nB10w[j] = 1.;
   
    // Multiply Rates
    RB1 = nB2[j]  * nB10[j];
    RB2 = nB5w[j] * nB7w[j];
    RB3 = nB5[j]  * nB7[j];
    RB4 = nB2w[j] * nB10w[j];
  
    RB1 = ( RB1 != 0.  ) ? RB1 : 1.;
    RB2 = ( RB2 != 0.  ) ? RB2 : 1.;
    RB3 = ( RB3 != 0.  ) ? RB3 : 1.;
    RB4 = ( RB4 != 0.  ) ? RB4 : 1.;
    
    erB1 = (nB2[j]   != 1.) ? ErrorMultiRate(nB2[j] ,nB2b[j] ,tB2 ,tB2b) : 0;
    erB2 = (nB5[j]   != 1.) ? ErrorMultiRate(nB5[j] ,nB5b[j] ,tB5 ,tB5b) : 0;
    erB3 = (nB7[j]   != 1.) ? ErrorMultiRate(nB7[j] ,nB7b[j] ,tB7 ,tB7b) : 0;
    erB4 = (nB10[j]  != 1.) ? ErrorMultiRate(nB10[j],nB10b[j],tB10,tB10b): 0;
    
    erB5 = (nB2w[j]  != 1.) ? ErrorMultiRate(nB2w[j] ,nB2bw[j] ,tB2w ,tB2bw) : 0;
    erB6 = (nB5w[j]  != 1.) ? ErrorMultiRate(nB5w[j] ,nB5bw[j] ,tB5w ,tB5bw) : 0;
    erB7 = (nB7w[j]  != 1.) ? ErrorMultiRate(nB7w[j] ,nB7bw[j] ,tB7w ,tB7bw) : 0;
    erB8 = (nB10w[j] != 1.) ? ErrorMultiRate(nB10w[j],nB10bw[j],tB10w,tB10bw) : 0;
   
    pow = 0.;
    
    if(nB2[j]  != 1)pow++;
    if(nB5[j]  != 1)pow++;
    if(nB7[j]  != 1)pow++;
    if(nB10[j] != 1)pow++;
    
    ersumB = Power(erB1/nB2[j],2) + Power(erB2/nB5[j],2) + Power(erB3/nB7[j],2) + Power(erB4/nB10[j],2) 
           + Power(erB5/nB2w[j],2) + Power(erB6/nB5w[j],2) + Power(erB7/nB7w[j],2) + Power(erB8/nB10w[j],2); 
    
    RR = Power((RB1*RB2)/(RB3*RB4),1./pow);
  
    A_multi_B[j]   = (1. - RR)/(1. + RR);
    A_multi_Ber[j] = (2./pow)*RR/Power((1. + RR),2) * sqrt(ersumB);
  
    // Calculate full Octet answer.....
    
    RR  = (RA1*RA2*RB3*RB4)/(RA3*RA4*RB1*RB2);
    
    pow = 0.; 
    
    if(nA2[j]  != 1)pow++;
    if(nA5[j]  != 1)pow++;
    if(nA7[j]  != 1)pow++;
    if(nA10[j] != 1)pow++;
    if(nB2[j]  != 1)pow++;
    if(nB5[j]  != 1)pow++;
    if(nB7[j]  != 1)pow++;
    if(nB10[j] != 1)pow++;
    
    RR = Power(RR,1./pow);
    
    A_multi[j] = (1. - RR)/(1. + RR);
    A_multier[j] = (2./pow)*(RR/Power((1.+ RR),2))*sqrt(ersumB + ersumA);
    
    pow = 0.;
  }

  
  return 0.; 
}

Double_t Octet::Calc_A_sum(){
  
  using namespace TMath;
 
  Double_t RA1,RA2,RA3,RA4;
  Double_t RB1,RB2,RB3,RB4;
  Double_t er1,er2,er3,er4,eArsum,eBrsum,RR,ersum;
  
  // ----------------------------------------------------
  // Super ratio will = (R1*R2)/(R3*R4)
  // fill terms based on the assigned runs
  for(Int_t j = 0 ; j < 9 ; j++){
    
    RA1 = Average_Quartet_Rate(nA2[j] ,tA2 ,nA2b[j] ,tA2b ,nA10[j] ,tA10 ,nA10b[j] ,tA10b,0);
    RA2 = Average_Quartet_Rate(nA5w[j],tA5w,nA5bw[j],tA5bw,nA7w[j] ,tA7w ,nA7bw[j] ,tA7bw,1);
    RA3 = Average_Quartet_Rate(nA5[j] ,tA5 ,nA5b[j] ,tA5b ,nA7[j]  ,tA7  ,nA7b[j]  ,tA7b,0);
    RA4 = Average_Quartet_Rate(nA2w[j],tA2w,nA2bw[j],tA2bw,nA10w[j],tA10w,nA10bw[j],tA10bw,1);
  
    // for broken octet some R's will be zero.  Set to 1 so that
    // they don't effect the super ratio.  Should be in combinations.
    // thus if R1 = 0 then R4 should also be 0.  
    RA1 = ( RA1 > 0 ) ? RA1 : 1.;  RA2 = ( RA2 > 0 ) ? RA2 : 1.;
    RA3 = ( RA3 > 0 ) ? RA3 : 1.;  RA4 = ( RA4 > 0 ) ? RA4 : 1.;
  
    // Calculate Errors --------------------------------------------  
    er1 = ErrorSum(nA2[j] ,nA10[j] ,tA2 ,tA10 ,nA2b[j] ,nA10b[j] ,tA2b ,tA10b,0);
    er2 = ErrorSum(nA5[j] ,nA7[j]  ,tA5 ,tA7  ,nA5b[j] ,nA7b[j]  ,tA5b ,tA7b,0);
    er3 = ErrorSum(nA5w[j],nA7w[j] ,tA5w,tA7w ,nA5bw[j],nA7bw[j] ,tA5bw,tA7bw,1);
    er4 = ErrorSum(nA2w[j],nA10w[j],tA2w,tA10w,nA2bw[j],nA10bw[j],tA2bw,tA10bw,1);
    
    // Sum the errors   
    eArsum = Power(er1/RA1,2) + Power(er2/RA3,2) + Power(er3/RA2,2) + Power(er4/RA4,2);
    // Calculate the super-ratio...
    RR = (RA1*RA2)/(RA3*RA4);
    // Check if we messed up and RR is a NaN
    if(IsNaN(RR)){
      // if so set it to A =0 and its error to 1,
      // it will be throw out in the collection function
      A_sum_A[j]   = 0.;
      A_sum_Aer[j] = 1.;
    } else {
      // Else calculate the quartet A
      A_sum_A[j]   = (1. - sqrt(RR))/(1. + sqrt(RR));
      A_sum_Aer[j] = (sqrt(RR))/Power((1. + sqrt(RR)),2) * sqrt(eArsum);
    }
    
    // Now repeat for the B quartet
    
    RB1 = Average_Quartet_Rate(nB2[j] ,tB2 ,nB2b[j] ,tB2b ,nB10[j] ,tB10 ,nB10b[j] ,tB10b,0);
    RB2 = Average_Quartet_Rate(nB5w[j],tB5w,nB5bw[j],tB5bw,nB7w[j] ,tB7w ,nB7bw[j] ,tB7bw,1);
    RB3 = Average_Quartet_Rate(nB5[j] ,tB5 ,nB5b[j] ,tB5b ,nB7[j]  ,tB7  ,nB7b[j]  ,tB7b,0);
    RB4 = Average_Quartet_Rate(nB2w[j],tB2w,nB2bw[j],tB2bw,nB10w[j],tB10w,nB10bw[j],tB10bw,1);
    
    RB1 = ( RB1 > 0 ) ? RB1 : 1.; RB2 = ( RB2 > 0 ) ? RB2 : 1.;
    RB3 = ( RB3 > 0 ) ? RB3 : 1.; RB4 = ( RB4 > 0 ) ? RB4 : 1.;
  
    er1 = ErrorSum(nB2[j] ,nB10[j] ,tB2 ,tB10,nB2b[j] ,nB10b[j] ,tB2b ,tB10b,1);
    er2 = ErrorSum(nB5[j] ,nB7[j]  ,tB5 ,tB7 ,nB5b[j] ,nB7b[j]  ,tB5b ,tB7b,1);
    er3 = ErrorSum(nB5w[j],nB7w[j] ,tB5w,tB7w,nB5bw[j],nB7bw[j] ,tB5bw,tB7bw,0);
    er4 = ErrorSum(nB2w[j],nB10w[j],tB2w,tB10w,nB2bw[j],nB10bw[j],tB2bw,tB10bw,0);
    
    RR = (RB1*RB2)/(RB3*RB4); 
    eBrsum = Power(er1/RB1,2) + Power(er2/RB3,2) + Power(er3/RB2,2) + Power(er4/RB4,2);
  
    if(IsNaN(RR)){
      A_sum_B[j] = 0.;
      A_sum_Ber[j] = 1.;
    } else {
      A_sum_B[j]   = (1.-sqrt(RR))/(1.+sqrt(RR));
      A_sum_Ber[j] = (sqrt(RR))/Power((1. + sqrt(RR)),2) * sqrt(eBrsum);
    }
    //======================================================================================
    RA1 = Average_Octet_Rate(nA2[j],tA2,nA2b[j],tA2b,nA10[j],tA10,nA10b[j],tA10b,
			     nB5[j],tB5,nB5b[j],tB5b,nB7[j] ,tB7 ,nB7b[j] ,tB7b,0);
			       
    RA2 = Average_Octet_Rate(nA5w[j],tA5w,nA5bw[j],tA5bw,nA7w[j],tA7w,nA7bw[j],tA7bw,
			     nB2w[j],tB2w,nB2bw[j],tB2bw,nB10w[j] ,tB10w ,nB10bw[j] ,tB10bw,1);

    RA3 = Average_Octet_Rate(nA2w[j],tA2w,nA2bw[j],tA2bw,nA10w[j],tA10w,nA10bw[j],tA10bw,
			     nB5w[j],tB5w,nB5bw[j],tB5bw,nB7w[j] ,tB7w ,nB7bw[j] ,tB7bw,1);

    RA4 = Average_Octet_Rate(nA5[j],tA5,nA5b[j],tA5b,nA7[j],tA7,nA7b[j],tA7b,
			     nB2[j],tB2,nB2b[j],tB2b,nB10[j] ,tB10 ,nB10b[j] ,tB10b,0);
			       
  // Finished with Quartet Asymmetries.....
    RA1 = ( RA1 != 1. ) ? RA1 : 1.; RA2 = ( RA2 != 1. ) ? RA2 : 1.;
    RA3 = ( RA3 != 1. ) ? RA3 : 1.; RA4 = ( RA4 != 1. ) ? RA4 : 1.;

    RR  = (RA1*RA2)/(RA3*RA4);
    
    er1 = ErrorSum(nA2[j] ,nA10[j] ,nB5[j] ,nB7[j]  ,tA2 ,tA10 ,tB5 ,tB7,
		   nA2b[j] ,nA10b[j] ,nB5b[j] ,nB7b[j]  ,tA2b ,tA10b ,tB5b ,tB7b,0);
    
    er2 = ErrorSum(nA5w[j],nA7w[j] ,nB2w[j],nB10w[j],tA5w,tA7w ,tB2w,tB10w,
		   nA5bw[j],nA7bw[j] ,nB2bw[j],nB10bw[j],tA5bw,tA7bw ,tB2bw,tB10bw,1);
    
    er3 = ErrorSum(nA2w[j],nA10w[j],nB5w[j],nB7w[j] ,tA2w,tA10w,tB5w,tB7w,
		   nA2bw[j],nA10bw[j],nB5bw[j],nB7bw[j] ,tA2bw,tA10bw,tB5bw,tB7bw,1);
    
    er4 = ErrorSum(nA5[j] ,nA7[j]  ,nB2[j] ,nB10[j] ,tA5 ,tA7  ,tB2 ,tB10,
		   nA5b[j] ,nA7b[j]  ,nB2b[j] ,nB10b[j] ,tA5b ,tA7b  ,tB2b ,tB10b,0);
  
    ersum = Power(er1/RA1,2) + Power(er2/RA2,2) + Power(er3/RA3,2) + Power(er4/RA4,2);
    
    
    //if(j == 2 ) cout << RA1 << "  " << RA2 << "   "  << RA3 << "   " << RA4 << endl;
 
    if( !(IsNaN(RR))){
      A_sum[j]    = (1. - sqrt(RR)) / ( 1. + sqrt(RR));
      A_sumer[j]  = sqrt(RR)/Power((1. + sqrt(RR)),2) * sqrt(ersum);
    } else {
      A_sum[j] = 0.;
      A_sumer[j] =1.;
    }
    
    //if(j==2)cout << "Octet asymmetry " << A_sum[j] <<  "  " << A_sumer[j] <<endl;
  }
  
  return 0.;
}

void Octet::Get_Rad_A()
{
using namespace TMath; 
  Double_t R1,R2,R3,R4;
  Double_t er1,er2,er3,er4,ersum;
  Double_t RR;
  
  for(Int_t i = 0 ; i < 12 ; i++){
   
    R1 = (A2rad[i]!=0) ? A2rad[i]/tA2 : 0;
    R1 += (A10rad[i]!=0) ? A10rad[i]/tA10 : 0;
    if(R1 == 0)R1=1;

    R2 = (A5wrad[i]!=0) ? A5wrad[i]/tA5w : 0;
    R2 += (A7wrad[i]!=0) ? A7wrad[i]/tA7w : 0;
    if(R2 == 0) R2 = 1.;

    R3 = (A5rad[i]!=0) ? A5rad[i]/tA5 : 0;
    R3 += (A7rad[i]!=0) ? A7rad[i]/tA7 : 0;
    if(R3 == 0)R3 = 1.;
    R4 = (A2wrad[i]!=0) ? A2wrad[i]/tA2w : 0;
    RR += (A10wrad[i]!=0) ?  A10wrad[i]/tA10w : 0;
    if(R4 == 0)R4 = 1.;

    TotRadCounts[i] = A2rad[i] + A10rad[i] + A5rad[i] + A7rad[i] 
                   + A5wrad[i] + A7wrad[i] + A2wrad[i] + A10wrad[i];
  
    er1 = A2rade[i]*A2rade[i]   + A10rade[i]*A10rade[i];
    er2 = A5rade[i]*A5rade[i]   + A7rade[i]*A7rade[i];
    er3 = A5wrade[i]*A5wrade[i] + A7wrade[i]*A7wrade[i];
    er4 = A2wrade[i]*A2wrade[i] + A10wrade[i]*A10wrade[i];

    // cout << R1 << "\t" << R2 << "\t" <<R3 <<"\t"<< R4 << endl;

    ersum = er1/(R1*R1) + er2/(R2*R2) + er3/(R3*R3) + er4/(R4*R4);
  
    RR = Abs((R1*R2)/(R3*R4));
  
    if(IsNaN(RR) || IsNaN(ersum)){
       RR    = 0.;     
       ersum = 1.;
       A_rad[i]   =  0.;
       A_rader[i] = 0.1;
    } else {
      A_rad[i]   = (1. - sqrt(RR))/(1. + sqrt(RR));
      A_rader[i] =  sqrt(RR)/Power((1. + sqrt(RR)),2) * sqrt(ersum);
    }
  }
  
  for(Int_t i = 0 ; i < 12; i++){

    R1 = (B2rad[i]!=0) ? B2rad[i]/tB2 : 0;
    R1 += (B10rad[i]!=0) ? B10rad[i]/tB10 : 0;
    if(R1 == 0)R1=1;

    R2 = (B5wrad[i]!=0) ? B5wrad[i]/tB5w : 0;
    R2 += (B7wrad[i]!=0) ? B7wrad[i]/tB7w : 0;
    if(R2 ==0) R2 = 1.;

    R3 = (B5rad[i]!=0) ? B5rad[i]/tB5 : 0;
    R3 += (B7rad[i]!=0) ? B7rad[i]/tB7 : 0;

    if(R3 == 0)R3 = 1.;
    R4 = (B2wrad[i]!=0) ? B2wrad[i]/tB2w : 0;
    R4 += (B10wrad[i]!=0) ?  B10wrad[i]/tB10w : 0;
    if(R4 == 0)R4 = 1.;
 
    TotRadCounts[i] += B2rad[i] + B10rad[i] + B5rad[i] + B7rad[i] 
                     + B5wrad[i] + B7wrad[i] + B2wrad[i] + B10wrad[i];
   
    er1 = B2rade[i]*B2rade[i]   + B10rade[i]*B10rade[i];
    er2 = B5rade[i]*B5rade[i]   + B7rade[i]*B7rade[i];
    er3 = B5wrade[i]*B5wrade[i] + B7wrade[i]*B7wrade[i];
    er4 = B2wrade[i]*B2wrade[i] + B10wrade[i]*B10wrade[i];
    //cout << R1 << "\t" << R2 << "\t" <<R3 <<"\t"<< R4 << endl;
    ersum = er1/(R1*R1) + er2/(R2*R2) + er3/(R3*R3) + er4/(R4*R4);
  
    RR = Abs((R1*R2)/(R3*R4));
  
    if(IsNaN(RR) ){
      RR = 0.;
      A_rad[i]   = 0.;
      A_rader[i] = 1.;
    } else {
      A_rad[i]   = (A_rad[i] + Abs((1. - sqrt(RR))/(1. + sqrt(RR))))/2.;
      A_rader[i] = sqrt(1./(1./(A_rader[i]*A_rader[i]) +
                       1./(RR/Power((1. + sqrt(RR)),4) * ersum)));
    }
  }

}
//---------------------------------------------------------------------------------
void Octet::Calc_A_sum_Bin()
{
  using namespace TMath; 
  Double_t R1,R2,R3,R4;
  Double_t er1,er2,er3,er4,ersum;
  Double_t RR;
  for(Int_t nn = 0 ; nn < 9 ; nn++){ 
  if(nA2[nn] != 1 || nA5[nn] != 1)  {
    for(Int_t i = 1 ; i < hAsyA[nn]->GetNbinsX() ; i++){
     
      R1 = 1.;R2 = 1.;R3 = 1.;R4 = 1.;
    
      R1 = Average_Quartet_Bin(hA10[nn]->GetBinContent(i),hA10[nn]->GetBinError(i)
			      ,hA2[nn]->GetBinContent(i),hA2[nn]->GetBinError(i),er1);
      R2 = Average_Quartet_Bin(hA7w[nn]->GetBinContent(i),hA7w[nn]->GetBinError(i)
			      ,hA5w[nn]->GetBinContent(i),hA5w[nn]->GetBinError(i),er2);
      R3 = Average_Quartet_Bin(hA7[nn]->GetBinContent(i),hA7[nn]->GetBinError(i)
			      ,hA5[nn]->GetBinContent(i),hA5[nn]->GetBinError(i),er3);
      R4 = Average_Quartet_Bin(hA10w[nn]->GetBinContent(i),hA10w[nn]->GetBinError(i)
			      ,hA2w[nn]->GetBinContent(i),hA2w[nn]->GetBinError(i),er4);
			    
      ersum = Power(er1/R1,2) + Power(er2/R2,2) + Power(er3/R3,2) + Power(er4/R4,2);
      
      RR = Power((R1*R2)/(R3*R4),1./2.);
  
      if(IsNaN(RR)){
	RR = 1.;
	ersum = 1.;
      } else {
	hAsyA[nn]->SetBinContent(i,(1. - RR)/(1. + RR));
	hAsyA[nn]->SetBinError(i, (2./2.)*RR/Power((1. + RR),2) * sqrt(ersum));
      }
    }
   }
  // now loop over B octets
  for(Int_t i = 1 ; i <= hAsyB[nn]->GetNbinsX() ; i++){
    
    R1 = 1.;
    R2 = 1.;
    R3 = 1.;
    R4 = 1.;
    
    R1 = Average_Quartet_Bin(hB10[nn]->GetBinContent(i),hB10[nn]->GetBinError(i)
			    ,hB2[nn]->GetBinContent(i),hB2[nn]->GetBinError(i),er1);
    R2 = Average_Quartet_Bin(hB7w[nn]->GetBinContent(i),hB7w[nn]->GetBinError(i)
			    ,hB5w[nn]->GetBinContent(i),hB5w[nn]->GetBinError(i),er2);
    R3 = Average_Quartet_Bin(hB7[nn]->GetBinContent(i),hB7[nn]->GetBinError(i)
			    ,hB5[nn]->GetBinContent(i),hB5[nn]->GetBinError(i),er3);
    R4 = Average_Quartet_Bin(hB10w[nn]->GetBinContent(i),hB10w[nn]->GetBinError(i)
			    ,hB2w[nn]->GetBinContent(i),hB2w[nn]->GetBinError(i),er4);
			    
    ersum = Power(er1/R1,2) + Power(er2/R2,2) + Power(er3/R3,2) + Power(er4/R4,2);
  
    RR = Power((R1*R2)/(R3*R4),1./2.);
    
    if(IsNaN(RR)){
  	RR = 1.;
	ersum = 1.;
    } else {
        hAsyB[nn]->SetBinContent(i, (1.- RR)/(1. + RR));
        hAsyB[nn]->SetBinError(i,RR/Power((1. + RR),2) * sqrt(ersum));
    }
  }

  for(Int_t bin = 1 ; bin <= hAsyTot[nn]->GetNbinsX() ; bin++){
    
    R1 = 1.;
    R2 = 1.;
    R3 = 1.;
    R4 = 1.;
    
    R1 = Average_Octet_Cts_Bin(hA2[nn]->GetBinContent(bin),tA2,
			       hA1[nn]->GetBinContent(bin),tA2b,
			       hA10[nn]->GetBinContent(bin),tA10,
			       hA12[nn]->GetBinContent(bin),tA10b,
			       hB5[nn]->GetBinContent(bin),tB5,
			       hB4[nn]->GetBinContent(bin),tB5b,
			       hB7[nn]->GetBinContent(bin),tB7,
			       hB9[nn]->GetBinContent(bin),tB7b,er1);
			   
    R2 = Average_Octet_Cts_Bin(hA5w[nn]->GetBinContent(bin),tA5w,
			       hA4w[nn]->GetBinContent(bin),tA5bw,
			       hA7w[nn]->GetBinContent(bin),tA7w,
			       hA9w[nn]->GetBinContent(bin),tA7bw,
			       hB2w[nn]->GetBinContent(bin),tB2w,
			       hB1w[nn]->GetBinContent(bin),tB2bw,
			       hB10w[nn]->GetBinContent(bin),tB10w,
			       hB12w[nn]->GetBinContent(bin),tB10bw,
			       er2);
			   
    R3 = Average_Octet_Cts_Bin(hA2w[nn]->GetBinContent(bin) ,tA2w,
			       hA1w[nn]->GetBinContent(bin) ,tA2bw,
			       hA10w[nn]->GetBinContent(bin),tA10w,
			       hA12w[nn]->GetBinContent(bin),tA10bw,
			       hB5w[nn]->GetBinContent(bin) ,tB5w,
			       hB4w[nn]->GetBinContent(bin) ,tB5bw,
			       hB7w[nn]->GetBinContent(bin) ,tB7w,
			       hB9w[nn]->GetBinContent(bin) ,tB7bw,
			       er3);
			   
    R4 = Average_Octet_Cts_Bin(hA5[nn]->GetBinContent(bin),tA5,
			       hA4[nn]->GetBinContent(bin),tA5b,
			       hA7[nn]->GetBinContent(bin),tA7,
			       hA9[nn]->GetBinContent(bin),tA7b,
			       hB2[nn]->GetBinContent(bin),tB2,
			       hB1[nn]->GetBinContent(bin),tB2b,
			       hB10[nn]->GetBinContent(bin),tB10,
			       hB12[nn]->GetBinContent(bin),tB10b,
			       er4);
			   
    ersum = Power(er1/R1,2) + Power(er2/R2,2) + Power(er3/R3,2) + Power(er4/R4,2);
    RR = Power((R1*R2)/(R3*R4),1./2.);
			   
    if(!(IsNaN(RR))){
      // Average the Errors
      hAsyTot[nn]->SetBinError(bin,RR/(Power((1. + RR),2))*sqrt(ersum));
      // Average the Asymmetries
      hAsyTot[nn]->SetBinContent(bin,(1.-RR)/(1.+RR));
    }
  
  }
 }
}
//=================================================================================================
void Octet::Find_Runs(std::vector<Beta_Run *>bta,std::vector<Bck_Run *>bck,Int_t nMax)
{
  //----------------------------------------------------------------------
  // Sort through the runs to fill the octet runs.  This uses the
  // MySql Database to determine the run type.  Multiple runs of the
  // same type are allowed and handled in the Fill() function for the
  // octet runs.
  //----------------------------------------------------------------------

  for(Int_t i = 0 ; i < nMax ; i++){
   // Check to see if the run is within the bounds of the Octet in question
   if(bta[i]->GetRunNumber() >= First && bta[i]->GetRunNumber() <= Last){
     // Start filling the run type histograms by comparing the octet type
     bta[i]->Load_Histograms(bck[bta[i]->Bkg_index],kTRUE);
     Int_t bkg_index = bta[i]->Bkg_index;
     // Use TObjArrayIter's to loop over the analysis histograms
     TObjArrayIter FgEast(bta[i]->HEastAn);
     TObjArrayIter FgWest(bta[i]->HWestAn);
     TObjArrayIter BkEast(bck[bkg_index]->HEastAn);
     TObjArrayIter BkWest(bck[bkg_index]->HWestAn);
     
     if(!strcmp(bta[i]->octtype,"a2") || !strcmp(bta[i]->octtype,"A2")){
     //------------------------------------------------------------------------
       // Add the run numbers to a vector for the run types
       A2.push_back( bta[i]->GetRunNumber());
       A1.push_back( bta[i]->GetBackgroundRun());
       //loop through the analysis choices and fill the histograms.
       for(Int_t j = 0 ; j < 9 ; j++){
	  Int_t nlast = 7;
	  if(j == 8)nlast = 8;
	  Fill(hA2[j] ,(TH1F*)FgEast.Next(),bta[i]->rtime_e, &tA2,nlast);
	  Fill(hA2w[j],(TH1F*)FgWest.Next(),bta[i]->rtime_w, &tA2w,nlast);
	  Fill(hA1[j] ,(TH1F*)BkEast.Next(),bta[i]->btime_e, &tA2b,nlast);
	  Fill(hA1w[j],(TH1F*)BkWest.Next(),bta[i]->btime_w, &tA2bw,nlast);
       }
       tA2ef = bta[i]->CountTimeEBeta;
       tA2ei = bta[i]->CountTimeEFirstBeta;
       tA2wf = bta[i]->CountTimeWBeta;
       tA2wi = bta[i]->CountTimeWFirstBeta;
       tA2efa = bta[i]->CountTimeEAll;
       tA2eia = bta[i]->CountTimeEFirst;
       tA2wfa = bta[i]->CountTimeWAll;
       tA2wia = bta[i]->CountTimeWFirst;
       for(Int_t j = 0 ; j < 9 ; j++){
	 nA2b[j]  = bta[i]->BkgRtE[j];
	 nA2bw[j] = bta[i]->BkgRtW[j];
       }
       
       Fill_Rad(A2rad,A2rade,bta[i]->erad,bta[i]->erader,bta[i]->rtime_e);
       Fill_Rad(A2wrad,A2wrade,bta[i]->wrad,bta[i]->wrader,bta[i]->rtime_w);
     //------------------------------------------------------------------------------   
     } else if(!strcmp(bta[i]->octtype,"a5") || !strcmp(bta[i]->octtype,"A5")){
       A5.push_back(bta[i]->GetRunNumber());
       A4.push_back(bta[i]->GetBackgroundRun());
       for(Int_t j = 0 ; j < 9 ; j++){
	  Int_t nlast = 7;
	  if(j == 8)nlast = 8;
	  Fill(hA5[j] ,(TH1F*)FgEast.Next(),bta[i]->rtime_e, &tA5,nlast);
	  Fill(hA5w[j],(TH1F*)FgWest.Next(),bta[i]->rtime_w, &tA5w,nlast);
	  Fill(hA4[j] ,(TH1F*)BkEast.Next(),bta[i]->btime_e, &tA5b,nlast);
	  Fill(hA4w[j],(TH1F*)BkWest.Next(),bta[i]->btime_w, &tA5bw,nlast);
       }
   //    cout << "About to fill a type A5" << endl;     
       tA5ef = bta[i]->CountTimeEBeta;
       tA5ei = bta[i]->CountTimeEFirstBeta;
       tA5wf = bta[i]->CountTimeWBeta;
       tA5wi = bta[i]->CountTimeWFirstBeta;
       tA5efa = bta[i]->CountTimeEAll;
       tA5eia = bta[i]->CountTimeEFirst;
       tA5wfa = bta[i]->CountTimeWAll;
       tA5wia = bta[i]->CountTimeWFirst;     

       for(Int_t j = 0 ; j < 9 ; j++){
	 nA5b[j] = bta[i]->BkgRtE[j];
	 nA5bw[j] = bta[i]->BkgRtW[j];
       }
	
       Fill_Rad(A5rad,A5rade,bta[i]->erad,bta[i]->erader,bta[i]->rtime_e);
       Fill_Rad(A5wrad,A5wrade,bta[i]->wrad,bta[i]->wrader,bta[i]->rtime_w);
     //----------------------------------------------------------------------------------  
     }  else if(!strcmp(bta[i]->octtype,"a7") || !strcmp(bta[i]->octtype,"A7")){
     //----------------------------------------------------------------------------------
       A7.push_back(bta[i]->GetRunNumber());
       A9.push_back(bta[i]->GetBackgroundRun());

       tA7ef = bta[i]->CountTimeEBeta;
       tA7ei = bta[i]->CountTimeEFirstBeta;
       tA7wf = bta[i]->CountTimeWBeta;
       tA7wi = bta[i]->CountTimeWFirstBeta;
       tA7efa = bta[i]->CountTimeEAll;
       tA7eia = bta[i]->CountTimeEFirst;
       tA7wfa = bta[i]->CountTimeWAll;
       tA7wia = bta[i]->CountTimeWFirst;
       
       for(Int_t j = 0 ; j < 9 ; j++){
	  Int_t nlast = 7;
	  if(j == 8)nlast = 8;
	  Fill(hA7[j] ,(TH1F*)FgEast.Next(),bta[i]->rtime_e, &tA7,nlast);
	  Fill(hA7w[j],(TH1F*)FgWest.Next(),bta[i]->rtime_w, &tA7w,nlast);
	  Fill(hA9[j] ,(TH1F*)BkEast.Next(),bta[i]->btime_e, &tA7b,nlast);
	  Fill(hA9w[j],(TH1F*)BkWest.Next(),bta[i]->btime_w, &tA7bw,nlast);
       }
        for(Int_t j = 0 ; j < 9 ; j++){
	 nA7b[j] = bta[i]->BkgRtE[j];
	 nA7bw[j] = bta[i]->BkgRtW[j];
       }
       
       Fill_Rad(A7rad,A7rade,bta[i]->erad,bta[i]->erader,bta[i]->rtime_e);
       Fill_Rad(A7wrad,A7wrade,bta[i]->wrad,bta[i]->wrader,bta[i]->rtime_w);
       //---------------------------------------------------------------------------------
     }  else if(strcmp(bta[i]->octtype,"a10")==0 || strcmp(bta[i]->octtype,"A10")==0){
       //---------------------------------------------------------------------------------
       A10.push_back(bta[i]->GetRunNumber());
       A12.push_back(bta[i]->GetBackgroundRun());
        
       for(Int_t j = 0 ; j < 9 ; j++){
	  Int_t nlast = 7;
	  if(j == 8)nlast = 8;
	  Fill(hA10[j] ,(TH1F*)FgEast.Next(),bta[i]->rtime_e, &tA10,nlast);
	  Fill(hA10w[j],(TH1F*)FgWest.Next(),bta[i]->rtime_w, &tA10w,nlast);
	  Fill(hA12[j] ,(TH1F*)BkEast.Next(),bta[i]->btime_e, &tA10b,nlast);
	  Fill(hA12w[j],(TH1F*)BkWest.Next(),bta[i]->btime_w, &tA10bw,nlast);
       }
       
       tA10ef = bta[i]->CountTimeEBeta;
       tA10ei = bta[i]->CountTimeEFirstBeta;
       tA10wf = bta[i]->CountTimeWBeta;
       tA10wi = bta[i]->CountTimeWFirstBeta;
       tA10efa = bta[i]->CountTimeEAll;
       tA10eia = bta[i]->CountTimeEFirst;
       tA10wfa = bta[i]->CountTimeWAll;
       tA10wia = bta[i]->CountTimeWFirst;

       for(Int_t j = 0 ; j < 9 ; j++){
	 nA10b[j] = bta[i]->BkgRtE[j];
	 nA10bw[j] = bta[i]->BkgRtW[j];
       }
       
       Fill_Rad(A10rad,A10rade,bta[i]->erad,bta[i]->erader,bta[i]->rtime_e);
       Fill_Rad(A10wrad,A10wrade,bta[i]->wrad,bta[i]->wrader,bta[i]->rtime_w);
       //----------------------------------------------------------------------------
     }  else if(!strcmp(bta[i]->octtype,"b2") || !strcmp(bta[i]->octtype,"B2")){
       //----------------------------------------------------------------------------
       B2.push_back(bta[i]->GetRunNumber());
       B1.push_back(bta[i]->GetBackgroundRun());
        for(Int_t j = 0 ; j < 9 ; j++){
	  Int_t nlast = 7;
	  if(j == 8)nlast = 8;
	  Fill(hB2[j] ,(TH1F*)FgEast.Next(),bta[i]->rtime_e, &tB2,nlast);
	  Fill(hB2w[j],(TH1F*)FgWest.Next(),bta[i]->rtime_w, &tB2w,nlast);
	  Fill(hB1[j] ,(TH1F*)BkEast.Next(),bta[i]->btime_e, &tB2b,nlast);
	  Fill(hB1w[j],(TH1F*)BkWest.Next(),bta[i]->btime_w, &tB2bw,nlast);
       }
     
       tB2ef = bta[i]->CountTimeEBeta;
       tB2ei = bta[i]->CountTimeEFirstBeta;
       tB2wf = bta[i]->CountTimeWBeta;
       tB2wi = bta[i]->CountTimeWFirstBeta;
       tB2efa = bta[i]->CountTimeEAll;
       tB2eia = bta[i]->CountTimeEFirst;
       tB2wfa = bta[i]->CountTimeWAll;
       tB2wia = bta[i]->CountTimeWFirst;

        for(Int_t j = 0 ; j < 9 ; j++){
	 nB2b[j]  = bta[i]->BkgRtE[j];
	 nB2bw[j] = bta[i]->BkgRtW[j];
       }
       
       Fill_Rad(B2rad,B2rade,bta[i]->erad,bta[i]->erader,bta[i]->rtime_e);
       Fill_Rad(B2wrad,B2wrade,bta[i]->wrad,bta[i]->wrader,bta[i]->rtime_w);
       
     }  else if(!strcmp(bta[i]->octtype,"b5") || !strcmp(bta[i]->octtype,"B5")){
       B5.push_back(bta[i]->GetRunNumber());
       B4.push_back(bta[i]->GetBackgroundRun());
       
       for(Int_t j = 0 ; j < 9 ; j++){
	  Int_t nlast = 7;
	  if(j == 8)nlast = 8;
	  Fill(hB5[j] ,(TH1F*)FgEast.Next(),bta[i]->rtime_e, &tB5,nlast);
	  Fill(hB5w[j],(TH1F*)FgWest.Next(),bta[i]->rtime_w, &tB5w,nlast);
	  Fill(hB4[j] ,(TH1F*)BkEast.Next(),bta[i]->btime_e, &tB5b,nlast);
	  Fill(hB4w[j],(TH1F*)BkWest.Next(),bta[i]->btime_w, &tB5bw,nlast);
       }
 
       tB5ef = bta[i]->CountTimeEBeta;
       tB5ei = bta[i]->CountTimeEFirstBeta;
       tB5wf = bta[i]->CountTimeWBeta;
       tB5wi = bta[i]->CountTimeWFirstBeta;
       tB5efa = bta[i]->CountTimeEAll;
       tB5eia = bta[i]->CountTimeEFirst;
       tB5wfa = bta[i]->CountTimeWAll;
       tB5wia = bta[i]->CountTimeWFirst;

       for(Int_t j = 0 ; j < 9 ; j++){
	 nB5b[j]  = bta[i]->BkgRtE[j];
	 nB5bw[j] = bta[i]->BkgRtW[j];
       }
       
       Fill_Rad(B5rad,B5rade,bta[i]->erad,bta[i]->erader,bta[i]->rtime_e);
       Fill_Rad(B5wrad,B5wrade,bta[i]->wrad,bta[i]->wrader,bta[i]->rtime_w);
       
     }  else if(!strcmp(bta[i]->octtype,"b7") || !strcmp(bta[i]->octtype,"B7")){
       B7.push_back(bta[i]->GetRunNumber());
       B9.push_back(bta[i]->GetBackgroundRun());
       for(Int_t j = 0 ; j < 9 ; j++){
	  Int_t nlast = 7;
	  if(j == 8)nlast = 8;
	  Fill(hB7[j] ,(TH1F*)FgEast.Next(),bta[i]->rtime_e, &tB7,nlast);
	  Fill(hB7w[j],(TH1F*)FgWest.Next(),bta[i]->rtime_w, &tB7w,nlast);
	  Fill(hB9[j] ,(TH1F*)BkEast.Next(),bta[i]->btime_e, &tB7b,nlast);
	  Fill(hB9w[j],(TH1F*)BkWest.Next(),bta[i]->btime_w, &tB7bw,nlast);
       }
       
       tB7ef = bta[i]->CountTimeEBeta;
       tB7ei = bta[i]->CountTimeEFirstBeta;
       tB7wf = bta[i]->CountTimeWBeta;
       tB7wi = bta[i]->CountTimeWFirstBeta;
       tB7efa = bta[i]->CountTimeEAll;
       tB7eia = bta[i]->CountTimeEFirst;
       tB7wfa = bta[i]->CountTimeWAll;
       tB7wia = bta[i]->CountTimeWFirst;

       for(Int_t j = 0 ; j < 9 ; j++){
	 nB7b[j]  = bta[i]->BkgRtE[j];
	 nB7bw[j] = bta[i]->BkgRtW[j];
       }
       
       Fill_Rad(B7rad,B7rade,bta[i]->erad,bta[i]->erader,bta[i]->rtime_e);
       Fill_Rad(B7wrad,B7wrade,bta[i]->wrad,bta[i]->wrader,bta[i]->rtime_w);

     }  else if(!strcmp(bta[i]->octtype,"b10") || !strcmp(bta[i]->octtype,"B10")){

       B10.push_back( bta[i]->GetRunNumber());
       B12.push_back( bta[i]->GetBackgroundRun());
       
       for(Int_t j = 0 ; j < 9 ; j++){
	  Int_t nlast = 7;
	  if(j == 8)nlast = 8;
	  Fill(hB10[j] ,(TH1F*)FgEast.Next(),bta[i]->rtime_e, &tB10,nlast);
	  Fill(hB10w[j],(TH1F*)FgWest.Next(),bta[i]->rtime_w, &tB10w,nlast);
	  Fill(hB12[j] ,(TH1F*)BkEast.Next(),bta[i]->btime_e, &tB10b,nlast);
	  Fill(hB12w[j],(TH1F*)BkWest.Next(),bta[i]->btime_w, &tB10bw,nlast);
       }
    
       tB10ef = bta[i]->CountTimeEBeta;
       tB10ei = bta[i]->CountTimeEFirstBeta;
       tB10wf = bta[i]->CountTimeWBeta;
       tB10wi = bta[i]->CountTimeWFirstBeta;
       tB10efa = bta[i]->CountTimeEAll;
       tB10eia = bta[i]->CountTimeEFirst;
       tB10wfa = bta[i]->CountTimeWAll;
       tB10wia = bta[i]->CountTimeWFirst;

       for(Int_t j = 0 ; j < 9 ; j++){
	 nB10b[j]  = bta[i]->BkgRtE[j];
	 nB10bw[j] = bta[i]->BkgRtW[j];
       }
       
       Fill_Rad(B10rad,B10rade,bta[i]->erad,bta[i]->erader,bta[i]->rtime_e);
       Fill_Rad(B10wrad,B10wrade,bta[i]->wrad,bta[i]->wrader,bta[i]->rtime_w);
     }
     bta[i]->Remove_Histograms(bck[bta[i]->Bkg_index]);
   }
  }
  
  QuartetAFirst = ((int)A1.size() > 0)  ? A1[0]  : 0; 
  QuartetALast  = ((int)A12.size() >0)  ? A12[0] : 0;
  QuartetBFirst = ((int)B1.size() > 0)  ? B1[0]  : 0;
  QuartetBLast  = ((int)B12.size() >0)  ? B12[0] : 0;
  SuperA1First  = ((int)A1.size() > 0)  ? A1[0]  : 0;
  SuperA1Last   = ((int)A5.size() > 0)  ? A5[0]  : 0;
  SuperA2First  = ((int)A7.size() > 0)  ? A7[0]  : 0;
  SuperA2Last   = ((int)A12.size()> 0)  ? A12[0] : 0;
  SuperB1First  = ((int)B1.size() > 0)  ? B1[0]  : 0;
  SuperB1Last   = ((int)B5.size() > 0)  ? B5[0]  : 0;
  SuperB2First  = ((int)B7.size() > 0)  ? B7[0]  : 0;
  SuperB2Last   = ((int)B12.size() > 0) ? B12[0] : 0;
  
  // Simple check to avoid t == 0 for error calculations
  tA2   = ( tA2  !=0) ? tA2  : 1;
  tA2w  = ( tA2w !=0) ? tA2w : 1;
  tA5   = ( tA5  !=0) ? tA5  : 1;
  tA5w  = ( tA5w !=0) ? tA5w : 1;
  tA7   = ( tA7  !=0) ? tA7  : 1;
  tA7w  = ( tA7w !=0) ? tA7w : 1;
  tA10  = ( tA10  !=0) ? tA10  : 1;
  tA10w = ( tA10w !=0) ? tA10w : 1;
  tB2   = ( tB2  !=0) ? tB2  : 1;
  tB2w  = ( tB2w !=0) ? tB2w : 1;
  tB5   = ( tB5  !=0) ? tB5  : 1;
  tB5w  = ( tB5w !=0) ? tB5w : 1;
  tB7   = ( tB7  !=0) ? tB7  : 1;
  tB7w  = ( tB7w !=0) ? tB7w : 1;
  tB10  = ( tB10  !=0) ? tB10  : 1;
  tB10w = ( tB10w !=0) ? tB10w : 1;
  
  GetIntegralCounts();
}
//-----------------------------------------------------------------------------------------
void Octet::Fill(TH1F *h,TH1F *h2,Double_t t1,Double_t* t2,Int_t inc)
{

  // Fill the Octet Type Histograms. If there are multiple runs counts
  // from the extra runs are added in. 
  using namespace TMath;
  
  Double_t t2val = *t2;
  // Check if the Octet histogram is empty.  
  if(h->Integral() == 0.){
    // If empty fill with the run histogram
    for(int i = 1 ; i<= h->GetNbinsX() ; i++){
      h->SetBinContent(i,h2->GetBinContent(i));
      h->SetBinError(i,h2->GetBinError(i));
    }
    t2val = t1;
  } else {
    // If not-empty then added the runs together.
    for(int i = 1 ; i<= h->GetNbinsX() ; i++){
      // t1 is the individual run time
      // t2val is the holder for t2 which is the total octet type run time. 
      h->SetBinContent(i,(h2->GetBinContent(i)*t1+h->GetBinContent(i)*t2val)/(t1+t2val));
      h->SetBinError(i,sqrt(Power(h->GetBinError(i),2)+Power(h2->GetBinError(i),2)));
    }
     if(inc == 8) t2val += t1;
  }
  *t2 = t2val;
  
}

void Octet::Initialize_Histo(Int_t n)
{
  
  // Energy spectra for type A analysis of Octets
  Int_t nxbins = 200;
  
  for(Int_t i = 0 ; i < 20 ; i++){
  
    hA2[i]  = new TH1F(Form("octet[%d]->hA2[%d]" ,n,i),"A2 spectrum",nxbins,0,2000);
    hA5[i]  = new TH1F(Form("octet[%d]->hA5[%d]" ,n,i),"A5 spectrum",nxbins,0,2000);
    hA7[i]  = new TH1F(Form("octet[%d]->hA7[%d]" ,n,i),"A7 spectrum",nxbins,0,2000);
    hA10[i] = new TH1F(Form("octet[%d]->hA10[%d]",n,i),"A10 spectrum",nxbins,0,2000);
    hB2[i]  = new TH1F(Form("octet[%d]->hB2[%d]" ,n,i),"B2 spectrum",nxbins,0,2000);
    hB5[i]  = new TH1F(Form("octet[%d]->hB5[%d]" ,n,i),"B5 spectrum",nxbins,0,2000);
    hB7[i]  = new TH1F(Form("octet[%d]->hB7[%d]" ,n,i),"B7 spectrum",nxbins,0,2000);
    hB10[i] = new TH1F(Form("octet[%d]->hB10[%d]",n,i),"B10 spectrum",nxbins,0,2000);
  
    hA2w[i]   = new TH1F(Form("octet[%d]->hA2w[%d]" ,n,i),"A2w spectrum",nxbins,0,2000);
    hA5w[i]  = new TH1F(Form("octet[%d]->hA5w[%d]" ,n,i),"A5w spectrum",nxbins,0,2000);
    hA7w[i]  = new TH1F(Form("octet[%d]->hA7w[%d]" ,n,i),"A7w spectrum",nxbins,0,2000);
    hA10w[i] = new TH1F(Form("octet[%d]->hA10w[%d]",n,i),"A10w spectrum",nxbins,0,2000);
    hB2w[i]  = new TH1F(Form("octet[%d]->hB2w[%d]" ,n,i),"B2w spectrum",nxbins,0,2000);
    hB5w[i]  = new TH1F(Form("octet[%d]->hB5w[%d]" ,n,i),"B5w spectrum",nxbins,0,2000);
    hB7w[i]  = new TH1F(Form("octet[%d]->hB7w[%d]" ,n,i),"B7w spectrum",nxbins,0,2000);
    hB10w[i] = new TH1F(Form("octet[%d]->hB10w[%d]",n,i),"B10w spectrum",nxbins,0,2000);
    
    hA1[i]  = new TH1F(Form("octet[%d]->hA1[%d]" ,n,i),"A2 spectrum",nxbins,0,2000);
    hA4[i]  = new TH1F(Form("octet[%d]->hA4[%d]" ,n,i),"A5 spectrum",nxbins,0,2000);
    hA9[i]  = new TH1F(Form("octet[%d]->hA9[%d]" ,n,i),"A7 spectrum",nxbins,0,2000);
    hA12[i] = new TH1F(Form("octet[%d]->hA12[%d]",n,i),"A10 spectrum",nxbins,0,2000);
    hB1[i]  = new TH1F(Form("octet[%d]->hB1[%d]" ,n,i),"B2 spectrum",nxbins,0,2000);
    hB4[i]  = new TH1F(Form("octet[%d]->hB4[%d]" ,n,i),"B5 spectrum",nxbins,0,2000);
    hB9[i]  = new TH1F(Form("octet[%d]->hB9[%d]" ,n,i),"B7 spectrum",nxbins,0,2000);
    hB12[i] = new TH1F(Form("octet[%d]->hB12[%d]",n,i),"B10 spectrum",nxbins,0,2000);
  
    hA1w[i]   = new TH1F(Form("octet[%d]->hA1w[%d]" ,n,i),"A2w spectrum",nxbins,0,2000);
    hA4w[i]  = new TH1F(Form("octet[%d]->hA4w[%d]" ,n,i),"A5w spectrum",nxbins,0,2000);
    hA9w[i]  = new TH1F(Form("octet[%d]->hA9w[%d]" ,n,i),"A7w spectrum",nxbins,0,2000);
    hA12w[i] = new TH1F(Form("octet[%d]->hA12w[%d]",n,i),"A10w spectrum",nxbins,0,2000);
    hB1w[i]  = new TH1F(Form("octet[%d]->hB1w[%d]" ,n,i),"B2w spectrum",nxbins,0,2000);
    hB4w[i]  = new TH1F(Form("octet[%d]->hB4w[%d]" ,n,i),"B5w spectrum",nxbins,0,2000);
    hB9w[i]  = new TH1F(Form("octet[%d]->hB9w[%d]" ,n,i),"B7w spectrum",nxbins,0,2000);
    hB12w[i] = new TH1F(Form("octet[%d]->hB12w[%d]",n,i),"B10w spectrum",nxbins,0,2000);

    hAsyA[i]   = new TH1F(Form("octet[%d]->hAsyA[%d]"  ,n,i),"A Octet Asymmetry"  ,nxbins,0,2000);
    hAsyB[i]   = new TH1F(Form("octet[%d]->hAsyB[%d]"  ,n,i),"B Octet Asymmetry"  ,nxbins,0,2000);
    hAsyTot[i] = new TH1F(Form("octet[%d]->hAsyTot[%d]",n,i),"Sum Octet Asymmetry",nxbins,0,2000);
 }
  // Zero the rad counters.
  
  for(Int_t i = 0; i < 12 ; i++){
    A2rad[i]    = 0.;
    A2wrad[i]   = 0.;
    A5rad[i]    = 0.;
    A5wrad[i]   = 0.;
    A7rad[i]    = 0.;
    A7wrad[i]   = 0.;
    A10rad[i]   = 0.;
    A10wrad[i ] = 0.;
    B2rad[i]    = 0.;
    B2wrad[i]   = 0.;
    B5rad[i]    = 0.;
    B5wrad[i]   = 0.;
    B7rad[i]    = 0.;
    B7wrad[i]   = 0.;
    B10rad[i]   = 0.;
    B10wrad[i]  = 0.;
    A2rade[i]   = 0.;
    A2wrade[i]  = 0.;
    A5rade[i]   = 0.;
    A5wrade[i]  = 0.;
    A7rade[i]   = 0.;
    A7wrade[i]  = 0.;
    A10rade[i]  = 0.;
    A10wrade[i] = 0.;
    B2rade[i]   = 0.;
    B2wrade[i]  = 0.;
    B5rade[i]   = 0.;
    B5wrade[i]  = 0.;
    B7rade[i]   = 0.;
    B7wrade[i]  = 0.;
    B10rade[i]  = 0.;
    B10wrade[i] = 0.;
  }
  
  Load_Background();
  
  
}

void Octet::Fill_Rad(Double_t *x1,Double_t *x1e,Float_t *x2,Float_t *x2e,Float_t t1)
{

  using namespace TMath;
  
  for(Int_t j = 0 ; j < 12 ; j++){
    x1[j]  += x2[j]*t1;
    x1e[j]  = (!IsNaN(x2e[j])) ? x2e[j] : 0; 
  }
 
}

void Octet::GetIntegralCounts()
{
  // Get the Integral Numbers for the A quartet
  
  // For the time being I've hard coded Jianglai's rates from octet 8278-8310 into then
  // analysis
  for(Int_t i = 0 ; i < nchoices ; i++){
    nA2[i]  = hA2[i]->Integral(nlow,nhigh);   nA5[i]   = hA5[i]->Integral(nlow,nhigh);
    nA7[i]  = hA7[i]->Integral(nlow,nhigh);   nA10[i]  = hA10[i]->Integral(nlow,nhigh);
    nA2w[i] = hA2w[i]->Integral(nlow,nhigh);  nA5w[i]  = hA5w[i]->Integral(nlow,nhigh);
    nA7w[i] = hA7w[i]->Integral(nlow,nhigh);  nA10w[i] = hA10w[i]->Integral(nlow,nhigh);  
  // Get the Integral Numbers for the A quartet
    nB2[i]  = hB2[i]->Integral(nlow,nhigh);   nB5[i]   = hB5[i]->Integral(nlow,nhigh);
    nB7[i]  = hB7[i]->Integral(nlow,nhigh);   nB10[i]  = hB10[i]->Integral(nlow,nhigh);
    nB2w[i] = hB2w[i]->Integral(nlow,nhigh);  nB5w[i]  = hB5w[i]->Integral(nlow,nhigh);
    nB7w[i] = hB7w[i]->Integral(nlow,nhigh);  nB10w[i] = hB10w[i]->Integral(nlow,nhigh);
    
  }

  TotCounts = nA2[0]*tA2  + nA5[0]*tA5  + nA7[0]*tA7  + nA10[0]*tA10
            +nA2w[0]*tA2w +nA5w[0]*tA5w +nA7w[0]*tA7w +nA10w[0]*tA10w
            + nB2[0]*tB2  + nB5[0]*tB5  + nB7[0]*tB7  + nB10[0]*tB10
            +nB2w[0]*tB2w +nB5w[0]*tB5w +nB7w[0]*tB7w +nB10w[0]*tB10w;
};
 
Double_t Octet::ErrorMultiRate(Double_t n1,Double_t n2,Double_t t1,Double_t t2)
{
  using namespace TMath;
  
  Double_t er = (n1+n2)/t1 + (n2/t2);
  er = sqrt(er);
  
  return er;
}
  
Double_t Octet::ErrorSum(Double_t n1,Double_t n2,Double_t t1,Double_t t2,
			 Double_t n3,Double_t n4,Double_t t3,Double_t t4,Int_t side)
{
  
  //Very simple intermediate step calculator for the super-ratio errors....
 using namespace TMath;
 vector<Double_t>rate,bck,t,tb;
 
 Double_t sig= 0.;
  if(n1>0){
    rate.push_back(n1);
    bck.push_back(n3);
    t.push_back(t1);
    tb.push_back(t3);
  }
  if(n2>0){
    rate.push_back(n2);
    bck.push_back(n4);
    t.push_back(t2);
    tb.push_back(t4);
  }

  
  for(Int_t i = 0 ; i < (int)rate.size() ; i++){
    sig += 1./((rate[i] + bck[i])/t[i] + bck[i]/tb[i]);
  };
  
  Double_t er = 1./TMath::Sqrt(sig);
   
  return er;
}

Double_t Octet::ErrorSum(Double_t n1,Double_t n2,Double_t n3,Double_t n4,
			 Double_t t1,Double_t t2,Double_t t3,Double_t t4,
			 Double_t n5,Double_t n6,Double_t n7,Double_t n8,
			 Double_t t5,Double_t t6,Double_t t7,Double_t t8,Int_t side)
{
  //Very simple intermediate step calculator for the super-ratio errors....
  using namespace TMath;
   
  
  Double_t er = 0.,sig=0.;
  vector<Double_t>rate,bck,t,tb;
  
  if(n1 > 0){
    rate.push_back(n1);
    bck.push_back(n5);
    t.push_back(t1);
    tb.push_back(t5);
  }
  if(n2 > 0){
    rate.push_back(n2);
    bck.push_back(n6);
    t.push_back(t2);
    tb.push_back(t6);
  }
  if(n3 > 0){
    rate.push_back(n3);
    bck.push_back(n7);
    t.push_back(t3);
    tb.push_back(t7);
  }
  if(n4 > 0){
    rate.push_back(n4);
    bck.push_back(n8);
    t.push_back(t4);
    tb.push_back(t8);
  }
  
  for(Int_t i = 0 ; i < (int)rate.size() ; i++){
    sig += 1./((rate[i] + bck[i])/t[i] + bck[i]/tb[i]);
  };
  
  er = 1./sqrt(sig);
  
  return er;
};
 
 void Octet::Debugger(Int_t nOct)
{
 
  Int_t ach = 2;  
  fstream fdebug;
  Int_t n = 0 ;
  fdebug.open(Form("output_files/debuging_octet_%d.txt",GetFirst()),fstream::out);
  
  fdebug << "Results for Octet starting with run " << GetFirst() << " and ending with " << GetLast() << endl;
  fdebug << "========================================================================================================="<<endl;
  if(A1.size() >0 ){
  fdebug << "A1 bkg run is = ";
  //for(Int_t n = 0 ; n < 5 ; n++)if(A1[n]!=0)
  fdebug << A1[n] << ",";
  fdebug <<  " East Rate in A1 = " << nA2b[ach] << " West Rate in A1 = " << nA2bw[ach] << endl;
  fdebug << "                           Run time is te = " << tA2b << "  tw = " << tA2bw << endl;

  fdebug << "A2 beta run is = ";
  //for(Int_t n = 0 ; n < 5 ; n++)if(A2[n]!=0)
  fdebug << A2[n] << ",";
  fdebug << " East Rate in A2 = " << nA2[ach] << " West Rate in A2 = " << nA2w[ach] << endl;
  fdebug << "                           Run time is te = " << tA2 << "  tw = " << tA2w << endl;

  fdebug << "A4 bkg run is = ";
  //for(Int_t n = 0 ; n < 5 ; n++)if(A4[n]!=0)
  fdebug << A4[n] << ",";
  fdebug << " East Rate in A4 = " << nA5b[ach] << " West Rate in A4 = " << nA5bw[ach] << endl;
  fdebug << "                           Run time is te = " << tA5b << "  tw = " << tA5bw << endl;

  fdebug << "A5 beta run is = ";
  //for(Int_t n = 0 ; n < 5 ; n++)if(A5[n]!=0)
  fdebug << A5[n] << ",";
  fdebug << " East Rate in A5 = " << nA5[ach] << " West Rate in A5 = " << nA5w[ach] << endl;
  fdebug << "                           Run time is te = " << tA5 << "  tw = " << tA5w << endl;
  }
  if(A9.size() >0){
  fdebug << "A9 bkg run is = ";  
  // for(Int_t n = 0 ; n < 5 ; n++)if(A9[n]!=0)
  fdebug << A9[n] << ",";
  fdebug << " East Rate in A9 = " << nA7b[ach] << " West Rate in A9 = " << nA7bw[ach] << endl;
  fdebug << "                           Run time is te = " << tA7b << "  tw = " << tA7bw << endl;

  fdebug << "A7 beta run is = "; 
  //for(Int_t n = 0 ; n < 5 ; n++)if(A7[n]!=0)
  fdebug << A7[n] << ",";
  fdebug << " East Rate in A7 = " << nA7[ach] << " West Rate in A7 = " << nA7w[ach] << endl;
  fdebug << "                           Run time is te = " << tA7 << "  tw = " << tA7w << endl;

  fdebug << "A12 bkg run is = "; 
  //for(Int_t n = 0 ; n < 5 ; n++)if(A12[n]!=0)
  fdebug << A12[n] << ",";
  fdebug << " East Rate in A12 = " << nA10b[ach] << " West Rate in A12 = " << nA10bw[ach] << endl;
  fdebug << "                           Run time is te = " << tA10b << "  tw = " << tA10bw << endl;
  
  fdebug << "A10 beta run is = ";
  //for(Int_t n = 0 ; n < 5 ; n++)if(A10[n]!=0)
  fdebug << A10[n] << ",";
  fdebug << " East Rate in A10 = " << nA10[ach] << " West Rate in A10 = " << nA10w[ach] << endl;
  fdebug << "                           Run time is te = " << tA10 << "  tw = " << tA10w << endl;
  }
  if(B1.size() > 0){
  fdebug << "B1 bkg run is = ";
  //for(Int_t n = 0 ; n < 5 ; n++)if(B1[n]!=0)
  fdebug << B1[n] << ",";
  fdebug << " East Rate in B1 = " << nB2b[ach] << " West Rate in B1 = " << nB2bw[ach] << endl;
  fdebug << "                           Run time is te = " << tB2b << "  tw = " << tB2bw << endl;

  fdebug << "B2 beta run is = "; 
  //for(Int_t n = 0 ; n < 5 ; n++)if(B2[n]!=0)
  fdebug << B2[n] << ",";
  fdebug << " East Rate in B2 = " << nB2[ach] << " West Rate in B2 = " << nB2w[ach] << endl;
  fdebug << "                           Run time is te = " << tB2 << "  tw = " << tB2w << endl;
  
  fdebug << "B4 bkg run is = "; 
  //for(Int_t n = 0 ; n < 5 ; n++)if(B4[n]!=0)
  fdebug << B4[n] << ",";
  fdebug << " East Rate in B4 = " << nB5b[ach] << " West Rate in B4 = " << nB5bw[ach] << endl;
  fdebug << "                           Run time is te = " << tB5b << "  tw = " << tB5bw << endl;

  fdebug << "B5 beta run is = ";
  //for(Int_t n = 0 ; n < 5 ; n++)if(B5[n]!=0)
  fdebug << B5[n] << ",";
  fdebug << " East Rate in B5 = " << nB5[ach] << " West Rate in B5 = " << nB5w[ach] << endl;
  fdebug << "                           Run time is te = " << tB5 << "  tw = " << tB5w << endl;
  }
  if(B9.size() > 0){
  fdebug << "B9 bkg run is = ";
  //for(Int_t n = 0 ; n < 5 ; n++)if(B9[n]!=0)
  fdebug << B9[n] << ",";
  fdebug << " East Rate in B9 = " << nB7b[ach] << " West Rate in B9 = " << nB7bw[ach] << endl;
  fdebug << "                           Run time is te = " << tB7b << "  tw = " << tB7bw << endl;

  fdebug << "B7 beta run is = "; 
  //for(Int_t n = 0 ; n < 5 ; n++)if(B7[n]!=0)
  fdebug << B7[n] << ",";
  fdebug << " East Rate in B7 = " << nB7[ach] << " West Rate in B7 = " << nB7w[ach] << endl;
  fdebug << "                           Run time is te = " << tB7 << "  tw = " << tB7w << endl;

  fdebug << "B12 bkg run is = ";
  //for(Int_t n = 0 ; n < 5 ; n++)if(B12[n]!=0)
  fdebug << B12[n] << ","; 
  fdebug << " East Rate in B12 = " << nB10b[ach] << " West Rate in B12 = " << nB10bw[ach] << endl;
  fdebug << "                           Run time is te = " << tB10b << "  tw = " << tB10bw << endl;

  fdebug << "B10 beta run is = ";
  //for(Int_t n = 0 ; n < 5 ; n++)if(B10[n]!=0)
  fdebug << B10[n] << ",";
  fdebug << " East Rate in B10 = " << nB10[ach] << " West Rate in B10 = " << nB10w[ach] << endl;
  fdebug << "                           Run time is te = " << tB10 << "  tw = " << tB10w << endl;
     }
  fdebug << "========================================================================================================" << endl;
  fdebug << "Super Ratios " << endl;
  fdebug << "A1-A6            A = " << Asuper1[ach] << " +/- " << Asuper1er[ach] << endl;
  fdebug << "A7-A12           A = " << Asuper2[ach] << " +/- " << Asuper2er[ach] << endl;
  fdebug << "B1-B6            A = " << Bsuper1[ach] << " +/- " << Bsuper1er[ach] << endl;
  fdebug << "B7-B12           A = " << Bsuper2[ach] << " +/- " << Bsuper2er[ach] << endl << endl;
  fdebug << "Quartets " << endl;
  fdebug << "A1-A12  Product  A = " << A_multi_A[ach] << " +/- " << A_multi_Aer[ach] << endl;
  fdebug << "A1-A12  Sum      A = " << A_sum_A[ach] << " +/- " << A_sum_Aer[ach] << endl;
  fdebug << "B1-B12  Product  A = " << A_multi_B[ach] << " +/- " << A_multi_Ber[ach] << endl;
  fdebug << "B1-B12  Sum      A = " << A_sum_B[ach] << " +/- " << A_sum_Ber[ach] << endl << endl;
  fdebug << "Octets " << endl;
  fdebug << "Product          A = " << A_multi[ach] << " +/- " << A_multier[ach] << endl;
  fdebug << "Sum              A = " << A_sum[ach] << " +/- " << A_sumer[ach] << endl;
  fdebug << "========================================================================================================" << endl;
  
  
};
 //=====================================================
Double_t Octet::Average_Quartet_Bin(Double_t r1,Double_t e1,Double_t r2,Double_t e2,Double_t &ett)
{
  Double_t ave = 0.;
  ett = 0.;
  
  if(r1 != 0 && r2 != 0 ){
    ett = 1./sqrt(1./(e1*e1) + 1./(e2*e2));
    ave = (r1/(e1*e1) + r2/(e2*e2))*(ett*ett);
  } else if( r1 !=0 && r2 == 0){
    ave = r1;
    ett = e1;
  } else if( r1 == 0 && r2 != 0){
    ave = r2;
    ett = e2;
  } else {
    ett = 1.;
  }
  
  return ave;
};
//======================================================================================

Double_t Octet::Average_Octet_Bin(Double_t r1,Double_t e1,Double_t r2,Double_t e2,
				  Double_t r3,Double_t e3,Double_t r4,Double_t e4, Double_t &ett)
{
  Double_t ave = 0.;
  Double_t r1e = 1.;
  Double_t r2e = 1.;
  Double_t r3e = 1.; 
  Double_t r4e = 1.;
  
  ett = 0.;
  
  if(r1 == 0){
    r1e = 0.;
    e1  = 1.;
  }
  if(r2 == 0){
    r2e = 0.;
    e2  = 1.;
  }
  if(r3 == 0){
    r3e = 0.;
    e3  = 1.;
  }
  if(r4 == 0){
    r4e = 0.;
    e4  = 1.;
  }
  
  ett = 1./sqrt(r1e/(e1*e1) + r2e/(e2*e2) + r3e/(e3*e3) + r4e/(e4*e4) );
  ave = (r1/(e1*e1) + r2/(e2*e2) + r3/(e3*e3) + r4/(e4*e4) )*(ett*ett);
 
  return ave; 
};

Double_t Octet::Average_Octet_Bin2(Double_t r1,Double_t t1,Double_t b1,Double_t bt1,
				   Double_t r2,Double_t t2,Double_t b2,Double_t bt2, 
				   Double_t r3,Double_t t3,Double_t b3,Double_t bt3,
				   Double_t r4,Double_t t4,Double_t b4,Double_t bt4, 
				   Double_t &ett)
{
 
  
  Double_t ave = 0.;
  Double_t sig = 0.;
  
  vector<Double_t>rate,bck,t,tb;
  
  if(r1>0){
    rate.push_back(r1);
    bck.push_back(b1);
    t.push_back(t1);
    tb.push_back(bt1);
  }
  if(r2>0){
    rate.push_back(r2);
    bck.push_back(b2);
    t.push_back(t2);
    tb.push_back(bt2);
  }
    if(r3>0){
    rate.push_back(r3);
    bck.push_back(b3);
    t.push_back(t3);
    tb.push_back(bt3);
  }
    if(r4>0){
    rate.push_back(r4);
    bck.push_back(b4);
    t.push_back(t4);
    tb.push_back(bt4);
  }
  
  for(Int_t i = 0 ; i < (int)rate.size() ; i++){
     ave += rate[i]/((rate[i] + bck[i])/t[i] + bck[i]/tb[i]);
     sig += 1./((rate[i] + bck[i])/t[i] + bck[i]/tb[i]);
  }
  
  ave = ave /sig;
  ett = 1./sqrt(sig);
  return ave;
};



Double_t Octet::Average_Octet_Cts_Bin(Double_t r1,Double_t t1,Double_t b1,Double_t bt1,
				      Double_t r2,Double_t t2,Double_t b2,Double_t bt2, 
				      Double_t r3,Double_t t3,Double_t b3,Double_t bt3,
				      Double_t r4,Double_t t4,Double_t b4,Double_t bt4, 
				      Double_t &ett)
{
  //------------------------------------------------------------------
  // this function combine the measured rates in a direct sum
  // instead of the weighted average.
  // the arguements are :
  // r1,r2,r3,r4 -> Signal Run rates
  // t1,t2,t3,t4 -> Signal Run times
  // b1,b2,b3,b4 -> Background rates
  // bt1,bt2,bt3,bt4 -> background times
  //
  Double_t rate  = r1*t1 + r2*t2 + r3*t3 + r4*t4;
  Double_t bsub  = b1*t1 + b2*t2 + b3*t3 + b4*t4;
  Double_t time  = t1 + t2 + t3 + t4;
  Double_t btime = bt1 + bt2 + bt3 + bt4;
  Double_t brate = 0;
  
  if(bt1 > 0)    brate += b1*bt1;
  if(bt2 > 0)    brate += b2*bt2;
  if(bt3 > 0)    brate += b3*bt3;
  if(bt4 > 0)    brate += b4*bt4;
 
  ett = sqrt(rate +bsub + brate*TMath::Power(time/btime,2))/time;
  
  return rate/time; 
};

 
Double_t Octet::Average_Quartet_Rate(Double_t r1,Double_t t1,Double_t b1,Double_t tb1,
				     Double_t r2,Double_t t2,Double_t b2,Double_t tb2,Int_t side)
{
  
  
  // average the rates for an entry in the quartet super ratio
  // if t1 or t2 is 0 then that run is missing so return the rate
  // of the existing run as the correct average.
  // if neither exists return 0.
  
  // r1,r2   - background subtracted beta rates
  // b1,b2   - background rates
  // t1,t2   - live times for the beta runs
  // tb1,tb2 - live times for the background runs 
  
  Double_t ave = 0.,sig=0.;
  
  vector<Double_t>rate,bck,t,tb;
  if(r1>0){
    rate.push_back(r1);
    bck.push_back(b1);
    t.push_back(t1);
    tb.push_back(tb1);
  }
  if(r2>0){
    rate.push_back(r2);
    bck.push_back(b2);
    t.push_back(t2);
    tb.push_back(tb2);
  }
  for(Int_t i = 0 ; i < (int)rate.size() ; i++){
    ave += rate[i]*t[i];
    sig += t[i];
  };
  
  ave = ave /sig;
  
  return ave;
}

Double_t Octet::Average_Octet_Rate(Double_t r1,Double_t t1,Double_t b1,Double_t tb1,
				   Double_t r2,Double_t t2,Double_t b2,Double_t tb2,
				   Double_t r3,Double_t t3,Double_t b3,Double_t tb3,
				   Double_t r4,Double_t t4,Double_t b4,Double_t tb4,Int_t side)
{
  
  
  // average the rates for an entree in the quartet super ratio
  // if t1 or t2 is 0 then that run is missing so return the rate
  // of the existing run as the correct average.
  // if neither exists return 0.
  
  // r1,r2,r3,r4     - background subtracted beta rates
  // b1,b2,b3,b4     - background rates
  // t1,t2,t3,t4     - live times for the beta runs
  // tb1,tb2,tb3,tb4 - live times for the background runs 
  
  Double_t ave = 0.;
  Double_t sig = 0.;
  
  vector<Double_t>rate,bck,t,tb;
  
  if(r1>0){
    rate.push_back(r1);
    bck.push_back(b1);
    t.push_back(t1);
    tb.push_back(tb1);
  }
  if(r2>0){
    rate.push_back(r2);
    bck.push_back(b2);
    t.push_back(t2);
    tb.push_back(tb2);
  }
    if(r3>0){
    rate.push_back(r3);
    bck.push_back(b3);
    t.push_back(t3);
    tb.push_back(tb3);
  }
    if(r4>0){
    rate.push_back(r4);
    bck.push_back(b4);
    t.push_back(t4);
    tb.push_back(tb4);
  }
  
  for(Int_t i = 0 ; i < (int)rate.size() ; i++){
    ave += rate[i]*t[i];
    sig += t[i];
  };
  
  ave = ave /sig;
  
  return ave;
};

void Octet::Load_Background()
{
  
  fstream fparon;
  Int_t bin;
  
  if(GetFirst() > 9300 && GetFirst() < 9814){
    fparon.open("input_files/bkg_rad_gas.rate",fstream::in);
  } else {
    fparon.open("input_files/bkg_normal.rate",fstream::in);
  }
  
  for(Int_t i = 0 ; i< 53 ; i++){
    fparon >> bin >> parente[i][0] >> parente[i][1] >> parentw[i][0] >> parentw[i][1];
  }

  fparon.close();
  
  return;
}
//---------------------------------------------------------------------------
void Octet::OutPutToDB()
{
  //-----------------------------------------------------------------------
  // This set of functions send the event rates to the analysis 
  // database.  
  //-----------------------------------------------------------------------
  // A Octets
  // Connect to the analysis database  ------------------------------------
  AnalysisDB *adb = AnalysisDB::getADB();
  //
  OutPutRatesToMPMDB(A2,AFP_OFF,TYPE_0_EVENT ,hA2[3],hA2w[3],tA2,tA2w,adb);
  OutPutRatesToMPMDB(A2,AFP_OFF,TYPE_I_EVENT ,hA2[5],hA2w[5],tA2,tA2w,adb);
  OutPutRatesToMPMDB(A2,AFP_OFF,TYPE_II_EVENT,hA2[6],hA2w[6],tA2,tA2w,adb);
  OutPutRatesToMPMDB(A5,AFP_ON,TYPE_0_EVENT ,hA5[3],hA5w[3],tA5,tA5w,adb);
  OutPutRatesToMPMDB(A5,AFP_ON,TYPE_I_EVENT ,hA5[5],hA5w[5],tA5,tA5w,adb);
  OutPutRatesToMPMDB(A5,AFP_ON,TYPE_II_EVENT,hA5[6],hA5w[6],tA5,tA5w,adb);
  OutPutRatesToMPMDB(A7,AFP_ON,TYPE_0_EVENT ,hA7[3],hA7w[3],tA7,tA7w,adb);
  OutPutRatesToMPMDB(A7,AFP_ON,TYPE_I_EVENT ,hA7[5],hA7w[5],tA7,tA7w,adb);
  OutPutRatesToMPMDB(A7,AFP_ON,TYPE_II_EVENT,hA7[6],hA7w[6],tA7,tA7w,adb);
  OutPutRatesToMPMDB(A10,AFP_OFF,TYPE_0_EVENT ,hA10[3],hA10w[3],tA10,tA10w,adb);
  OutPutRatesToMPMDB(A10,AFP_OFF,TYPE_I_EVENT ,hA10[5],hA10w[5],tA10,tA10w,adb);
  OutPutRatesToMPMDB(A10,AFP_OFF,TYPE_II_EVENT,hA10[6],hA10w[6],tA10,tA10w,adb);
  // B Octets -----------------------------------------------------------------
  OutPutRatesToMPMDB(B2,AFP_ON,TYPE_0_EVENT ,hB2[3],hB2w[3],tB2,tB2w,adb);
  OutPutRatesToMPMDB(B2,AFP_ON,TYPE_I_EVENT ,hB2[5],hB2w[5],tB2,tB2w,adb);
  OutPutRatesToMPMDB(B2,AFP_ON,TYPE_II_EVENT,hB2[6],hB2w[6],tB2,tB2w,adb);
  OutPutRatesToMPMDB(B5,AFP_OFF,TYPE_0_EVENT ,hB5[3],hB5w[3],tB5,tB5w,adb);
  OutPutRatesToMPMDB(B5,AFP_OFF,TYPE_I_EVENT ,hB5[5],hB5w[5],tB5,tB5w,adb);
  OutPutRatesToMPMDB(B5,AFP_OFF,TYPE_II_EVENT,hB5[6],hB5w[6],tB5,tB5w,adb);
  OutPutRatesToMPMDB(B7,AFP_OFF,TYPE_0_EVENT ,hB7[3],hB7w[3],tB7,tB7w,adb);
  OutPutRatesToMPMDB(B7,AFP_OFF,TYPE_I_EVENT ,hB7[5],hB7w[5],tB7,tB7w,adb);
  OutPutRatesToMPMDB(B7,AFP_OFF,TYPE_II_EVENT,hB7[6],hB7w[6],tB7,tB7w,adb);
  OutPutRatesToMPMDB(B10,AFP_ON,TYPE_0_EVENT ,hB10[3],hB10w[3],tB10,tB10w,adb);
  OutPutRatesToMPMDB(B10,AFP_ON,TYPE_I_EVENT ,hB10[5],hB10w[5],tB10,tB10w,adb);
  OutPutRatesToMPMDB(B10,AFP_ON,TYPE_II_EVENT,hB10[6],hB10w[6],tB10,tB10w,adb);
};
//----------------------------------------------------------------------------
void Octet::OutPutRatesToMPMDB(vector<Int_t> nRuns,AFPState afp,
  EventType EType,TH1F *hEast,TH1F *hWest,Double_t tEast,Double_t tWest,
  AnalysisDB *adb)
{
        if((int)nRuns.size() == 0 ) return;
        AnaCutSpec cuts;
	
        cuts.emin   = hAsyA[2]->GetBinCenter(nlow) - hAsyA[2]->GetBinWidth(nlow)/2.; 
        cuts.emax   = hAsyA[2]->GetBinCenter(nhigh) + hAsyA[2]->GetBinWidth(nhigh)/2.;
        cuts.radius = GetRadCut();

       // AnalysisDB *adb = AnalysisDB::getADB();
        AnaResult Asym;
	//Set Analysis Description
        Asym.anatp     = AnaResult::ANA_COUNTS;
  	Asym.datp      = AnaResult::REAL_DATA;
	Asym.grouping  = AnaResult::GROUP_RUN;
	Asym.anach     = ANCHOICE_C;
	//Set Runs
  	Asym.startRun  = nRuns.front();
  	Asym.endRun    = nRuns.back();
	//Output East
  	Asym.etypes.insert(EType);
	if(EType == TYPE_II_EVENT)Asym.etypes.insert(TYPE_III_EVENT);
  	Asym.s         = EAST;
  	Asym.afp       = afp;
	
  	vector<AnaResult> adbm = adb->findMatching(Asym);
  	for(Int_t i = 0 ; i < (int)adbm.size(); i++)
	  adb->deleteAnaResult(adbm[i].arid);
  	Asym.value     = hEast->Integral(nlow,nhigh)*tEast;
  	Asym.err       = TMath::Sqrt(hEast->Integral(nlow,nhigh)*tEast);
  	Asym.csid      = adb->uploadCutSpec(cuts);
  	adb->uploadAnaResult(Asym);
	
	// OutPut West
  	Asym.s         = WEST;
        vector<AnaResult> adbmW = adb->findMatching(Asym);
  	for(Int_t i = 0 ; i < (int)adbmW.size(); i++)
	  adb->deleteAnaResult(adbmW[i].arid);
  	Asym.value     = hWest->Integral(nlow,nhigh)*tWest;
  	Asym.err       = TMath::Sqrt(hWest->Integral(nlow,nhigh)*tWest);
  	Asym.csid      = adb->uploadCutSpec(cuts);
  	adb->uploadAnaResult(Asym);
};
