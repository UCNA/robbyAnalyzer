#ifndef __Types_h__
#define __Types_h__

struct Asym_t{
    std::vector<Double_t> run_number;
    std::vector<Double_t> A_ave;
    std::vector<Double_t> A_error;
};
 
// struct AnaChoice_t{
//   std::vector<Double_t> Asuper1 (9,0);
//   std::vector<Double_t> Asuper1er (9,0);
//   std::vector<Double_t> Asuper2 (9,0);
//   std::vector<Double_t> Asuper2er (9,0);
//   std::vector<Double_t> Bsuper1 (9,0);
//   std::vector<Double_t> Bsuper1er (9,0);
//   std::vector<Double_t> Bsuper2 (9,0);
//   std::vector<Double_t> Bsuper2er (9,0);
//   std::vector<Double_t> A_multi_Aer (9,0);
//   std::vector<Double_t> A_multi_Ber (9,0);
//   std::vector<Double_t> A_multier (9,0);
//   std::vector<Double_t> A_multi_A (9,0);
//   std::vector<Double_t> A_multi_B (9,0);
//   std::vector<Double_t> A_multi (9,0);
//   std::vector<Double_t> A_sum_Aer (9,0);
//   std::vector<Double_t> A_sum_Ber (9,0);
//   std::vector<Double_t> A_sumer (9,0);
//   std::vector<Double_t> A_sum_A (9,0);
//   std::vector<Double_t> A_sum_B (9,0);
//   std::vector<Double_t> A_sum (9,0);
// };

struct Scint_t {
     Float_t q1;
     Float_t q2;
     Float_t q3;
     Float_t q4;
     Float_t e1;
     Float_t e2;
     Float_t e3;
     Float_t e4;
     Float_t nPE1;
     Float_t nPE2;
     Float_t nPE3;
     Float_t nPE4;
     Float_t de1;
     Float_t de2;
     Float_t de3;
     Float_t de4;
     Float_t energy;
     Float_t denergy;
  }; 

struct mpmmwpc_t {
     Float_t center;
     Float_t width;
     Float_t maxValue;
     Float_t cathSum;
     Int_t   maxWire;
     Int_t   nClipped;
     Int_t   mult;
     Int_t   err;
  };

struct backs_t {
      
      Float_t e_all;
      Float_t etype_0;
      Float_t etype_1;
      Float_t etype_23;
      
      Float_t e_alle;
      Float_t etype_0e;
      Float_t etype_1e;
      Float_t etype_23e;
      
      Float_t w_all;
      Float_t wtype_0;
      Float_t wtype_1;
      Float_t wtype_23;
      
      Float_t w_alle;
      Float_t wtype_0e;
      Float_t wtype_1e;
      Float_t wtype_23e;
      
      Float_t etype_1_bck;
      Float_t etype_23_bck;
      Float_t wtype_1_bck;
      Float_t wtype_23_bck;
      Float_t w_all_bck;
      Float_t e_all_bck;
     
  };

// struct AnaChoiceHists_t {
//     TH1F *heq,*hwq;               // Qadc Spectrum for each peak
//     TH1F *heqC,*hwqC;             // Analysis C option
//     TH1F *heqB,*hwqB;             // Analysis B option
//     TH1F *heqD,*hwqD;             // Analysis D option
//     TH1F *heqE,*hwqE;             // Analysis E option
//     TH1F *heqF,*hwqF;             // Analysis F option
//     TH1F *heqG,*hwqG;             // Analysis G option
//     TH1F *heqH,*hwqH;             // Analysis H option
//     TH1F *heqI,*hwqI;             // Analysis I option
// }
  
#endif
