{

  fstream f1;
  fstream f2;
  fstream f3;

  f1.open("histogram_geoA.txt",fstream::in);
  f2.open("histogram_geoB.txt",fstream::in);
  f3.open("histogram_geoC.txt",fstream::in);


  TH1F *hgeoA = new TH1F("hgeoA","Geometry A",100,0,1000);
  TH1F *hgeoB = new TH1F("hgeoB","Geometry B",100,0,1000);
  TH1F *hgeoC = new TH1F("hgeoC","Geometry C",100,0,1000);


  Double_t x,y,t1,t2,t3,t4;
 
  for(Int_t i = 1 ; i < 100 ; i++){
    f1 >> x >> y >> t1 >> t2;
    hgeoA->SetBinContent(i,y);
    f2 >> x >> y >> t1 >> t2;
    hgeoB->SetBinContent(i,y);
    f3 >> x >> y >> t1 >> t2;
    hgeoC->SetBinContent(i,y);
  }

  TCanvas *c1 = new TCanvas("c1","c1");

  hgeoA->DrawNormalized();
  hgeoB->SetLineColor(2);
  hgeoB->DrawNormalized("same");
  hgeoC->SetLineColor(4);
  hgeoC->DrawNormalized("same");

}
