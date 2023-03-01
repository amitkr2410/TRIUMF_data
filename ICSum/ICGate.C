#include "TGraph.h"
#include "TH1F.h"
void ICGate()
{
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(2,2);
  c1->cd(1);
  TFile *f = new TFile("~/Amit/treeIris/root_files/output02327.root");
  TH1F *h = (TH1F *h)f->FindObjectAny("adc15");
  
  TH1F *h10B = (TH1F*)h->Clone("h10B");
  TH1F *h10C = (TH1F*)h->Clone("h10C"); 
  h10B->GetXaxis()->SetRangeUser(1000,1700);
  h10C->GetXaxis()->SetRangeUser(1500,3000);
 
  TH1F *h1 = (TH1F*)h10B->Clone("h1");
  TH1F *h2 = (TH1F*)h10C->Clone("h2");
  //  h1->Draw();
  c1->cd(2);  //10C
  // h2->Draw();
 
  c1->cd(1);  //10B
  Int_t npeaks = 1;
  //Use TSpectrum to find the peak candidates
  TSpectrum *s10B = new TSpectrum(npeaks);
  TSpectrum *s10C = new TSpectrum(npeaks);

  Float_t *Triangles10C;
  Float_t *Triangles10B;
  Double_t inSig = 1.5;
  Int_t nfound = s10B->Search(h1,inSig,"noBackground");//new
  Triangles10B = s10B->GetPositionX();
  printf("Found %d candidate peaks to fit\t %f\n",nfound,Triangles10B[0]);
  
  //Fitting Using gaussian function                                                           
  Float_t LB=50.0;
  Float_t UB=50.0;
  Float_t xxmin=Triangles10B[0]-LB;
  Float_t xxmax=Triangles10B[0]+UB;
  TF1 *myfunc10B = new TF1("myfunc10B","gaus",xxmin,xxmax);//Define a function "myfun" to u\
se range option in fitting function                                                         
  h1->Fit("myfunc10B","R");//R tell Fit function to use user defined range in "myfun" funct\
ion                                                                                         
 Double_t chi2_10B = myfunc10B->GetChisquare();
 Double_t ndf_10B = myfunc10B->GetNDF();//The number of degrees of freedom corresponds to the number of points used in the fit minus the number of free parameters.                         
 Double_t p0_10B = myfunc10B->GetParameter(0);  //Height of gaussian function                  
 Double_t p1_10B = myfunc10B->GetParameter(1);  //mean of gaussian function                    
 Double_t p2_10B = myfunc10B->GetParameter(2);  //sigma of gaussian function                   
 Double_t ep1_10B = myfunc10B->GetParError(1);  //Error in peak                                
 cout<<"Chi2="<<chi2_10B;
 cout<<"mean ="<<p1_10B<<"\tSigma="<<p2_10B<<"\n";
 cout<<"ChiSquarePerNo_of_deg_Freedom = \t"<<chi2_10B/ndf_10B<<"\n";
 cout<<"Error in Peak channel number = "<<ep1_10B<<"\n";
 //h1->Draw();
 
 // myfunc10B->Draw("same");
 gStyle->SetOptFit(1112);
 gStyle->SetOptStat(0000000);
 c1->cd(2);
//Fitting Using gaussian function
 Int_t nfound = s10C->Search(h2,inSig,"noBackground");//new                                                                                                   
 Triangles10C = s10C->GetPositionX();
 printf("Found %d candidate peaks to fit\t %f\n",nfound,Triangles10C[0]);
  Float_t LB=70.0;
  Float_t UB=60.0;
  Float_t xxmin=Triangles10C[0]-LB;
  Float_t xxmax=Triangles10C[0]+UB;
  TF1 *myfunc10C = new TF1("myfunc10C","gaus",xxmin,xxmax);//Define a function "myfun" to use range option in fitting function
  h2->Fit("myfunc10C","R");//R tell Fit function to use user defined range in "myfun" function
  //TF1 *myfun = h2->GetFunction("gaus");
  Double_t chi2_10C = myfunc10C->GetChisquare();
  Double_t ndf_10C=myfunc10C->GetNDF();//The number of degrees of freedom corresponds to the number of points used in the fit minus the number of free parameters.
  Double_t p0_10C = myfunc10C->GetParameter(0);  //Height of gaussian function
  Double_t p1_10C = myfunc10C->GetParameter(1);  //mean of gaussian function
  Double_t p2_10C = myfunc10C->GetParameter(2);  //sigma of gaussian function
  Double_t ep1_10C = myfunc10C->GetParError(1);  //Error in peak
  cout<<"Chi2="<<chi2_10C;
  cout<<"mean ="<<p1_10C<<"\tSigma="<<p2_10C<<"\n";
  cout<<"ChiSquarePerNo_of_deg_Freedom = \t"<<chi2_10C/ndf_10C<<"\n";
  cout<<"Error in Peak channel number = "<<ep1_10C<<"\n";
  // h2->Draw();

  // myfunc10C->Draw("same");
  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(0000000);

  c1->cd(3);
  TCanvas *c2 = new TCanvas("c2");
  xxmin = p1_10B-4*p2_10B;
  xxmax = p1_10B+8*p2_10B;

  TF1 *sgf = new TF1("sgf","[0]*exp(-0.5*pow((x-[1])/([2]+(x<[1])*[3]*(x-[1])),2))");//SKewed gaussians with 4 parameters: for gaussian fits set par#3 as zero
  TF1 *POL4 = new TF1("POL4","[4] +[5]*x + [6]*x*x + [7]*x*x*x"); 
 TF1 *MYFUN10B = new TF1("MYFUN10B","sgf+POL4",xxmin,xxmax);

   MYFUN10B->FixParameter(1,p1_10B); //Peak position
  //  MYFUN10B->FixParameter(2,p2_10B); //Sigam of gaussian
 MYFUN10B->FixParameter(0,p0_10B);   //Height
  h->Draw();
  c2->SetLogy();
  h->Fit("MYFUN10B","R");
  //MYFUN10B->Draw("same");
  MYFUN10B->SetLineColor(kRed);

  xxmin = p1_10C-8*p2_10C;
  xxmax = p1_10C+8*p2_10C;
  TF1 *MYFUN10C = new TF1("MYFUN10C","sgf+POL4",xxmin,xxmax);
    MYFUN10C->FixParameter(1,p1_10C);
    //MYFUN10C->FixParameter(2,p2_10C);  
    MYFUN10C->FixParameter(0,p0_10C);
    MYFUN10C->SetLineColor(8);
    h->Fit("MYFUN10C","R");
    /*
    TLegend *leg = new TLegend(0.05,0.75,0.35,0.2);
    leg->AddEntry(MYFUN10B,"Gaussian plot with 3 Sigam ","l");   // h1 and h2 are histogram pointers
    leg->AddEntry(MYFUN10C,"Gaussian plot with 3 Sigma","l");
    leg->Draw();
    
    c1->cd(4);
    h->Draw();
    xxmin = p1_10C-3*p2_10C;
    xxmax = p1_10C+3*p2_10C;
    TF1 *MYFUN10C = new TF1("MYFUN10C","gaus",xxmin,xxmax);
    MYFUN10C->FixParameter(1,p1_10C);
    MYFUN10C->FixParameter(2,p2_10C);
    MYFUN10C->FixParameter(0,p0_10C);
    MYFUN10C->SetLineColor(8);
    MYFUN10C->Draw("same");
    TLine *l1 = new TLine(xxmin,0,xxmin,650);
    l1->Draw("same");  			
    l1->SetLineStyle(3);
    TLine *l2 = new TLine(xxmax,0,xxmax,650);
    l2->Draw("same");
    l2->SetLineStyle(3);
    */
} 
