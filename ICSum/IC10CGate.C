//Make 10C selection from IC spectrum and choose how much sigma one should choose
//Only for decided the 10C gate selection for further analysis
#include "iostream"
#include "fstream"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLine.h"

using namespace std;

void IC10CGate()
{
  TCanvas *c1 = new TCanvas("c1","c1",1700,900);
    
  TFile *f = new TFile("~/Amit/treeIris/root_files/output02274.root"); //Read IC spectrum from this file
  TH1F *h = (TH1F *) f->FindObjectAny("adc15")->Clone("h");
  
  TH1F *h10C = (TH1F *) h->Clone("h10C"); 
  h10C->GetXaxis()->SetRangeUser(600,1000);
  TH1F *h2 = (TH1F *) h10C->Clone("h2");

  Int_t npeaks = 1;
  //Use TSpectrum to find the peak candidates
  TSpectrum *s10C = new TSpectrum(npeaks);

  Float_t *Triangles10C;
  Double_t inSig = 1.5;

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
  c1->SaveAs("10CICPeak8torr.png")

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  xxmin = p1_10C-4*p2_10C;
  xxmax = p1_10C+4*p2_10C;
  TF1 *sgf = new TF1("sgf","[0]*exp(-0.5*pow((x-[1])/([2]),2))");//SKewed gaussians with 4 parameters: for gaussian fits set par#3 as zero
  //  TF1 *POL4 = new TF1("POL4","[4] +[5]*x + [6]*x*x + [7]*x*x*x"); 
  TF1 *MYFUN10C = new TF1("MYFUN10C","sgf",xxmin,xxmax);

   MYFUN10C->FixParameter(1,p1_10C); //Peak position
  //  MYFUN10C->FixParameter(2,p2_10C); //Sigam of gaussian
  MYFUN10C->FixParameter(0,p0_10C);   //Height
  h10C->Draw();
  c2->SetLogy();
  h10C->Fit("MYFUN10C","R");
  MYFUN10C->SetLineColor(kRed);
  Double_t p2_Sigma = MYFUN10C->GetParameter(2);
  xxmin = p1_10C - 2.3*p2_Sigma;
  xxmax = p1_10C + 3.0*p2_Sigma;
  TLine *l1 = new TLine(xxmin,0,xxmin,650);
  TLine *l2 = new TLine(xxmax,0,xxmax,650);
  l1->Draw("same");
  l2->Draw("same");
  cout<<"Lower Channel number for 10C = "<<xxmin<<endl;
  cout<<"Upper Channel number for 10C = "<<xxmax<<endl;
  //  c2->SaveAs("10CGateSelectionFinal8torr.png");
} 
