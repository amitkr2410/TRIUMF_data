#include"iostream"
#include"fstream"
#include"TGraph.h"
#include"TCanvas.h"
#include"TChain.h"
#include"TGraph.h"
#include"TSpectrum.h"
#include"TH1F.h"
#include"TStyle.h"
using namespace std;

void Plot()
{
  TCanvas *c1 = new TCanvas("c1","c1",900,700);
  TGraphErrors *g1 = new TGraphErrors("MeanEachRingQValue8torr.txt","%*lf %lf %lf %lf");
  //TGraphErrors *GErrorStatisticalCM = new TGraphErrors(n,theta_degCM,CrosssectionCM,Null,ErrorStatisticalCM);
  g1->SetTitle("Peak position of ground state in QValue spectrum Vs Laboratory angle");
  g1->GetYaxis()->SetTitle(" Peak Position (in KeV)"); 
  g1->GetXaxis()->SetTitle("Laboratory angle (in degree)");
  g1->SetMarkerColor(2);g1->SetMarkerStyle(8);
  g1->GetXaxis()->SetRangeUser(70.,130.);
  g1->GetYaxis()->SetRangeUser(-400.,100.);
  c1->SetFillColor(42);  
  c1->GetFrame()->SetFillColor(42);
  g1->Draw("AP");
  
  c1->SaveAs("PeakPosition8Torr.jpeg");

  TCanvas *c2 = new TCanvas("c2","c2",900,700);
  TGraphErrors *g2 = new TGraphErrors("MeanEachRingQValue8torr.txt","%*lf %lf %*lf %*lf %lf %lf");
  //TGraphErrors *GErrorStatisticalCM = new TGraphErrors(n,theta_degCM,CrosssectionCM,Null,ErrorStatisticalCM);             
  g2->SetTitle("FWHM of ground state in QValue spectrum Vs Laboratory angle");
  g2->GetYaxis()->SetTitle(" FWHM (in KeV)");
  g2->GetXaxis()->SetTitle("Laboratory angle (in degree)");
  g2->SetMarkerColor(2);
  g2->SetMarkerStyle(8);
  g2->GetXaxis()->SetRangeUser(25.,55.);
  g2->GetYaxis()->SetRangeUser(300.,1000.);
  c2->SetFillColor(42);
  c2->GetFrame()->SetFillColor(42);
  g2->Draw("AP");

  c2->SaveAs("FWHM8Torr.jpeg"); 

  TCanvas *c3 = new TCanvas("c3","c3",900,700);
  TGraph *g3 = new TGraph("ProtonsEachRingQValue8torr.txt","%*lf %lf %lf ");
  TGraph *g4 = new TGraph("ProtonsEachRingQValue8torr.txt","%*lf %lf %*lf %lf ");
  g3->SetTitle("Protons counts for ground state Vs Rings of Yd detector");
  g3->GetYaxis()->SetTitle("Number of protons ");
  g3->GetXaxis()->SetTitle("Laboratory angle (in degree)");
  g3->SetMarkerColor(3);
  g3->SetMarkerStyle(8);
  g3->GetXaxis()->SetRangeUser(25.,55.);
  g3->GetYaxis()->SetRangeUser(0.,2000.);
  // c3->SetFillColor(42);
  //c3->GetFrame()->SetFillColor(42);
  g3->Draw("AP");g4->SetMarkerColor(2);
  g4->Draw("*P"); c3->SetLogy();c3->Update();
  
  c3->SaveAs("ProtonNumber.jpeg");

  TCanvas *c5 = new TCanvas("c5","c5",900,700);
  TGraphErrors *g5 = new TGraphErrors("ProtonsEachRingQValue8torr.txt","%*lf %lf %*lf %*lf %lf");
  g5->SetTitle("Ratio of total protons to the background events for each ring ");
  g5->GetYaxis()->SetTitle(" Ratio of total protons to the background events ");
  g5->GetXaxis()->SetTitle("Laboratory angle (in degree)");
  g5->SetMarkerColor(2);
  g5->SetMarkerStyle(8);
  g5->GetXaxis()->SetRangeUser(25.,55.);
  g5->GetYaxis()->SetRangeUser(0.,54.);
  c5->SetFillColor(42);
  c5->GetFrame()->SetFillColor(42);
  g5->Draw("AP");
  c5->SaveAs("ProtonRatio.jpeg");
}

