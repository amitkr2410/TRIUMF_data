#include "TGraph"
#include "fstream"
void ICScalarPlot()
{
  TCanvas *c1 = new TCanvas("c1");
  c1->Divide(1,2);

  //  TGraph *g1 = new TGraph("~/Amit/Analysis_10C/Plot/ICSum/Scalar10CEntries195torr.txt","%lf %lf %*lf %*lf");
  //TGraph *g2 = new TGraph("~/Amit/Analysis_10C/Plot/ICSum/Scalar10CEntries195torr.txt","%lf %*lf  %lf %*lf");
  //TGraph *g3 = new TGraph("~/Amit/Analysis_10C/Plot/ICSum/Scalar10CEntries195torr.txt","%lf %*lf  %*lf %lf");

  TGraph *g1 = new TGraph("~/Amit/Analysis_10C_8torr/Plot/ICSum/Scalar10CEntriesWithICPedestalSubtraction8torr.txt","%lf %lf  %*lf %*lf %*lf");
  TGraph *g2 = new TGraph("~/Amit/Analysis_10C_8torr/Plot/ICSum/Scalar10CEntriesWithICPedestalSubtraction8torr.txt","%lf %*lf  %lf %*lf %*lf");
  TGraph *g3 = new TGraph("~/Amit/Analysis_10C_8torr/Plot/ICSum/Scalar10CEntriesWithICPedestalSubtraction8torr.txt","%lf %*lf  %*lf %*lf %lf");
  
  g1->GetXaxis()->SetTitle("Run number");
  g1->GetYaxis()->SetTitle("Total ICScaler counts(With deadtime correction) for each run");
  g1->GetYaxis()->SetRangeUser(0,44000000);
  g1->SetMarkerColor(2);
  g1->SetMarkerStyle(3);
  g1->SetMarkerSize(1);  

  g2->GetXaxis()->SetTitle("Run number");
  g2->GetYaxis()->SetTitle("10C particle counts for each run");
  g2->GetYaxis()->SetRangeUser(0,44000000);
  g2->SetMarkerColor(3);
  g2->SetMarkerStyle(3);
  g2->SetMarkerSize(1);

    g3->GetXaxis()->SetTitle("Run number");
    g3->GetYaxis()->SetTitle("Ratio of Total AcceptedTrigger to Total FreeTrigger");
    g3->GetYaxis()->SetRangeUser(0.,1.0);
    g3->SetMarkerColor(2);
    g3->SetMarkerStyle(3);
    g3->SetMarkerSize(1);

  c1->cd(1);
  g1->Draw("ALP");
  g2->Draw("LP");
  
  c1->cd(2);
  g3->Draw("ALP");
  c1->SaveAs("/home/iris/Amit/Analysis_10C_8torr/Plot/ICSum/Ratio_TotalIC_10CeventsICWithICPedestalSubtraction8torrdata.png");
}
