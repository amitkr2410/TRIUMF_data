//Evaluate full expression of differential cross section
//Use text file "ScatteringDataStandardThicknessWithICPedestalSubtraction.txt" or "ScatteringDataStandardThickness.txt"
//and "Efficiency.txt" and "Kinematics_10C_pp_10CGroundState19_5torr.root"
//Generate png files displaying the differential cross section in Lab-frame and CM-frame
//Generate png files displaying Lab and CM relation
//Generate png files displaying efficiency vs Lab-angle, Jacobian of dOmega for lab-to-cm transformation
#include "iostream"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TFrame.h"
#include "TMath.h"
#include "cstdlib"
#include "fstream"
using namespace std;

void CrossSectionFinal()
{
/*
  TGraph *gp1 = new TGraph("ScatteringDataStandardThickness.txt","%*lf %lf %lf %*lf"); // YdNo : NScattered: NIncident: ErrorNScattered
  TGraph *gp1error = new TGraph("ScatteringDataStandardThickness.txt","%lf %*lf %*lf %lf"); // YdNo : NScattered: NIncident: ErrorNScattered
*/
TGraph *gp1 = new TGraph("ScatteringDataStandardThicknessWithICPedestalSubtraction.txt","%*lf %lf %lf %*lf %*lf"); //NScaledProton: NIncident 
TGraph *gp1error = new TGraph("ScatteringDataStandardThicknessWithICPedestalSubtraction.txt","%lf %*lf %*lf %lf %*lf"); //YdNo:SigmaSquareProton
 
  TGraph *gp2 = new TGraph("Efficiency.txt","%lf %lf");           //YdNo : Efficiency
  TGraphErrors *gpNN = new TGraphErrors("TheoryCrossSection.txt","%lf %lf %lf"); //NN crossSectionCM
  TGraphErrors *gpNN3N = new TGraphErrors("TheoryCrossSection.txt","%lf %*lf %*lf %lf %lf"); //NN+3N crossSectionCM
  TGraphErrors *gpNN3NInd = new TGraphErrors("TheoryCrossSection.txt","%lf %*lf %*lf %*lf %*lf %lf %lf"); //NN+3NInduced crossSectionCM

  TFile *f = new TFile("Kinematics_10C_pp_10CGroundState19_5torr.root"); // ThetadegLab(x): ThetadegCM(y)
  f->ls();

  TGraph * GLabCM; 
  GLabCM = (TGraph *)f->FindObjectAny("g3")->Clone(); //thetaLab(x) :thetacm(y)

  
  Double_t Ns[16],Ni[16],SigmaSquareProton[16],SigmaSquareSHT[16];
  Double_t CrosssectionLab[16],theta_degLab[16],dtheta_radLab[16],theta_radLab;
  Double_t CrosssectionCM[16],theta_degCM[16],dtheta_radCM[16], theta_radCM,garbage, CM1rad,CM2rad;
  Double_t Eff[16]={0.9,0.9,0.9,0.9,0.8,0.8,0.6,0.6,0.5,0.5,0.5,0.5,0.45,0.45,0.4,0.4};
  Double_t CrosssectionLabError[16];
  Double_t CrosssectionCMError[16];
  Double_t CrosssectionCMWEff[16],CrosssectionCMWEffError[16],CrosssectionLabWEff[16],CrosssectionLabWEffError[16]; //Optional Ritu( for efficiency=1 all rings)
  Double_t ErrorSystematicLab[16],ErrorStatisticalLab[16],Ns2Error[16];
  Double_t ErrorSystematicCM[16],ErrorStatisticalCM[16];
  Double_t Jacobian[16]={0.};
  Double_t Null[16] = {0.};  
  Int_t n = gp1->GetN();

  Double_t Yd1r=50., Yd2r = 129.; // inner and outer radii in mm
  Double_t YdDistance = 100.; //distance from source to YY1 for run 2627 (140 mm from target position) tk check this

  Double_t MH = 2*1.008;  //in g/mol ; =1.008g/mol;
  Double_t Density = 0.086*pow(10,6) ; // in g/m3;  =0.086 g/cm3
  Double_t NA = 2*6.023*pow(10,23); //
  Double_t StandardSHTThickness = 100.*pow(10,-6); //in meter; = 100 micron
  ofstream OutFile,OutFileRitu;
  OutFile.open("CrossSectionResult195torr.txt",ios::out);
  OutFileRitu.open("CrossSectionResultRitu195torr.txt",ios::out);
  OutFile<<"YdNo"<<"DCSLab+/-ErrorDCSLab"<<"\t"<<"%ErrorDCSLab"<<"\t"<<"ErrorStatitcalLab"<<"\t"<<"%ErrorStatitcalLab"<<"\t"<<"ErrorSystematicLab"<<"\t"<<"%ErrorSystematicLab"<<"\t"<<"DCSCM+/-ErrorDCSCM"<<"\t"<<"%ErrorDCSCM"<<"\t"<<"ErrorStatitcalCM"<<"\t"<<"%ErrorStatitcalCM"<<"\t"<<"ErrorSystematicCM"<<"\t"<<"%ErrorSystematicCM"<<endl;
  OutFileRitu<<"#YdNo\t"<<"LabAngle(deg)\t"<<"CMangle(deg)\t"<<"dSigmadOmegaLab\t"<<"dSigmdOmegaLabError\t"<<"dSigmadOmegaCM\t"<<"dSigmdOmegaCMError\t"<<endl;
  for(Int_t i=0; i<n;i++)
    {
      gp1->GetPoint(i,Ns[i],Ni[i]);
      gp1error->GetPoint(i,garbage,SigmaSquareProton[i]);
      gp2->GetPoint(i,garbage,Eff[i]);
      Ns2Error[i] = SigmaSquareProton[i];
      SigmaSquareSHT[i]= pow((0.05*StandardSHTThickness),2);

      dtheta_radLab[i] = atan((Yd1r+((i+1)*(Yd2r-Yd1r)/16))/YdDistance)-atan((Yd1r+((i)*(Yd2r-Yd1r)/16))/YdDistance) ;
      cout<<"dthetaLab="<<dtheta_radLab[i]<<endl;
      CM1rad = (3.1415926/180.)*(GLabCM->Eval((180./3.1415926)*atan((Yd1r+((i+1)*(Yd2r-Yd1r)/16))/YdDistance)));
      CM2rad = (3.1415926/180.)*(GLabCM->Eval((180./3.1415926)*atan((Yd1r+((i)*(Yd2r-Yd1r)/16))/YdDistance)));

      cout<<"CMRad = "<<CM1rad<<"\t"<<CM2rad<<endl;
      dtheta_radCM[i] = CM2rad - CM1rad;
      cout<<"dthetaCM="<<dtheta_radCM[i]<<endl;
      
      //theta_radLab = atan((Yd1r*(16.-(i)%16-0.5)+Yd2r*( (i)%16 + 0.5))/16./YdDistance);
      Double_t thetaSum = atan((Yd1r+((i+1)*(Yd2r-Yd1r)/16))/YdDistance) + atan((Yd1r+((i)*(Yd2r-Yd1r)/16))/YdDistance) ;
      theta_radLab = thetaSum/2.;
      theta_degLab[i] = theta_radLab*180./3.1415926;
      theta_degCM[i] = GLabCM->Eval(theta_degLab[i]);
      theta_radCM = theta_degCM[i]*(3.1415926/180.);cout<<"\n thetadegCM = "<<theta_degCM[i]<<"\t"<<theta_radCM<<endl;

      CrosssectionLab[i]= pow(10,+31)*(Ns[i]/Ni[i])*MH/(NA*Density*StandardSHTThickness*sin(theta_radLab)*2.*3.1415926*Eff[i]*dtheta_radLab[i]);  
      CrosssectionLabError[i]= sqrt(pow(CrosssectionLab[i],2)*( (Ns2Error[i]/pow(Ns[i],2)) + (SigmaSquareSHT[i]/pow(StandardSHTThickness,2))  + 0.0025 ));
      ErrorStatisticalLab[i] =  sqrt(pow(CrosssectionLab[i],2)*((Ns2Error[i]/pow(Ns[i],2)) ));
      ErrorSystematicLab[i] = sqrt(pow(CrosssectionLab[i],2)*( (SigmaSquareSHT[i]/pow(StandardSHTThickness,2)) + 0.0025 ) );
      cout<<"CrossSection in Lab = "<<CrosssectionLab[i]<<"+/-"<<CrosssectionLabError[i]<<endl;
      cout<<"Error in CrossSection in Lab from Statistics = "<<ErrorStatisticalLab[i]<<"\t due Systematics ="<<ErrorSystematicLab[i]<<endl;     

      CrosssectionLabWEff[i] = CrosssectionLab[i]*Eff[i]; //Optional Ritu (For efficiency=1 for all rings)
      CrosssectionLabWEffError[i] = sqrt(pow(CrosssectionLabWEff[i],2)*( (Ns2Error[i]/pow(Ns[i],2)) + (SigmaSquareSHT[i]/pow(StandardSHTThickness,2)) )); //Optional Ritu (For efficiency=1 for all rings)
      cout<<"CrossSection in Lab for eff as one = "<<CrosssectionLabWEff[i]<<"+/-"<<CrosssectionLabWEffError[i]<<endl;

      Jacobian[i] = sin(theta_radLab)*dtheta_radLab[i]/(sin(theta_radCM)*dtheta_radCM[i]);cout<<"Jacobian = "<<Jacobian[i]<<endl;
      cout<<"Crosssection in CM frame found using Jacobian= "<<CrosssectionLab[i]*Jacobian[i]<<endl;

      CrosssectionCM[i]= pow(10,+31)*(Ns[i]/Ni[i])*MH/(NA*Density*StandardSHTThickness*sin(theta_radCM)*2.*3.1415926*Eff[i]*dtheta_radCM[i]);
      CrosssectionCMError[i] = sqrt(pow(CrosssectionCM[i],2)*( (Ns2Error[i]/pow(Ns[i],2)) + (SigmaSquareSHT[i]/pow(StandardSHTThickness,2)) + 0.0025 ));
      ErrorStatisticalCM[i] = sqrt(pow(CrosssectionCM[i],2)*( (Ns2Error[i]/pow(Ns[i],2)) )); 
      ErrorSystematicCM[i] = sqrt(pow(CrosssectionCM[i],2)*( (SigmaSquareSHT[i]/pow(StandardSHTThickness,2)) + 0.0025));      
  
      cout<<"Crosssection in CM frame directly = "<<CrosssectionCM[i]<<"+/- "<<CrosssectionCMError[i]<<endl;

      CrosssectionCMWEff[i] = CrosssectionCM[i]*Eff[i]; //Optional Ritu (For efficiency=1 for all rings)
      CrosssectionCMWEffError[i] = sqrt(pow(CrosssectionCMWEff[i],2)*( (Ns2Error[i]/pow(Ns[i],2)) + (SigmaSquareSHT[i]/pow(StandardSHTThickness,2)) )); //Optional Ritu (For efficiency=1 for all rings)
      cout<<"For (efficiency =1) leads DCS error = "<<CrosssectionCMWEffError[i]<<endl;
      //Double_t mm=10.016853/1.007825;
      // cout<<"Jacobian from paper = "<<(1.+ mm*cos(3.1415926-theta_radCM))/pow((1+mm*mm + 2*mm*cos(3.1415926-theta_radCM)),1.5)<<endl;
       OutFile<<i<<"\t"<<CrosssectionLab[i]<<"+/-"<<CrosssectionLabError[i]<<"\t"<<CrosssectionLabError[i]*100./CrosssectionLab[i]<<"\t"<<ErrorStatisticalLab[i]<<"\t"<<ErrorStatisticalLab[i]*100./CrosssectionLab[i]<<"\t"<<ErrorSystematicLab[i]<<"\t"<<ErrorSystematicLab[i]*100./CrosssectionLab[i]<<"\t"<<CrosssectionCM[i]<<"+/-"<<CrosssectionCMError[i]<<"\t"<<CrosssectionCMError[i]*100./CrosssectionCM[i]<<"\t"<<ErrorStatisticalCM[i]<<"\t"<<ErrorStatisticalCM[i]*100./CrosssectionCM[i]<<"\t"<<ErrorSystematicCM[i]<<"\t"<<ErrorSystematicCM[i]*100./CrosssectionCM[i]<<endl;

      OutFileRitu<<i<<"\t"<<theta_degLab[i]<<"\t"<<theta_degCM[i]<<"\t"<<CrosssectionLab[i]<<"\t"<<CrosssectionLabError[i]<<"\t"<<CrosssectionCM[i]<<"\t"<<CrosssectionCMError[i]<<endl;
 
}
 
  OutFile.close(); 
  OutFileRitu.close();
  TCanvas *c1 = new TCanvas("c1","c1",1600,700);
  c1->Divide(2,1);

  c1->cd(1); //Cross-section in Lab frame
  TGraphErrors *Glab = new TGraphErrors(n,theta_degLab,CrosssectionLab,Null,CrosssectionLabError);
  Glab->SetTitle("Differential CrossSection in Lab-frame and CM-frame");
  Glab->GetYaxis()->SetTitle("d#sigma/d#Omega(in mb/sr)"); Glab->GetXaxis()->SetTitle("Laboratory angle (in degree)");
  Glab->SetMarkerColor(2);  Glab->SetMarkerStyle(8);  Glab->SetMarkerSize(1);
  c1->SetFillColor(42);  c1->GetFrame()->SetFillColor(42);
  Glab->Draw("AP");

  c1->cd(2); //Cross-section in CM frame
  TGraphErrors *Gcm = new TGraphErrors(n,theta_degCM,CrosssectionCM,Null,CrosssectionCMError);
  Gcm->GetYaxis()->SetTitle("d#sigma/d#Omega(in mb/sr)cm");  Gcm->GetXaxis()->SetTitle("CM angle (in degree)");
  Gcm->SetMarkerColor(2);  Gcm->SetMarkerStyle(8);  Gcm->SetMarkerSize(1);
  c1->SetFillColor(42);  c1->GetFrame()->SetFillColor(42);
  Gcm->GetXaxis()->SetRangeUser(70.,130.);Gcm->GetYaxis()->SetRangeUser(0.,60.);
  Gcm->Draw("AP");
  gpNN->Draw("Same"); gpNN->SetLineColor(9);
  gpNN3N->Draw("Same");gpNN3N->SetLineColor(3);
  gpNN3NInd->Draw("Same");gpNN3NInd->SetLineColor(1);
                                                                                                                                                
  TLegend *leg12 = new TLegend(0.1,0.75,0.3,0.89);
  leg12->AddEntry(Gcm,"Data","p");
  leg12->AddEntry(gpNN3NInd,"NN+3NInd","l");
  leg12->AddEntry(gpNN3N,"NN+3N","l");
  leg12->AddEntry(gpNN,"NN","l");
  leg12->Draw();

  c1->SaveAs("CrosssectionInLabCMFrame.png");
 
  TCanvas *c2 = new TCanvas("c2","c2",1600,700);

  TGraph *Geff = new TGraph(n,theta_degLab,Eff);
  Geff->GetYaxis()->SetTitle("Efficiency of Yd detector");  Geff->GetXaxis()->SetTitle("Lab angle (in degree)");
  Geff->SetMarkerColor(2);  Geff->SetMarkerStyle(8);  Geff->SetMarkerSize(1);
  c2->SetFillColor(42);  c2->GetFrame()->SetFillColor(42);
  Geff->Draw("AP");
  c2->SaveAs("Efficiency.png"); 
 
  TCanvas *c3 = new TCanvas("c3","c3",1600,700);
  c3->Divide(2,1);
  c3->cd(1); //Lab-angle and CM-angle relation 

  GLabCM->GetYaxis()->SetTitle("CM angle ( in degree )");  GLabCM->GetXaxis()->SetTitle("Laboratory angle (in degree)");
  GLabCM->SetMarkerColor(2);  GLabCM->SetMarkerStyle(8);  GLabCM->SetMarkerSize(1);
  c3->SetFillColor(42);  c3->GetFrame()->SetFillColor(42);
  GLabCM->Draw("AL");
  TLine *line1 = new TLine(theta_degLab[0],0.,theta_degLab[0],180.);
  TLine *line2 = new TLine(theta_degLab[15],0.,theta_degLab[15],180.); 
  line1->SetLineColor(30); line1->Draw(); 
  line2->SetLineColor(30); line2->Draw();

  c3->cd(2); //Jacobian                                            
  TGraph *GJacobian = new TGraph(n,theta_degLab,Jacobian);
  GJacobian->GetYaxis()->SetTitle("Jacobian = d#Omega/d#Omega'");  GJacobian->GetXaxis()->SetTitle("Laboratory angle (in degree)");
  GJacobian->SetMarkerColor(2);  GJacobian->SetMarkerStyle(8);  GJacobian->SetMarkerSize(1);
  c3->SetFillColor(42);  c3->GetFrame()->SetFillColor(42);
  GJacobian->GetYaxis()->SetRangeUser(0.2,0.5);
  GJacobian->Draw("AP");
  c3->SaveAs("LabCMRelationandJacobian.png"); 
 
  //Witout taking Efficiency into account "Ritu"
  TCanvas *c4 = new TCanvas("c4","c4",1600,700);
  c4->Divide(2,1);

  c4->cd(1); //Cross-section in Lab frame                                                                                                             
  TGraphErrors *Glabweff = new TGraphErrors(n,theta_degLab,CrosssectionLabWEff,Null,CrosssectionLabWEffError);
  Glabweff->SetTitle("Differential CrossSection in Lab-frame for  Efficiency =1");
  Glabweff->GetYaxis()->SetTitle("(For efficiency = 1) d#sigma/d#Omega(in mb/sr)"); Glabweff->GetXaxis()->SetTitle("Laboratory angle (in degree)");
  Glabweff->SetMarkerColor(2);  Glabweff->SetMarkerStyle(8);  Glabweff->SetMarkerSize(1);
  c4->SetFillColor(42);  c4->GetFrame()->SetFillColor(42);
  Glabweff->Draw("AP");

  c4->cd(2); //Cross-section in CM frame                                                                                                              
  TGraphErrors *Gcmweff = new TGraphErrors(n,theta_degCM,CrosssectionCMWEff,Null,CrosssectionCMWEffError);
  Gcmweff->SetTitle("Differential CrossSection in CM-frame for  Efficiency =1");
  Gcmweff->GetYaxis()->SetTitle("(For efficiency = 1 ) d#sigma/d#Omega(in mb/sr)cm");  Gcmweff->GetXaxis()->SetTitle("CM angle (in degree)");
  Gcmweff->SetMarkerColor(2);  Gcmweff->SetMarkerStyle(8);  Gcmweff->SetMarkerSize(1);
  Gcmweff->GetXaxis()->SetRangeUser(70.,130.);Gcmweff->GetYaxis()->SetRangeUser(0.,60.);
  c4->SetFillColor(42);  c4->GetFrame()->SetFillColor(42);
  Gcmweff->Draw("AP");
  gpNN->Draw("Same"); gpNN->SetLineColor(9);
  gpNN3N->Draw("Same");gpNN3N->SetLineColor(3);//c5->Update();                                                                                                                                             
  gpNN3NInd->Draw("Same");gpNN3NInd->SetLineColor(1);

  TLegend *leg42 = new TLegend(0.1,0.75,0.3,0.89);
  leg42->AddEntry(Gcmweff,"Data","p");
  leg42->AddEntry(gpNN3NInd,"NN+3NInduced","l");
  leg42->AddEntry(gpNN3N,"NN+3N","l");
  leg42->AddEntry(gpNN,"NN","l");
  leg42->Draw();

  c4->SaveAs("ConstantEfficiencyCrosssectionInLabCMFrame.png");

  TCanvas *c5 = new TCanvas("c5","c5",1600,700);
  c5->Divide(2,1);
  
  c5->cd(1); //Cross-section in Lab frame                                                                                                                                     
  TGraphErrors *GErrorSystematicLab = new TGraphErrors(n,theta_degLab,CrosssectionLab,Null,ErrorSystematicLab);
  GErrorSystematicLab->SetTitle("Differential CrossSection with systematic error in Lab-frame");
  GErrorSystematicLab->GetYaxis()->SetTitle("( d#sigma/d#Omega(in mb/sr)"); GErrorSystematicLab->GetXaxis()->SetTitle("Laboratory angle (in degree)");
  GErrorSystematicLab->SetMarkerColor(2);  GErrorSystematicLab->SetMarkerStyle(8);  GErrorSystematicLab->SetMarkerSize(1);
  GErrorSystematicLab->GetXaxis()->SetRangeUser(25.,55.); GErrorSystematicLab->GetYaxis()->SetRangeUser(0.,180.);
  c5->SetFillColor(42);  c5->GetFrame()->SetFillColor(42); 
  GErrorSystematicLab->Draw("AP");
  //gpNN->Draw("LP"); gpNN->->SetMarkerColor(2);
  //gpNN3N->Draw("LP"); gpNN3N->SetMarkerColor(3);
  
  c5->cd(2);
  TGraphErrors *GErrorStatisticalLab = new TGraphErrors(n,theta_degLab,CrosssectionLab,Null,ErrorStatisticalLab);
  GErrorStatisticalLab->SetTitle("Differential CrossSection with statistical error in Lab-frame");
  GErrorStatisticalLab->GetYaxis()->SetTitle("( d#sigma/d#Omega(in mb/sr)"); GErrorStatisticalLab->GetXaxis()->SetTitle("Laboratory angle (in degree)");
  GErrorStatisticalLab->SetMarkerColor(2);GErrorStatisticalLab->SetMarkerStyle(8);  
  GErrorStatisticalLab->Draw("AP");
  
  c5->SaveAs("CrosssectionInLabFrameErrorStatisticalSystematics.png"); 
  
  TCanvas *c6 = new TCanvas("c6","c6",1600,700);
  c6->Divide(2,1);

  c6->cd(1); //Cross-section in CM frame                                                                                                                                       
  TGraphErrors *GErrorSystematicCM = new TGraphErrors(n,theta_degCM,CrosssectionCM,Null,ErrorSystematicCM);
  GErrorSystematicCM->SetTitle("Differential CrossSection with systematic error in CM-frame");
  GErrorSystematicCM->GetYaxis()->SetTitle(" d#sigma/d#Omega(in mb/sr)cm");  GErrorSystematicCM->GetXaxis()->SetTitle("CM angle (in degree)");
  GErrorSystematicCM->SetMarkerColor(2); GErrorSystematicCM->SetMarkerStyle(8);  GErrorSystematicCM->SetMarkerSize(1);
  GErrorSystematicCM->GetXaxis()->SetRangeUser(70.,130.); GErrorSystematicCM->GetYaxis()->SetRangeUser(0.,60.);
  c6->SetFillColor(42);  c6->GetFrame()->SetFillColor(42);
  GErrorSystematicCM->Draw("AP");
  gpNN->Draw("Same"); gpNN->SetLineColor(9); 
  gpNN3N->Draw("Same");gpNN3N->SetLineColor(3);//c5->Update();
  gpNN3NInd->Draw("Same");gpNN3NInd->SetLineColor(1);  

  TLegend *leg = new TLegend(0.1,0.75,0.3,0.89);
  leg->AddEntry(GErrorSystematicCM,"Data","p");
  leg->AddEntry(gpNN3NInd,"NN+3NInduced","l");
  leg->AddEntry(gpNN3N,"NN+3N","l");
  leg->AddEntry(gpNN,"NN","l");
  leg->Draw();
  
  c6->cd(2);
  TGraphErrors *GErrorStatisticalCM = new TGraphErrors(n,theta_degCM,CrosssectionCM,Null,ErrorStatisticalCM);
  GErrorStatisticalCM->SetTitle("Differential CrossSection with statistical error in CM-frame");
  GErrorStatisticalCM->GetYaxis()->SetTitle("( d#sigma/d#Omega(in mb/sr)cm"); GErrorStatisticalCM->GetXaxis()->SetTitle("CM angle (in degree)");
  GErrorStatisticalCM->SetMarkerColor(2);GErrorStatisticalCM->SetMarkerStyle(8);
  GErrorStatisticalCM->GetXaxis()->SetRangeUser(70.,130.);GErrorStatisticalCM->GetYaxis()->SetRangeUser(0.,60.);
  c6->SetFillColor(42);  c6->GetFrame()->SetFillColor(42);
  GErrorStatisticalCM->Draw("AP"); 
  gpNN->Draw("Same"); gpNN->SetLineColor(9);
  gpNN3N->Draw("Same");gpNN3N->SetLineColor(3);//c5->Update();
  gpNN3NInd->Draw("Same");gpNN3NInd->SetLineColor(1);

  TLegend *leg54 = new TLegend(0.1,0.75,0.3,0.89);
  leg54->AddEntry(GErrorSystematicCM,"Data","p");
  leg54->AddEntry(gpNN3NInd,"NN+3NInduced","l");
  leg54->AddEntry(gpNN3N,"NN+3N","l");
  leg54->AddEntry(gpNN,"NN","l");
  leg54->Draw();

  c6->SaveAs("CrosssectionInCMFrameErrorStatisticalSystematics.png");
  
}
