//Fit the QValue spectrum for each spectrum using "Gaussian+LinearFit" and find the fit parameter
//Using Fit parameters draw gaussian and linear function over the QValue_Histogram
//Find the 3*Sigma range and hence count the no of proton in QValue_Histogram and also find the background usign linear fitting curve under 3*Sigma range. 
//Using the Fitting parameter we find following quantity of intereset:
// Mean of ground state, Sigma of ground state and hence FWHM and 3*Sigma range

#include"iostream"
#include"fstream"
#include"TGraph.h"
#include"TCanvas.h"
#include"TChain.h"
#include"TGraph.h"
#include"TSpectrum.h"
#include"TH1F.h"
#include"TStyle.h"
#inclue"TCutG.h"

using namespcae std;

void QValueEachRing(Int_t YdNo = 1)
	{

  TFile *f = new TFile("~/Amit/Analysis_10C/Plot/QValueEachRing195torr.root");
  f->ls();  
  TCanvas *c1 = new TCanvas("c1","c1",1700,900);
  //  gStyle->SetOptFit(0000);
  //  gStyle->SetOptStat(0);
 
  Char_t var[50];
  sprintf(var,"h%i",YdNo);
  TH1F *h = (TH1F *)f->FindObjectAny(var);

  if (YdNo >=10) //==>10
    {
      h->Rebin(2);
      TH1F *H = (TH1F* )h->Clone("H");
      h->GetXaxis()->SetRangeUser(-2,+2);
      H->GetXaxis()->SetRangeUser(-2,+2);
      h->GetYaxis()->SetTitle("Counts/150KeV");
      h->GetXaxis()->SetTitle("QValue (in MeV)");
      H->GetYaxis()->SetTitle("Counts/150KeV");
      H->GetXaxis()->SetTitle("QValue (in MeV)");
    }
  else
    {
      TH1F *H = (TH1F* )h->Clone("H");
      h->GetXaxis()->SetRangeUser(-2,+2);
      H->GetXaxis()->SetRangeUser(-2,+2);
      h->GetYaxis()->SetTitle("Counts/75KeV");
      h->GetXaxis()->SetTitle("QValue (in MeV)");  
      H->GetYaxis()->SetTitle("Counts/75KeV");
      H->GetXaxis()->SetTitle("QValue (in MeV)");
    }

  c1->Divide(2,1);
  c1->cd(1);
  h->Draw();
  Double_t Pedestal = 0.;  
  TSpectrum *s = new TSpectrum(1);
  Int_t nfound = s->Search(h,1);
  Float_t * Triangle = s->GetPositionX();
  cout<<"Peak is at "<<Triangle[0]<<endl; 


  Double_t xxmin = Triangle[0]-1.2;
  Double_t xxmax = Triangle[0]+1.2;  
  
  TF1 *myfunc_gaus = new TF1("myfunc_gaus","gaus",xxmin,xxmax);//Define a function "myfun" to use range option in fitting function
  myfunc_gaus->SetParameter(1,Triangle[0]);
  h->Fit("myfunc_gaus","R");//R tell Fit function to use user defined range in "myfun" function
  //TF1 *myfun = h2->GetFunction("gaus");
  Double_t chi2_gaus = myfunc_gaus->GetChisquare();
  Double_t ndf_gaus=myfunc_gaus->GetNDF();//The number of degrees of freedom corresponds to the number of points used in the fit minus the number of free parameters.
  Double_t p0_gaus = myfunc_gaus->GetParameter(0);  //Height of gaussian function
  Double_t p1_gaus = myfunc_gaus->GetParameter(1);  //mean of gaussian function
  Double_t p2_gaus = myfunc_gaus->GetParameter(2);  //sigma of gaussian function
  Double_t ep1_gaus = myfunc_gaus->GetParError(1);  //Error in peak
  cout<<"Chi2 using Gaussian Fit alone ="<<chi2_gaus;
  cout<<"mean ="<<p1_gaus<<"\tSigma="<<p2_gaus<<"\n";
  cout<<"ChiSquarePerNo_of_deg_Freedom for gaussian fitting = \t"<<chi2_gaus/ndf_gaus<<"\n";
  cout<<"Error in Peak channel number = "<<ep1_gaus<<"\n";

  //myfunc_gaus->Draw("same");

  //  gStyle->SetOptFit(0000);
  h->SetStats(kFALSE);

  c1->cd(2);
  H->Draw();
  H->SetStats(kFALSE);

  xxmin = p1_gaus-6.0*p2_gaus; //3.2*Sigma 6.0*p2_gaus
  xxmax = 2.1;//p1_gaus+; // 3.2*Sigma 6.0*p2_gaus
  Double_t par[5]={10.};
  TF1 *gf  = new TF1("gf","[0]*exp(-0.5*pow((x-[1])/([2]),2))");
  // TF1 *sgf = new TF1("sgf","[0]*exp(-0.5*pow((x-[1])/([2]+(x>[1])*[3]*(x-[1])),2))",xxmin,xxmax,4);//SKewed gaussians with 4 parameters: for gaussian fits set par#3 as zero
  TF1 *POL5 = new TF1("POL5","[0] +[1]*x");
  //TF1 *LORENTZ = new TF1("LORENTZ","[10]/(pow(x-[11],2)+[10]*[10])");
  TF1 *MyfunTotal = new TF1("MyfunTotal","gf+POL5",xxmin,xxmax);
  //par[3]=p2_gaus; //Skewed Factor
   par[2]=p2_gaus; //Sigam of gaussian
   par[1]=p1_gaus; //Mean of gaussian
   par[0]=p0_gaus;   //Height  of gaussian

  par[3]=0.1;
  par[4]=0.2;   
  // par[5]=0.1; //Poly4 C0
  //par[6]=0.1; //Poly4 C1
  //par[7]=0.1; //Poly4 C2
  //par[8]=0.1; //Poly4 C3
  //par[9]=0.1; //LORENTZ
  //par[10]=0.1; //LORENTZ
  MyfunTotal->SetParameters(par);
  //  MyfunTotal->SetParameters(p0_gaus, p1_gaus, p2_gaus,5.10,3.);
  if (YdNo ==15 || YdNo ==14 || YdNo == 12 ) 
    {
      MyfunTotal->FixParameter(4,0.);    
    }

   MyfunTotal->SetParNames("Height(p0)","Mean(p1)","Sigma(p2)","PolyC0(p3)","PolyC1(p4)");  
  //c1->SetLogy();
   H->Fit("MyfunTotal","R");
   MyfunTotal->Draw("Same");  
   MyfunTotal->GetParameters(par);MyfunTotal->SetLineColor(3);

  TF1 *gf1 = new TF1("gf1","gf",xxmin,xxmax);
  TF1 *POL = new TF1("POL","POL5",xxmin,xxmax);
  //TF1 *MF = new TF1("MF","MyfunTotal",xxmin-0.3,xxmax+0.3);
  gf1->SetParameters(par[0],par[1],par[2]);
  POL->SetParameters(par[3],par[4]);
  //MF->SetParameters(par); 
  gf1->Draw("Same");gf1->SetLineColor(4);
  POL->Draw("Same");POL->SetLineColor(28);
  //  MF->Draw("Same");MF->SetLineColor(3); 
  Double_t Sigma = MyfunTotal->GetParameter(2);
  Double_t Mean  = MyfunTotal->GetParameter(1);
  Double_t ErrorSigma = MyfunTotal->GetParError(2);
  Double_t ErrorMean  = MyfunTotal->GetParError(1);
 
 //( 3 ) Analysis under 3Sigma region
  xxmin = Mean - 3*Sigma;
  xxmax = Mean + 3*Sigma;
  TLine *l1 = new TLine(xxmin,0,xxmin,MyfunTotal->GetParameter(0));
  TLine *l2 = new TLine(xxmax,0,xxmax,MyfunTotal->GetParameter(0));
  l1->Draw("same");
  l2->Draw("same");
  cout<<"Lower Qvalue (in MeV)  = "<<xxmin<<endl;
  cout<<"Uppre Qvalue (in MeV)  = "<<xxmax<<endl;

    
  TAxis *axis = H->GetXaxis();
  Int_t bmin = axis->FindBin(xxmin); cout<<"bmin = "<<bmin<<endl;
  Int_t bmax = axis->FindBin(xxmax); cout<<"bmax = "<<bmax<<endl;
  Double_t integral = h->Integral(bmin,bmax);cout<<"H entries under 3Sigma are = "<<integral<<endl;
  Double_t BinWidth = h->GetBinWidth(1);cout<<"BinWidth = "<<BinWidth<<endl;
  Int_t TotalBin = (xxmax-xxmin)/BinWidth;cout<<"Total Bin under 3Sigma ="<<TotalBin<<endl; 
  Double_t integral_noise = TotalBin*(par[3]+par[4]*0.5*(xxmin+xxmax));cout<<"Noise Integral under 3Sigma  = "<<integral_noise<<endl;
  Double_t Resolution = Sigma*2.35*100/(Mean);
  Double_t ErrorResolution = Resolution*( (ErrorSigma/Sigma) + (ErrorMean/Mean) );

  cout<<"Mean of QValue for 19.5 torr data is "<<Mean<<"+/-"<<ErrorMean <<"(in MeV)"<<endl;
  cout<<"Sigma of QValue for 19.5 torr data is "<<Sigma<<"+/-"<<ErrorSigma<<"(in MeV)"<<endl;
  cout<<"FWHM for QValue for 19.5 torr data is "<<Sigma*2.35<<"+/-"<<ErrorSigma<<"(in MeV)"<<endl;   
  cout<<"ChiSquarePerNo_of_deg_Freedom for SkewedGauss+POL4 fitting = \t"<<MyfunTotal->GetChisquare()/MyfunTotal->GetNDF()<<"\n";
   
  Double_t theta;
  Double_t Yd1r=50., Yd2r = 129.; // inner and outer radii in mm
  Double_t YdDistance = 100.; //distance from source to YY1 for run 2627 (140 mm from target position) tk check this
  theta = atan((Yd1r*(16.-(YdNo-1)%16-0.5)+Yd2r*( (YdNo-1)%16 + 0.5))/16./YdDistance);
  //...................(4)............StatusBox....................//
  //..............................................................//
  TPaveText * pt = new TPaveText(0.62,0.58,0.89,0.88,"NDC");
  pt->SetFillColor(0);
  pt->SetTextSize(0.02);
  pt->SetTextAlign(12);
  sprintf(var,"Mean = %4.3f   #pm %4.3f MeV",Mean, ErrorMean);pt->AddText(var); 
  sprintf(var,"Sigma = %4.3f  #pm %4.3f MeV",Sigma, ErrorSigma);pt->AddText(var);
  sprintf(var,"FWHM = %4.3f   #pm %4.3f MeV",2.35*Sigma,ErrorSigma);pt->AddText(var);
  sprintf(var,"DF : p0*e^{-0.5*(#frac{x-p1}{p2})^2} + c0+c1*x");pt->AddText(var);
  sprintf(var,"p0= %4.1f, p1=%4.2f, p2= %4.3f",par[0],par[1],par[2]);pt->AddText(var);
  sprintf(var,"c0=%4.3f, c1= %4.3f",par[3],par[4]);pt->AddText(var);
  sprintf(var,"Chi2/NDF = %4.2f ", MyfunTotal->GetChisquare()/MyfunTotal->GetNDF());pt->AddText(var);
  sprintf(var,"Total Protons = %i ", integral);pt->AddText(var);
  sprintf(var,"Background events =  %i",integral_noise);pt->AddText(var);
  
  pt->Draw();
  TLegend *newbox = new TLegend(0.62,0.48,0.89,0.58);
  newbox->SetFillColor(0);
  newbox->SetTextSize(0.02);
  newbox->AddEntry(MyfunTotal,"Gauss+LinearFit","l");
  newbox->AddEntry(gf1,"Gauss","l");
  newbox->AddEntry(POL,"LinearFit","l");
  newbox->Draw();
  ofstream OutFile1,OutFile2,OutFile3;
  OutFile1.open("MeanEachRingQValue195torr.txt",ios::app);
  OutFile2.open("FittingParEachRingQValue195torr.txt",ios::app);
  OutFile3.open("ProtonsEachRingQValue195torr.txt",ios::app);
  if(YdNo==1)
    {
      OutFile1<<"#YdNo"<<"\t"<<"Lab Angle (in degree) \t"<<"Mean(KeV) \t"<<"ErrorMean \t"<<"FWHM( KeV)\t"<<"ErrorFWHM \t"<<"Sigma(KeV)"<<"ErrorSigma"<<endl;
    }
  OutFile1<<YdNo<<"\t"<<theta*180/3.1415926<<"\t"<<Mean*1000<<"\t"<<ErrorMean*1000<<"\t"<<Sigma*2.35*1000<<"\t"<<ErrorSigma*1000<<"\t"<<Sigma*1000<<"\t"<<ErrorSigma*1000<<endl;

  if(YdNo==1)
    {
  OutFile2<<"#YdNo"<<"\t"<<"par[0]:Height"<<"\t"<<"par[1]:Mean"<<"\t"<<"par[2]:Sigma"<<"\t"<<"par[3]:Constant"<<"\t"<<"par[4]:Slope"<<endl;
    }
  OutFile2<<YdNo<<"\t"<<par[0]<<"\t"<<par[1]<<"\t"<<par[2]<<"\t"<<par[3]<<"\t"<<par[4]<<endl;

  if(YdNo==1)
    {
      OutFile3<<"#YdNo"<<"\t"<<"theta(degree)"<<"\t"<<"NoOfProton"<<"\t"<<"BackgroundProton"<<"\t"<<"NoOfProton/BackgroundProton"<<"\t"<<"NetProton"<<endl;
    }
  OutFile3<<YdNo<<"\t"<<theta*180/3.1415926<<"\t"<<integral<<"\t"<<integral_noise<<"\t"<<integral/integral_noise<<"\t"<<integral-integral_noise<<endl;

  OutFile1.close();
  OutFile2.close(); 
  OutFile3.close(); 

  sprintf(var,"Image%i.jpeg",YdNo);
  c1->SaveAs(var);  
}
