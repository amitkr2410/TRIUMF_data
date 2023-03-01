//To find differential Cross-section of 10C(p,p)10C for each ring,
//we are going to find the ratio of No. of Scattered proton(scaled to standard thickness after background subtraction) to the 10C incident beam particle
//We will not evaluate full expression of differential crosssection
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"
#include "fstream"
#include "TCutG.h"
#include "TChain.h"
#include "iostream"
using namespace std;

void CrossSectionPreliminary()
{

  TCanvas *c1 = new TCanvas("c1","c1",1000,900);
  TGraph *g1 = new TGraph("ProtonsEachRingQValue195torr.txt","%lf %*lf %*lf %*lf %lf");        //YdNo : Ratio=Proton/Background)
  TGraph *g2 = new TGraph("MeanEachRingQValue195torr.txt","%lf %*lf %lf %*lf %*lf %*lf %*lf"); // YdNo : Mean(KeV)
  TGraph *g3 = new TGraph("MeanEachRingQValue195torr.txt","%lf %*lf %*lf %*lf %*lf %*lf %lf"); // YdNo : Sigma(KeV)
  TGraph *g4 = new TGraph("~/Amit/treeIris/scripts/H2thickness.txt","%lf %lf %*lf"); // run : SHT thickness(um) 
  TGraph *g4error = new TGraph("~/Amit/treeIris/scripts/H2thickness.txt","%lf %*lf %lf"); // run : SHT thicknessError(um)
//  TGraph *g5 = new TGraph("~/Amit/Analysis_10C/Plot/ICSum/Scalar10CEntries195torr.txt","%lf %*lf %lf"); //(Includes IC pedestal events)run :10C beam particles 
  TGraph *g5 = new TGraph("~/Amit/Analysis_10C/Plot/ICSum/Scalar10CEntriesWithICPedestalSubtraction195torr.txt","%lf %*lf %lf %*lf"); //(excluding the IC pedestal events in counting)run :10C beam particles

  Int_t run;
  Double_t NScatterFinal[16] ={0.}, NIncidentFinal[16]={0.}, TotalProtons[16]={0.};
  Double_t SigmaSquareProton[16]={0.};
  Double_t SigmaSquareSHT[16]={0.};
  Double_t StandardSHTThickness = 100.; //We scale protons from all run to 100um SHT Thickness 

  TFile *f = new TFile("~/Amit/treeIris/gates/protons.root");
  //f->ls();
  //  protonsNewCsI1;
TCutG *protons1;
protons1 = (TCutG *)f->FindObjectAny("protonsNewCsI1");
protons1->SetName("protons1");
   f->Close();

  ifstream MyFile;
  MyFile.open("/home/iris/Amit/Analysis_10C/data195torr.txt",ios::in); //All the 10C+SHT run for 19.5torr IC pressure
  if(MyFile.is_open())
    {
      cout<<"File is good"<<endl;
      MyFile>>run;
      cout<<"processing run number#"<<run<<endl;
       while(!MyFile.eof())
	{
	  Char_t var1[70];
	  sprintf(var1,"~/Amit/treeIris/root_files/tree0%d.root",run);
	  TChain *ch = new TChain("Iris");
	  ch->AddFile(var1);

	  for(Int_t YdNo=1 ; YdNo<17; YdNo++)
	    { 

	      Char_t var2[70];
              Char_t varHist[70];
              Char_t varPlot[70];
	      sprintf(varPlot,"QValue1>>H%d(160,-6.,6.)",YdNo);
              sprintf(varHist,"H%d",YdNo);
	      sprintf(var2,"protons1 && TYdRing == %d && TIC15 > 1786 && TIC15 <2189",YdNo-1);cout<<endl;
	      cout<<var2 <<endl;
              ch->Draw(varPlot,var2);

	      TH1F * H = (TH1F*)gDirectory->FindObjectAny(varHist);	     
	      cout<< "H = " << H << " YdNo = "<< YdNo <<endl;
              Double_t Sigma = g3->Eval(YdNo)/1000.;// in MeV
              Double_t Mean = g2->Eval(YdNo)/1000.; //in MeV
              Double_t SHTTh = g4->Eval(run);
              Double_t SHTThError = g4error->Eval(run);
              Double_t xxmin = Mean - 3*Sigma;cout<<"xxmin = "<<xxmin<<endl;
              Double_t xxmax = Mean + 3*Sigma;cout<<"xxmax = "<<xxmax<<endl;
	      //Num = H->GetEntries();cout<<"Num = "<<Num<<endl;
	      Int_t bmin = H->GetXaxis()->FindBin(xxmin); cout<<"bmin = "<<bmin<<endl;
	      Int_t bmax = H->GetXaxis()->FindBin(xxmax); cout<<"bmax = "<<bmax<<endl;
	      Double_t integral=0.;integral = H->Integral(bmin,bmax);cout<<"H entries under 3Sigma are = "<<integral<<endl;
	      Double_t BinWidth = H->GetBinWidth(1);cout<<"BinWidth = "<<BinWidth<<endl;
	      Int_t TotalBin = (xxmax-xxmin)/BinWidth;cout<<"Total Bin under 3Sigma ="<<TotalBin<<endl;
	      Double_t integral_noise = 0.;integral_noise = integral/g1->Eval(YdNo);cout<<"Noise Integral under 3Sigma  = "<<integral_noise<<"\t NetProton afterbackground subt ="<<integral - integral_noise<<endl;
              
              TotalProtons[YdNo-1] = TotalProtons[YdNo-1] + (integral - integral_noise);  //NetProtons after background subtraction;
	      cout<<"Cummulative number of protons for YdNo = "<<YdNo<<"\t "<<TotalProtons[YdNo-1]<<endl;
	      
	      NScatterFinal[YdNo-1] = NScatterFinal[YdNo-1] + ((integral - integral_noise)*StandardSHTThickness/SHTTh); //Net proton scaled to the Standard SHT thickness =100.um  or 1.um
	      cout<<"Total Scaled Scattered proton for YdNo = "<<YdNo<<"\t "<<NScatterFinal[YdNo-1]<<endl;

	      //SigmaSquareSHT[YdNo-1] = SigmaSquareSHT[YdNo-1] + pow(((integral - integral_noise)*SHTThError/(SHTTh*SHTTh)),2); //Amit
              //SigmaSquareProton[YdNo-1] = SigmaSquareProton[YdNo-1] + ((integral - integral_noise)/(SHTTh*SHTTh)); ///Amit
              SigmaSquareProton[YdNo-1] = NScatterFinal[YdNo-1];	     
	      NIncidentFinal[YdNo-1] = NIncidentFinal[YdNo-1] + g5->Eval(run); //No of 10C in each run */
	      cout<<"Incident 10C particle in the run#"<<run<<"\t is \t"<<g5->Eval(run)<<endl;
	      cout<<"Processed run number = "<<run<<"and Yd ring = "<<YdNo<<endl;
	      //H->Reset();
	    }//end of forLoop YdNo1-16
	  delete ch;
	  MyFile>>run;cout<<"Started processing run number = "<<run<<endl;
	} //end of While.eof()Loop
       
    } //end of if Myfile.openLoop
  
 
  //Write Scattering particle data on File
  ofstream OutFile;
//  OutFile.open("/home/iris/Amit/Analysis_10C/CrossSection/ScatteringDataStandardThickness.txt",ios::app);
    OutFile.open("/home/iris/Amit/Analysis_10C/CrossSection/ScatteringDataStandardThicknessWithICPedestalSubtraction.txt",ios::out);
  //OutFile.open("test.txt",ios::app);
    OutFile<<"#YdNo"<<"\t"<<"NScaledScatterFinal"<<"\t"<<"NIncidentFinal"<<"\t"<<"SigmaSquareNScattered"<<"\t"<<"TotalProtons"<<endl;
  for(Int_t YdNo =1 ; YdNo<17 ; YdNo++)
    {
      OutFile<<YdNo<<"\t"<<NScatterFinal[YdNo-1]<<"\t"<<NIncidentFinal[YdNo-1]<<"\t"<<SigmaSquareProton[YdNo-1]<<"\t"<<TotalProtons[YdNo-1]<<endl;
    }
    OutFile.close(); 
  MyFile.close();	
}


