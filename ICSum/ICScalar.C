//Counting the 10C beam particle for each run  at 19.5 pressure events
//This will be used in finding the differential cross section

//...............This needed before doing current task...../

//........................................................./
#include"TGraph.h"
#include"TFile.h"
#include "fstream"

void ICScalar()
{

 Int_t run;
 Char_t var1[50],var2[50];
 TCanvas *c1 = new TCanvas("c1");
 c1->Divide(1,2);

 ifstream MyFile;
 // ofstream OutFile;OutFile.open("Scalar10CEntries195torr.txt",ios::out);
 ofstream OutFile;OutFile.open("/home/iris/Amit/Analysis_10C_8torr/Plot/ICSum/Scalar10CEntriesWithICPedestalSubtraction8torr.txt",ios::out);
 OutFile<<"#run \t"<<"NetICScalarCounts \t"<<"C10Entries \t"<<"NetICScalarCounts/C10Entries \t"<<"TotalAcceptedTrigger/TotalFreeTrigger"<<endl;

 TGraph *g1 = new TGraph("/home/iris/Amit/Analysis_10C_8torr/Plot/ICSum/Scalar.txt","%lf %lf %*lf %*lf");  //run : ICScalar
 TGraph *g2 = new TGraph("/home/iris/Amit/Analysis_10C_8torr/Plot/ICSum/Scalar.txt","%lf %*lf %lf %*lf");  //run: FreeTrigger Total
 TGraph *g3 = new TGraph("/home/iris/Amit/Analysis_10C_8torr/Plot/ICSum/Scalar.txt","%lf %*lf %*lf %lf");  //run: AcceptedTrigger Total
 MyFile.open("/home/iris/Amit/Analysis_10C_8torr/data8torr.txt",ios::in);
 if(MyFile.is_open())
   {
   cout<<"File is good"<<endl;
   MyFile>>run;
   cout<<run<<endl;
     while(!MyFile.eof())
     {
       TChain *ch = new TChain("Iris");
       sprintf(var1,"~/Amit/treeIris/root_files/tree0%d.root",run);

       ch->AddFile(var1);

       c1->cd(1);
       //ch->Draw("TIC15>>h0","","");        //Taking enrire IC adc spectrum
       ch->Draw("TIC15>>h0","TIC15 >200. ",""); //Taking entire IC adc spectrum without pedestal
       c1->cd(2);
       ch->Draw("TIC15>>h1","TIC15 >680. && TIC15 < 871 ","");
              
       Double_t Integral0,Integral1;
       Integral0 = h0->GetEntries(); // TOtal no of particles present in IC adc spectrum
       Integral1 = h1->GetEntries(); //No of C10 present in IC adc spectrum
       h0->GetXaxis()->SetTitle("Channel Number");
       h0->GetYaxis()->SetTitle("Counts");
       h1->GetXaxis()->SetTitle("Channel Number");
       h1->GetYaxis()->SetTitle("Counts");
       gStyle->SetOptStat(1100001);

       Double_t ICScalarCounts = g1->Eval(run); //IC Scalar Counts in given run#
       Double_t TotalFreeTrigger = g2->Eval(run); //FreeTrigger
       Double_t TotalAcceptedTrigger = g3->Eval(run); //Accepted Trigger
       Double_t NetICScalarCounts = ICScalarCounts*(TotalAcceptedTrigger/TotalFreeTrigger); //ICScalar counts including deadtime effect
       Double_t Ratio = NetICScalarCounts/Integral0; //Ratio of ICScalar and total IC adc spectrum events including deadtime effect  
       Double_t C10Entries = Ratio*Integral1; //Total No of 10C in beam   
       cout<<"\n ICScalarCounts for run# "<<run<<"\t is \t " <<ICScalarCounts;
       cout<<"\n Ratio of AcceptedTrigger to the FreeTrigger is \t"<<TotalAcceptedTrigger/TotalFreeTrigger;   
       cout<<"\n NetICScalarCounts i.e. icluding deadtime effect is \t" <<NetICScalarCounts;
       cout<<"\n Total Intergal in IC spectrum = "<<Integral0;
       cout<<"\n Total Integral under C10 spectrum = "<<Integral1;
       cout<<"\n Ratio of ICScalar to ICSpectrumCounts including deadtime effect is  "<<Ratio;
       cout<<"\n 10C Beam particles during this is \t"<<C10Entries<<endl;
       OutFile<<run<<"\t"<<NetICScalarCounts<<"\t"<<C10Entries<<"\t"<<NetICScalarCounts/C10Entries<<"\t"<<TotalAcceptedTrigger/TotalFreeTrigger<<endl;

       delete ch;
     MyFile>>run;
     cout<<run<<endl;
    }
   }
 MyFile.close(); 
 OutFile.close();

}
