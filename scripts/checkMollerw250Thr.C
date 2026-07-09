/**
 * Macro that gets estimated collected e-e and e-p events based on the yield histograms in the supplied root file.
 *
 * @author Erik Wrightson
 * @version 05.14.2025
 * @creation 05.14.2025
 */
 
 //ROOT Includes that may be handy to have.
 #include <TROOT.h>
 #include <TChain.h>
 #include <TFile.h>
 #include <TStyle.h>
 #include <TCanvas.h>
 #include <TObject.h>
 #include <TGraphAsymmErrors.h>
 #include <TGraphErrors.h>
 #include <TTree.h>
 #include <TLeaf.h>
 #include <TFitResult.h>
 #include <TH1.h>
 #include <TH1F.h>
 #include <TH2.h>
 #include <TF1.h>
 #include <TH3.h>
 #include <TH3F.h>
 #include <TLegend.h>
 #include <TMath.h>
 #include <TColor.h>
 #include <TRandom3.h>
 //Needed for reading in vector types from root files.
 #include <TInterpreter.h>
 //#include <RDataFrame.h>
 //Allows for 4D vector handling and Boosting between frames in ROOT 
 //(remember to include -lGenVector in the make file since it is not included in the default `root-config --glibs`)
 #include <Math/GenVector/LorentzVector.h>
 #include <Math/PxPyPzE4D.h>
 #include <Math/Vector4Dfwd.h>
 #include <Math/Boost.h>
 
 typedef ROOT::Math::PxPyPzEVector p4vec;
 typedef ROOT::Math::Boost boost;
 
 //Generally useful includes for non-ROOT file handling and in/output stream management
 #include <iostream>
 #include <fstream>
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <tuple>
 void checkMollerw250Thr(TString inFile, TString inFile2){
    
    TString base = "0p7_TypeA_Thr200OverThr300";
    TString pdf = base+".pdf";
    TString rootName = base + ".root";

    TFile *f = new TFile(inFile);

    TH1F* h_ee_holder = (TH1F*) f->Get("h_ee_Yield_type_a_elast");
    TH1F* h_ee_200Thr_over_300Thr = (TH1F*) h_ee_holder->Clone();
    h_ee_200Thr_over_300Thr->SetName("h_ee_200Thr_over_300Thr");
    h_ee_200Thr_over_300Thr->SetTitle("e-e Yield 200Thr/300Thr ;#theta (#circ);Count");

    TFile *f2 = new TFile(inFile2);
    TH1F* h_ee_300Thr = (TH1F*) f2->Get("h_ee_Yield_type_a_elast");

    h_ee_200Thr_over_300Thr->Divide(h_ee_300Thr);
    
    TCanvas *c = new TCanvas("c", "Estimated Counts Canvas",1200,1000);
    c->cd(1);
    h_ee_200Thr_over_300Thr->SetMarkerStyle(20);
    h_ee_200Thr_over_300Thr->SetAxisRange(0.9,1.4, "Y");
    h_ee_200Thr_over_300Thr->Draw();
    c->Print(pdf);
    c->Clear();

    TObjArray* arr = new TObjArray(0,0);
    (*arr).Add(h_ee_200Thr_over_300Thr);

    TFile file(rootName,"RECREATE");
    (*arr).Write();
	file.Close();

 }