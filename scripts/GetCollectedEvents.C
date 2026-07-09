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
 void GetCollectedEvents(TString inFile, Double_t time){
    
    TString base = "0p7_counts";
    TString pdf = base+".pdf";
    TString rootName = base + ".root";

    TFile *f = new TFile(inFile);

    TH1F* h_ee_holder = (TH1F*) f->Get("h_ee_Yield_type_a_elast");
    TH1F* h_ee_Count = (TH1F*) h_ee_holder->Clone();
    h_ee_Count->SetName("h_ee_Count");
    h_ee_Count->SetTitle("e-e Collected Count;#theta (#circ);Count");

    TH1F* h_percentH2_ee = (TH1F*) f->Get("h_percentH2_ee");

    h_ee_Count->Scale(80000);
    h_ee_Count->Multiply(h_percentH2_ee);
    //h_ee_Count->Scale(0.95);
    h_ee_Count->Scale(1.0/1000); 
    h_ee_Count->Scale(time);

    Float_t countSum = 0;
    for(Int_t i = 0; i < 8; i++){
        countSum += h_ee_Count->GetBinContent(38 + i);
    }
    cout<< h_ee_Count->GetXaxis()->GetBinCenter(38) << " "<<countSum << endl;

    
    TH1F* h_ep_holder = (TH1F*) f->Get("h_ep_Yield_type_a_numBlocks");
    TH1F* h_ep_Count = (TH1F*) h_ep_holder->Clone();
    h_ep_Count->SetName("h_ee_Count");
    h_ep_Count->SetTitle("e-p Collected Count;#theta (#circ);Count");

    TH1F* h_percentH2_ep = (TH1F*) f->Get("h_percentH2_ep");

    
    h_ep_Count->Scale(80000);
    h_ep_Count->Multiply(h_percentH2_ep);
    //h_ep_Count->Scale(0.95);
    h_ep_Count->Scale(1.0/1000);
    h_ep_Count->Scale(time);

    Float_t countSum_ep = 0;
    for(Int_t x = 0; x < 8; x++){
        countSum_ep += h_ep_Count->GetBinContent(38 + x);
    }
    cout<< h_ep_Count->GetXaxis()->GetBinCenter(38) << " "<<countSum_ep << endl;

    TH1F* h_rG_ee = (TH1F*) f->Get("h_resGas_ee");
    TH1F* h_resGas_ee = (TH1F*) h_rG_ee->Clone();
    h_resGas_ee->SetName("h_resGas_ee_coll");

    TH1F* h_tR_ee = (TH1F*) f->Get("h_theRest_ee");
    TH1F* h_theRest_ee = (TH1F*) h_tR_ee->Clone();
    h_theRest_ee->SetName("h_theRest_ee_coll");

    TH1F* h_rG_ep = (TH1F*) f->Get("h_resGas_ep");
    TH1F* h_resGas_ep = (TH1F*) h_rG_ep->Clone();
    h_resGas_ep->SetName("h_resGas_ep_coll");

    TH1F* h_tR_ep = (TH1F*) f->Get("h_theRest_ep");
    TH1F* h_theRest_ep = (TH1F*) h_tR_ep->Clone();
    h_theRest_ee->SetName("h_theRest_ep_coll");
    
    TCanvas *c = new TCanvas("c", "Estimated Counts Canvas",1200,1000);
    auto legend = new TLegend(0.1,0.8,0.4,0.9);
    c->cd(1);
    h_ee_Count->SetMarkerStyle(20);
    gPad->SetLogy(1);
    h_ee_Count->Draw();
    c->Print(pdf + "(");
    c->Clear();

    c->cd(1);
    gPad->SetLogy(0);
    h_resGas_ee->SetAxisRange(0,0.3, "Y");
    h_resGas_ee->SetStats(0);
    h_resGas_ee->Draw("E P");
    h_theRest_ee->Draw("SAME E P");
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_resGas_ee,"(b)-(c) Res. Gas","p");
	legend->AddEntry(h_theRest_ee,"(c) Cell & Coll.","p");
    legend->Draw();
    c->Print(pdf);
    legend->Clear();
    c->Clear();

    c->cd(1);
    h_ep_Count->SetMarkerStyle(20);
    //gPad->SetLogy(1);
    h_ep_Count->SetAxisRange(3.7, 4.6, "X");
    h_ep_Count->Draw();
    c->Print(pdf);
    c->Clear();

    c->cd(1);
    gPad->SetLogy(0);
    h_resGas_ep->SetAxisRange(0,0.6, "Y");
    h_resGas_ep->SetStats(0);
    h_resGas_ep->Draw("E P");
    h_theRest_ep->Draw("SAME E P");
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_resGas_ep,"(b)-(c) Res. Gas","p");
	legend->AddEntry(h_theRest_ep,"(c) Cell & Coll.","p");
    legend->Draw();
    c->Print(pdf + ")");
    legend->Clear();
    c->Clear();

    TObjArray* arr = new TObjArray(0,0);
    (*arr).Add(h_ee_Count);
    (*arr).Add(h_ep_Count);
    (*arr).Add(h_ee_holder);
    (*arr).Add(h_percentH2_ee);

    TFile file(rootName,"RECREATE");
    (*arr).Write();
	file.Close();

 }