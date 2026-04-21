/**
 * Background Contamination Processor
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 04.21.2026
 * @creation 04.15.2026
 */

#include <iostream>
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <sstream>

//ROOT Includes that may be handy to have.
#include <TROOT.h>
#include <TSystem.h>
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
#include <TLegend.h>
#include <TMath.h>
#include <TColor.h>
#include <TString.h>
#include "ROOT/TThreadedObject.hxx"
//Needed for reading in vector types from root files.
#include <TInterpreter.h>

#include "includes/Yields.h"

using namespace std;

/**
 * Processes the file with the input file name and reads out the root file names contained within and adds them to a vector.
 *
 * @param fileListFileName - the name of the file with the list of ROOT file names.
 */
vector<TString> processFileList(string fileListFileName){
        ifstream file(fileListFileName);
        string line;
        
        vector<TString> list;
        while(getline(file,line)){
            TString l(line);
            list.push_back(l);
        }

        file.close();

        return list;
}

/**
 * Makes a TChain of the entries in a TString vector assuming they are valid paths to ROOT files.
 *
 * @param names - the vector of ROOT file names to be linked in the chain.
 */
TChain* makeChain(vector<TString> names){
        TChain* chain = new TChain("recon");
        for(unsigned int i = 0; i < names.size(); i++){
            chain->Add(names.at(i));
        }

        return chain;
}


/**
 * Prints out the proper usage directions for this program and what each flag means.
 *
 * @param prog - the progam name that is being currently run.
 */
static void printUsage(const char *prog)
{
    cerr << "Usage: " << prog << " [options]\n"
              << "\t-a <fileName_a> Marks which file pertains to the filled gas target with residual gas in the chamber.\n"
              << "\t-b <fileName_b> Marks which file pertains to only gas pumped into the vaccuum chamber.\n"
              << "\t-c <fileName_c> Marks which file pertains to an empty gas target and chamber with cell in place.\n"
              << "\t-d <fileName_d> Marks which file pertains to an empty target and chamber with the cell lifted.\n"
              << "\t-L Indicates that FileList files are being provided instead of individual ROOT files.\n"
              << "\t-f <filename (NO EXTENSION)> Custom name for both the output pdf and the ROOT file.\n"
              << "\t-D <liveChargeDb> Different liveCharge database file location and name.\n"
              << "\t-v verbose: Option to plot and fill histograms before and after every cut.\n"
              << "\t-m <nThreads> Multithreading that will use the input amount of CPU cores \n"
              << "\t-h Show this help\n"
              << "\tNOTE: Either option -f or -L are REQUIRED for running properly.";
}

/**
 * Makes a map where the keys are run numbers and the values are the livecharge values.
 *
 * @param dbName - the name of the .dat file that contains the liveCharge values.
 * @param m - the map to fill with keys and values passed by reference.
 */
void makeMap(string dbName, map<Int_t, Double_t>& m){
    ifstream file(dbName);
    string line;
    
    getline(file, line); //Skip the header line.

    while(getline(file,line)){
        stringstream ss(line);
        string word;

        Int_t wordNum = 0;
        Int_t key;
        Double_t val;

        while (ss >> word) {
            if(wordNum == 0){
                key = stoi(word);
            }
            if(wordNum == 3){
                val = stod(word);
                wordNum = -1;
            }
            wordNum++;
        }
        m[key] = val;
    }

    cout<<"Made the quick reference livecharge map.\n";
}

/**
 * The main function that launches the trigger analysis.
 *
 * @param argc - the number of input arguments
 * @param argv - an array of the different arguments as an array of char* (strings).
 */
int main (int argc, char **argv){

    bool a = false;
    bool b = false;
    bool c = false;
    bool d = false;
    bool list = false;
    bool allPDF = false;
    bool f = false;

    string fileName_a;
    string fileName_b;
    string fileName_c;
    string fileName_d;

    string cpus = "1";

    TString outfile = "./outfiles/defaultOutput.pdf";
    TString root_outfile = "./rootFiles/defaultRootOutput.root";
    TString rootfileLocation = "./rootFiles/";
    TString outfileLocation = "./outfiles/";
    string liveChargeDb = "./database/beam_charge.dat";

	if (argc<1) {
		cout<<"ERR: Incorrect amount of arguments."<<endl;
        printUsage(argv[0]);
		return -1;
	}

    // ── Parse command-line ───────────────────────────────────────────────
    int opt;
    while ((opt = getopt(argc, argv, "a:b:c:d:Lf:D:Am:")) != -1) {
        switch (opt) {
            case 'a': a = true; fileName_a = optarg; break;
            case 'b': b = true; fileName_b = optarg; break;
            case 'c': c = true; fileName_c = optarg;break;
            case 'd': d = true; fileName_d = optarg;break;
            case 'L': list = true; break;
            case 'f': outfile = optarg; f = true; break;
            case 'D': liveChargeDb = optarg; break;
            case 'v': allPDF = true; break;
            case 'm': cpus = optarg; break;
            case 'h':
            default: printUsage(argv[0]); return (opt == 'h') ? 0 : 1;
        }
    }

    if(f){
        root_outfile = rootfileLocation + outfile + ".root";
        outfile = outfileLocation + outfile + ".pdf";
    }

    if(!(a || b || c || d)){
        cerr<<"No file with a valid configuration label was provided. At least one option [-a] [-b] [-c] [-d] must be used.";
        return -2;
    }

    map<Int_t, Double_t> lcMap;
    makeMap(liveChargeDb, lcMap);

    struct stat buffer;   
    bool existA = (stat(fileName_a.c_str(), &buffer) == 0);
    bool existB = (stat(fileName_b.c_str(), &buffer) == 0);
    bool existC = (stat(fileName_c.c_str(), &buffer) == 0);
    bool existD = (stat(fileName_d.c_str(), &buffer) == 0);

    if(!existA &&  !existB && !existC && !existD){
        cerr<<"A single valid input file or a filelist txt file was not provided.\n";
        return -3;
    }

    ROOT::EnableThreadSafety();
    ROOT::EnableImplicitMT(stoi(cpus));//fileNameVec_a.size()*10);

    vector<TString> fileNameVec_a;
    if(list && existA){
        fileNameVec_a = processFileList(fileName_a);
    }
    if(!list && existA){
        fileNameVec_a.push_back((TString) fileName_a);
    }

    vector<TString> fileNameVec_b;
    if(list && existB){
        fileNameVec_b = processFileList(fileName_b);
    }
    if(!list && existB){
        fileNameVec_b.push_back((TString) fileName_b);
    }

    vector<TString> fileNameVec_c;
    if(list && existC){
        fileNameVec_c = processFileList(fileName_c);
    }
    if(!list && existC){
        fileNameVec_c.push_back((TString) fileName_c);
    }

    vector<TString> fileNameVec_d;
    if(list && existD){
        fileNameVec_d = processFileList(fileName_d);
    }
    if(!list && existD){
        fileNameVec_d.push_back((TString) fileName_d);
    }

    
    TChain* fChain_a;
    Yields* a_obj;
    TH1F* h_a_ee_Yield;
    TH1F* h_a_ep_Yield;
    if(a){
        fChain_a = makeChain(fileNameVec_a);
        a_obj = new Yields(fChain_a, 0, lcMap, allPDF);
        a_obj->Evaluate();
        if(!b && !c && !d){
            a_obj->printPDF(outfile, true, true);
        }
        else{
            a_obj->printPDF(outfile, true, false);
        }
        
        h_a_ee_Yield = a_obj->get_ee_YieldHisto();
        h_a_ee_Yield->SetName("h_a_ee_Yield");
        h_a_ep_Yield = a_obj->get_ep_YieldHisto();
        h_a_ep_Yield->SetName("h_a_ep_Yield");
        a_obj->save_Histos(root_outfile, true);
    }

    TChain* fChain_b;
    Yields* b_obj;
    TH1F* h_b_ee_Yield;
    TH1F* h_b_ep_Yield;
    if(b){
        fChain_b = makeChain(fileNameVec_b);
        b_obj = new Yields(fChain_b, 1, lcMap, allPDF);
        b_obj->Evaluate();
        if(!a){
            if(!c && !d){
                b_obj->printPDF(outfile, true, true);
            }
            else{
                b_obj->printPDF(outfile, true, false);
            }
            
        }
        else{
            if(!c && !d){
                b_obj->printPDF(outfile, false, true);
            }
            else{
                b_obj->printPDF(outfile, false, false);
            }
        }
        
        h_b_ee_Yield = b_obj->get_ee_YieldHisto();
        h_b_ee_Yield->SetName("h_b_ee_Yield");
        h_b_ep_Yield = b_obj->get_ep_YieldHisto();
        h_b_ep_Yield->SetName("h_b_ep_Yield");
        b_obj->save_Histos(root_outfile, false);
    }

    TChain* fChain_c;
    Yields* c_obj;
    TH1F* h_c_ee_Yield;
    TH1F* h_c_ep_Yield;
    if(c){
        fChain_c = makeChain(fileNameVec_c);
        c_obj = new Yields(fChain_c, 2, lcMap, allPDF);
        c_obj->Evaluate();
        if(!a && !b){
            if(!d){
                c_obj->printPDF(outfile, true, true);
            }
            else{
                c_obj->printPDF(outfile,true, false);
            }
        }
        else{
            if(!d){
                c_obj->printPDF(outfile, false, true);
            }
            else{
                c_obj->printPDF(outfile, false, false);
            }
        }
        
        h_c_ee_Yield = c_obj->get_ee_YieldHisto();
        h_c_ee_Yield->SetName("h_c_ee_Yield");
        h_c_ep_Yield = c_obj->get_ep_YieldHisto();
        h_c_ep_Yield->SetName("h_c_ep_Yield");
        c_obj->save_Histos(root_outfile, false);
    }

    TChain* fChain_d;
    Yields* d_obj;
    TH1F* h_d_ee_Yield;
    TH1F* h_d_ep_Yield;
    if(d){
        fChain_d = makeChain(fileNameVec_d);
        d_obj = new Yields(fChain_d, 3, lcMap, allPDF);
        d_obj->Evaluate();
        if(!a && !b && !c){
            d_obj->printPDF(outfile, true, true);
        }
        else{
            if(a && b && c){
                d_obj->printPDF(outfile, false, false);
            }
            else{
                d_obj->printPDF(outfile, false, true);
            }
        }
        
        h_d_ee_Yield = d_obj->get_ee_YieldHisto();
        h_d_ee_Yield->SetName("h_d_ee_Yield");
        h_d_ep_Yield = d_obj->get_ep_YieldHisto();
        h_d_ep_Yield->SetName("h_d_ep_Yield");
        d_obj->save_Histos(root_outfile, false);
    }

    TCanvas *c1 = new TCanvas("c1", "BackgroundSubtraction_Canvas",1000,1000);
    auto legend = new TLegend(0.1,0.8,0.4,0.9);

    //Declare the Background contamination histograms for the e-e events
    TH1F* h_resGas_ee = new TH1F("h_resGas_ee", "e-e Background Contamination of H2 Signal;#theta (#circ);N_Cont/N_H2", 60, 0, 6);
    TH1F* h_collimators_ee = new TH1F("h_collimators_ee", "e-e Background Contamination of H2 Signal;#theta (#circ);N_Cont/N_H2", 60, 0, 6);
    TH1F* h_cell_ee = new TH1F("h_cell_ee", "e-e Background Contamination of H2 Signal;#theta (#circ);N_Cont/N_H2", 60, 0, 6);
    TH1F* h_H2Gas_ee = new TH1F("h_H2Gas_ee", "e-e Background Contamination of H2 Signal;#theta (#circ);N_Cont/N_H2", 60, 0, 6);

    //Declare the background contamination histograms for the e-p events
    TH1F* h_resGas_ep = new TH1F("h_resGas_ep", "e-p Background Contamination of H2 Signal;#theta (#circ);N_Cont/N_H2", 60, 0, 6);
    TH1F* h_collimators_ep = new TH1F("h_collimators_ep", "e-p Background Contamination of H2 Signal;#theta (#circ);N_Cont/N_H2", 60, 0, 6);
    TH1F* h_cell_ep = new TH1F("h_cell_ep", "e-p Background Contamination of H2 Signal;#theta (#circ);N_Cont/N_H2", 60, 0, 6);
    TH1F* h_H2Gas_ep = new TH1F("h_H2Gas_ep", "e-p Background Contamination of H2 Signal;#theta (#circ);N_Cont/N_H2", 60, 0, 6);

    if(a && b && c && d){

        //e-e Yield Background Contamintation
        h_H2Gas_ee->Add(h_a_ee_Yield, h_b_ee_Yield,1, -1);

        /*c1->Divide(2);
        c1->cd(1);
        h_b_ee_Yield->Draw("E P");
        c1->cd(2);
        h_c_ee_Yield->Draw("E P");
        c1->Print(outfile);
        c1->Clear();*/

        h_resGas_ee->Add(h_b_ee_Yield, h_c_ee_Yield, 1, -1);
        h_resGas_ee->Divide(h_H2Gas_ee);
        h_resGas_ee->SetMarkerStyle(20);
        h_resGas_ee->SetMarkerColor(kBlack);
        h_resGas_ee->SetLineColor(kBlack);

        h_cell_ee->Add(h_c_ee_Yield, h_d_ee_Yield, 1, -1);
        h_cell_ee->Divide(h_H2Gas_ee);
        h_cell_ee->SetMarkerStyle(21);
        h_cell_ee->SetMarkerColor(kRed);
        h_cell_ee->SetLineColor(kRed);

        h_collimators_ee->Divide(h_d_ee_Yield,h_H2Gas_ee);
        h_collimators_ee->SetMarkerStyle(23);
        h_collimators_ee->SetMarkerColor(kBlue);
        h_collimators_ee->SetLineColor(kBlue);

        c1->cd(1);
        h_resGas_ee->SetAxisRange(0,0.10, "Y");
        h_resGas_ee->SetStats(0);
        h_resGas_ee->Draw("P E");
        h_cell_ee->Draw("P E SAME");
        h_collimators_ee->Draw("P E SAME");


        legend->SetHeader("Legend","C"); // option "C" allows to center the header
	    legend->AddEntry(h_resGas_ee,"(b)-(c) Res. Gas","p");
	    legend->AddEntry(h_cell_ee,"(c)-(d) Cell","p");
        legend->AddEntry(h_collimators_ee,"(d) Coll.","p");
	    legend->Draw();
        c1->Print(outfile);
        c1->Clear();
        legend->Clear();

        //e-p Yield Background Contamination
        h_H2Gas_ep->Add(h_a_ep_Yield, h_b_ep_Yield,1, -1);

        h_resGas_ep->Add(h_b_ep_Yield, h_c_ep_Yield, 1, -1);
        h_resGas_ep->Divide(h_H2Gas_ep);
        h_resGas_ep->SetMarkerStyle(20);
        h_resGas_ep->SetMarkerColor(kBlack);
        h_resGas_ep->SetLineColor(kBlack);

        h_cell_ep->Add(h_c_ep_Yield, h_d_ep_Yield, 1, -1);
        h_cell_ep->Divide(h_H2Gas_ep);
        h_cell_ep->SetMarkerStyle(21);
        h_cell_ep->SetMarkerColor(kRed);
        h_cell_ep->SetLineColor(kRed);

        h_collimators_ep->Divide(h_d_ep_Yield,h_H2Gas_ep);
        h_collimators_ep->SetMarkerStyle(23);
        h_collimators_ep->SetMarkerColor(kBlue);
        h_collimators_ep->SetLineColor(kBlue);

        c1->cd(1);
        h_resGas_ep->SetAxisRange(0,0.2, "Y");
        h_resGas_ep->SetStats(0);
        h_resGas_ep->Draw("P E");
        h_cell_ep->Draw("P E SAME");
        h_collimators_ep->Draw("P E SAME");
        legend->SetHeader("Legend","C"); // option "C" allows to center the header
	    legend->AddEntry(h_resGas_ee,"(b)-(c) Res. Gas","p");
	    legend->AddEntry(h_cell_ee,"(c)-(d) Cell","p");
        legend->AddEntry(h_collimators_ee,"(d) Coll.","p");
        legend->Draw();
        c1->Print(outfile + ")");
        c1->Clear();

        TObjArray* arr = new TObjArray(0,0);
        (*arr).Add(h_H2Gas_ee);
        (*arr).Add(h_resGas_ee);
        (*arr).Add(h_cell_ee);
        (*arr).Add(h_collimators_ee);

        (*arr).Add(h_H2Gas_ep);
        (*arr).Add(h_resGas_ep);
        (*arr).Add(h_cell_ep);
        (*arr).Add(h_collimators_ep);

         TFile file1(root_outfile,"UPDATE");
        (*arr).Write();
	    file1.Close();
    }


    return 0;
}