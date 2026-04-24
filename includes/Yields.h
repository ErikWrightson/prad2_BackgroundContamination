/**
 * Header file for Yields object which evaluates reconstructed HyCal Events for e-e and e-p yields for background subtraction.
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 04.21.2026
 * @creation 04.16.2026
 */

#ifndef Yields_H
#define Yields_H

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
#include <TLegend.h>
#include <TMath.h>
#include <TColor.h>
#include <TString.h>
//Needed for reading in vector types from root files.
#include <TInterpreter.h>
#include "ROOT/TThreadedObject.hxx"

//Generally useful includes
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>

#include <sqlite3.h>
#include <iomanip>

using namespace std;

class Yields{

    public:
        TChain* chain;

        static constexpr Double_t rad2Deg = 180/TMath::Pi(); //Conversion from radians to degrees

        static constexpr Int_t  MAX_CLUSTERS = 400; //Maximum number of clusters.
        static constexpr Int_t EE_CUT_NUM = 4; //Number of cuts to apply for e-e
        static constexpr Int_t EP_CUT_NUM = 2; //Number of cuts to apply for e-p

        static constexpr Double_t M_p = 938.272; //Mass of Proton MeV/c^2
        static constexpr Double_t M_e = 0.511; //Mass of Electron MeV/c^2

        //Constructor that ensures the chain tree is set up.
        Yields(TChain* c, Int_t type, map<Int_t, Double_t>& m, bool a, bool g);

        void Evaluate();

        //Histogram handlers.
        void delete_Histos();
        void save_Histos(TString rootFile, bool first);
        void printPDF(TString pdfName,bool begin, bool end);
        TH1F* get_ee_YieldHisto();
        TH1F* get_ep_YieldHisto();

    private:
        
        Int_t curFileNum;
        Int_t type;
        Long64_t entries;
        bool all;
        bool gems;

        map<Int_t, Double_t> lcMap;
        Double_t lc;
        vector<TString> runlist;

        TString ee_cutNames[EE_CUT_NUM] = {"_noCut", "_expectedE","_coplanarity" ,"_elast"};
        TString ee_cut[EE_CUT_NUM] = {"None","Energy", "Coplanarity", "Elasticity"};
        
        TString ep_cutNames[EP_CUT_NUM] = {"_expectedE", "_numBlocks"};
        TString ep_cut[EP_CUT_NUM] = {"Energy", "Number of Blocks"};
        
        TString typeArr[4] = {"_a", "_b", "_c","_d"};

        Float_t prev_x[2];
        Float_t prev_y[2];
        Float_t prev_z;
        
        Float_t prev_x1[2];
        Float_t prev_y1[2];
        Float_t prev_z1;

        //Event information from root tree.
        Int_t evNum;
        Int_t totalE;
        Float_t EBeam;
        Int_t nClust;
        Float_t cl_x[MAX_CLUSTERS];
        Float_t cl_y[MAX_CLUSTERS];
        Float_t cl_z[MAX_CLUSTERS];
        Float_t cl_E[MAX_CLUSTERS];
        UChar_t cl_nblocks[MAX_CLUSTERS];
        UShort_t cl_center[MAX_CLUSTERS];
        UInt_t cl_flag[MAX_CLUSTERS];
        
        //GEM Matching Information
        UInt_t match_flag[MAX_CLUSTERS];
        Float_t matchGEMx[MAX_CLUSTERS][2];
        Float_t matchGEMy[MAX_CLUSTERS][2];
        Float_t matchGEMz[MAX_CLUSTERS][2];

        //Declare arrays of histograms that will be made before and after each cut.
        TH2F* h_ee_HC_XY[EE_CUT_NUM];
        TH2F* h_ee_EvTheta[EE_CUT_NUM];
        TH1F* h_ee_Yield[EE_CUT_NUM];
        TH1F* h_ee_YieldPerLC[EE_CUT_NUM];

        TH2F* h_ep_HC_XY[EP_CUT_NUM];
        TH2F* h_ep_EvTheta[EP_CUT_NUM];
        TH1F* h_ep_Yield[EP_CUT_NUM];
        TH1F* h_ep_YieldPerLC[EE_CUT_NUM];

        TH2F* h_eeCenters;

        Double_t ee_ExpectedE(Double_t theta);
        Double_t ep_ExpectedE(Double_t theta);

        vector<Double_t> findCenter(Float_t* pr_x, Float_t* pr_y, Int_t in1, Int_t in2);

        Double_t EnergyRes(Double_t energy);

        void settup_Histos(Int_t en);

        void fill_ee_Histos(Int_t c, Double_t* theta, Int_t index);
        void fill_ep_Histos(Int_t c, Double_t* theta, Int_t index);

        void find_Events();

        Float_t projToZPlane(Float_t nonZ, Float_t ogZ, Float_t newZ);
};

#endif