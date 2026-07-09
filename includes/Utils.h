/**
 * Header file for Utils class which provides useful utilities functions to this project.
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 04.21.2026
 * @creation 04.16.2026
 */

#ifndef Utils_H
#define Utils_H

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
#include <TPrincipal.h>
#include <TVectorD.h>
#include "ROOT/TThreadedObject.hxx"

//Generally useful includes
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>

#include <iomanip>

using namespace std;

class Utils{

    public:

        struct LineOfBestFit{
            Float_t x0;
            Float_t y0;
            Float_t z0;

            Float_t vx;
            Float_t vy;
            Float_t vz;
        };

        struct Point{
            Float_t x;
            Float_t y;
            Float_t z;
        };
        

        static LineOfBestFit FitLine(vector<Float_t> x, vector<Float_t> y, vector<Float_t> z);
        static Point ClosestApproachToZAxis(Utils::LineOfBestFit line);
        static Double_t Q_2(Float_t Ep, Float_t E_0, Double_t theta);
};

#endif