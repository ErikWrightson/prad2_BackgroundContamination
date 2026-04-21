/**
 * Yields class that defines the functionality for separating e-e and e-p and getting live charge normalized yields.
 *
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 04.21.2026 
 * @creation 04.16.2026
 */
#include "Yields.h"

/**
 * Creates a Yields object and properly sets up the needed reconstructed HyCal branches for the tree.
 *
 * @param c - the TChain to set up the branching and trees for.
 * @param t - the type of files this will be processing [a=0, b=1, c=2, d=3]
 * @param m - the LiveCharge map to use
 * @param a - Flag for weather or not all histograms should be filled or just the final ones.
 */
Yields::Yields(TChain* c, Int_t t, map<Int_t, Double_t>& m, bool a){

    curFileNum = -1;
    all = a;
    
    type = t;
    lcMap = m;
    
    chain = (TChain*) c;

    chain->SetMakeClass(1);

    //General Event Data
    chain->SetBranchAddress("event_num",       &evNum);
    chain->SetBranchAddress("total_energy",    &totalE);
    chain->SetBranchAddress("EBeam",           &EBeam);

    //HyCal Information
    chain->SetBranchAddress("n_clusters", &nClust);
    chain->SetBranchAddress("cl_x",       cl_x);       //Cluster x position
    chain->SetBranchAddress("cl_y",       cl_y);       //Cluster y position
    chain->SetBranchAddress("cl_z",       cl_z);       //Cluster z position
    chain->SetBranchAddress("cl_energy",  cl_E);       //Cluster energy
    chain->SetBranchAddress("cl_nblocks", cl_nblocks); //Number of blocks in the cluster
    chain->SetBranchAddress("cl_center",  cl_center);  //center module id for this cluster
    chain->SetBranchAddress("cl_flag",    cl_flag);    //Cluster flags

    entries = chain->GetEntries();

    chain->GetEntry(0);
    Int_t en = (Int_t) EBeam;

    cout<<"\t\t\t\t\t\tSet all Branch Addresses. Entries in chain = " << entries;

    settup_Histos(en);
        
}

/**
 * Calculates the expected energy of an electron in Moller (e-e) Scattering at a given angle and beam energy.
 * 
 * @param theta - the theta position of the electron on the calorimeter.
 *
 * @return - the expected energy if this is a Moller event
 */
Double_t Yields::ee_ExpectedE(Double_t theta){
    
    Double_t cosTheta = TMath::Cos(theta);
    Double_t cosTheta_2 = cosTheta*cosTheta; 

    Double_t num = M_e*(EBeam + M_e + (EBeam-M_e)*cosTheta_2);
    Double_t denom  = EBeam + M_e - (EBeam-M_e)*cosTheta_2;

    return num/denom;
}

/**
 * Calculates the expected energy of an electron in Mott (e-p) Scattering at a given angle and beam energy.
 * 
 * @param theta - the theta position of the electron on the calorimeter.
 *
 * @return - the expected energy if this is a Mott event
 */
Double_t Yields::ep_ExpectedE(Double_t theta){
    Double_t sinTheta = TMath::Sin(theta);
    Double_t sinTheta_2 = sinTheta*sinTheta;

    Double_t cosTheta = TMath::Cos(theta);
    Double_t cosTheta_2 = cosTheta*cosTheta;

    Double_t num = (EBeam+M_p)*(M_p*EBeam + M_e*M_e) + TMath::Sqrt(M_p*M_p - M_e*M_e*sinTheta_2)*(EBeam*EBeam - M_e*M_e)*cosTheta;
    Double_t denom = (EBeam+M_p)*(EBeam+M_p) - (EBeam*EBeam - M_e*M_e)*cosTheta_2;

    return num/denom;
}

/**
 * Gets the energy resolution at at the provided energy for the HyCal detector.
 *
 * @param energy - the energy to get the resolution of in MeV.
 *
 * @return - the energy resolution.
 */
Double_t Yields::EnergyRes(Double_t energy){
    return energy*0.026/TMath::Sqrt(energy/1000.0);
}

/**
 * Sets up the histograms for this object.
 *
 * @param en - energy to use for setting the upper limit of the histograms.
 */
void Yields::settup_Histos(Int_t en){

    Int_t numFiles = chain->GetListOfFiles()->GetEntries();

    h_eeCenters = new TH2F("h_eeCenters_type"+typeArr[type], "e-e Centers Type"+typeArr[type]+";x(mm);y(mm)", 120, -60, 60, 120, -60, 60);

    for(int i = 0; i < EE_CUT_NUM; i++){
        /*h_ee_HC_XY.emplace_back("h_ee_HC_XY_type"+typeArr[type]+ee_cutNames[i], "e-e HyCal XY Type"+typeArr[type]+" Cut: "+ee_cut[i]+";x(mm);y(mm)", 700, -700, 700, 700, -700, 700);
        h_ee_EvTheta.emplace_back("h_ee_EvTheta_type"+typeArr[type]+ee_cutNames[i], "e-e E vs. #theta Type"+typeArr[type]+" Cut: "+ee_cut[i]+";#theta (#circ);E (MeV)", 80, 0, 8, 1500, 0, 1500);
        h_ee_Yield.emplace_back("h_ee_Yield_type"+typeArr[type]+ee_cutNames[i], "e-e LiveCharge Normalized Yield vs. #theta Type"+typeArr[type]+" Cut: "+ee_cut[i]+ ";#theta (#circ);Counts", 80, 0, 8);*/
        h_ee_HC_XY[i] = new TH2F("h_ee_HC_XY_type"+typeArr[type]+ee_cutNames[i], "e-e HyCal XY Type"+typeArr[type]+" Cut: "+ee_cut[i]+";x(mm);y(mm)", 700, -700, 700, 700, -700, 700);
        h_ee_EvTheta[i] = new TH2F("h_ee_EvTheta_type"+typeArr[type]+ee_cutNames[i], "e-e E vs. #theta Type"+typeArr[type]+" Cut: "+ee_cut[i]+";#theta (#circ);E (MeV)", 60, 0, 6, en+400, 0, en+400);
        h_ee_Yield[i] = new TH1F("h_ee_Yield_type"+typeArr[type]+ee_cutNames[i], "e-e LiveCharge Normalized Yield vs. #theta Type"+typeArr[type]+" Cut: "+ee_cut[i]+ ";#theta (#circ);Counts", 60, 0, 6);
        h_ee_YieldPerLC[i] = new TH1F("h_ee_YieldPerLC_type"+typeArr[type]+ee_cutNames[i],"e-e Yield per LiveCharge Per File Type"+typeArr[type]+" Cut: "+ee_cut[i]+ ";File Number; Counts",numFiles+1, 0, numFiles+1);
    }

    for(int j; j < EP_CUT_NUM; j++){
        /*h_ep_HC_XY.emplace_back("h_ep_HC_XY_type"+typeArr[type]+ep_cutNames[j], "e-p HyCal XY Type"+typeArr[type]+" Cut: "+ep_cut[j]+";x(mm);y(mm)", 700, -700, 700, 700, -700, 700);
        h_ep_EvTheta.emplace_back("h_ep_EvTheta_type"+typeArr[type]+ep_cutNames[j], "e-p E vs. #theta Type"+typeArr[type]+" Cut: "+ep_cut[j]+ ";#theta (#circ);E (MeV)", 80, 0, 8, 1500, 0, 1500);  
        h_ep_Yield.emplace_back("h_ep_Yield_type"+typeArr[type]+ep_cutNames[j], "e-p LiveCharge Normalized Yield vs. #theta Type"+typeArr[type]+" Cut: "+ep_cut[j]+ ";#theta (#circ);Counts", 80, 0, 8);*/
        h_ep_HC_XY[j] = new TH2F("h_ep_HC_XY_type"+typeArr[type]+ep_cutNames[j], "e-p HyCal XY Type"+typeArr[type]+" Cut: "+ep_cut[j]+";x(mm);y(mm)", 700, -700, 700, 700, -700, 700);
        h_ep_EvTheta[j] = new TH2F("h_ep_EvTheta_type"+typeArr[type]+ep_cutNames[j], "e-p E vs. #theta Type"+typeArr[type]+" Cut: "+ep_cut[j]+ ";#theta (#circ);E (MeV)", 60, 0, 6, en+400, 0, en+400);  
        h_ep_Yield[j] = new TH1F("h_ep_Yield_type"+typeArr[type]+ep_cutNames[j], "e-p LiveCharge Normalized Yield vs. #theta Type"+typeArr[type]+" Cut: "+ep_cut[j]+ ";#theta (#circ);Counts", 60, 0, 6);
        h_ep_YieldPerLC[j] = new TH1F("h_ep_YieldPerLC_type"+typeArr[type]+ep_cutNames[j],"e-p Yield per LiveCharge Per File Type"+typeArr[type]+" Cut: "+ep_cut[j]+ ";File Number; Counts",numFiles+1, 0, numFiles+1);
    }
}

/**
 * Deletes all histograms in this object.
 */
void Yields::delete_Histos(){

    delete h_eeCenters;
    for(int i = 0; i < EE_CUT_NUM; i++){
        delete h_ee_HC_XY[i];
        delete h_ee_EvTheta[i];
        delete h_ee_Yield[i];
        delete h_ee_YieldPerLC[i];
    }

    for(int j = 0; j < EP_CUT_NUM; j++){
        delete h_ep_HC_XY[j];
        delete h_ep_EvTheta[j];
        delete h_ep_Yield[j];
        delete h_ep_YieldPerLC[j];
    }

}

/**
 * Saves all histograms in this object to the given ROOT file.
 *
 * @param rootFile - the ROOT File to save to.
 * @param first - notes if the file needs to be recreated or just updated.
 */
void Yields::save_Histos(TString rootFile, bool first){
    TObjArray* arr = new TObjArray(0,0);

    if(h_eeCenters->GetEntries()){(*arr).Add(h_eeCenters);}

    for(int i = 0; i < EE_CUT_NUM; i++){
        if(h_ee_HC_XY[i]->GetEntries()){(*arr).Add(h_ee_HC_XY[i]);}
        if(h_ee_EvTheta[i]->GetEntries()){(*arr).Add(h_ee_EvTheta[i]);}
        if(h_ee_Yield[i]->GetEntries()){(*arr).Add(h_ee_Yield[i]);}
        if(h_ee_YieldPerLC[i]->GetEntries()){(*arr).Add(h_ee_YieldPerLC[i]);}
    }

    for(int j = 0; j < EP_CUT_NUM; j++){
        if(h_ep_HC_XY[j]->GetEntries()){(*arr).Add(h_ep_HC_XY[j]);}
        if(h_ep_EvTheta[j]->GetEntries()){(*arr).Add(h_ep_EvTheta[j]);}
        if(h_ep_Yield[j]->GetEntries()){(*arr).Add(h_ep_Yield[j]);}
        if(h_ep_YieldPerLC[j]->GetEntries()){(*arr).Add(h_ep_YieldPerLC[j]);}
    }

    if(first){
        TFile file1(rootFile,"RECREATE");
        (*arr).Write();
	    file1.Close();
    }
    else{
        TFile file2(rootFile,"UPDATE");
        (*arr).Write();
	    file2.Close();
    }
	
    delete_Histos();
}

/**
 * Evaluate the TChain to separate the e-e and e-p events.
 */
void Yields::Evaluate(){

    TString curFileName = "";
    Int_t runNum = 0;

    prev_x[0] = -100000;
    prev_x[1]= -100000;
    prev_y[0] = -100000;
    prev_y[1] = -100000;


    for(Long64_t i = 0; i < entries; i++){
        //chain->LoadTree(); //Load the current tree
        if((i)%10000 == 0 || (entries-i)<10000){cout<<"\rType" << typeArr[type] << " Events: " << i << "/" << entries << flush;}

        //Check if we are on a new file and if so get the liveCharge for this run.
        TString temp = (chain->GetCurrentFile())->GetName();
        if(curFileName != temp){
            curFileNum++;
            curFileName = temp;

            TString s = curFileName;
            s.Remove(0,s.Last('/')+1);
            s.ReplaceAll("prad2Replay_", "");
            s.ReplaceAll(".root","");

            runNum = s.Atoi();
            lc = lcMap.at(runNum);
        }

        chain->GetEntry(i);

        find_Events();

    }

}

/**
 * Fills the e-e histograms for the current cut level with events that survived the cut.
 *
 * @param c - the cut number
 * @param theta - the theta of each cluster in this event.
 * @param index - the index of the hit to put in the histograms.
 */
void Yields::fill_ee_Histos(Int_t c, Double_t* theta, Int_t index){
    h_ee_HC_XY[c]->Fill(cl_x[index],cl_y[index]);
    h_ee_EvTheta[c]->Fill(theta[index]*rad2Deg, cl_E[index]);
    h_ee_Yield[c]->Fill(theta[index]*rad2Deg,1/lc*1000);
    h_ee_YieldPerLC[c]->Fill(curFileNum, 1/lc*1000);
}

/**
 * Fills the e-p histograms for the current cut level with events that survived the cut.
 *
 * @param c - the cut number
 * @param theta - the theta of each cluster in this event.
 * @param index - the index of the event to add to the histograms.
 */
void Yields::fill_ep_Histos(Int_t c, Double_t* theta, Int_t index){
    h_ep_HC_XY[c]->Fill(cl_x[index],cl_y[index]);
    h_ep_Yield[c]->Fill(theta[index]*rad2Deg,1/lc*1000);
    h_ep_EvTheta[c]->Fill(theta[index]*rad2Deg, cl_E[index]);
    h_ep_YieldPerLC[c]->Fill(curFileNum, 1/lc*1000);

}

/**
 * Finds the intersection of the current moller pair and the previous one.
 *
 * @param pr_x - the two x positions of the previous double arm moller pair.
 * @param pr_y - the two y positions of the previous double arm moller pair.
 */
vector<Double_t> Yields::findCenter(Float_t* pr_x, Float_t* pr_y, Int_t in0, Int_t in1){
    
    Float_t m_1 = (cl_y[in1]-cl_y[in0])/(cl_x[in1]-cl_x[in0]);
    Float_t b_1 = -1.0*(cl_x[in0]*m_1) +(cl_y[in0]);

    Float_t m_2 = (pr_y[1]-pr_y[0])/(pr_x[1]-pr_x[0]);
    Float_t b_2 = -1.0*(pr_x[0]*m_2) + (pr_y[0]);

    Double_t A[2][2] = {{-1.0*m_1, 1}, {-1.0*m_2, 1}};
    Double_t B[2] = {b_1, b_2};

    Double_t oneOvdet = 1/(A[0][0]*A[1][1] - A[0][1]*A[1][0]);
    Double_t A_inv[2][2] = {{A[1][1], -1*A[0][1]},{-1*A[1][0], A[0][0]}};

    Double_t x_cen = oneOvdet*(A_inv[0][0]*B[0] + A_inv[0][1]*B[1]);
    Double_t y_cen = oneOvdet*(A_inv[1][0]*B[0] + A_inv[1][1]*B[1]);

    vector<Double_t> ans;
    ans.push_back(x_cen);
    ans.push_back(y_cen);

    return ans;
}

/**
 * Helper method for looping through each event and finding if it qualifies as an e-e or e-p event.
 */
void Yields::find_Events(){
    
    Double_t theta[nClust];
    Double_t phi[nClust];
    Double_t expE[nClust];
    Double_t ep_expE[nClust];
    vector<Int_t> ee_passedEHits;
    vector<Int_t> ee_passedCopHits;
    vector<Int_t> ee_passedElastHits;

    for(Int_t j = 0; j < nClust; j++){

        //Geometric cut on the Absorber
        if(TMath::Abs(cl_x[j]) > 20.77 * 2.25 || TMath::Abs(cl_y[j]) > 20.75*2.25){
            
            theta[j] = TMath::ATan2(TMath::Sqrt(cl_x[j]*cl_x[j]+cl_y[j]*cl_y[j]),cl_z[j]);
            
            //Geometric cut on high angles. 
            if(theta[j]*rad2Deg<6){
            
                phi[j] = TMath::ATan2(cl_y[j],cl_x[j]);
                if(phi[j]<0){
                    phi[j] += 2*TMath::Pi();
                }
                phi[j] = phi[j] * rad2Deg;
                expE[j] = ee_ExpectedE(theta[j]);
                ep_expE[j] = ep_ExpectedE(theta[j]);

                //Find e-e events
                //No cut
                if(all){fill_ee_Histos(0, theta, j);}

                //Expected ee Energy Cut
                if((TMath::Abs(cl_E[j] - expE[j]) <  3.0*EnergyRes(expE[j]))){
                
                    if(all){fill_ee_Histos(1, theta, j);}
                    ee_passedEHits.push_back(j);

                    for(Int_t k = 0; k < (Int_t) ee_passedEHits.size() && k != j; k++){

                        //Coplanarity Cut for double arm moller
                        if((TMath::Abs(TMath::Abs(phi[ee_passedEHits.at(k)]-phi[j])-180) < 10)){ //Check that the double arm mollers are coplanar
                    
                            //If this was the first pair to pass the coplanarity cut, make sure to put both hits in the histogram.
                            if(ee_passedCopHits.size()==0){
                                ee_passedCopHits.push_back(ee_passedEHits.at(k));
                                if(all){fill_ee_Histos(2, theta, ee_passedEHits.at(k));}
                            }
                            ee_passedCopHits.push_back(j);
                            if(all){fill_ee_Histos(2, theta, j);}

                            //Cut for elasticity
                            if(TMath::Abs(cl_E[ee_passedEHits.at(k)] + cl_E[j] - EBeam - M_e) < 3*EnergyRes(EBeam)){
                        
                                //If this was the first pair to pass the elasticity cut, make sure to put both hits in the histogram.
                                if(ee_passedElastHits.size()==0){
                                    ee_passedElastHits.push_back(ee_passedEHits.at(k));
                                    fill_ee_Histos(3, theta, ee_passedEHits.at(k));
                                }
                        
                                fill_ee_Histos(3, theta, j);
                                ee_passedElastHits.push_back(ee_passedEHits.at(k));

                                if(prev_x[0] > -10000 && prev_x[1] > -10000 && prev_y[0] > -10000 && prev_y[1] > -10000){
                                    vector<Double_t> center = findCenter(prev_x, prev_y, j, ee_passedEHits.at(k));

                                    h_eeCenters->Fill(center.at(0), center.at(1));

                                    prev_x[0] = -100000;
                                    prev_y[0] = -100000;
                                    prev_x[1] = -100000;
                                    prev_y[1] = -100000;
                                }
                                else{
                                    prev_x[0] = cl_x[ee_passedEHits.at(k)];
                                    prev_x[1] = cl_x[j];

                                    prev_y[0] = cl_y[ee_passedEHits.at(k)];
                                    prev_y[1] = cl_y[j];
                                }
                            }
                        }
                    }
                }

                //Find e-p events
                if(TMath::Abs(cl_E[j]-ep_expE[j]) < 3.0*EnergyRes(ep_expE[j])){
                    if(all){fill_ep_Histos(0, theta, j);}

                    //Number of blocks cut
                    if(cl_nblocks[j]>=5){
                        fill_ep_Histos(1, theta, j);
                    }
                }
            }
        }
    }
}

/**
 * Print the relevant histograms to a pdf.
 * 
 * @param pdfName - the name of the pdf to save to
 * @param begin - Notes if this is the first thing being added to this pdf.
 * @param end - Notes is this is the last thing being added to this pdf.
 */
void Yields::printPDF(TString pdfName,bool begin, bool end){
    TCanvas *c = new TCanvas("c"+typeArr[type], "Type" + typeArr[type] + "_Yield_Canvas",1000,1000);

    if(all){
        for(int i = 0; i < EE_CUT_NUM; i++){

            c->Divide(2,2);
            c->cd(1);
            h_ee_HC_XY[i]->Draw("COLZ");
            c->cd(2);
            h_ee_EvTheta[i]->Draw("COLZ");
            c->cd(3);
            gPad->SetLogy(1);
            h_ee_Yield[i]->Draw("E");
            c->cd(4);
            h_ee_YieldPerLC[i]->Draw("E");
            if(i==0 && begin){
                c->Print(pdfName+"(");
            }
            else{
                c->Print(pdfName);
            }
		    c->Clear();

            if(i == 3){
                c->cd(1);
                h_eeCenters->Draw("COLZ");
                c->Print(pdfName);
                c->Clear();

            }
        }
        for(int j = 0; j < EP_CUT_NUM; j++){
            c->Divide(2,2);
            c->cd(1);
            h_ep_HC_XY[j]->Draw("COLZ");
            c->cd(2);
            h_ep_EvTheta[j]->Draw("COLZ");
            c->cd(3);
            gPad->SetLogy(1);
            h_ep_Yield[j]->Draw("E");
            c->cd(4);
            h_ep_YieldPerLC[j]->Draw("E");
            if(j==1 && end){
                c->Print(pdfName + ")");
            }
            else{
                c->Print(pdfName);
            }
            c->Clear();
        }
    }
    else{
        c->Divide(2,2);
        c->cd(1);
        h_ee_HC_XY[3]->Draw("COLZ");
        c->cd(2);
        h_ee_EvTheta[3]->Draw("COLZ");
        c->cd(3);
        gPad->SetLogy(1);
        h_ee_Yield[3]->Draw("E");
        c->cd(4);
        h_ee_YieldPerLC[3]->Draw("E");
        if(begin){
                c->Print(pdfName+"(");
        }
        else{
            c->Print(pdfName);
        }
        c->Clear();


        c->cd(1);
        h_eeCenters->Draw("COLZ");
        c->Print(pdfName);
        c->Clear();

        c->Divide(2,2);
        c->cd(1);
        h_ep_HC_XY[1]->Draw("COLZ");
        c->cd(2);
        h_ep_EvTheta[1]->Draw("COLZ");
        c->cd(3);
        gPad->SetLogy(1);
        h_ep_Yield[1]->Draw("E");
        c->cd(4);
        h_ep_YieldPerLC[1]->Draw("E");
        if(end){
            c->Print(pdfName + ")");
        }
        else{
            c->Print(pdfName);
        }
        c->Clear();


    }
    
}

/**
 * Gets a clone of the e-e Yield Histogram after all cuts.
 *
 * @return - the clone of the e-e Yield Histogram after all cuts.
 */
TH1F* Yields::get_ee_YieldHisto(){
    return (TH1F*) h_ee_Yield[EE_CUT_NUM-1]->Clone();
}

/**
 * Gets a clone of the e-p Yield Histogram after all cuts.
 *
 * @return - the clone of the e-p Yield Histogram after all cuts.
 */
TH1F* Yields::get_ep_YieldHisto(){
    return (TH1F*) h_ep_Yield[EP_CUT_NUM-1]->Clone();
}