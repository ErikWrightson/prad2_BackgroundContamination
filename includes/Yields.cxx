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
 * @param a - Flag for whether or not all histograms should be filled or just the final ones.
 * @param g - Flag indicating whether or not to use the GEM positions for the Moller center calculation.
 * @param h - Flag indicating if history should be searched or not.
 */
Yields::Yields(TChain* c, Int_t t, map<Int_t, Double_t>& m, bool a, bool g, bool h, Float_t EB, bool vert){

    gErrorIgnoreLevel = 3000;

    curFileNum = -1;
    all = a;
    gems = g;
    hist = h;
    z = vert;
    
    type = t;
    lcMap = m;

    EBeam = EB;
    
    chain = (TChain*) c;

    chain->SetMakeClass(1);

    //General Event Data
    chain->SetBranchAddress("event_num",       &evNum);
    chain->SetBranchAddress("total_energy",    &totalE);
    //chain->SetBranchAddress("EBeam",           &EBeam);

    //HyCal Information
    chain->SetBranchAddress("n_clusters", &nClust);
    chain->SetBranchAddress("cl_x",       cl_x);       //Cluster x position
    chain->SetBranchAddress("cl_y",       cl_y);       //Cluster y position
    chain->SetBranchAddress("cl_z",       cl_z);       //Cluster z position
    chain->SetBranchAddress("cl_energy",  cl_E);       //Cluster energy
    chain->SetBranchAddress("cl_nblocks", cl_nblocks); //Number of blocks in the cluster
    chain->SetBranchAddress("cl_center",  cl_center);  //center module id for this cluster
    chain->SetBranchAddress("cl_flag",    cl_flag);    //Cluster flags

    if(g){
        chain->SetBranchAddress("matchFlag", match_flag); //Matching Flag bit 0 for GEM0, bit 1 for GEM1 etc.
        chain->SetBranchAddress("mHit_gx", matchGEMx);   //The x-coordinate of the Matches found on each GEM plane.
        chain->SetBranchAddress("mHit_gy", matchGEMy);   //The y-coordinate of the Matches found on each GEM plane.
        chain->SetBranchAddress("mHit_gz", matchGEMz);   //The z-coordinate of the Matches found on each GEM plane.
    }

    entries = chain->GetEntries();

    Int_t en = (Int_t) EBeam;

    cout<<"Set all Branch Addresses. Entries in chain = " << entries;

    setup_Histos(en);
        
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
void Yields::setup_Histos(Int_t en){

    Int_t numFiles = chain->GetListOfFiles()->GetEntries();

    h_eeCenters = new TH2F("h_eeCenters_type"+typeArr[type], "e-e Centers Type"+typeArr[type]+";x(mm);y(mm)", 240, -60, 60, 240, -60, 60);

    if(gems){
        h_eeCenters_GEM[0]= new TH2F("h_eeCenters_GEM0_type"+typeArr[type], "e-e Centers Type"+typeArr[type]+" GEM Index 0"+";x(mm);y(mm)", 240, -60, 60, 240, -60, 60);
        h_eeCenters_GEM[1]= new TH2F("h_eeCenters_GEM1_type"+typeArr[type], "e-e Centers Type"+typeArr[type]+" GEM Index 1"+";x(mm);y(mm)", 240, -60, 60, 240, -60, 60);
        h_ee_zVert_DoubleArmMoller[0] = new TH1F("h_ee_zVertz_DoubleArmMoller_type"+typeArr[type]+"_coplanarity","e-e Reconstructed Distance from Z Vertex Type"+typeArr[type]+" Cut: "+ee_cut[EE_CUT_NUM-2],1000,5000,6000);
        h_ee_zVert_DoubleArmMoller[1] = new TH1F("h_ee_zVertz_DoubleArmMoller_type"+typeArr[type]+"_elast","e-e Reconstructed Distance from Z Vertex Type"+typeArr[type]+" Cut: "+ee_cut[EE_CUT_NUM-1],1000,5000,6000);
    }

    for(int i = 0; i < EE_CUT_NUM; i++){
        /*h_ee_HC_XY.emplace_back("h_ee_HC_XY_type"+typeArr[type]+ee_cutNames[i], "e-e HyCal XY Type"+typeArr[type]+" Cut: "+ee_cut[i]+";x(mm);y(mm)", 700, -700, 700, 700, -700, 700);
        h_ee_EvTheta.emplace_back("h_ee_EvTheta_type"+typeArr[type]+ee_cutNames[i], "e-e E vs. #theta Type"+typeArr[type]+" Cut: "+ee_cut[i]+";#theta (#circ);E (MeV)", 80, 0, 8, 1500, 0, 1500);
        h_ee_Yield.emplace_back("h_ee_Yield_type"+typeArr[type]+ee_cutNames[i], "e-e LiveCharge Normalized Yield vs. #theta Type"+typeArr[type]+" Cut: "+ee_cut[i]+ ";#theta (#circ);Counts", 80, 0, 8);*/
        h_ee_HC_XY[i] = new TH2F("h_ee_HC_XY_type"+typeArr[type]+ee_cutNames[i], "e-e HyCal XY Type"+typeArr[type]+" Cut: "+ee_cut[i]+";x(mm);y(mm)", 700, -700, 700, 700, -700, 700);
        h_ee_EvTheta[i] = new TH2F("h_ee_EvTheta_type"+typeArr[type]+ee_cutNames[i], "e-e E vs. #theta Type"+typeArr[type]+" Cut: "+ee_cut[i]+";#theta (#circ);E (MeV)", 60, 0, 6, en+400, 0, en+400);
        h_ee_Yield[i] = new TH1F("h_ee_Yield_type"+typeArr[type]+ee_cutNames[i], "e-e LiveCharge Normalized Yield vs. #theta Type"+typeArr[type]+" Cut: "+ee_cut[i]+ ";#theta (#circ);Counts", 60, 0, 6);
        h_ee_YieldPerLC[i] = new TH1F("h_ee_YieldPerLC_type"+typeArr[type]+ee_cutNames[i],"e-e Yield per LiveCharge Per File Type"+typeArr[type]+" Cut: "+ee_cut[i]+ ";File Number; Counts",numFiles+1, 0, numFiles+1);

        if(gems && z){
            h_ee_zVert[i] = new TH1F("h_ee_zVert_type"+typeArr[type]+ee_cutNames[i],"e-e Reconstructed Z Vertex Type"+typeArr[type]+" Cut: "+ee_cut[i]+";z (mm);Counts",1000,-2000,7000);
        }
    }

    for(int j; j < EP_CUT_NUM; j++){
        /*h_ep_HC_XY.emplace_back("h_ep_HC_XY_type"+typeArr[type]+ep_cutNames[j], "e-p HyCal XY Type"+typeArr[type]+" Cut: "+ep_cut[j]+";x(mm);y(mm)", 700, -700, 700, 700, -700, 700);
        h_ep_EvTheta.emplace_back("h_ep_EvTheta_type"+typeArr[type]+ep_cutNames[j], "e-p E vs. #theta Type"+typeArr[type]+" Cut: "+ep_cut[j]+ ";#theta (#circ);E (MeV)", 80, 0, 8, 1500, 0, 1500);  
        h_ep_Yield.emplace_back("h_ep_Yield_type"+typeArr[type]+ep_cutNames[j], "e-p LiveCharge Normalized Yield vs. #theta Type"+typeArr[type]+" Cut: "+ep_cut[j]+ ";#theta (#circ);Counts", 80, 0, 8);*/
        h_ep_HC_XY[j] = new TH2F("h_ep_HC_XY_type"+typeArr[type]+ep_cutNames[j], "e-p HyCal XY Type"+typeArr[type]+" Cut: "+ep_cut[j]+";x(mm);y(mm)", 700, -700, 700, 700, -700, 700);
        h_ep_EvTheta[j] = new TH2F("h_ep_EvTheta_type"+typeArr[type]+ep_cutNames[j], "e-p E vs. #theta Type"+typeArr[type]+" Cut: "+ep_cut[j]+ ";#theta (#circ);E (MeV)", 60, 0, 6, en+400, 0, en+400);  
        h_ep_Yield[j] = new TH1F("h_ep_Yield_type"+typeArr[type]+ep_cutNames[j], "e-p LiveCharge Normalized Yield vs. #theta Type"+typeArr[type]+" Cut: "+ep_cut[j]+ ";#theta (#circ);Counts", 60, 0, 6);
        h_ep_YieldPerLC[j] = new TH1F("h_ep_YieldPerLC_type"+typeArr[type]+ep_cutNames[j],"e-p Yield per LiveCharge Per File Type"+typeArr[type]+" Cut: "+ep_cut[j]+ ";File Number; Counts",numFiles+1, 0, numFiles+1);

        if(gems && z){
            h_ep_zVert[j] = new TH1F("h_ep_zVert_type"+typeArr[type]+ep_cutNames[j],"e-p Reconstructed Z Vertex Type"+typeArr[type]+" Cut: "+ep_cut[j]+";z (mm);Counts",1000,-2000,7000);
        }
    }
}

/**
 * Deletes all histograms in this object.
 */
void Yields::delete_Histos(){

    delete h_eeCenters;
    delete h_epToee_ratio;

    if(gems){
        delete h_eeCenters_GEM[0];
        delete h_eeCenters_GEM[1];
    }
    for(int i = 0; i < EE_CUT_NUM; i++){
        delete h_ee_HC_XY[i];
        delete h_ee_EvTheta[i];
        delete h_ee_Yield[i];
        delete h_ee_YieldPerLC[i];

        if(gems && z){
            delete h_ee_zVert[i];
        }
    }

    for(int j = 0; j < EP_CUT_NUM; j++){
        delete h_ep_HC_XY[j];
        delete h_ep_EvTheta[j];
        delete h_ep_Yield[j];
        delete h_ep_YieldPerLC[j];

        if(gems && z){
            delete h_ep_zVert[j];
        }
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
    if(h_epToee_ratio->GetEntries()){(*arr).Add(h_epToee_ratio);}

    if(gems){
        if(h_eeCenters_GEM[0]->GetEntries()){(*arr).Add(h_eeCenters_GEM[0]);}
        if(h_eeCenters_GEM[1]->GetEntries()){(*arr).Add(h_eeCenters_GEM[1]);}
        if(h_ee_zVert_DoubleArmMoller[0]->GetEntries()){(*arr).Add(h_ee_zVert_DoubleArmMoller[0]);}
        if(h_ee_zVert_DoubleArmMoller[1]->GetEntries()){(*arr).Add(h_ee_zVert_DoubleArmMoller[1]);}
    }

    for(int i = 0; i < EE_CUT_NUM; i++){
        if(h_ee_HC_XY[i]->GetEntries()){(*arr).Add(h_ee_HC_XY[i]);}
        if(h_ee_EvTheta[i]->GetEntries()){(*arr).Add(h_ee_EvTheta[i]);}
        if(h_ee_Yield[i]->GetEntries()){(*arr).Add(h_ee_Yield[i]);}
        if(h_ee_YieldPerLC[i]->GetEntries()){(*arr).Add(h_ee_YieldPerLC[i]);}

        if(gems && z && h_ee_zVert[i]->GetEntries()){(*arr).Add(h_ee_zVert[i]);}
    }

    for(int j = 0; j < EP_CUT_NUM; j++){
        if(h_ep_HC_XY[j]->GetEntries()){(*arr).Add(h_ep_HC_XY[j]);}
        if(h_ep_EvTheta[j]->GetEntries()){(*arr).Add(h_ep_EvTheta[j]);}
        if(h_ep_Yield[j]->GetEntries()){(*arr).Add(h_ep_Yield[j]);}
        if(h_ep_YieldPerLC[j]->GetEntries()){(*arr).Add(h_ep_YieldPerLC[j]);}

        if(gems && z && h_ep_zVert[j]->GetEntries()){(*arr).Add(h_ep_zVert[j]);}
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

    prev_x[0] =  -100000;
    prev_x[1] =  -100000;
    prev_y[0] =  -100000;
    prev_y[1] =  -100000;
    prev_z    =  -100000;
    
    prev_x1[0] = -100000;
    prev_x1[1] = -100000; 
    prev_y1[0] = -100000;
    prev_y1[1] = -100000;
    prev_z1    = -100000;

    cout<<endl;
    for(Long64_t i = 0; i < entries; i++){
        //chain->LoadTree(); //Load the current tree
        if((i+1)%10000 == 0 || (entries-i+1)<10000){
            cout<<"\rType" << typeArr[type] << " Events: " << i+1 << "/" << entries << flush;
            if(i+1 == entries){
                cout<<endl;
            }
        }

        //Check if we are on a new file and if so get the liveCharge for this run.
        TString temp = (chain->GetCurrentFile())->GetName();
        if(curFileName != temp){
            curFileNum++;
            curFileName = temp;

            TString s = curFileName;
            s.Remove(0,s.Last('/')+1);
            //s.ReplaceAll("prad2Replay_", "");
            Int_t start = s.Index("_") + 1;
            Int_t end = s.Index(".");
            s = s(start, end-start);
            //s.ReplaceAll("prad_", "");
            //s.ReplaceAll(".root","");

            runlist.push_back(s);

            runNum = s.Atoi();
            lc = lcMap.at(runNum);
        }

        chain->GetEntry(i);

        if(gems){
            find_Events_wGEMs();
        }
        else{
            find_Events_OnlyHyCal();
        }

    }
    
    for(Int_t q = 0; q < EE_CUT_NUM; q++){
        h_ee_Yield[q]->Scale(1.0/runlist.size());
    }
    for(Int_t w = 0; w < EP_CUT_NUM; w++){
        h_ep_Yield[w]->Scale(1.0/runlist.size());
    }

}

/**
 * Fills the e-e histograms for the current cut level with events that survived the cut.
 *
 * @param c - the cut number
 * @param theta - the theta array of each cluster in this event.
 * @param index - the index of the hit to put in the histograms.
 * @param v - the reconstructed z vertex assuming it came from the beamline
 */
void Yields::fill_ee_Histos(Int_t c, Double_t* theta, Int_t index, Float_t v){
    h_ee_HC_XY[c]->Fill(cl_x[index],cl_y[index]);
    h_ee_EvTheta[c]->Fill(theta[index]*rad2Deg, cl_E[index]);
    h_ee_Yield[c]->Fill(theta[index]*rad2Deg,1/lc*1000);
    h_ee_YieldPerLC[c]->Fill(curFileNum, 1/lc*1000);

    if(gems && z){
        h_ee_zVert[c]->Fill(v, 1/lc*1000);
    }
}

/**
 * Fills the e-p histograms for the current cut level with events that survived the cut.
 *
 * @param c - the cut number
 * @param theta - the theta array of each cluster in this event.
 * @param index - the index of the event to add to the histograms.
 * @param v - the reconstructed z vertex assuming it came from the beamline
 */
void Yields::fill_ep_Histos(Int_t c, Double_t* theta, Int_t index, Float_t v){
    h_ep_HC_XY[c]->Fill(cl_x[index],cl_y[index]);
    h_ep_Yield[c]->Fill(theta[index]*rad2Deg,1/lc*1000);
    h_ep_EvTheta[c]->Fill(theta[index]*rad2Deg, cl_E[index]);
    h_ep_YieldPerLC[c]->Fill(curFileNum, 1/lc*1000);

    if(gems && z){
        if(c > EP_CUT_NUM-2 && theta[index]*rad2Deg>2){
            h_ep_zVert[c]->Fill(v,1/lc*1000);
        }
        else{
            h_ep_zVert[c]->Fill(v,1/lc*1000);
        }
    }

}

/**
 * Finds the intersection of the current moller pair and the previous one.
 *
 * @param pr_x - the two x positions of the previous double arm moller pair.
 * @param pr_y - the two y positions of the previous double arm moller pair.
 */
vector<Double_t> Yields::findCenter(Float_t* pr_x, Float_t* pr_y, Float_t x0, Float_t y0, Float_t x1, Float_t y1){
    
    Float_t m_1 = (y1-y0)/(x1-x0);
    Float_t b_1 = -1.0*(x0*m_1) +(y0);

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
 * Helper method for looping through each event and finding if it qualifies as an e-e or e-p event only using information from HyCal.
 */
void Yields::find_Events_OnlyHyCal(){
    
    Double_t theta[nClust];
    Double_t phi[nClust];
    Double_t expE[nClust];
    Double_t ep_expE[nClust];

    Float_t vert[nClust];

    vector<Int_t> ee_passedEHits;
    vector<Int_t> ee_passedCopHits;
    vector<Int_t> ee_passedElastHits;

    for(Int_t j = 0; j < nClust; j++){

        //Geometric cut on the Absorber
        if(TMath::Abs(cl_x[j]) > 20.77 * 2.25 || TMath::Abs(cl_y[j]) > 20.75*2.25){
            
            theta[j] = TMath::ATan2(TMath::Sqrt(cl_x[j]*cl_x[j]+cl_y[j]*cl_y[j]),cl_z[j]);
            
            //Geometric cut on high angles. 
            if(theta[j]*rad2Deg<6){

                vert[j] = find_VertZ_beamline(j);
            
                phi[j] = TMath::ATan2(cl_y[j],cl_x[j]);
                if(phi[j]<0){
                    phi[j] += 2*TMath::Pi();
                }
                phi[j] = phi[j] * rad2Deg;
                expE[j] = ee_ExpectedE(theta[j]);
                ep_expE[j] = ep_ExpectedE(theta[j]);

                //Find e-e events
                //No cut
                if(all){fill_ee_Histos(0, theta, j, vert[j]);}

                //Expected ee Energy Cut
                if((TMath::Abs(cl_E[j] - expE[j]) <  3.0*EnergyRes(expE[j]))){
                
                    if(all){fill_ee_Histos(1, theta, j, vert[j]);}
                    ee_passedEHits.push_back(j);

                    for(Int_t k = 0; k < (Int_t) ee_passedEHits.size() && k != j; k++){

                        //Coplanarity Cut for double arm moller
                        if((TMath::Abs(TMath::Abs(phi[ee_passedEHits.at(k)]-phi[j])-180) < 10)){ //Check that the double arm mollers are coplanar
                    
                            //If this was the first pair to pass the coplanarity cut, make sure to put both hits in the histogram.
                            if(ee_passedCopHits.size()==0){
                                ee_passedCopHits.push_back(ee_passedEHits.at(k));
                                if(all){fill_ee_Histos(2, theta, ee_passedEHits.at(k), vert[ee_passedEHits.at(k)]);}
                            }
                            ee_passedCopHits.push_back(j);
                            if(all){fill_ee_Histos(2, theta, j, vert[j]);}

                            //Cut for elasticity
                            if(TMath::Abs(cl_E[ee_passedEHits.at(k)] + cl_E[j] - EBeam - M_e) < 3*EnergyRes(EBeam)){
                        
                                //If this was the first pair to pass the elasticity cut, make sure to put both hits in the histogram.
                                if(ee_passedElastHits.size()==0){
                                    ee_passedElastHits.push_back(ee_passedEHits.at(k));
                                    fill_ee_Histos(3, theta, ee_passedEHits.at(k), vert[ee_passedEHits.at(k)]);
                                }
                        
                                fill_ee_Histos(3, theta, j, vert[j]);
                                ee_passedElastHits.push_back(ee_passedEHits.at(k));

                                if(prev_x[0] > -10000 && prev_x[1] > -10000 && prev_y[0] > -10000 && prev_y[1] > -10000){
                                    vector<Double_t> centerHC = findCenter(prev_x, prev_y, cl_x[j], cl_y[j], cl_x[ee_passedEHits.at(k)], cl_y[ee_passedEHits.at(k)]);

                                    h_eeCenters->Fill(centerHC.at(0), centerHC.at(1));

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
                    if(all){fill_ep_Histos(0, theta, j, vert[j]);}

                    //Number of blocks cut
                    if(cl_nblocks[j]>=5){
                        fill_ep_Histos(1, theta, j, vert[j]);
                    }
                }
            }
        }
    }
}

/**
 * Helper method for looping through each event and finding if it qualifies as an e-e or e-p event while using information from the GEMs.
 */
void Yields::find_Events_wGEMs(){
    
    Double_t theta[nClust];
    Double_t phi[nClust];
    //Double_t phi_GEM[nClust][2];
    Double_t expE[nClust];
    Double_t ep_expE[nClust];

    Double_t vert[nClust];
    Double_t vert_DoubleArmMoller[nClust];

    vector<Int_t> ee_passedEHits;
    vector<Int_t> ee_passedCopHits;
    vector<Int_t> ee_passedElastHits;

    for(Int_t j = 0; j < nClust; j++){

        //Geometric cut on the Absorber
        if(TMath::Abs(cl_x[j]) > 20.77 * 2.25 || TMath::Abs(cl_y[j]) > 20.75*2.25){
            
            theta[j] = TMath::ATan2(TMath::Sqrt(cl_x[j]*cl_x[j]+cl_y[j]*cl_y[j]),cl_z[j]);
            
            //Geometric cut on high angles. 
            if(theta[j]*rad2Deg<6){

                vert[j] = find_VertZ_beamline(j); //find the z vertex of this event assuming it came from the beamline.

                phi[j] = TMath::ATan2(cl_y[j],cl_x[j]);//(matchGEMx[j][1],matchGEMy[j][1]);
                if(phi[j]<0){
                    phi[j] += 2*TMath::Pi();
                }
                phi[j] = phi[j] * rad2Deg;

                /*phi_GEM[j][0] = TMath::ATan2(matchGEMy[j][0],matchGEMx[j][0]);
                if(phi_GEM[j][0]<0){
                    phi_GEM[j][0] += 2*TMath::Pi();
                }
                phi_GEM[j][0] = phi_GEM[j][0] * rad2Deg;

                phi_GEM[j][1] = TMath::ATan2(matchGEMy[j][1],matchGEMx[j][1]);
                if(phi_GEM[j][1]<0){
                    phi_GEM[j][1] += 2*TMath::Pi();
                }
                phi_GEM[j][1] = phi_GEM[j][1] * rad2Deg;*/
                 
                expE[j] = ee_ExpectedE(theta[j]);
                ep_expE[j] = ep_ExpectedE(theta[j]);

                //Find e-e events
                //No cut
                if(all){fill_ee_Histos(0, theta, j, vert[j]);}

                //Expected ee Energy Cut and ensure that this hit has a match on both GEMs
                if((TMath::Abs(cl_E[j] - expE[j]) <  3.0*EnergyRes(expE[j])) && ((match_flag[j] & (1<<0)) || (match_flag[j] & (1<<1))) && ((match_flag[j] & (1<<2)) || (match_flag[j] & (1<<3)))){
                
                    if(all){fill_ee_Histos(1, theta, j, vert[j]);}
                    ee_passedEHits.push_back(j);

                    for(Int_t k = 0; k < (Int_t) ee_passedEHits.size() && k != j; k++){

                        //Coplanarity Cut for double arm moller
                        if((TMath::Abs(TMath::Abs(phi[ee_passedEHits.at(k)]-phi[j])-180) < 10)){ //Check that the double arm mollers are coplanar
                    
                            //If this was the first pair to pass the coplanarity cut, make sure to put both hits in the histogram.
                            if(ee_passedCopHits.size()==0){
                                ee_passedCopHits.push_back(ee_passedEHits.at(k));
                                if(TMath::Abs(matchGEMz[j][1]-matchGEMz[ee_passedEHits.at(k)][1]) < 30){
                                    vert_DoubleArmMoller[ee_passedEHits.at(k)] = find_DoubleArm_ee_VertZ(ee_passedEHits.at(k), j);
                                }
                                else{
                                    vert_DoubleArmMoller[ee_passedEHits.at(k)] = 10000;
                                }
                                if(all){
                                    fill_ee_Histos(2, theta, ee_passedEHits.at(k), vert[ee_passedEHits.at(k)]);
                                    if(vert_DoubleArmMoller[ee_passedEHits.at(k)] < 10000){
                                        h_ee_zVert_DoubleArmMoller[0]->Fill(vert_DoubleArmMoller[ee_passedEHits.at(k)]);
                                    }
                                }
                            }
                            ee_passedCopHits.push_back(j);
                            if(TMath::Abs(matchGEMz[j][1]-matchGEMz[ee_passedEHits.at(k)][1])<30){
                                vert_DoubleArmMoller[j] = find_DoubleArm_ee_VertZ(j, ee_passedEHits.at(k));
                            }
                            else{
                                vert_DoubleArmMoller[j] = 100000;
                            }
                            if(all){
                                fill_ee_Histos(2, theta, j, vert[j]);
                                if(vert_DoubleArmMoller[j] < 10000){
                                    h_ee_zVert_DoubleArmMoller[0]->Fill(vert_DoubleArmMoller[j]);
                                }
                            }

                            //Cut for elasticity
                            if(TMath::Abs(cl_E[ee_passedEHits.at(k)] + cl_E[j] - EBeam - M_e) < 3*EnergyRes(EBeam)){
                        
                                //If this was the first pair to pass the elasticity cut, make sure to put both hits in the histogram.
                                if(ee_passedElastHits.size()==0){
                                    ee_passedElastHits.push_back(ee_passedEHits.at(k));
                                    fill_ee_Histos(3, theta, ee_passedEHits.at(k), vert[ee_passedEHits.at(k)]);

                                    if(vert_DoubleArmMoller[ee_passedEHits.at(k)]<10000){
                                        h_ee_zVert_DoubleArmMoller[1]->Fill(vert_DoubleArmMoller[ee_passedEHits.at(k)]);
                                    }
                                }
                        
                                fill_ee_Histos(3, theta, j, vert[j]);
                                if(vert_DoubleArmMoller[j]<10000){
                                    h_ee_zVert_DoubleArmMoller[1]->Fill(vert_DoubleArmMoller[j]);
                                }
                                ee_passedElastHits.push_back(ee_passedEHits.at(k));

                                Int_t p_ind = ee_passedEHits.at(k);

                                if(prev_x[0] > -10000 && prev_x[1] > -10000 && prev_y[0] > -10000 && prev_y[1] > -10000 && prev_z > -10000 && prev_x1[0] > -10000 && prev_x1[1] > -10000 && prev_y1[0] > -10000 && prev_y1[1] > -10000 && prev_z1 > -10000){
                                    
                                    //Get center from the first GEM plane coordinates.
                                    //Project to a common z-plane to find the center so there is no skewing.
                                    Float_t GEM0_hit0proj[2] = {projToZPlane(matchGEMx[j][0],matchGEMz[j][0], prev_z), projToZPlane(matchGEMy[j][0],matchGEMz[j][0], prev_z)};
                                    //cout<<GEM0_hit0proj[0] << " " << matchGEMx[j][0] << endl;
                                    Float_t GEM0_hit1proj[2] = {projToZPlane(matchGEMx[p_ind][0],matchGEMz[p_ind][0], prev_z), projToZPlane(matchGEMy[p_ind][0],matchGEMz[p_ind][0], prev_z)};
                                    
                                    vector<Double_t> centerG0 = findCenter(prev_x, prev_y, GEM0_hit0proj[0], GEM0_hit0proj[1], GEM0_hit1proj[0], GEM0_hit1proj[1]);
                                    h_eeCenters_GEM[0]->Fill(centerG0.at(0), centerG0.at(1));

                                    prev_x[0] = -100000;
                                    prev_y[0] = -100000;
                                    prev_x[1] = -100000;
                                    prev_y[1] = -100000;
                                    prev_z    = -100000;

                                    //Get center from the second GEM plane coordinates.
                                    //Project to a common z-plane to find the center so there is no skewing.
                                    Float_t GEM1_hit0proj[2] = {projToZPlane(matchGEMx[j][1],matchGEMz[j][1], prev_z1), projToZPlane(matchGEMy[j][1],matchGEMz[j][1], prev_z1)};
                                    Float_t GEM1_hit1proj[2] = {projToZPlane(matchGEMx[p_ind][1],matchGEMz[p_ind][1], prev_z1), projToZPlane(matchGEMy[p_ind][1],matchGEMz[p_ind][1], prev_z1)};
                                    
                                    vector<Double_t> centerG1 = findCenter(prev_x1, prev_y1, GEM1_hit0proj[0], GEM1_hit0proj[1], GEM1_hit1proj[0], GEM1_hit1proj[1]);
                                    h_eeCenters_GEM[1]->Fill(centerG1.at(0), centerG1.at(1));

                                    prev_x1[0] = -100000;
                                    prev_y1[0] = -100000;
                                    prev_x1[1] = -100000;
                                    prev_y1[1] = -100000;
                                    prev_z1    = -100000;
                                }
                                else{
                                    
                                    prev_z    = matchGEMz[j][0];
                                    prev_x[0] = projToZPlane(matchGEMx[p_ind][0], matchGEMz[p_ind][0], prev_z);
                                    prev_x[1] = matchGEMx[j][0];
                                    prev_y[0] = projToZPlane(matchGEMy[p_ind][0], matchGEMz[p_ind][0], prev_z);
                                    prev_y[1] = matchGEMy[j][0];

                                    prev_z1    = matchGEMz[j][1];
                                    prev_x1[0] = projToZPlane(matchGEMx[p_ind][1], matchGEMz[p_ind][1], prev_z1);
                                    prev_x1[1] = matchGEMx[j][1];
                                    prev_y1[0] = projToZPlane(matchGEMx[p_ind][1], matchGEMz[p_ind][1], prev_z1);
                                    prev_y1[1] = matchGEMy[j][1];

                                }
                            }
                        }
                    }
                }

                //Find e-p events
                if(TMath::Abs(cl_E[j]-ep_expE[j]) < 3.0*EnergyRes(ep_expE[j]) && ((match_flag[j] & (1<<0)) || (match_flag[j] & (1<<1))) && ((match_flag[j] & (1<<2)) || (match_flag[j] & (1<<3)))){
                    if(all){fill_ep_Histos(0, theta, j, vert[j]);}

                    //Number of blocks cut
                    if(cl_nblocks[j]>=5){
                        fill_ep_Histos(1, theta, j, vert[j]);
                    }
                }
            }
        }
    }
}

/**
 * Project the non-Z component of a vector to a new Z plane.
 *
 * @param nonZ - the non-Z component to be projected
 * @param ogZ - the original Z component
 * @param newZ - the new Z plane value to project to
 *
 * @return - the projected value of the nonZ component where it is projected
 *           in a straight line to the new Z plane provided
 */
Float_t Yields::projToZPlane(Float_t nonZ, Float_t ogZ, Float_t newZ){
	Float_t factor = newZ/ogZ;
	return factor*nonZ;
}

/**
 * Print the relevant histograms to a pdf.
 * 
 * @param pdfName - the name of the pdf to save to
 * @param begin - Notes if this is the first thing being added to this pdf.
 * @param end - Notes is this is the last thing being added to this pdf.
 */
void Yields::printPDF(TString pdfName,bool begin, bool end){
    h_epToee_ratio = (TH1F*) h_ep_Yield[EP_CUT_NUM-1]->Clone();
    h_epToee_ratio->Divide(h_ee_Yield[EE_CUT_NUM-1]);
    h_epToee_ratio->SetName("h_epToee_ratio"+typeArr[type]);
    h_epToee_ratio->SetTitle("e-p/e-e Ratio After All Cuts");
    if(gems){
        h_epToee_ratio->SetAxisRange(0,20,"Y");
    }
    else{
        h_epToee_ratio->SetAxisRange(0,1,"Y");
    }
    
    

    TCanvas *c = new TCanvas("c"+typeArr[type], "Type" + typeArr[type] + "_Yield_Canvas",1000,1000);

    if(z){
        h_ep_zVert[EP_CUT_NUM-1]->SetTitle("e-p Reconstructed Z Vertex w/ e-p over 2#circ Type"+typeArr[type]+" Cut: "+ep_cut[EP_CUT_NUM-1]);
        h_ep_zVert[EP_CUT_NUM-1]->SetAxisRange(-1000,1000, "X");

        /*h_ee_zVert[EE_CUT_NUM-1]->SetTitle("e-e Reconstructed Distance from Z Vertex Type"+typeArr[type]+" Cut: "+ee_cut[EE_CUT_NUM-1]);
        h_ee_zVert[EE_CUT_NUM-1]->SetAxisRange(5000,6000, "X");
        h_ee_zVert[EE_CUT_NUM-2]->SetTitle("e-e Reconstructed Distance from Z Vertex Type"+typeArr[type]+" Cut: "+ee_cut[EE_CUT_NUM-2]);
        h_ee_zVert[EE_CUT_NUM-2]->SetAxisRange(5000,6000, "X");*/
    }
    

    for(Int_t l = 0; l < (Int_t) runlist.size(); l++){
        for(int m = 0; m < EE_CUT_NUM; m++){
            h_ee_YieldPerLC[m]->GetXaxis()->SetBinLabel(l+1,runlist.at(l));
        }
        for(int k = 0; k < EP_CUT_NUM; k++){
            h_ep_YieldPerLC[k]->GetXaxis()->SetBinLabel(l+1,runlist.at(l));
        }
    }

    if(all){                                 //Verbose Option
        for(int i = 0; i < EE_CUT_NUM; i++){

            c->Divide(2,2);
            c->cd(1);
            gPad->SetLogz(1);
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
            
            if(gems && z){
                if(i<EE_CUT_NUM-2){
                    c->cd(1);
                    h_ee_zVert[i]->Draw("HIST");
                    c->Print(pdfName);
                    c->Clear();
                }
                else{
                    c->Divide(1,2);
                    c->cd(1);
                    h_ee_zVert[i]->Draw("HIST");
                    c->cd(2);
                    h_ee_zVert_DoubleArmMoller[i-2]->Draw("HIST");
                    c->Print(pdfName);
                    c->Clear();
                    
                }
                
            }

            if(i == 3 && !gems){
                gStyle->SetOptFit(1011);
                c->Divide(2,2);
                c->cd(1);
                h_eeCenters->Draw("COLZ");
                c->cd(3);
                TH1D* temp_x = h_eeCenters->ProjectionX();
                temp_x->SetTitle("e-e Type"+typeArr[type] +" Center X-Projection");
                temp_x->Fit("gaus","Q","",-10,10);
                temp_x->Draw("P L");
                c->cd(4);
                TH1D* temp_y = h_eeCenters->ProjectionY();
                temp_y->SetTitle("e-e Type"+typeArr[type] +" Center Y-Projection");
                temp_y->Fit("gaus", "Q", "", -10, 10);
                temp_y->Draw("P L");
                c->Print(pdfName);
                c->Clear();
                gStyle->SetOptFit(0);

            }
            else if(i == 3){ 
                gStyle->SetOptFit(1011);
                c->Divide(2,2);
                c->cd(1);
                h_eeCenters_GEM[0]->Draw("COLZ");
                c->cd(3);
                TH1D* temp_x = h_eeCenters_GEM[0]->ProjectionX();
                temp_x->SetTitle("e-e Type"+typeArr[type] +" GEM Index 0 Center X-Projection");
                temp_x->SetMarkerStyle(20);
                temp_x->Fit("gaus","Q","",-10,10);
                temp_x->Draw("P L");
                c->cd(4);
                TH1D* temp_y = h_eeCenters_GEM[0]->ProjectionY();
                temp_y->SetTitle("e-e Type"+typeArr[type] +" GEM Index 0 Center Y-Projection");
                temp_y->SetMarkerStyle(20);
                temp_y->Fit("gaus", "Q", "", -10, 10);
                temp_y->Draw("P L");
                c->Print(pdfName);
                c->Clear();

                c->Divide(2,2);
                c->cd(1);
                h_eeCenters_GEM[1]->Draw("COLZ");
                c->cd(3);
                TH1D* temp_x1 = h_eeCenters_GEM[1]->ProjectionX();
                temp_x1->SetTitle("e-e Type"+typeArr[type] +" GEM Index 1 Center X-Projection");
                temp_x1->SetMarkerStyle(20);
                temp_x1->Fit("gaus","Q","",-10,10);
                temp_x1->Draw("P L");
                c->cd(4);
                TH1D* temp_y1 = h_eeCenters_GEM[1]->ProjectionY();
                temp_y1->SetTitle("e-e Type"+typeArr[type] +" GEM Index 1 Center Y-Projection");
                temp_y1->SetMarkerStyle(20);
                temp_y1->Fit("gaus", "Q", "", -10, 10);
                temp_y1->Draw("P L");
                c->Print(pdfName);
                c->Clear();
                gStyle->SetOptFit(0);
            }
        }
        for(int j = 0; j < EP_CUT_NUM; j++){
            c->Divide(2,2);
            c->cd(1);
            gPad->SetLogz(1);
            h_ep_HC_XY[j]->Draw("COLZ");
            c->cd(2);
            h_ep_EvTheta[j]->Draw("COLZ");
            c->cd(3);
            gPad->SetLogy(1);
            h_ep_Yield[j]->Draw("E");
            c->cd(4);
            h_ep_YieldPerLC[j]->Draw("E");

            if(gems && z){
                c->Print(pdfName);
                c->Clear();

                c->cd(1);
                h_ep_zVert[j]->Draw("HIST");
            }
            
            c->Print(pdfName);
            c->Clear();
        }
        
        c->cd(1);
        //h_epToee_ratio->Draw("E");
        if(end){
            c->Print(pdfName + ")");
        }
        else{
            c->Print(pdfName);
        }
        c->Clear();

        
    } 
    else{
        c->Divide(2,2);
        c->cd(1);
        gPad->SetLogz(1);
        h_ee_HC_XY[EE_CUT_NUM-1]->Draw("COLZ");
        c->cd(2);
        h_ee_EvTheta[EE_CUT_NUM-1]->Draw("COLZ");
        c->cd(3);
        gPad->SetLogy(1);
        h_ee_Yield[EE_CUT_NUM-1]->Draw("E");
        c->cd(4);
        h_ee_YieldPerLC[EE_CUT_NUM-1]->Draw("E");
        if(begin){
            c->Print(pdfName+"(");
        }
        else{
            c->Print(pdfName);
        }
        c->Clear();

        if(gems && z){
            c->cd(1);
            h_ee_zVert[EE_CUT_NUM-1]->Draw("HIST");
            c->Print(pdfName);
            c->Clear();
        }

        c->Divide(2,2);
        c->cd(1);
        h_eeCenters->Draw("COLZ");
        c->cd(3);
        TH1D* temp_x2 = h_eeCenters->ProjectionX();
        temp_x2->SetTitle("e-e Type"+typeArr[type] +" Center X-Projection");
        //temp_x2->Fit("gaus","Q","",-60,60);
        temp_x2->Draw("HIST");
        c->cd(4);
        TH1D* temp_y2 = h_eeCenters->ProjectionY();
        temp_y2->SetTitle("e-e Type"+typeArr[type] +" Center Y-Projection");
        //temp_y2->Fit("gaus", "Q", "", -60, 60);
        temp_y2->Draw("HIST");
        c->Print(pdfName);
        c->Clear();

        c->Divide(2,2);
        c->cd(1);
        gPad->SetLogz(1);
        h_ep_HC_XY[EP_CUT_NUM-1]->Draw("COLZ");
        c->cd(2);
        h_ep_EvTheta[EP_CUT_NUM-1]->Draw("COLZ");
        c->cd(3);
        gPad->SetLogy(1);
        h_ep_Yield[EP_CUT_NUM-1]->Draw("E");
        c->cd(4);
        h_ep_YieldPerLC[EP_CUT_NUM-1]->Draw("E");
        
        if(gems && z){
            c->Print(pdfName);
            c->Clear();

            c->cd(1);
            h_ep_zVert[EP_CUT_NUM-1]->Draw("HIST");
        }
        c->Print(pdfName);

        c->cd(1);
        //h_epToee_ratio->Draw("E");
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

/**
 * Helper method that finds the vertex Z for the current Mott event or any event assuming it comes from beamline.
 *
 * @param j - index of the hit being investigated.
 *
 * @return - the reconstructed vertex z of all particles assuming origin of path along the beamline.
 */
Float_t Yields::find_VertZ_beamline(Int_t j){
    //Float_t beamline[3] = {0, 0, 1};

    Float_t eventVec[3] = {matchGEMx[j][0]-matchGEMx[j][1], matchGEMy[j][0]-matchGEMy[j][1], matchGEMz[j][0]-matchGEMz[j][1]};

    Float_t mag = TMath::Sqrt(eventVec[0]*eventVec[0] + eventVec[1]*eventVec[1] + eventVec[2]*eventVec[2]);

    Float_t u[3] = {eventVec[0]/mag, eventVec[1]/mag, eventVec[2]/mag};

    return matchGEMz[j][1] - (matchGEMx[j][1]/u[0])*u[2];
}

/**
 * Helper Method that finds the vertez Z for the current Double arm e-e event.
 *
 * @param j - the index of the hit being investigated
 * @param k - the index of the identified Moller partner
 *
 * @return - the reconstructed z vertex position of the double arm Moller pairs.
 */
Float_t Yields::find_DoubleArm_ee_VertZ(Int_t j, Int_t k){

    Float_t x1 = matchGEMx[j][1];
    Float_t y1 = matchGEMy[j][1];
    Float_t r_1 = TMath::Sqrt(x1*x1 + y1*y1);

    Float_t x2 = matchGEMx[k][1];
    Float_t y2 = matchGEMy[k][1];
    Float_t r_2 = TMath::Sqrt(x2*x2 + y2*y2);

    return TMath::Sqrt(((M_e+EBeam)*r_1*r_2)/(2*M_e));
}