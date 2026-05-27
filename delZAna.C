#define delZAna_cxx

#include "delZAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void delZAna::Begin()
{
   //Initialize global variables here.
  GEV = 1000.;
  MEV2GEV = .001;
  nEvtTotal = 0;
  njet=0;
  //Create the histogram file
  _HstFile = new TFile(_HstFileName,"recreate");
  //Call the function to book the histograms we declared in Hists.
  BookHistograms();
}

void delZAna::End()
{
  //Write histograms and close histogram file
  _HstFile->Write();
  _HstFile->Close();
  //Output to screen.
  cout<<"Total events = "<<nEvtTotal<<endl;
  //Open the text output file
  ofstream fout(_SumFileName);
  //Put text output in the summary file.
  fout<<"Total events = "<<nEvtTotal<<endl;
}

void delZAna::Loop()
{
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
  
  if(fChain==0) return;
  Long64_t numberOfEntries = treeReader->GetEntries();
  //This prints total entries before processing begins
  cout<<"Total entries = "<<numberOfEntries<<endl;
  nEvtTotal=0;

  //This is the main loop over the number of events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){

    treeReader->ReadEntry(entry);
    nEvtTotal++;
    //Output processing information to screen based on verbosity level.
    if(_verbosity>5000)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;          
    if(_verbosity>1000 && nEvtTotal%5000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
    else if(_verbosity>0 && nEvtTotal%50000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
    
    
    //Your CODE starts here
    //###############################
    //Construction of the arrays:
    // 1.jetgen
    // 2.goodEle
    //#################################
    



    //MissingET *met = (MissingET*)brMissingET->At(0); //not using this.
    TLorentzVector mymet(0,0,0,0); // build our own met
    
    //goodMu.clear(); //clear array from previous event
    goodEle.clear();
    jetgen.clear();
    
    int nele=0; int nmuons=0; int nparts=0;
    
    //GENJET
    //##################################################################
    // These are genJets. Here we use the Delphes Collection.
    for(int ij=0; ij< brGenJet->GetEntries(); ij++){
      Jet *jt = (Jet*)brGenJet->At(ij);
      Lepton temp; temp.v.SetPtEtaPhiM(jt->PT,jt->Eta,jt->Phi,jt->Mass); temp.ind = ij;
      jetgen.push_back(temp);
    }
    Sort(2);//sorting the genjets
    
    
    
    TLorentzVector q1,q2;
    q1.SetPtEtaPhiM(0,0,0,0); q2.SetPtEtaPhiM(0,0,0,0);
    
    //Now we loop over truth particles and keep what we need.
    for(int no=0; no<brParticle->GetEntries(); no++){
      GenParticle *par = (GenParticle*)brParticle->At(no);
      
      /*
      // Example of getting daughters of a particle.
      //if particle is Z, save its daughters.
      if(par->PID==23){
	GenParticle *dau1 = (GenParticle*)brParticle->At(par->D1);
	GenParticle *dau2 = (GenParticle*)brParticle->At(par->D2);
	if(dau1->PID==23 || dau2->PID==23) continue;
	q1.SetPtEtaPhiM(dau1->PT,dau1->Eta,dau1->Phi,dau1->Mass);
	q2.SetPtEtaPhiM(dau2->PT,dau2->Eta,dau2->Phi,dau2->Mass);	
	if(q1.Pt()>0 && q2.Pt()>0){	
	  float massqq = (q1+q2).M();
	  h.massjj[0]->Fill(massqq);
	}
      }
      */
      
      // Building my missing et from neutrinos
      if(par->Status==1){ //its a stable particle
	if(abs(par->PID)==12 || abs(par->PID)==14 || abs(par->PID)==16){
	  TLorentzVector neut; neut.SetPtEtaPhiM(par->PT,par->Eta,par->Phi,0);
	  mymet += neut;
	}
      }

      //Electrons
      if(par->Status==1 && abs(par->PID)==11){//its a stable electron
	nele++;
	Lepton temp; temp.v.SetPtEtaPhiM(par->PT,par->Eta,par->Phi,par->Mass);
	temp.id = par->PID; temp.ind=no;
	float isol = calc_iso(0.4,temp.v,no);
	if(temp.v.Pt()>20 && fabs(temp.v.Eta())<2.4 && isol/temp.v.Pt()<0.15){//check if this candidate passes selections.
	  goodEle.push_back(temp); //store it in the goodele array.
	}
      }
      
      
      
     
      Sort(1); //sort
      
    }//end of truth particles loop


    //#######################################################
    //Kinematic Plots
    //#######################################################

    //Electrons
    /*
    h.ele[0]->Fill((int)goodEle.size());
       
    for(int i=0; i<(int)goodEle.size(); i++){
      h.ele[1]->Fill(goodEle.at(i).v.E());
      h.ele[2]->Fill(goodEle.at(i).v.Eta());
      h.ele[3]->Fill(goodEle.at(i).v.Phi());
    }
    */
    
    //jets

    /*
    h.jet[0]->Fill((int)jetgen.size());
    
    for(int i=0; i<(int)jetgen.size(); i++){
      h.jet[1]->Fill(jetgen.at(i).v.E());
      h.jet[2]->Fill(jetgen.at(i).v.Theta());
      h.jet[3]->Fill(jetgen.at(i).v.Phi());
    }
    */
    
    if((int)jetgen.size()>1){//atlest 2 jets
      h.jet[4]->Fill(jetgen.at(0).v.E());
      h.jet[5]->Fill(jetgen.at(0).v.Theta());
      h.jet[6]->Fill(jetgen.at(0).v.Phi());
      
      h.jet[7]->Fill(jetgen.at(1).v.E());
      h.jet[8]->Fill(jetgen.at(1).v.Theta());
      h.jet[9]->Fill(jetgen.at(1).v.Phi());

      //invar mass
      float  invar_mass= (jetgen.at(0).v + jetgen.at(1).v).M() ;
      h.jet[10]->Fill(invar_mass);
    }

    //lund jet var
    
   
    //MET
    float missinget = mymet.Pt();
    float metphi = mymet.Phi();
    
    // Call the Lund Plane function for the leading jet (index 0)
    AnalyzeLundPlane(0);
    
    
  }//end loop over events.
}//end of delzana


//for fastjets

void delZAna::AnalyzeLundPlane(int iJet)
{
  //the function at least started
  if (nEvtTotal <= 3) {
    cout << "\n--- Starting Lund Analysis for Event " << nEvtTotal << ", Target Jet " << iJet << " ---" << endl;
  }
  
  // Safety Check 1
  if (jetgen.size() <= 1) {
    if (nEvtTotal <= 3) cout << "  -> ABORTED: jetgen.size() is " << jetgen.size() << " (Needs to be > 1)" << endl;
    return;
  }
  
  //Safety Check 2
  if (brParticle->GetEntries() == 0) {
    if (nEvtTotal <= 3) cout << "  -> ABORTED: brParticle is empty!" << endl;
    return;
  }
  
  // if (jetgen.size() <= 1 || brParticle->GetEntries() == 0) return;
  std::vector<fastjet::PseudoJet> subjets;
  int status1_count = 0;

  //#############################################################
  //CONSTRUCTING MY JET
  //#############################################################

  
  //Collect Constituents
  for (int no = 0; no < brParticle->GetEntries(); no++) {
    GenParticle *par = (GenParticle*)brParticle->At(no);
    if (par->Status == 1) {
      fastjet::PseudoJet constituent(par->Px, par->Py, par->Pz, par->E);
      subjets.push_back(constituent);
      status1_count++; // Count how many stable (=1) particles we actually find
    }
  }
  
  if (nEvtTotal <= 3) {
    cout << "  -> Scanned " << brParticle->GetEntries() << " total particles." << endl;
    cout << "  -> Found " << status1_count << " Status=1 particles." << endl;
  }
  
  //Safety Check 3
  if (subjets.empty()) {
    if (nEvtTotal <= 3) cout << "  -> ABORTED: subjets array is empty! No Status 1 particles found." << endl;
    return;
  }

  //PRIMARY JET CLUSTERING (Generalized e+e- k_t algorithm)
  //since its a e+e- collider so the radius of the jers can be infinite
  //fastjet::JetDefinition jet_def(fastjet::ee_genkt_algorithm, fastjet::JetDefinition::max_allowable_R, -1);
  fastjet::JetDefinition jet_def(fastjet::ee_genkt_algorithm, 1.0, -1); //the deltaR for jets is 1.0
  fastjet::ClusterSequence cluster_seq(subjets, jet_def);

  
  //the clustering also looked at the entire event, kept merging particles together because they were all  within an "infinite" distance of each other
  //created one single massive super-jet that ate the entire Z boson ans peaks at 91
  
  std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(cluster_seq.inclusive_jets());
  if (inclusive_jets.empty()) return;
  fastjet::PseudoJet leadJet = inclusive_jets[0];
  
  
  //EXCLUSIVE CLUSTERING
  //explicitly forces the FastJet math engine to stop clustering the moment it has narrowed the event down to exactly 2 jets- exclusive_jets[0] and exclusive_jets[1]
  //in a perfect detector they both will have the sane energy - 45.5 GeV
  //but that's not true in the real world so we sort them by Energy
    //std::vector<fastjet::PseudoJet> exclusive_jets = sorted_by_E(cluster_seq.exclusive_jets(2));
  
  // Safety check to ensure we actually got 2 jets
  if (inclusive_jets.size() < 2) return;
  
  // Grab the leading jet
  fastjet::PseudoJet recoJet = inclusive_jets[0];

  
  h.recojet[1]->Fill(recoJet.E());
  h.recojet[2]->Fill(recoJet.theta());               
  //  h.jet[15]->Fill(recoJet.phi()); //this gives me a distribution from 0 to 2pi
  h.recojet[3]->Fill(recoJet.phi_std()); //this gives me a distribution from - pi to pi
  h.recojet[4]->Fill(recoJet.m());                 
  h.recojet[5]->Fill(recoJet.constituents().size());

  //#############################################################
  //RECLUSTERING VIA EE-CAMBRIDGE PLUGIN (The Lund Declustering Tree)
  //#############################################################

  //Create the array to hold our splits so .size() works later
  std::vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet>> lund_subjets;
  // memory allocation via pointers 
  fastjet::JetDefinition::Plugin* ee_plugin = new fastjet::EECambridgePlugin(1.0);
  fastjet::JetDefinition jd(ee_plugin);
  fastjet::Recluster rc(jd);
    
  fastjet::PseudoJet j = rc.result(recoJet);
  fastjet::PseudoJet jj, j1, j2;
    
  jj = j;
  // Follow the tracking parent splits all the way down the decay tree
  //We have a jet jj and we want to pull it apart into the two pieces j1 and j2 that originally merged to create it
  while (jj.has_parents(j1, j2)) {
    // Enforce a strict hierarchy where
    //j1 is always the harder core component
    //j2 is the branch that was emitted
    
    if (j1.modp2() < j2.modp2()) {
      fastjet::PseudoJet jTemp = j1;
      j1 = j2;
      j2 = jTemp;
    }
    //save the split to the array
    lund_subjets.push_back(make_pair(j1, j2));

        
    //#############################################################
    //getting the variab;es
    //#############################################################

    // Call the function and get the pair back
    std::pair<float, float> lund_vars = GetKinematics(j1, j2);
    
    // Unpack the results
    float log_1_over_theta = lund_vars.first;
    float log_kt           = lund_vars.second;
  
    h.recojet[6]->Fill(log_kt);           
    h.recojet[7]->Fill(log_1_over_theta); 


 
    // LUND PLANE 2D HISTOGRAM
    h.lundPlane[0]->Fill(log_1_over_theta, log_kt);
    
    //store the harder particle j1, throw away the emitted particle j2
    jj = j1;
  }// while has parenys loop ends
  //The loop stops when you hit a particle (jj = j1) that cannot be un-merged any further

  h.jet[8]->Fill(lund_subjets.size()); // Fill total number of splits
  //everytime we unmerge a jet into 2 peices - thats 1 split 

  //  Free the plugin resource explicitly to completely prevent memory leaks
  delete ee_plugin;

 
}



//func to get the variables
std::pair<float, float> delZAna::GetKinematics(fastjet::PseudoJet j1, fastjet::PseudoJet j2) 
{
  //Calculate Splitting Angle Theta Matrix
  TVector3 J1(j1.px(), j1.py(), j1.pz());
  TVector3 J2(j2.px(), j2.py(), j2.pz());
  double omc = 1.0 - cos(J1.Angle(J2));

  double theta;
  if (omc > 0.058) {
    theta = acos(1.0 - omc); // Large angle case use cos inverse formula
  } else {
    theta = sqrt(2.0 * omc); // small angle case use taylor expansion
  }
       
  float log_1_over_theta = std::log(1.0 / theta);
  
  
  // Calculate Relative Transverse Momentum kt
  double sin_theta = sin(theta);
  double kt = j2.modp() * sin_theta;
  float log_kt = std::log(kt);
  
  return std::make_pair(log_1_over_theta, log_kt);
}
    
  




  
void delZAna::Sort(int opt)
{
  //Sort selected objects by pT (always descending).
  //option 1 sorts the gooMu array
  if(opt==1){
    for(int i=0; i<(int)goodMu.size()-1; i++){
      for(int j=i+1; j<(int)goodMu.size(); j++){
	if( goodMu[i].v.Pt() < goodMu[j].v.Pt() ) swap(goodMu.at(i),goodMu.at(j)); }}
    for(int i=0; i<(int)goodTau.size()-1; i++){
      for(int j=i+1; j<(int)goodTau.size(); j++){
	if( goodTau[i].v.Pt() < goodTau[j].v.Pt() ) swap(goodTau.at(i),goodTau.at(j)); }}
  }
  if(opt==2){
    for(int i=0; i<(int)jetgen.size()-1; i++){
      for(int j=i+1; j<(int)jetgen.size(); j++){
	if( jetgen[i].v.Pt() < jetgen[j].v.Pt() ) swap(jetgen.at(i),jetgen.at(j)); }}
  }
    //Once you have other arrays (goodEle, goodPho), write code here to sort them.
}
float delZAna::delta_phi(float phi1, float phi2)
{
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}


Double_t delZAna::delR(TLorentzVector v1,  TLorentzVector v2)
{ 
  double del_eta = abs(v1.Eta() - v2.Eta());
  return sqrt( pow(del_eta,2) + pow(delta_phi(v1.Phi(),v2.Phi()),2) );
}
float delZAna::calc_iso(float dR, TLorentzVector v1, int skip)
{
  // Sum the pT of all stable particles with deltaR<dR of particle v1
  // Skip the particle itself (skip) and skip neutrinos.
  float isol = 0;
  for(int no=0; no<brParticle->GetEntries(); no++){
    GenParticle *par = (GenParticle*)brParticle->At(no);
    if(no==skip) continue;
    if(par->Status==1){
      if(abs(par->PID)==12 || abs(par->PID)==14 || abs(par->PID)==16) continue;
      TLorentzVector v2; v2.SetPtEtaPhiM(par->PT,par->Eta,par->Phi,par->Mass);
      if(v1.DeltaR(v2)<dR) isol+=v2.Pt();
    }
  }
  return isol;

}
void delZAna::TauBuilder()
{
  //First we get the list of taus and keep them in tauind
  // we only keep those taus which dont decay to leptons
  // i.e. we want only those taus which decay hadronically
  vector<int> tauind;
  for(int no=0; no<brParticle->GetEntries(); no++){
    GenParticle *par = (GenParticle*)brParticle->At(no);
    if( abs(par->PID)==15 ){
      bool gottau=true;
      for(int da=par->D1;da<par->D2+1;da++){
	GenParticle *parda = (GenParticle*)brParticle->At(da);
	if(abs(parda->PID)==15 || abs(parda->PID)==11 || abs(parda->PID)==13){
	  gottau=false; break;
	}
      }
      if(gottau) tauind.push_back(no);
    }
  }
  //cout<<(int)tauind.size()<<endl;

  // Now for the taus which decay hadronically, we
  // loop over the visible daughters (dav) and build the visible 4-vector (vistau)
  // We only keep those taus with pT>20 GeV in goodTau.
  for(int i=0; i<(int)tauind.size(); i++){
    GenParticle *tau = (GenParticle*)brParticle->At(tauind.at(i));
    TLorentzVector vistau(0,0,0,0);
    for(int da=tau->D1;da<tau->D2+1;da++){
      GenParticle *dt = (GenParticle*)brParticle->At(da);
      if(abs(dt->PID)>10 && abs(dt->PID)<17) continue;
      //cout<<"NonLeptonic Daughter"<<endl;
      TLorentzVector dav; dav.SetPtEtaPhiM(dt->PT,dt->Eta,dt->Phi,dt->Mass);
      vistau += dav;
    }
    Lepton temp; temp.v = vistau; temp.id=tau->PID; temp.ind=tauind.at(i);
    if(temp.v.Pt()>20 && fabs(temp.v.Eta())<2.4)goodTau.push_back(temp);
      
  }
      
  
}

void delZAna::BookHistograms()
{
  // Booking syntax for histogram pointers (that we declared in struct Hists in header)
  // obj = new CLASS("localname","Histogram title",NumberOfBins,LowerEdge,HigherEdge)

  //ele
  /*
  h.ele[0]         =   new TH1F("ngoodele","Number of Electron",5,0,5);
  h.ele[1]         =   new TH1F("ele_E","Electron Energy",1000,0,2000); 
  h.ele[2]         =   new TH1F("ele_Eta","Electron Eta",64,-3.2,3.2); 
  h.ele[3]         =   new TH1F("ele_Phi","Electron Phi",64,-3.2,3.2); 
  */
  
  //jet
  /*
  h.jet[0]         =   new TH1F("genjet_n",          "No of GenJets",        6,0,6);
  h.jet[1]         =   new TH1F("genjet_E",          "GenJet Energy",        50,0,100); 
  h.jet[2]         =   new TH1F("genjet_Theta",      "GenJet Theta",         64,0,3.2); 
  h.jet[3]         =   new TH1F("genjet_Phi",        "GenJet Phi",           64,-3.2,3.2); 
  */
  
  h.jet[4]         =   new TH1F("genjet0_E",         "Leading GenJet Energy",          50,0,100); 
  h.jet[5]         =   new TH1F("genjet0_theta",     "Leading GenJet theta",           64,0,3.2); 
  h.jet[6]         =   new TH1F("genjet0_Phi",       "Leading GenJet Phi",             64,-3.2,3.2); 


  h.jet[7]         =   new TH1F("genjet1_E",          "Sub Leading GenJet Energy",     50,0,100); 
  h.jet[8]         =   new TH1F("genjet1_theta",      "Sub Leading GenJet theta",      64,0,3.2); 
  h.jet[9]         =   new TH1F("genjet1_Phi",        "Sub Leading GenJet Phi",        64,-3.2,3.2); 

  h.jet[10]        =   new TH1F("genjet_invar_mass",  "Invariant Mass of Leading and subleading Jets",100,0,100); 

  //LundJetPlane Var
  h.recojet[1] = new TH1F("recojet0_energy",             "Reco Jet Energy; Energy [GeV]; Events",                                      50, 0, 100);
  h.recojet[2] = new TH1F("recojet0_theta",              "Reco Jet #theta; #theta [rad]; Events",                                            64, 0, 3.2);
  h.recojet[3] = new TH1F("recojet0_phi",                "Reco Jet #phi; #phi [rad]; Events",                                          64, -3.2, 3.2);
  h.recojet[4] = new TH1F("recojet0_mass",               "Reco Jet Mass; Mass [GeV]; Events",                                          50, 0, 50);
  h.recojet[5] = new TH1F("recojet0_constituents",       "Reco Jet Constituents; N_{particles}; Events",                               50, 0, 50);


  h.recojet[6] = new TH1F("lund_kt",                     "Lund Split Momentum; ln(k_{t}) [GeV]; Number of Splits",                     50, -8, 6);
  h.recojet[7] = new TH1F("lund_theta",                  "Lund Split Angle #theta; ln(1.0 / #theta); Number of Splits",                50, 0, 7);
  h.recojet[8] = new TH1F("n_splits",                    "Number of Lund Splits per Jet; Splits; Events",                              30, 0, 30);   


  //Lund Plane
  
  // X-axis: ln(1/theta) from 0 to 7
  // Y-axis: ln(kt) from -8 to 6
  h.lundPlane[0] = new TH2F("lund_plane", "Primary Lund Jet Plane; ln(1.0 / #theta); ln(k_{t}) [GeV]", 50, 0, 7, 50, -8, 6);
  
  // Note the calling of Sumw2() for each histogram we declare.
  // This is needed so that if scale the histograms later, the
  // errors on the points are still correct after scaling.
  // This is relevant.. so ask someone if you don't understand.
  
}















/*
void delZAna::AnalyzeLundPlane(int iJet)
{
  // Protect against empty arrays or single-constituent clusters
  if (jetgen.size() <= 1 || brParticle->GetEntries() == 0) return;
  
  std::vector<fastjet::PseudoJet> subjets;
  //#############################################################
  //CONSTRUCTING MY JET------------------------------------------
  //#############################################################
  
  //COLLECT CONSTITUENTS
  // Using the stable truth particles loop to find tracking components
  for (int no = 0; no < brParticle->GetEntries(); no++) {
    GenParticle *par = (GenParticle*)brParticle->At(no);
    if (par->Status == 1) {
      // Reconstruct the 4-momentum coordinates natively
      fastjet::PseudoJet constituent(par->Px, par->Py, par->Pz, par->E);
      subjets.push_back(constituent);
    }
  }
  
  if (subjets.empty()) return;

  //PRIMARY JET CLUSTERING (Generalized e+e- k_t algorithm)
  fastjet::JetDefinition jet_def(fastjet::ee_genkt_algorithm, fastjet::JetDefinition::max_allowable_R, -1);
  fastjet::ClusterSequence cluster_seq(subjets, jet_def);
    
  std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(cluster_seq.inclusive_jets());
  if (inclusive_jets.empty()) return;
    
  fastjet::PseudoJet leadJet = inclusive_jets[0];


  //checking
  if (nEvtTotal <= 3) {
    cout << "\n==========================================" << endl;
    cout << "EVENT " << nEvtTotal << " | MACRO JET BUILT!" << endl;
    cout << "Target Jet Index: " << iJet << endl;
    cout << "  -> Jet pT:   " << leadJet.pt() << " GeV" << endl;
    cout << "  -> Jet Eta:  " << leadJet.eta() << endl;
    cout << "  -> Jet Phi:  " << leadJet.phi() << endl;
    cout << "  -> Jet Mass: " << leadJet.m() << " GeV" << endl;
    cout << "  -> Number of Particles Inside: " << leadJet.constituents().size() << endl;
    cout << "==========================================\n" << endl;
  }



}
*/
