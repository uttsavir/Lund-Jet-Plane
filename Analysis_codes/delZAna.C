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
    // 3.jetreco
    // 4.bjetreco
    //#################################

    //MissingET *met = (MissingET*)brMissingET->At(0); //not using this.
    TLorentzVector mymet(0,0,0,0); // build our own met
    
    goodEle.clear();
    jetgen.clear();
    jetreco.clear();
    bjetreco.clear();

    
    int nele=0; int nmuons=0; int nparts=0;


    //##################################################################
    //GENJET (truth level) - using delphes collection
    //##################################################################
    /* 
    //Delphes uses anti kt (R=0.4) to tag jet.
    for(int ij=0; ij< brGenJet->GetEntries(); ij++){
      Jet *jt = (Jet*)brGenJet->At(ij);
      Lepton temp;
      temp.v.SetPtEtaPhiM(jt->PT,jt->Eta,jt->Phi,jt->Mass); temp.ind = ij;
      jetgen.push_back(temp);
    }
    Sort(2);//sorting the genjets
    */

    //##################################################################
    // GENJET (truth level) - Reclustered using ee_genkt
    //##################################################################

    /*
    for(int ij = 0; ij < brGenJet->GetEntries(); ij++){
      Jet *jt = (Jet*)brGenJet->At(ij);
      
      // 1. Gather all constituents inside this specific Delphes GenJet
      std::vector<fastjet::PseudoJet> jet_constituents;
      for (int c = 0; c < jt->Constituents.GetEntriesFast(); c++) {
        TObject *obj = jt->Constituents.At(c);
        if (!obj) continue;
        GenParticle *par = (GenParticle*)obj;
        
        // Convert Delphes GenParticles to FastJet PseudoJets
        jet_constituents.push_back(fastjet::PseudoJet(par->Px, par->Py, par->Pz, par->E));
      }
      
      // Safety check: skip empty jets
      if (jet_constituents.size() < 1) continue;


      // 2. Recluster those exact constituents using ee_genkt
      fastjet::JetDefinition reclust_jd(fastjet::ee_genkt_algorithm, 1.0, -1);
      fastjet::ClusterSequence reclust_cs(jet_constituents, reclust_jd);
      
      // Get the resulting reclustered jet(s)
      std::vector<fastjet::PseudoJet> reclustered_jets = fastjet::sorted_by_E(reclust_cs.inclusive_jets());
      
      if (reclustered_jets.size() > 0) {
        // Grab the leading reclustered jet
        fastjet::PseudoJet myNewJet = reclustered_jets[0];
        
        // 3. Save the RECLUSTERED kinematics to your array
        Lepton temp;
        // Notice we are asking FastJet for the kinematics now, not Delphes!
        temp.v.SetPtEtaPhiM(myNewJet.pt(), myNewJet.eta(), myNewJet.phi(), myNewJet.m()); 
        temp.ind = ij; // Save the original Delphes index just in case
        
        jetgen.push_back(temp);
	AnalyzeLundPlane(myNewJet, 0); //DRAW THE TRUTH-LEVEL LUND PLANE HERE
      }
    }
    
    //sort(2); // Sort the newly reclustered genjets by pT
    */


    //##################################################################
    // GENJET (truth level) - Using direct Delphes Kinematics
    //not clustering it using gen kt
    //##################################################################
    
    for(int ij = 0; ij < brGenJet->GetEntries(); ij++){
      Jet *jt = (Jet*)brGenJet->At(ij);

      Lepton temp;
      temp.v.SetPtEtaPhiM(jt->PT, jt->Eta, jt->Phi, jt->Mass); 
      temp.ind = ij; 
      jetgen.push_back(temp);//Saving the truth level jets in my jetgen array

      
      //Prepare the jet for the Lund Plane function
      // we still need to bundle its constituents together so the Lund Plane 
      // function can see the internal particles to build its C/A tree.
      
      std::vector<fastjet::PseudoJet> jet_constituents;
      for (int c = 0; c < jt->Constituents.GetEntriesFast(); c++) {
        TObject *obj = jt->Constituents.At(c);
        if (!obj) continue;
        GenParticle *par = (GenParticle*)obj;
        
        jet_constituents.push_back(fastjet::PseudoJet(par->Px, par->Py, par->Pz, par->E));
      }
      
      // If the jet has internal structure (at least 2 particles), map it!
      if (jet_constituents.size() >= 2) {
        
        // 'join' packages the particles into one PseudoJet WITHOUT reclustering
        fastjet::PseudoJet trueJet = fastjet::join(jet_constituents);
        
        // Draw the Truth-Level Lund Plane
        AnalyzeLundPlane(trueJet, 0); 
      }
    }


    
    //##################################################################
    // RecoJets (detector level)
    //##################################################################
    
    for(int ij = 0; ij < brJet->GetEntries(); ij++){
      Jet *jt = (Jet*)brJet->At(ij);
      
      // 1. Save the original Delphes kinematics directly to your array!
      Lepton temp;
      temp.v.SetPtEtaPhiM(jt->PT, jt->Eta, jt->Phi, jt->Mass); 
      temp.ind = ij; 
      jetreco.push_back(temp);

      /*
      //making the bjetreco array
      h.btag_score[0]->Fill(Jet_btagDeepFlavB[ij]);
      //this discriminator goes from 0 to 1 
      //where 1 means its def a bjet
      */
      h.btag_flag->Fill(jt->BTag);
      
      bool is_bjet = false;
      if(is_bjet = (jt->BTag == 1))
	//BTah ==1 -> its a bjet
	//BTag == 0 -> its not a bjet
	is_bjet = true;

      if(is_bjet){
	h.btag_flag_aftercut->Fill(jt->BTag);
	bjetreco.push_back(temp);
      }

      //Prepare the jet for the Lund Plane function
      std::vector<fastjet::PseudoJet> jet_constituents;
      
      for (int c = 0; c < jt->Constituents.GetEntriesFast(); c++) {
        TObject *obj = jt->Constituents.At(c);
        if (!obj) continue;
	GenParticle *par = (GenParticle*)obj;
	
	jet_constituents.push_back(fastjet::PseudoJet(par->Px, par->Py, par->Pz, par->E));
      }
      
      // If the jet has internal structure (at least 2 particles), map it!
      if (jet_constituents.size() >= 2) {
        
        // 'join' packages the particles into one PseudoJet WITHOUT reclustering
        fastjet::PseudoJet myJet = fastjet::join(jet_constituents);
        
        // Draw the Reco-Level Lund Plane (histIndex = 1)
        AnalyzeLundPlane(myJet, 1);
	
	// Draw the b jet Reco-Level Lund Plane (histIndex = )
	if(is_bjet)
        AnalyzeLundPlane(myJet, 4); 
      }
    }












    

    /*

    // ===============================================================
    // 2. Prepare the jet for the Lund Plane function
    // ===============================================================

    for(int ij = 0; ij < brJet->GetEntries(); ij++){
      Jet *jn = (Jet*)brJet->At(ij);//pointing to the jet at the above index

      Lepton temp;
      temp.v.SetPtEtaPhiM(jn->PT, jn->Eta, jn->Phi, jn->Mass); 
      temp.ind = ij; 
      jetreco.push_back(temp);
      
      std::vector<fastjet::PseudoJet> jet_constituents;
      
      // Loop over every single internal constituent inside this specific detector jet
      for (int c = 0; c < jn->Constituents.GetEntriesFast(); c++) {
	TObject *obj = jn->Constituents.At(c);
	if (!obj) continue;
	
	TLorentzVector v;
          
	// SMART CASTING: Detector hits are Tracks or Towers, NOT GenParticles!
	if (obj->IsA() == Track::Class()) {
	  Track *trk = (Track*)obj;
	  v.SetPtEtaPhiM(trk->PT, trk->Eta, trk->Phi, trk->Mass);
	}
	
	else if (obj->IsA() == Tower::Class()) {
	  Tower *tow = (Tower*)obj;
	  v.SetPtEtaPhiE(tow->ET, tow->Eta, tow->Phi, tow->E); 
	  }
	else {
	  
	  continue; // Skip unknown objects to prevent crashing
	}
          
	// Translate the correctly read physics properties into a FastJet PseudoJet
	jet_constituents.push_back(fastjet::PseudoJet(v.Px(), v.Py(), v.Pz(), v.E()));
      }
	
        
      // If the jet has internal structure (at least 2 particles), map it!
      if (jet_constituents.size() >= 2) {
	  
	// 'join' packages the particles into one PseudoJet WITHOUT reclustering
	fastjet::PseudoJet myJet = fastjet::join(jet_constituents);
	  
	// Draw the Reco-Level Lund Plane (histIndex = 1)
	AnalyzeLundPlane(myJet, 1); 
      }
    }
    
    */




    







        
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
      
    // SORT ARRAYS BEFORE PLOTTING KINEMATICS!
    Sort(1); // Sorts Leptons
    Sort(2); // Sorts GenJets and RecoJets by pT  

    }//end of truth particles loop
    
    
    //#######################################################
    //Kinematic Plots
    //#######################################################
    
    //genjets
    
    h.jet[0]->Fill((int)jetgen.size());
    
    for(int i=0; i<(int)jetgen.size(); i++){
      h.jet[1]->Fill(jetgen.at(i).v.E());
      h.jet[2]->Fill(jetgen.at(i).v.Theta());
      h.jet[3]->Fill(jetgen.at(i).v.Phi());
    }
    
    
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

    //recojets
    h.recojet[0]->Fill((int)jetreco.size());
    
    for(int i=0; i<(int)jetreco.size(); i++){
      h.recojet[1]->Fill(jetreco.at(i).v.E());
      h.recojet[2]->Fill(jetreco.at(i).v.Theta());
      h.recojet[3]->Fill(jetreco.at(i).v.Phi());
    }
    
    if((int)jetreco.size() > 1){ //at least 2 jets
      h.recojet[4]->Fill(jetreco.at(0).v.E());
      h.recojet[5]->Fill(jetreco.at(0).v.Theta());
      h.recojet[6]->Fill(jetreco.at(0).v.Phi());
      
      h.recojet[7]->Fill(jetreco.at(1).v.E());
      h.recojet[8]->Fill(jetreco.at(1).v.Theta());
      h.recojet[9]->Fill(jetreco.at(1).v.Phi());

      float invar_mass = (jetreco.at(0).v + jetreco.at(1).v).M();
      h.recojet[10]->Fill(invar_mass);
    }

    
    //recobjets
    h.recobjet[0]->Fill((int)bjetreco.size());
    
    for(int i=0; i<(int)bjetreco.size(); i++){
      h.recobjet[1]->Fill(bjetreco.at(i).v.E());
      h.recobjet[2]->Fill(bjetreco.at(i).v.Theta());
      h.recobjet[3]->Fill(bjetreco.at(i).v.Phi());
    }
    
    if((int)bjetreco.size() > 1){ //at least 2 jets
      h.recobjet[4]->Fill(bjetreco.at(0).v.E());
      h.recobjet[5]->Fill(bjetreco.at(0).v.Theta());
      h.recobjet[6]->Fill(bjetreco.at(0).v.Phi());
      
      h.recobjet[7]->Fill(bjetreco.at(1).v.E());
      h.recobjet[8]->Fill(bjetreco.at(1).v.Theta());
      h.recobjet[9]->Fill(bjetreco.at(1).v.Phi());

      float invar_mass = (bjetreco.at(0).v + bjetreco.at(1).v).M();
      h.recobjet[10]->Fill(invar_mass);
    }
    
    
    





    //---------------------------------------------------------------------------
    //dR Matching b/w Leading Gen Jet and Leading Reco Jet 
    //---------------------------------------------------------------------------
    
    // Safety check: Ensure both arrays actually have at least one jet!
    if (jetreco.size() > 0 && jetgen.size() > 0) {

      float dR = delR(jetreco.at(0).v, jetgen.at(0).v);
      h.jet[11] -> Fill(dR);
      
      // If the leading jets are an exact geometric match (dR < 0.3)
      if (dR < 0.3) {
	
	//ploting vars for verification
	h.jet[12] -> Fill(dR);
	
	h.jet[13]->Fill(jetgen.at(0).v.E());
	
	h.recojet[11]->Fill(jetreco.at(0).v.E()); 
		
        //-----------------------------------------------------
	//PLOT THE LUND PLANE FOR GENJETS (hist index = 2)
        //-----------------------------------------------------
        int orig_gen_idx = jetgen.at(0).ind; //index of our leading jet
	
	//using a pointer instead of a var to save RAM
        Jet *g_jet = (Jet*)brGenJet->At(orig_gen_idx);//pointing to the truth level jet at the above index

	//creating my pseudojets array to give as an input to the LJP func later
	std::vector<fastjet::PseudoJet> gen_constituents;
        
        // Loop over every single internal constituent inside this specific truth jet
	//take the jet g_jet -> get its contituents
        for (int c = 0; c < g_jet->Constituents.GetEntriesFast(); c++) {
          TObject *obj = g_jet->Constituents.At(c);//ROOT stores everything as a generic base-class object called a TObject
          if (!obj) continue; // if the pointer in empty then continue with the loop

	  //DOWNCASTING-
	  //rn the compiler has no idea what specific data is inside obj
	  // Because this is a GenJet,we logically know that its constituents are ALWAYS GenParticles.
	  //we are telling the code that this specific piece of the jet is a Truth-Level Particle
          GenParticle *par = (GenParticle*)obj;
	  
	  //FastJet has no idea what a GenParticle is. It uses its own custom class called PseudoJet.
	  //So,we read the basic properties of the particle and translate them into a PseudoJet, which is the only language FastJet understands.
          gen_constituents.push_back(fastjet::PseudoJet(par->Px, par->Py, par->Pz, par->E));
        }

	//LJP measures how a jet splits apart. if the jet has only 1 particle in it, it can't split
	if (gen_constituents.size() >= 2) {
	  // Package the particles and give them as input to the func
          fastjet::PseudoJet genjet_aftermatching = fastjet::join(gen_constituents);
	  //calling our func
          AnalyzeLundPlane(genjet_aftermatching, 2); 
        }
	
	
        //-----------------------------------------------------
	//PLOT THE LUND PLANE FOR RECOJETS (hist index = 3)
        //-----------------------------------------------------
        int orig_reco_idx = jetreco.at(0).ind; //index of our leading jet
	
        Jet *r_jet = (Jet*)brJet->At(orig_reco_idx);//pointing to the truth level jet at the above index

	//creating my pseudojets array to give as an input to the LJP func later
	std::vector<fastjet::PseudoJet> reco_constituents;
        
        // Loop over every single internal constituent inside this specific truth jet
	//take the jet g_jet -> get its contituents
        for (int c = 0; c < r_jet->Constituents.GetEntriesFast(); c++) {
          TObject *obj = r_jet->Constituents.At(c);//ROOT stores everything as a generic base-class object called a TObject
          if (!obj) continue; // if the pointer in empty then continue with the loop

	  //we are telling the code that this specific piece of the jet is a Truth-Level Particle
          GenParticle *par = (GenParticle*)obj;
	  
	  //FastJet has no idea what a GenParticle is. It uses its own custom class called PseudoJet.
	  //So,we read the basic properties of the particle and translate them into a PseudoJet, which is the only language FastJet understands.
	  reco_constituents.push_back(fastjet::PseudoJet(par->Px, par->Py, par->Pz, par->E));
        }

	//LJP measures how a jet splits apart. if the jet has only 1 particle in it, it can't split
	if (reco_constituents.size() >= 2) {
          // Package the particles and give them as input to the func
          fastjet::PseudoJet recojet_aftermatching = fastjet::join(reco_constituents);
	  //callint our func
          AnalyzeLundPlane(recojet_aftermatching, 3); 
        }


      }//if dr< 0.3 ends
      
    }//jet size



    
    //---------------------------------------------------------------------------
    //dR Matching b/w Leading Gen Jet and Leading Reco Jet 
    //---------------------------------------------------------------------------

    if (bjetreco.size() > 0 && jetgen.size() > 0) {
      
      // Calculate the distance STRICTLY between Leading Reco b-Jet and Leading Gen Jet
      float dR_bjet = delR(bjetreco.at(0).v, jetgen.at(0).v);
      h.recobjet[14] -> Fill(dR_bjet);
      
      // If they are an exact geometric match (dR < 0.3)
      if (dR_bjet < 0.3) {
	h.recobjet[15] -> Fill(dR_bjet);
	
	h.jet[17]->Fill(jetgen.at(0).v.E());
	
	h.recobjet[16]->Fill(jetreco.at(0).v.E()); 


	//-----------------------------------------------------
	//PLOT THE LUND PLANE FOR GENJETS (hist index =5)
        //-----------------------------------------------------
        int orig_gen_idx = jetgen.at(0).ind; //index of our leading jet
	
	//using a pointer instead of a var to save RAM
        Jet *g_jet = (Jet*)brGenJet->At(orig_gen_idx);//pointing to the truth level jet at the above index

	
	//creating my pseudojets array to give as an input to the LJP func later
	std::vector<fastjet::PseudoJet> gen_constituents;
        
        // Loop over every single internal constituent inside this specific truth jet
	//take the jet g_jet -> get its contituents
        for (int c = 0; c < g_jet->Constituents.GetEntriesFast(); c++) {
          TObject *obj = g_jet->Constituents.At(c);//ROOT stores everything as a generic base-class object called a TObject
          if (!obj) continue; // if the pointer in empty then continue with the loop

	  //DOWNCASTING-
	  //rn the compiler has no idea what specific data is inside obj
	  // Because this is a GenJet,we logically know that its constituents are ALWAYS GenParticles.
	  //we are telling the code that this specific piece of the jet is a Truth-Level Particle
          GenParticle *par = (GenParticle*)obj;
	  
	  //FastJet has no idea what a GenParticle is. It uses its own custom class called PseudoJet.
	  //So,we read the basic properties of the particle and translate them into a PseudoJet, which is the only language FastJet understands.
          gen_constituents.push_back(fastjet::PseudoJet(par->Px, par->Py, par->Pz, par->E));
        }

	//LJP measures how a jet splits apart. if the jet has only 1 particle in it, it can't split
	if (gen_constituents.size() >= 2) {
	  // Package the particles and give them as input to the func
          fastjet::PseudoJet genjet_aftermatching = fastjet::join(gen_constituents);
	  //calling our func
          AnalyzeLundPlane(genjet_aftermatching, 5); 
        }
	

	
        //-------------------------------------------------------------------------
        // PLOT THE LUND PLANE FOR THE MATCHED RECO B-JET (hist index = 6)
        //-------------------------------------------------------------------------
	int orig_reco_idx = bjetreco.at(0).ind; //index of our leading jet
	
        Jet *r_jet = (Jet*)brJet->At(orig_reco_idx);//pointing to the truth level jet at the above index

	//creating my pseudojets array to give as an input to the LJP func later
	std::vector<fastjet::PseudoJet> reco_constituents;
        
        // Loop over every single internal constituent inside this specific truth jet
	//take the jet g_jet -> get its contituents
        for (int c = 0; c < r_jet->Constituents.GetEntriesFast(); c++) {
          TObject *obj = r_jet->Constituents.At(c);//ROOT stores everything as a generic base-class object called a TObject
          if (!obj) continue; // if the pointer in empty then continue with the loop

	  //we are telling the code that this specific piece of the jet is a Truth-Level Particle
          GenParticle *par = (GenParticle*)obj;
	  
	  //FastJet has no idea what a GenParticle is. It uses its own custom class called PseudoJet.
	  //So,we read the basic properties of the particle and translate them into a PseudoJet, which is the only language FastJet understands.
	  reco_constituents.push_back(fastjet::PseudoJet(par->Px, par->Py, par->Pz, par->E));
        }

	//LJP measures how a jet splits apart. if the jet has only 1 particle in it, it can't split
	if (reco_constituents.size() >= 2) {
          // Package the particles and give them as input to the func
          fastjet::PseudoJet recojet_aftermatching = fastjet::join(reco_constituents);
	  //callint our func
          AnalyzeLundPlane(recojet_aftermatching, 6); 
        }

       

      } // if dR_bjet < 0.3 ends
      
    } // bjet size safety check ends
  
   
    //MET
    float missinget = mymet.Pt();
    float metphi = mymet.Phi();
    

  
    
  }//end loop over events.
}//end of delzana


//Function to get the Lund Jet Plane

//To call the func - AnalyzeLundPlane(Pseudo Jet Name, hist no to fill);
void delZAna::AnalyzeLundPlane(fastjet::PseudoJet inputJet, int histIndex)
  
{
  //EXTRACT CONSTITUENTS
  std::vector<fastjet::PseudoJet> constituents = inputJet.constituents();

  if (inputJet.constituents().size() < 2) return;

  int split_counter = 0; //to count no of splits
  //everytime we unmerge a jet into 2 peices - thats 1 split 


  //####################################################################
  //CA Reclustering (The Lund Declustering Tree)
  //####################################################################
  
  fastjet::JetDefinition::Plugin* ee_plugin = new fastjet::EECambridgePlugin(1.0);
  fastjet::JetDefinition ca_jd(ee_plugin);
  
  // fastjet::Recluster reclustering(ca_jd);

  // This creates a perfect memory history locally!
  fastjet::ClusterSequence ca_cs(constituents, ca_jd);
  std::vector<fastjet::PseudoJet> ca_jets = fastjet::sorted_by_E(ca_cs.inclusive_jets());

  if (ca_jets.size() == 0) {
      delete ee_plugin;
      return;
  }


  // Unroll the jet
  // We grab the leading jet from our local CA tree to start the unrolling
  fastjet::PseudoJet jj = ca_jets[0]; 
  fastjet::PseudoJet j1, j2;
  
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

    
    //####################################################################
    //Calling the func to get the variables to make the lund jet plane
    //####################################################################
    
    std::pair<float, float> lund_vars = GetKinematics(j1, j2);

    float log_1_over_theta = lund_vars.first;
    float log_kt           = lund_vars.second;

    // Fill the correct histograms based on Gen (0) or Reco (1)
    //for genjets
    if(histIndex == 0){
      h.jet[15]->Fill(log_kt);
      h.jet[16]->Fill(log_1_over_theta);                
    }
    //for recojets
    else if(histIndex == 1){
      h.recojet[12]->Fill(log_kt);
      h.recojet[13]->Fill(log_1_over_theta);                
    }
    //for recobjets
    else if(histIndex == 4){
      h.recobjet[12]->Fill(log_kt);
      h.recobjet[13]->Fill(log_1_over_theta);                
    }
    

    // Fill the correct Lund Jet Plane
    h.lundPlane[histIndex]->Fill(log_1_over_theta, log_kt);
    
    split_counter++;

    
    //store the harder particle j1, throw away the emitted particle j2
    jj = j1;

  }// while has parenys loop ends
  //The loop stops when you hit a particle (jj = j1) that cannot be un-merged any further


  // Save the split count
  if (histIndex == 0) h.jet[14]->Fill(split_counter); //no of splits for genjets 
  else if (histIndex == 1) h.recojet[14]->Fill(split_counter); //no of splits for recojets 
  else if (histIndex == 4) h.recobjet[11]->Fill(split_counter); //no of splits for recobjets
  
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
    









// SORT FUNCTION
void delZAna::Sort(int opt)
{
  //arrange by descending pT
  auto sortByPt = [](const Lepton &a, const Lepton &b) { 
      return a.v.Pt() > b.v.Pt(); 
  };
  
  // Option 1: Sort Leptons
  if (opt == 1) {
    std::sort(goodMu.begin(), goodMu.end(), sortByPt);
    std::sort(goodTau.begin(), goodTau.end(), sortByPt);
    std::sort(goodEle.begin(), goodEle.end(), sortByPt); 
  }
  
  // Option 2: Sort Jets
  else if (opt == 2) {
    std::sort(jetgen.begin(), jetgen.end(), sortByPt);
    std::sort(jetreco.begin(), jetreco.end(), sortByPt); 
  }
}

/*
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
*/


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

 
  //Genjet
  h.jet[0]         =   new TH1F("genjet_n",          "No of GenJets",        6,0,6);
  h.jet[1]         =   new TH1F("genjet_E",          "GenJet Energy",        50,0,100); 
  h.jet[2]         =   new TH1F("genjet_Theta",      "GenJet Theta",         64,0,3.2); 
  h.jet[3]         =   new TH1F("genjet_Phi",        "GenJet Phi",           64,-3.2,3.2); 
   
  h.jet[4]         =   new TH1F("genjet0_E",         "Leading GenJet Energy",          50,0,100); 
  h.jet[5]         =   new TH1F("genjet0_theta",     "Leading GenJet theta",           64,0,3.2); 
  h.jet[6]         =   new TH1F("genjet0_Phi",       "Leading GenJet Phi",             64,-3.2,3.2); 

  h.jet[7]         =   new TH1F("genjet1_E",          "Sub Leading GenJet Energy",     50,0,100); 
  h.jet[8]         =   new TH1F("genjet1_theta",      "Sub Leading GenJet theta",      64,0,3.2); 
  h.jet[9]         =   new TH1F("genjet1_Phi",        "Sub Leading GenJet Phi",        64,-3.2,3.2); 

  h.jet[10]        =   new TH1F("genjet_invar_mass",  "Invariant Mass of Leading and subleading Jets",100,0,100); 

 //dr b/w genjet0 and recojet0
  h.jet[11]         =   new TH1F("gen0_reco0_dR",        "Delta R between Leading GenJet & RecoJet",        40,0,8); 

  //dr b/w genjet0 and recojet0 after cut
  h.jet[12]         =   new TH1F("gen0_reco0_dR_aftermatching",        "Delta R between Leading GenJet & RecoJet- After Matching",        40,0,0.4); 

   
   //E of the jet after dr cut
  h.jet[13]         =   new TH1F("genjet0_E_after_dr_matching_w_recojet",          "GenJet Energy - After Matching",        50,0,100); 
     
  h.jet[17]         =   new TH1F("genjet0_E_after_dr_matching_w_recobjet",          "GenJet Energy - After Matching",        50,0,100); 



  

  //RecoJet
  h.recojet[0]         =   new TH1F("recojet_n",          "No of RecoJets",        6,0,6);
  h.recojet[1]         =   new TH1F("recojet_E",          "RecoJet Energy",        50,0,100); 
  h.recojet[2]         =   new TH1F("recojet_Theta",      "RecoJet Theta",         64,0,3.2); 
  h.recojet[3]         =   new TH1F("recojet_Phi",        "RecoJet Phi",           64,-3.2,3.2); 
  
  
  h.recojet[4]         =   new TH1F("recojet0_E",         "Leading RecoJet Energy",          50,0,100); 
  h.recojet[5]         =   new TH1F("recojet0_theta",     "Leading RecoJet theta",           64,0,3.2); 
  h.recojet[6]         =   new TH1F("recojet0_Phi",       "Leading RecoJet Phi",             64,-3.2,3.2); 


  h.recojet[7]         =   new TH1F("recojet1_E",          "Sub Leading RecoJet Energy",     50,0,100); 
  h.recojet[8]         =   new TH1F("recojet1_theta",      "Sub Leading RecoJet theta",      64,0,3.2); 
  h.recojet[9]         =   new TH1F("recojet1_Phi",        "Sub Leading RecoJet Phi",        64,-3.2,3.2); 

  h.recojet[10]        =   new TH1F("recojet_invar_mass",  "Invariant Mass of Leading and subleading RecoJet",100,0,100); 

  h.recojet[11]         =   new TH1F("recojet0_E_after_dr_matching",          "RecoJet Energy - After Matching",        50,0,100);
 


  
 
 




  
  //RecobJet

  //Btag score

  //  h.btag_flag = new TH1F("btag_flag", "b-tag Flag", 2, -1, 2);//this is making 2 bins in the range -1 to  2 ie 3
  h.btag_flag = new TH1F("btag_flag", "b-tag Flag", 100, 0, 1);
  h.btag_flag_aftercut = new TH1F("btag_flag_aftercut", "b-tag Flag- After cut", 100, 0, 1);

  
  h.recobjet[0]         =   new TH1F("recobjet_n",          "No of RecoJets",        6,0,6);
  h.recobjet[1]         =   new TH1F("recobjet_E",          "RecoJet Energy",        50,0,100); 
  h.recobjet[2]         =   new TH1F("recobjet_Theta",      "RecoJet Theta",         64,0,3.2); 
  h.recobjet[3]         =   new TH1F("recobjet_Phi",        "RecoJet Phi",           64,-3.2,3.2); 
  
  
  h.recobjet[4]         =   new TH1F("recobjet0_E",         "Leading RecoJet Energy",          50,0,100); 
  h.recobjet[5]         =   new TH1F("recobjet0_theta",     "Leading RecoJet theta",           64,0,3.2); 
  h.recobjet[6]         =   new TH1F("recobjet0_Phi",       "Leading RecoJet Phi",             64,-3.2,3.2); 


  h.recobjet[7]         =   new TH1F("recobjet1_E",          "Sub Leading RecoJet Energy",     50,0,100); 
  h.recobjet[8]         =   new TH1F("recobjet1_theta",      "Sub Leading RecoJet theta",      64,0,3.2); 
  h.recobjet[9]         =   new TH1F("recobjet1_Phi",        "Sub Leading RecoJet Phi",        64,-3.2,3.2); 

  h.recobjet[10]        =   new TH1F("recobjet_invar_mass",  "Invariant Mass of Leading and subleading RecoJet",100,0,100);


  //dr b/w genjet0 and recobjet0
  h.recobjet[14]         =   new TH1F("recobjet0_genjet0_dR",        "Delta R between Leading GenJet & RecoJet",        40,0,8); 

  //dr b/w genjet0 and recobjet0 - after cut
   h.recobjet[15]         =   new TH1F("recobjet0_genjet0_dR_aftercut",        "Delta R between Leading GenJet & RecoJet - After cut",        40,0,0.4);

   //E of the recobjet after dr cut
  h.recobjet[16]         =   new TH1F("recobjet0_E_after_dr_matching",          "RecobJet Energy - After Matching",        50,0,100);

  





  //lund jet plane var
  
  h.jet[14] = new TH1F("n_splits_genjet",                    "Number of Lund Splits per genJet; Splits; Events",                              30, 0, 30);
  h.jet[15] = new TH1F("lund_kt_genjet",                     "Lund Split Momentum for GenJets; ln(k_{t}) [GeV]; Number of Splits",                     50, -8, 6);
  h.jet[16] = new TH1F("lund_theta_genjet",                  "Lund Split Angle for GenJets #theta; ln(1.0 / #theta); Number of Splits",                50, 0, 7);  

   
  h.recojet[14] = new TH1F("n_splits_recojet",                    "Number of Lund Splits per recoJet; Splits; Events",                              30, 0, 30);   
  
  h.recojet[12] = new TH1F("lund_kt_recojet",                     "Lund Split Momentum for RecoJets; ln(k_{t}) [GeV]; Number of Splits",                     50, -15, 5);
  h.recojet[13] = new TH1F("lund_theta_recojet",                  "Lund Split Angle for RecoJets #theta; ln(1.0 / #theta); Number of Splits",                50, -1, 14);
  
  
  h.recobjet[11] = new TH1F("n_splits_recobjet",                    "Number of Lund Splits per RecobJet; Splits; Events",                              30, 0, 30);   
  
  h.recobjet[12] = new TH1F("lund_kt_recobjet",                     "Lund Split Momentum for RecobJets; ln(k_{t}) [GeV]; Number of Splits",                     50, -15, 5);
  h.recobjet[13] = new TH1F("lund_theta_recobjet",                  "Lund Split Angle for RecobJets #theta; ln(1.0 / #theta); Number of Splits",                50, - 1, 14);

 

   
  //Lund Plane
  // X-axis: ln(1/theta) from 0 to 7
  // Y-axis: ln(kt) from -8 to 6
  
  h.lundPlane[0] = new TH2F("lund_plane_genjet", "Primary Lund Jet Plane for GenJets; ln(1.0 / #theta); ln(k_{t}) [GeV]", 50, 0, 10, 50, -6, 4);
  
  h.lundPlane[1] = new TH2F("lund_plane_recojet", "Primary Lund Jet Plane for RecoJets; ln(1.0 / #theta); ln(k_{t}) [GeV]", 50, 0, 10, 50, -8, 4);
  
  h.lundPlane[4] = new TH2F("lund_plane_recobjet", "Primary Lund Jet Plane for Reco bJets; ln(1.0 / #theta); ln(k_{t}) [GeV]", 50, 0, 10, 50, -8, 4);
  


  //After Matching 

  h.lundPlane[2] = new TH2F("lund_plane_genjet_aftermatching_with_recojet", "Primary Lund Jet Plane for GenJets - After Matching; ln(1.0 / #theta); ln(k_{t}) [GeV]", 50, 0, 7, 50, -6, 4);
  
  h.lundPlane[3] = new TH2F("lund_plane_recojet_aftermatching", "Primary Lund Jet Plane for RecoJets - After Matching; ln(1.0 / #theta); ln(k_{t}) [GeV]", 50, 0, 7, 50, -6, 4);
  
  h.lundPlane[5] = new TH2F("lund_plane_genjet_aftermatching_with_recobjet", "Primary Lund Jet Plane for RecobJets - After Matching; ln(1.0 / #theta); ln(k_{t}) [GeV]", 50, 0, 7, 50, -6, 4);
  
  h.lundPlane[6] = new TH2F("lund_plane_recobjet_aftermatching", "Primary Lund Jet Plane for Reco bJets - After Matching; ln(1.0 / #theta); ln(k_{t}) [GeV]", 50, 0, 7, 50, -8, 4);
  








  // Note the calling of Sumw2() for each histogram we declare.
  // This is needed so that if scale the histograms later, the
  // errors on the points are still correct after scaling.
  // This is relevant.. so ask someone if you don't understand.
  
}

/*
  //JET MATCHING
    for (int ir = 0; ir < (int)jetreco.size(); ir++) {
      
      //initiallizing variables   
      float min_dR = 999.0; //Before searching we assumes the distance to the nearest jet is impossibly huge
      int matched_index = -1;//an invalid array index signifing that no match has been found.
      
      for (int i = 0; i < (int)jetgen.size(); i++) {
	float current_dR = delR(jetreco.at(ir).v, jetgen.at(i).v); //delta r calculation b/w recojet ir and all the genjets
	if (current_dR < min_dR) {
          min_dR = current_dR;
          matched_index = i; // Remember which GenJet was closest
	}
      }
       
      // If the closest GenJet is within dR < 0.5, it is a match
      if (min_dR < 0.3 && matched_index != -1) {
	
      }
    }//matching loop ends 
    

 */


    /*
    //##################################################################
    // RECOJETS 
    //##################################################################
    std::vector<fastjet::PseudoJet> all_particles;

    // Declare this outside the if block so it survives for the matching later!
    std::vector<fastjet::PseudoJet> manual_jets;
    
    // Scan the whole event for stable particles
    for (int no = 0; no < brParticle->GetEntries(); no++) {
      GenParticle *par = (GenParticle*)brParticle->At(no);
      if (par->Status == 1) {
        all_particles.push_back(fastjet::PseudoJet(par->Px, par->Py, par->Pz, par->E));
      }
    }
    
    if (all_particles.size() >= 2) {
      // Cluster the raw event
      fastjet::JetDefinition temp_jd(fastjet::ee_genkt_algorithm, 1.0, -1);
      fastjet::ClusterSequence temp_cs(all_particles, temp_jd);
      // std::vector<fastjet::PseudoJet> manual_jets = fastjet::sorted_by_E(temp_cs.inclusive_jets());
      manual_jets = fastjet::sorted_by_E(temp_cs.inclusive_jets());
      
      for (int im = 0; im < (int)manual_jets.size(); im++) {
        fastjet::PseudoJet myJet = manual_jets[im];
        
        //Only keep jets with pT > 10 GeV
        if (myJet.pt() < 5.0) continue;
	
        // Save its kinematics to your array
        Lepton temp;
        temp.v.SetPtEtaPhiM(myJet.pt(), myJet.eta(), myJet.phi(), myJet.m());
        temp.ind = im;
        jetreco.push_back(temp);

	AnalyzeLundPlane(myJet, 1);
	
      }
    }
    */












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
