//////////////////////////////////////////////////
// Author: Sourabh Dube
// Date: May 30th, 2014
// Run on Delphes output and have same structure
// as makeclass code
//////////////////////////////////////////////////

#ifndef delZAna_h
#define delZAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>

#include <TH2.h>
#include <TMath.h>
#include "TLorentzVector.h"
#include "TVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TClonesArray.h"
// Header file for the classes stored in the TTree if any.
#include <vector>
#include <fstream>
#include <iostream>
#include <limits>
// external header files needed from Delphes
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

//Added for fastjet
//#include "external/fastjet/PseudoJet.hh" 
#include "external/fastjet/ClusterSequence.hh"
#include "external/fastjet/plugins/EECambridge/fastjet/EECambridgePlugin.hh"
//Recluster is not a member of fastjet. It is categorized as an add-on tool within the FastJet library structure. Thus it needs to be linked in the header file separately.
#include "external/fastjet/tools/Recluster.hh"

//declaration of the delphes classes we shall use
class Particle;
class ExRootTreeReader;
class GenParticle;
class Photon;
class Electron;
class Muon;
class Jet;
class Track;
class MissingET;


// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;

class delZAna {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // List of branches

  ExRootTreeReader *treeReader;// = new ExRootTreeReader(&chain);
  //declaring the branches which hold the classes
  TClonesArray  *brMissingET;
  //TClonesArray  *brScalarHT;
  TClonesArray  *brJet;
  TClonesArray  *brGenJet;
  //TClonesArray  *brFatJet;
  //TClonesArray  *brPhoton;
  //TClonesArray  *brTrack;
  //TClonesArray  *brMuon;
  //TClonesArray  *brEle;
  TClonesArray  *brParticle;

  delZAna(TTree *tree=0);
  virtual ~delZAna();
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  //User added functions
  //You add functions here
  void Begin();
  void End();
  void BookHistograms();
  void SetHstFileName(const char *HstFileName){ _HstFileName = HstFileName;} //const 
  void SetSumFileName(const char *SumFileName){ _SumFileName = SumFileName;} //const 
  void SetVerbose(int verbose){ _verbosity = verbose; }
  void Sort(int opt);

  float delta_phi(float phi1, float phi2);
  Double_t delR(TLorentzVector v1,  TLorentzVector v2);
  void TauBuilder();
  float calc_iso(float,TLorentzVector,int);


  //for fastjets
  void AnalyzeLundPlane(Jet* targetJet, int histIndex);
  void AnalyzeLundPlane(fastjet::PseudoJet manualJet, int histIndex);

  // Helper function to calculate the Lund kinematics of a single split
  std::pair<float, float> GetKinematics(fastjet::PseudoJet j1, fastjet::PseudoJet j2);


  
public:
  struct Hists {
    //Declare the histograms you want in here.
    TH1F *etmiss,*ngoodtau,*goodtaupt;
    TH1F *massjj[2];
    TH1F *nmuons;
    TH1F *nele;
    TH1F *njets;
    TH1F *jetpt , *jeteta, *jetphi ;
    
     
    TH1F *jetpt0 , *jeteta0, *jetphi0, *jet_invarmass;
    TH1F *jetpt1 , *jeteta1, *jetphi1 ;
    TH1F *ele[20], *jet[50], *recojet[50] ;
    TH1F *eventselect[20];

    TH2F *lundPlane[50];
     
  };
  struct Lepton {
    TLorentzVector v;
    int id;
    int ind;
    float wt;
    int flavor;
  };

   
 protected:
  Hists h;

private:
  TFile *_HstFile;
  const char *_HstFileName; //const 
  const char *_SumFileName; // const 
  int _verbosity;
  float GEV, MEV2GEV;
  int nEvtTotal, njet;
  vector<Lepton> goodMu;
  vector<Lepton> goodTau;
  vector<Lepton> jetgen;
  vector<Lepton>  goodbJets ;
  vector<Lepton>  goodEle ;
  vector<Lepton>  leptons ;
  vector<Lepton> jetreco;
};

#endif

#ifdef delZAna_cxx
delZAna::delZAna(TTree *tree) : fChain(0) 
{
  Init(tree);
}
delZAna::~delZAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}
void delZAna::Init(TTree *tree)
{
  if(!tree) return;
  fChain = tree;
  
  treeReader = new ExRootTreeReader(fChain);
  brMissingET = treeReader->UseBranch("MissingET");
  //brScalarHT = treeReader->UseBranch("ScalarHT");
  brJet = treeReader->UseBranch("Jet");
  brGenJet = treeReader->UseBranch("GenJet");
  //  brFatJet = treeReader->UseBranch("FatJet");
  //brPhoton = treeReader->UseBranch("Photon");
  //brTrack = treeReader->UseBranch("Track");
  //brMuon= treeReader->UseBranch("Muon");
  //brEle = treeReader->UseBranch("Electron");
  brParticle= treeReader->UseBranch("Particle");
}
#endif
