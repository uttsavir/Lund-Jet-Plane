#include <iostream>
#include <stdlib.h>
#include "TChain.h"
#include "delZAna.h"

//This is the driver script, which becomes our main program
//Here we set the options which we wish to use, which files
//we want to run over and so on.

//the argument decides what input sample we want to run over.
//we give separate names of output files for each set of
//input files.

int main(int argc, char *argv[])
{
  using std::cout;
  using std::endl;

  if(argc<2){
    cout<<" Please give one integer argument "<<endl;
    return 0;
  }
  int sample_id = atoi(argv[1]);

  
  const char *hstfilename, *sumfilename;

  TChain *chain = new TChain("Delphes");

  if(sample_id==1){
    chain->Add("input/ee_lightjetsNLO_50k_91GeV_op.root");
    hstfilename = "output_files/hst_ee_lightjetsNLO_25may.root";//output histogram file
    sumfilename = "output_files/sum_ee_lightjetsNLO_25may.txt";//output text file
  }


  if(sample_id==2){
    chain->Add("input/ee_lightjetNLO_27thmay.root");
    hstfilename = "output_files/hst_ee_lightjetsNLO_27may.root";//output histogram file
    sumfilename = "output_files/sum_ee_lightjetsNLO_27may.txt";//output text file
  }
  


  
  cout<<"Set Options... ";
  std::cout<<"Output files are "<<hstfilename<<" and "<<sumfilename<<std::endl;
  
  delZAna t(chain);//declared an instance of our class.
  t.SetHstFileName(hstfilename);
  t.SetSumFileName(sumfilename);
  t.SetVerbose(1500);//set verbosity level for output.
  cout<<" Processing... ";
  t.Begin();
  t.Loop();
  t.End();

  cout<<"Completed."<<endl;
  return 0;
}
