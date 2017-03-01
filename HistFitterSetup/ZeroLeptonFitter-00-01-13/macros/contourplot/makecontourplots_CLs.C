#include <vector>
#include "TPRegexp.h"
#include "contourmacros/SUSY_contourplots.C"

void makecontourplots_CLs(const TString Grid="GG_direct",TString dirname = "Outputs/3.2ifb/GG_direct",TString outdir="plots",bool showSR=true,bool unblind=false,bool showAllSRExp=false, bool useTGraph=false, TString SRLabel="", bool removeLowMLine=false) 
{
  int  discexcl(1); // 0=discovery, 1=exclusion
  bool showcms(false);
  bool drawOldLimits(true);
  vector<TString> infilelist;

  if(showAllSRExp) drawOldLimits = false;

  // XSec Nominal Up Down
  TString combined[3]={"combined_fixSigXSecNominal", "combined_fixSigXSecUp","combined_fixSigXSecDown"};

  TString listSuffix="__1_harvest_list";
  if(Grid=="GG_onestepCC"||Grid=="SS_onestepCC"){
     listSuffix="__mlspNE60_harvest_list";
  }
  if(Grid.EndsWith("_mlsp60")){
     listSuffix="__mlspEE60_harvest_list";
  }
  
  for(int i=0; i<3; i++){
    if(useTGraph==0){
       infilelist.push_back(dirname+"/"+Grid+"_"+combined[i]+listSuffix+".root");
     } else {
       infilelist.push_back(dirname+"/output.root");      
     }
     cout<< infilelist[i]<<" ";
  }
  cout<<endl;
  
  TPMERegexp tmp('/');
  int nsplit = tmp.Split(dirname);
  TString lumi="11.3";
  for( int i=0 ; i<nsplit ; i++ ){
    if( tmp[i].Contains("ifb") ){
      lumi = tmp[i].ReplaceAll("ifb","");
    }
  }
  cout<<"lumi="<<lumi<<" fb-1"<<endl;

  (void) SUSY_contourplots(
      infilelist.at(0), infilelist.at(1), infilelist.at(2),
      "0-leptons, 2-6 jets", 
      lumi,
      Grid,
      outdir,
      showSR,
      unblind,
      showAllSRExp,
      drawOldLimits,
      useTGraph,
      discexcl=1,
      SRLabel,
      removeLowMLine);
}

