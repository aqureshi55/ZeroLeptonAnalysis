
// cannot change the order of include
#include "summary_harvest_tree_description.h"
#include "contourmacros/m0_vs_m12_nofloat.C"
#include "contourmacros/mX_vs_mY_nofloat.C"

void makecontourhists(const TString& combo = "all", const TString& gridName = "msugra")  {

  std::cout<<" makecontourhists : Start! "<<std::endl;
  if(gridName=="GG_onestepCC"){
    const char* ehistfile = mX_vs_mY_nofloat(combo, 0, "mgluinomlsp_nofloat.root", "mgluino", "mlsp", 100, 100, 0, 1600, 0, 1600);
  }else if(gridName=="SS_onestepCC"){
    const char* ehistfile = mX_vs_mY_nofloat(combo, 0, "msquarkmlsp_nofloat.root", "msquark", "mlsp", 100, 100, 0, 1600, 0, 1600);
  }else if(gridName=="SS_onestepCC_mlsp60"){
    const char* ehistfile = mX_vs_mY_nofloat(combo, 0, "msquarkmlsp_nofloat.root", "msquark", "(mchargino-mlsp)/(msquark-mlsp)", 100, 100, 0, 1600, 0, 1.5);
  }else if(gridName=="GG_onestepCC_mlsp60"){
    const char* ehistfile = mX_vs_mY_nofloat(combo, 0, "msquarkmlsp_nofloat.root", "mgluino", "(mchargino-mlsp)/(mgluino-mlsp)", 100, 100, 0, 1600, 0, 1.5);
  }else if(gridName=="SS_onestepN2C_mlsp60"){
    const char* ehistfile = mX_vs_mY_nofloat(combo, 0, "msquarkmlsp_nofloat.root", "msquark", "(mneutralino2-mlsp)/(msquark-mlsp)", 100, 100, 0, 1600, 0, 1.5);
  }else if(gridName=="GG_onestepN2C_mlsp60"){
    const char* ehistfile = mX_vs_mY_nofloat(combo, 0, "msquarkmlsp_nofloat.root", "mgluino", "(mneutralino2-mlsp)/(mgluino-mlsp)", 100, 100, 0, 1600, 0, 1.5);
  }else if(gridName=="SM_GG_N2"){
    const char* ehistfile = mX_vs_mY_nofloat(combo, 0, "mgluinomlsp2_nofloat.root", "mgluino", "mlsp2", 100, 100, 0, 1600, 0, 1600);
  }else if(gridName=="long_lived_gl"){
    const char* ehistfile = mX_vs_mY_nofloat(combo, 0, "mgluinomlsp2_nofloat.root", "log10(lifetime)", "mgluino", 1000, 1000, 0, 10, 0, 3000);
  }else if(gridName.BeginsWith("SG_pMSSM_")){
    const char* ehistfile = mX_vs_mY_nofloat(combo, 0, "mgluinomsquark_nofloat.root", "mgluino", "msquark", 100, 100, 0, 1600, 0, 1600);
  }else{
    const char* ehistfile = mX_vs_mY_nofloat(combo, 0, "m0m12_nofloat.root", "m0", "m12", 100, 100, 0, 1600, 0, 1600);
  }

  return;
}

