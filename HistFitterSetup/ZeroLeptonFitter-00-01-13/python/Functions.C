#include <vector>
#include <algorithm>
#include <iostream>
//#include <map>
using namespace std;

#include "TLorentzVector.h"

auto kappaMapFunction = [](int region){
  // if (region == "Meff" ) return 1.55;
  //if (region == "SRG" )  return 1.55;
  // if (region == "SRC" )  return 1.55;
  // if (region == "SRS" )  return 2.00;


  // if (region == 0 )  return 1.55;
  // if (region == 1 )  return 1.55;
  //  if (region == 2 )  return 3.2;
  // if (region == 3 )  return 2.00;
  return 1.334;
};

float gammaCorWeight(int RunNumber, int regionEnum = 0){// std::string region = "Meff") {
  if (RunNumber >=361039 && RunNumber <= 361061+1){
    if (regionEnum == 1)
      return 1.668;
    else
      return 1.56;
  }
  else return 1.0;
}

