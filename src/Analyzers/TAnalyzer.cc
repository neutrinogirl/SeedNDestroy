//
// Created by Stephane Zsoldos on 7/6/22.
//

//
#include "../Readers/TData.hh"

//
#include "TAnalyzer.hh"

// #### #### #### #### #### #### #### #### #### #### #### #### //
void TAnalyzer::Do(void *Data) {

  //
  auto fData = static_cast<TData*>(Data);
  //
  fData->GetPosition().Print();
  fData->GetDirection().Print();
  std::cout << fData->GetEnergy() << std::endl;
  auto vHits = fData->GetVHits();
  vHits[0].Print();

}

#include <TFile.h>
void TAnalyzer::Export(const char *filename) {
  TFile f(filename, "RECREATE");
  f.cd();
  f.Close();
}
