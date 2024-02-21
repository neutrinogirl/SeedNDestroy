//
// Created by Stephane Zsoldos on 7/6/22.
//

//
#include "Templates/TData.hh"

//
#include "TAnalyzer.hh"

// #### #### #### #### #### #### #### #### #### #### #### #### //
void TAnalyzer::Do(void *Data) {

  //
  auto fData = static_cast<TData*>(Data);
  //
  std::cout << "E:" << fData->GetEnergy() << std::endl;
  auto vHits = fData->GetVHits();
  // Print hits
  for(auto& h : vHits)
	std::cout << h << std::endl;
  // Ask user to hit a key to continue
  std::cout << "Press any key to continue..." << std::endl;
  std::cin.ignore();
  //

}

#include <TFile.h>
void TAnalyzer::Export(const char *filename) {
  TFile f(filename, "RECREATE");
  f.cd();
  f.Close();
}
