//
// Created by Stephane Zsoldos on 7/6/22.
//

#include "TAppAnalysis.hh"

// #### #### #### #### #### #### #### #### #### #### #### #### //
void TAppAnalysis::Do(TData *Data) {
  for(const auto& hit : Data->GetVHits()) {
	hit.PMTPos.Print();
  }

}

#include <TFile.h>
void TAppAnalysis::Export(const char *filename) {
  TFile f(filename, "RECREATE");
  f.cd();
  f.Close();
}
