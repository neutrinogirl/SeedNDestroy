//
// Created by Stephane Zsoldos on 7/10/22.
//

#include "ERecAnalysis.hh"
ERecAnalysis::ERecAnalysis() {

  vEdges = {0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0};

  for(const auto& E: vEdges){
	v2D.push_back(new TH2D(Form("h2D_E%.1f", E), Form("E=%.1f ; N_{hits} ; Q ", E), 100, 0., 2000., 100, 0., 2000.));
  }


}
// 0 == underflow
// 1 == [0., 0.5[
// ...
// n == [n-1, n[
// size() == overflow
int ERecAnalysis::GetInterval(const double& E) {
  return std::upper_bound(vEdges.begin(), vEdges.end(), E) - vEdges.begin();
}

#include <SnD/RATData.hh>
void ERecAnalysis::Do(void *Data) {

  RATData *RData = static_cast<RATData*>(Data);

  const int EBin = GetInterval(RData->E);

  if(EBin > 0 && EBin <= vEdges.size()) {
	v2D[EBin-1]->Fill(RData->NHits, RData->Q);
  }

}
#include <TFile.h>
void ERecAnalysis::Export(const char *filename) {

  TFile *f = new TFile(filename, "RECREATE");
  for(const auto& h: v2D) {
	h->Write();
  }
  f->Close();

}