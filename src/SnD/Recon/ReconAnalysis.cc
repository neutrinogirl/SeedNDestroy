//
// Created by Stephane Zsoldos on 7/6/22.
//

#include "ReconAnalysis.hh"

#include <SnD/RATData.hh>
#include <SnD/Multilateration.hh>

ReconAnalysis::ReconAnalysis(TH1D* h, Cylinder* c, const std::string& treename) {
  hPDF = h;
  Cyl = c;
  Tree = new TTree(treename.c_str(), treename.c_str());
}

void ReconAnalysis::Do(void *Data) {

  // Get Data
  auto *RData = reinterpret_cast<RATData*>(Data);

  // // Get centroid seed
  // TVector3 Centroid = GetCentroid(RData->vHits);
  //
  // // Get time seed
  // double TSeed = Cyl->GetDWall(Centroid);
  //
  // // Get SnD seeds
  // std::vector<PosT> vSeeds = GetVPosTSeeds(RData->vHits, hPDF, Cyl);

  // Fill Tree
  Tree->Fill();

}