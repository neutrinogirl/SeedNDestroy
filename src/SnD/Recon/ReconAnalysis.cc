//
// Created by Stephane Zsoldos on 7/6/22.
//

#include "ReconAnalysis.hh"

#include <SnD/RATData.hh>
#include <SnD/Multilateration.hh>

#include <ROOT/Utils.hh>

ReconAnalysis::ReconAnalysis(const char *filename, const double &R, const double &HH, const std::string& treename){
  hPDF = GetROOTObj<TH1D>(filename, "");
  Cyl = new Cylinder(R, HH);
  Tree = new TTree(treename.c_str(), treename.c_str());
  Seed.SetTree(Tree);
}
ReconAnalysis::~ReconAnalysis(){
  delete Tree;
  delete Cyl;
  delete hPDF;
}

void ReconAnalysis::Do(void *Data) {

  // Get Data
  auto *RData = static_cast<RATData*>(Data);

  // Get centroid seed
  TVector3 Centroid = GetCentroid(RData->vHits);

  // Get time seed
  double TSeed = Cyl->GetDWall(Centroid);

  // Read
  Seed.Pos = Centroid;
  Seed.T = TSeed;

  // Get SnD seeds
  std::vector<PosT> vSeeds = GetVPosTSeeds(RData->vHits, hPDF, Cyl);

  // Fill Tree
  Tree->Fill();

}

#include <TFile.h>
void ReconAnalysis::Export(const char *filename) {
  TFile f(filename, "RECREATE");
  Tree->Write();
  f.Close();
}