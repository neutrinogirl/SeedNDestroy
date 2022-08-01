//
// Created by Stephane Zsoldos on 7/6/22.
//

#include <TH2D.h>

#include "ReconAnalysis.hh"

#include "SnD/RATData.hh"
#include "SnD/Multilateration.hh"

#include "ROOT/Utils.hh"

ReconAnalysis::ReconAnalysis(const char *pdfname, const char *histname,
							 const double &R, const double &HH,
							 const char *treename){
  hPDF = GetROOTObj<TH2D>(pdfname, histname)->ProjectionX("hPDF");
  Cyl = new Cylinder(R, HH);
  Tree = new TTree(treename, treename);
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

  Tree->Fill();

  // Get SnD seeds
  std::vector<PosT> vSeeds = GetVPosTSeeds(RData->vHits, hPDF, Cyl);
  Seed = *vSeeds.begin();

  Tree->Fill();

}

#include <TFile.h>
void ReconAnalysis::Export(const char *filename) {
  TFile f(filename, "RECREATE");
  Tree->Write();
  f.Close();
}