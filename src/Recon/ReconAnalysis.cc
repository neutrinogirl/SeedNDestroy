//
// Created by Stephane Zsoldos on 7/6/22.
//

#include <TH2D.h>

#include "ReconAnalysis.hh"

#include <SnD/RATData.hh>
#include <SnD/Multilateration.hh>
#include <SnD/Recon.hh>
#include <SnD/Map.hh>

#include <ROOT/Utils.hh>

ReconAnalysis::ReconAnalysis(const char *pdfname, const char *histname,
							 const double &R, const double &HH,
							 const char *treename){
  hPDF = GetROOTObj<TH2D>(pdfname, histname)->ProjectionX("hPDF");
  Cyl = new Cylinder(R, HH);
  Tree = new TTree(treename, treename);
  RT.SetTree(Tree);
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
  double TSeed = Cyl->GetTWall(Centroid);

  // Get SnD seeds
  std::vector<PosT> vSeeds = GetVPosTSeeds(RData->vHits, hPDF, Cyl, 0, 10);
  vSeeds.emplace_back(Centroid, TSeed);

  // Recon
  RT = Recon(RData->vHits, hPDF, Cyl, vSeeds);

  // Map
  std::vector<TCanvas*> vMap = GetMap(RData->vHits, hPDF, Cyl);
  TFile f("MAP.root", "UPDATE");
  for (auto &c : vMap) {
	c->SetName(Form("%s_%s",
					RData->tag.c_str(), c->GetName()));
	c->Write();
  }
  f.Close();
  for (auto &&obj: *gDirectory->GetList()) {
	if (!std::string(obj->GetName()).find("hGrid_")) {
	  delete obj;
	}
  }
  for(auto &c: vMap) {
	delete c;
  }
  vMap.clear();

  // Fill
  Tree->Fill();
}

#include <TFile.h>
void ReconAnalysis::Export(const char *filename) const {
  TFile f(filename, "RECREATE");
  Tree->Write();
  f.Close();
}