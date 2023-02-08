//
// Created by Stephane Zsoldos on 7/6/22.
//

#include <csignal>

#include <TH2D.h>

#include "Recon.hh"
#include "Templates/TData.hh"

#include <SnD/Multilateration.hh>
#include <SnD/Recon.hh>
#include <SnD/Map.hh>

#include <ROOT/Utils.hh>

ReconAnalysis::ReconAnalysis(const char *pdfname, const char *histname, const char* perpmthistname ,
							 const double &R, const double &HH,
							 int me, int a, int ms,
							 bool im, const char *mn,
							 bool iv,
							 bool ib, bool iu, bool ip,
							 bool itt,
               const char *filename,
							 const char *treename)
	: nMaxEvts(me), algo(a), max_seed(ms), ismap(im), mapname(mn), isverbose(iv), isbinned(ib), isunbinned(iu), isperpmt(ip), istrigtime(itt) {
  //
  hPDF = GetROOTObj<TH2D>(pdfname, histname)->ProjectionX("hPDF");
  std::cout << "Load PDF: " << hPDF->GetName() << std::endl;
  mPDF2D = GetROOTMObj<TH2D>(pdfname, perpmthistname, "TH2D");
  std::transform(
	  mPDF2D.begin(), mPDF2D.end(),
	  std::inserter(mPDF1D, mPDF1D.begin()),
	  [](const std::pair<int, TH2D*>& p){
		return std::make_pair(p.first, p.second->ProjectionX());
	  }
  );
  //
  Cyl = new Cylinder(R, HH);
  //
  OFile = new TFile(filename, "RECREATE");
  Tree = new TTree(treename, treename);
  RT.SetTree(Tree);
  //
  max_seed = max_seed < 0 ? std::numeric_limits<int>::max() : max_seed;
  nMaxEvts = nMaxEvts < 0 ? std::numeric_limits<int>::max() : nMaxEvts;
  //
  if(!isbinned && !isunbinned && !isperpmt)
	isbinned = true;
}
ReconAnalysis::~ReconAnalysis(){
  delete OFile;
  delete Cyl;
  delete hPDF;
  for(auto& p : mPDF2D)
	delete p.second;
}

void ReconAnalysis::Do(void *Data) {

  // Get Data
  auto wData = static_cast<TData*>(Data);
  auto vHits = wData->GetVHits();
  auto iEvt = wData->GetEventID();
  auto iTrig = wData->GetTriggerID();
  const char *tag = Form("Evt%d_Trigger%d", iEvt, iTrig);
  //
  if(iEvt > nMaxEvts)
	raise(SIGINT);

  // Get centroid seed
  TVector3 Centroid = GetCentroid(vHits);

  // Get time seed
  double TSeed = Cyl->GetTWall(Centroid);
  if(istrigtime)
	TSeed = 0;

  // Get SnD seeds
  std::vector<PosT> vSeeds = GetVPosTSeeds(vHits, hPDF, Cyl, max_seed, istrigtime);
  vSeeds.emplace_back(Centroid, TSeed);

  // Set bounds limit depending on trigtime
  auto fPosBounds = istrigtime ? SetPosBounds: SetBounds;

  // Recon
  if(isunbinned){
	FitStruct FS = {vHits, hPDF};
	RT = Recon(&FS, Cyl, vSeeds, GetAlgo(algo), fPosTU, {fPosBounds, SetPars});
	// Fill
	Tree->Fill();
  } else if(isperpmt) {
	FitMapStruct FMS = {vHits, mPDF1D};
	RT = Recon(&FMS, Cyl, vSeeds, GetAlgo(algo), fPosTPerPMT, {fPosBounds, SetPars});
	// Fill
	Tree->Fill();
  } else {
	FitStruct FS = {vHits, hPDF};
	RT = Recon(&FS, Cyl, vSeeds, GetAlgo(algo), fPosT, {fPosBounds, SetPars});
	// Fill
	Tree->Fill();
  }

  // Verbose
  if(isverbose)
	RT.Print();

  // Map
  if(ismap)
	SaveMap(vHits, hPDF, Cyl, tag, mapname.c_str());

}

void ReconAnalysis::Export() const {
  OFile->cd();
  Tree->Write();
  OFile->Close();
}