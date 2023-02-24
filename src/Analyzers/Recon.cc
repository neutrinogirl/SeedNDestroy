//
// Created by Stephane Zsoldos on 7/6/22.
//

#include <csignal>

#include "Recon.hh"

#include "Templates/TData.hh"
#include "Algo/WOpt.hh"
#include "Algo/VHits.hh"
#include "ROOT/Utils.hh"

ReconAnalysis::ReconAnalysis(const char *pdfname, const char *histname, const char* perpmthistname ,
							 const double &R, const double &HH,
							 int me, int a, int ms,
							 bool iv,
							 bool ib, bool iu, bool ip,
							 bool iat,
							 bool ijs,
							 const char *filename,
							 const char *treename)
	: nMaxEvts(me), algo(a), max_seed(ms),
	  isverbose(iv),
	  isbinned(ib), isunbinned(iu), isperpmt(ip),
	  isapplytrigger(iat),
	  isjustseed(ijs){
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
  std::vector<Hit> vHits = wData->GetVHits();
  //
  if(isapplytrigger){
	double T = GetFirstHitTime(vHits, 1.f);
	std::transform(
		vHits.begin(),
		vHits.end(),
		vHits.begin(),
		[T](const Hit& h){
		  return h-T;
		}
	);
  }
  //
  int iEvt = wData->GetEventID();
  //
  if(iEvt > nMaxEvts)
	raise(SIGINT);

  // Init vSeeds with Centroid
  std::vector<PosT> vSeeds = {
	  GetCentroidBasedSeed(vHits, Cyl)
  };

  // Get LS seed
  vSeeds.emplace_back(GetLSBasedSeed(vHits, Cyl, vSeeds));

  auto ConvertAndFill = [&](const TVector3& p){
	TVector3 v3mm = ConvertTVector3Unit<double>(p, SpaceUnit::dm, SpaceUnit::mm);
	RT.X = v3mm.x();
	RT.Y = v3mm.y();
	RT.Z = v3mm.z();
	Tree->Fill();
  };

  if(isjustseed){
	ConvertAndFill(vSeeds.back().GetTVector3());
	return;
  }


  std::vector<void (*)(nlopt::opt &opt, Bnd *c)> vfPars = {
	  SetBounds,
	  SetPars,
  };
  if(algo == 2)
	vfPars.emplace_back(SetInequalityConstraint);

  // Recon
  if(isunbinned){
	FitStruct FS = {vHits, hPDF};
	RT = Recon(&FS, Cyl, vSeeds, GetAlgo(algo), fPosTU, vfPars);
	// Fill
	// Convert back to mm
	ConvertAndFill(RT.GetTVector3());
  } else if(isperpmt) {
	FitMapStruct FMS = {vHits, mPDF1D};
	RT = Recon(&FMS, Cyl, vSeeds, GetAlgo(algo), fPosTPerPMT, vfPars);
	// Fill
	// Convert back to mm
	ConvertAndFill(RT.GetTVector3());
  } else {
	FitStruct FS = {vHits, hPDF};
	RT = Recon(&FS, Cyl, vSeeds, GetAlgo(algo), fPosT, vfPars);
	// Fill
	// Convert back to mm
	ConvertAndFill(RT.GetTVector3());
  }

  // Verbose
  if(isverbose)
	RT.Print();

}

void ReconAnalysis::Export() const {
  OFile->cd();
  Tree->Write();
  OFile->Close();
}