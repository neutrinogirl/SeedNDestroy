//
// Created by Stephane Zsoldos on 7/6/22.
//

#include <csignal>
#include <iomanip>
#include <numeric>
#include <algorithm>

#include "Recon.hh"

#include "Templates/TData.hh"
#include "Algo/VHits.hh"
#include "Algo/WOpt.hh"
#include "ROOT/Utils.hh"

ReconAnalysis::ReconAnalysis(const char *pdf_name, const char *hist_name, const char* per_pmt_hist_name,
							 float R, float HH,
							 int me, int a, int ms,
							 bool iv,
							 bool iu, bool ip,
							 const char *file_name,
							 const char *tree_name)
	: NMaxEvts_(me), Algo_(a), NMaxSeeds_(ms),
	  IsVerbose_(iv),
	  IsUnbinned_(iu), IsPerPMT_(ip)
	  {
  //
  hPDF = GetROOTObj<TH2D>(pdf_name, hist_name)->ProjectionX("hPDF");
  std::cout << "Load PDF: " << hPDF->GetName() << std::endl;
  mPDF2D = GetROOTMObj<TH2D>(pdf_name, per_pmt_hist_name, "TH2D");
  std::transform(
	  mPDF2D.begin(), mPDF2D.end(),
	  std::inserter(mPDF1D, mPDF1D.begin()),
	  [](const std::pair<int, TH2D*>& p){
		return std::make_pair(p.first, p.second->ProjectionX());
	  }
  );
  //
  OFile = new TFile(file_name, "RECREATE");
  Tree = new TTree(tree_name, tree_name);
  ReconCoord.SetTree(Tree);
  //
  DetEdges = new CylEdges(R, HH, SpaceUnit::mm);
  //
  NMaxSeeds_ = NMaxSeeds_ < 0 ? std::numeric_limits<int>::max() : NMaxSeeds_;
		NMaxEvts_ = NMaxEvts_ < 0 ? std::numeric_limits<int>::max() : NMaxEvts_;
  //
  fNLL = IsUnbinned_ ? GetUNLL : GetNLL;
}
ReconAnalysis::~ReconAnalysis(){
  delete OFile;
  delete hPDF;
  for(auto& p : mPDF2D)
	delete p.second;
}

void ReconAnalysis::Do(void *Data) {
  //
  // Get Data
  auto wData = static_cast<TData*>(Data);
  std::vector<Hit> vHits = wData->GetVHits();
  if(vHits.empty())
	return;

  //
  // Get Seeds
  std::vector<Coord> vSeeds = GetSeeds(vHits, DetEdges);

  //
  // Prune Seeds
  std::vector<RecCoord> vPruneSeeds = PruneSeeds(vSeeds, vHits, DetEdges);

  //
  // Define fit parameters
  std::vector<void (*)(nlopt::opt &opt, CylEdges *c)> vfPars = {
	  SetBounds,
	  SetPars,
  };
  if(Algo_ == 2)
	vfPars.emplace_back(SetInequalityConstraint);

  if (IsPerPMT_) {
	FitMapStruct FMS = {vHits, mPDF1D, GetMUNLL};
	ReconCoord = Recon(&FMS, DetEdges, vSeeds, GetAlgo(Algo_), fPosTPerPMT, vfPars);
  } else {
	FitStruct FS = {vHits, hPDF, !IsUnbinned_, fNLL};
	ReconCoord = Recon(&FS, DetEdges, vSeeds, GetAlgo(Algo_), fPosT, vfPars);
  }

  //
  // Convert
  ReconCoord.ConvertTo(SpaceUnit::mm);
  Tree->Fill();

  //
  // Verbose
  if(IsVerbose_)
	std::cout << ReconCoord << std::endl;

  //
  // Get Event ID
  int iEvt = wData->GetEventID();
  //
  // Raise sigint if max number of events is reached
  if(iEvt >= NMaxEvts_)
	raise(SIGINT);

}

void ReconAnalysis::Export() const {
  OFile->cd();
  Tree->Write();
  OFile->Close();
}
