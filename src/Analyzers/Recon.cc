//
// Created by Stephane Zsoldos on 7/6/22.
//

#include <csignal>
#include <iomanip>
#include <numeric>
#include <algorithm>

#include "Recon.hh"

#include "Templates/TData.hh"
#include "Algo/WOpt.hh"
#include "Algo/VHits.hh"
#include "ROOT/Utils.hh"

ReconAnalysis::ReconAnalysis(const char *pdfname, const char *histname, const char* perpmthistname ,
							 const double &R, const double &HH,
							 int me, int a, int ms,
							 bool iv,
							 bool iu, bool ip,
							 bool iat,
							 bool ijs,
							 const char *filename,
							 const char *treename)
	: nMaxEvts(me), algo(a), max_seed(ms),
	  isverbose(iv),
	  isunbinned(iu), isperpmt(ip),
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
  ReconCoord.SetTree(Tree);
  //
  if(istrack){
	Tree->Branch("IterX", &vIterX);
	Tree->Branch("IterY", &vIterY);
	Tree->Branch("IterZ", &vIterZ);
	Tree->Branch("IterT", &vIterT);
	Tree->Branch("Iterf", &vf);
	Tree->Branch("nIter", &nIter, "nIter/I");
  }
  //
  max_seed = max_seed < 0 ? std::numeric_limits<int>::max() : max_seed;
  nMaxEvts = nMaxEvts < 0 ? std::numeric_limits<int>::max() : nMaxEvts;
  //
  fNLL = isunbinned ? GetUNLL : GetNLL;
}
ReconAnalysis::~ReconAnalysis(){
  delete OFile;
  delete Cyl;
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
  // Apply trigger if necessary (files with True hit time for each PMTs)
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
  // Get Event ID
  int iEvt = wData->GetEventID();
  //
  // Interrupt if user-defined max number of events reached
  if(iEvt-1 == nMaxEvts)
	raise(SIGINT);

  //
  //
  std::vector<Coord> vSeeds = GetSeeds(vHits, Cyl);

  // Sort seeds from lowest to highest chi2
  std::sort(vSeeds.begin(), vSeeds.end(),
			[&](const PosT& a, const PosT& b){
			  return
				  fNLL(*hPDF, ConvertTVector3Unit<double>(a.GetTVector3(), SpaceUnit::dm, SpaceUnit::mm), a.T, vHits) <
				  fNLL(*hPDF, ConvertTVector3Unit<double>(b.GetTVector3(), SpaceUnit::dm, SpaceUnit::mm), b.T, vHits);
			}
  );

  //
  if(isverbose){
	auto PosSeed = ConvertTVector3Unit<double>(vSeeds.front().GetTVector3(), SpaceUnit::dm, SpaceUnit::mm);
	PosSeed.Print();
	std::cout << fNLL(*hPDF, PosSeed, vSeeds.front().T, vHits) << std::endl;
  }

  //
  auto ConvertAndFill = [&](const TVector3& p){
	TVector3 v3mm = ConvertTVector3Unit<double>(p, SpaceUnit::dm, SpaceUnit::mm);
	RT.X = v3mm.x();
	RT.Y = v3mm.y();
	RT.Z = v3mm.z();
	Tree->Fill();
  };

  //
  if(isjustseed){
	ConvertAndFill(vSeeds.back().GetTVector3());
	return;
  }

  //
  std::vector<void (*)(nlopt::opt &opt, Bnd *c)> vfPars = {
	  SetBounds,
	  SetPars,
  };
  if(algo == 2)
	vfPars.emplace_back(SetInequalityConstraint);

  // Recon
  FitStruct FS = {vHits, hPDF, !isunbinned, fNLL};
  FS.filldata = istrack;
  RT = Recon(&FS, Cyl, {vSeeds.front()}, GetAlgo(algo), Walk, {SetBounds, SetPars});
  if(isperpmt) {
	FitMapStruct FMS = {vHits, mPDF1D, GetMUNLL};
	RT = Recon(&FMS, Cyl, {RT}, GetAlgo(algo), fPosTPerPMT, vfPars);
  }
  // Fill
  // FS.FillSliceIterateData(&vIterX, &vIterY, &vIterZ, &vIterT, &vf, &nIter);
  // Convert back to mm
  ConvertAndFill(RT.GetTVector3());

  // Verbose
  if(isverbose)
	RT.Print();

}

void ReconAnalysis::Export() const {
  OFile->cd();
  Tree->Write();
  OFile->Close();
}

void Debug(void* Data, TH1D* hPDF, const std::map<int, TH1D*>& mPDF1D){
  // Get Data
  auto wData = static_cast<TData*>(Data);
  std::vector<Hit> vHits = wData->GetVHits();
  const int nHits = static_cast<int>(vHits.size());
  // Get true position of the event
  auto Pos = wData->GetPosition();
  auto T = wData->GetTime();
  // Get NLL event
  double TrueNLL   = GetNLL(*hPDF, Pos, T, vHits) / hPDF->GetNbinsX();
  double TrueUNLL  = GetUNLL(*hPDF, Pos, T, vHits) / nHits;
  double TrueMUNLL = GetMUNLL(mPDF1D, Pos, T, vHits) / nHits;
  // Sort all hits from closest to farthest from true position
  SortHitsFromPos(vHits, Pos);
  // Print print print
  std::cout << "EventID: " << wData->GetEventID() << " "
			<< "vHits.size(): " << vHits.size() << std::endl;
  std::cout << std::endl;
  // Print NLLs
  std::cout << std::setprecision(6) << "TrueNLL:" << TrueNLL << " "
			<< std::setprecision(6) << "TrueUNLL:" << TrueUNLL << " "
			<< std::setprecision(6) << "TrueMUNLL:" << TrueMUNLL << std::endl;
  std::cout << std::endl;
  // Print all hits using the sorted order, the distance from true position and hit.Print()
  for(const auto& h : vHits){
	std::cout << std::setprecision(4) << h.T << "ns "
			  << std::setprecision(4) << h.Q << "Q "
			  << std::round(h.GetD(Pos)) << "mm " << std::endl;
  }

  // Slice vHits
  std::vector<double> vT = GetTs(vHits), vdT(vT.size()), vD = GetDs(vHits, Pos), vdD(vD.size());
  // Differentiate to get dT and dD
  std::adjacent_difference(vT.begin(), vT.end(), vdT.begin());
  vdT[0] = 0;
  std::adjacent_difference(vD.begin(), vD.end(), vdD.begin());
  vdD[0] = 0;

  //
  std::cout << std::endl;
}