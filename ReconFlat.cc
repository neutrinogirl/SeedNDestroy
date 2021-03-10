//
// Created by zsoldos on 1/5/21.
//

#include <iostream>
#include <string>
#include <vector>

#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>

#include <Wrapper.hh>
#include <ProgressBar.hpp>

#include "include/Recon.hh"
#include "Centroid.hh"
#include "Multilateration.hh"
#include "TriggerTimeMap.hh"

#include "DebugRecon.hh"
#include "Output.hh"
#include "RATDSCopier.hh"

int main(int argc, char *argv[]){

  // ######################################## //
  // Create TApp
  TApplication theApp("App", &argc, argv);


  // ######################################## //
  // Parse arguments
  // Simple struct containing filename and verbosity level
  Args args;
  ProcessArgs(theApp, args);
  const bool isVerbose = args.isVerbose;
  const std::string pdf = args.pdfname;
  const std::string output = args.outname;


  // ######################################## //
  // Load PDF
  auto hPDF = GetRootHisto<TH2D>(pdf.c_str(), Form("hCTVSTResPDF_TTOF_QW%d", args.wPower));
  TH1D *hPDF_TRes = hPDF->ProjectionX();
  if(args.isVerbose)
    std::cout << "PDF LOADED: " << hPDF_TRes->GetName() << std::endl;


  // ######################################## //
  // Create wrapper object
  wRAT w_rat(args.filename);
  const unsigned long int nEvts = args.nEvts > 0 ? args.nEvts : w_rat.GetNEvts();
  const unsigned int wPower = args.wPower;
  if(args.isVerbose)
	std::cout << "INPUT FILE LOADED: " << args.filename.front() << std::endl;


  // ######################################## //
  // Create structure holding data
  DataStruct1D ds = {hPDF_TRes, wPower};


  // ######################################## //
  // DET Boundaries
  const std::vector<double> DetBnds =
	  args.isBox ?
	  std::vector<double>({args.bnds.x(), args.bnds.y(), args.bnds.z()}) :
	  std::vector<double>({args.bnds.Perp(), args.bnds.Perp(), args.bnds.z()});
  const std::vector<double> TBnds =
	  args.isBox ?
	  std::vector<double>({0., args.bnds.Mag() / SOL}) :
	  std::vector<double>({0., sqrt(2*std::pow(args.bnds.Perp(), 2) + std::pow(args.bnds.z(), 2)) / SOL});
  Bnds bnds = {DetBnds, TBnds};
  if(args.isVerbose)
	std::cout << "DET BNDS SET: Rho=" << bnds.GetTVector3().Perp() << " z=" << bnds.GetTVector3().z() << std::endl;


  // ######################################## //
  // LOAD TTPDF
  TrigTimePDF TimePDF;
  TimePDF.Load(pdf);
  if(args.isVerbose)
	std::cout << "TimePDF LOADED" << std::endl;

  // Save fit in ttree
  typedef struct PolFitResults {
	double A, AErr, B, BErr, Chi2NDF;
  } PolFitResults;
  PolFitResults pfr;
  TFile file(pdf.c_str(), "READ");
  auto T = file.Get<TTree>("pfr");
  T->SetBranchAddress("A", &pfr.A);
  T->SetBranchAddress("AErr", &pfr.AErr);
  T->SetBranchAddress("B", &pfr.B);
  T->SetBranchAddress("BErr", &pfr.BErr);
  T->GetEntry(0);
  delete T;
  file.Close();
  std::cout << pfr.A << "[ns/mm]" << std::endl;

  // ######################################## //
  // Create structure holding boundaries
  DetParams dp = {args.bnds.Perp(), args.bnds.z(), pfr.A};

  const double MaxDWall = *std::min_element(DetBnds.begin(), DetBnds.end());
  TF1 fDWallVSTTrig("fDWallVSTTrig", "pol1", 0, MaxDWall);
  fDWallVSTTrig.SetParameter(0, pfr.B); fDWallVSTTrig.SetParError(0, pfr.BErr);
  fDWallVSTTrig.SetParameter(1, pfr.A); fDWallVSTTrig.SetParError(1, pfr.AErr);


  // ######################################## //
  // Output file
  TFile fOut(output.c_str(), "RECREATE");
  fOut.Close();
  // ######################################## //
  // OUTPUT Tree
  TTree tree("T", "Recon");
  Event evt;
  SetTTree(tree, evt);


  // ######################################## //
  // Loop and get vector of NHits
  ProgressBar progress_bar(nEvts, 70);
  for(auto iEvt=0; iEvt<nEvts; iEvt++){

	// Record the tick
	++progress_bar;

	// Point to evt
	w_rat.SetEvt(iEvt);

	// Get number of trigger associated with an event
	// i.e, number of EV inside the rat DS
	auto nTriggers = w_rat.GetNTriggers();

	for(auto iTrig=0; iTrig<nTriggers; iTrig++){

	  // Get ID tag
	  const std::string tag = Form("Evt%dTrig%d", iEvt, iTrig);
	  evt.iEvt = iEvt;
	  evt.iTrig = iTrig;

	  // Load MC Truth
	  LoadMCInfo2Evt(w_rat, evt);

	  // Get vector of hits
	  std::vector<Hit> vHits = w_rat.GetVHits(iTrig);
	  if(vHits.empty())
	    continue;
	  std::sort(vHits.begin(), vHits.end());

	  LoadVHits(evt, vHits);

	  //
	  // DO STUFF
	  //

	  // Prep Recon
	  ds.vHits.clear();
	  ds.vHits = vHits;


	  //
	  // #### #### #### TIME SEEDING #### #### #### //
	  //

	  auto CentroidSeed = GetCentroidSeed(vHits, bnds, 2);
	  if(!bnds.IsIn(CentroidSeed)){
		std::cerr << "NASTY Centroid seed" << std::endl;
		continue;
	  }

	  const double DWallSeed = GetDWall(CentroidSeed, args.bnds.Perp(), args.bnds.z());
	  const double TDWallSeed = fDWallVSTTrig.Eval(DWallSeed);

	  // const double TSeed = TimePDF.GetTrigTime(CentroidSeed, TBounds);
	  const double TSeed = TDWallSeed;

	  //
	  // #### #### #### POS SEEDING #### #### #### //
	  //

	  const std::size_t MaxSeeds = 5;
	  std::vector<TVector3> vSeeds = GetVSeeds(vHits, hPDF_TRes, -TSeed, bnds, wPower, MaxSeeds);

	  //
	  // #### #### #### Time to fit some sinsemilia #### #### #### //
	  //

	  std::vector< std::vector<double> > vX;
	  vX.reserve(vSeeds.size());

	  for(auto &Seed: vSeeds){

		// Prep Recon
		ds.Reset();

		// Recon
		// X = {XRec, YRec, ZRec, TRec, NLL, NLOPT::Results}
		auto x = ReconPosTime(ds, bnds, dp, Seed, -TSeed);
		vX.emplace_back(x);

	  }

	  //
	  // #### #### #### Sort fits by NLL #### #### #### //
	  //

	  std::sort(vX.begin(), vX.end(), [](const std::vector<double>& v1, const std::vector<double>& v2){
		return v1[4] < v2[4];
	  });

	  evt.RecPos = Vec(vX.front());
	  evt.RecT   = vX.front()[3];
	  evt.Chi2   = vX.front()[4];
	  evt.NLOPT  = vX.front()[5];

	  tree.Fill();

	  //
	  // ...
	  //

	}

	if(isVerbose)
	  progress_bar.display();

  }

  if(isVerbose)
	progress_bar.done();

  fOut.OpenFile(output.c_str(), "UPDATE");
  tree.Write();
  fOut.Close();

  return EXIT_SUCCESS;
}
