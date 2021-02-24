//
// Created by zsoldos on 1/5/21.
//

#include <iostream>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>

#include <Wrapper.hh>
#include <ProgressBar.hpp>

#include <EventDisplay.hh>

#include "include/Recon.hh"
#include "Centroid.hh"
#include "Multilateration.hh"
#include "TriggerTimeMap.hh"

#include "DebugRecon.hh"
#include "FitPerfMonitor.hh"

#include "Output.hh"

int main(int argc, char *argv[]){

  // ######################################## //
  // Create TApp
  TApplication theApp("App", &argc, argv);
  gROOT->SetBatch(true);
  gStyle->SetPalette(kInvertedDarkBodyRadiator);


  // ######################################## //
  // Parse arguments
  // Simple struct containing filename and verbosity level
  Args args;
  ProcessArgs(theApp, args);
  const bool isVerbose = args.isVerbose;
  const std::string pdf = args.pdfname;
  const std::string output = args.outname;


  // ######################################## //
  // Output file
  TFile fOut(output.c_str(), "RECREATE");
  fOut.Close();


  // ######################################## //
  // Load PDF
  auto hPDF = GetRootHisto<TH2D>(pdf.c_str(), "hCTVSTResPDF_TTOF");
  TH1D *hPDF_TRes = hPDF->ProjectionX();
  SetBasicTStyle(hPDF_TRes, kBlack, 2, kDashed);


  // ######################################## //
  // Create wrapper object
  wRAT w_rat(args.filename);
  const unsigned long int nEvts = args.nEvts > 0 ? args.nEvts : w_rat.GetNEvts();
  const unsigned int wPower = args.wPower;


  // ######################################## //
  // Create structure holding data
  DataStruct1D ds = {hPDF_TRes, wPower};


  // ######################################## //
  // DET Boundaries
  const std::vector<double> DetBnds = {args.bnds.x(), args.bnds.y(), args.bnds.z()};
  const std::vector<double> TBnds = {-10., args.bnds.Mag() / SOL};
  Bnds bnds = {DetBnds, TBnds};


  // ######################################## //
  // LOAD TTPDF
  TrigTimePDF TimePDF;
  TimePDF.Load(pdf);


  // ######################################## //
  // FitterDebbuger and co
  std::vector<FitPerfMonitor> v_centroid_per_monitor = {
	  FitPerfMonitor("centroid_NoWeight",bnds.GetMin(), bnds.GetMin()/100.),
	  FitPerfMonitor("centroid_Q",bnds.GetMin(), bnds.GetMin()/100.),
	  FitPerfMonitor("centroid_Q2",bnds.GetMin(), bnds.GetMin()/100.),
	  FitPerfMonitor("centroid_Q3",bnds.GetMin(), bnds.GetMin()/100.),
	  FitPerfMonitor("centroid_Q4",bnds.GetMin(), bnds.GetMin()/100.)
  };
  FitPerfMonitor recon_perf_monitor("rec", bnds.GetMin(), bnds.GetMin()/100.);

  MapNLL map_nll(DetBnds, {args.bnds.x()/10, args.bnds.y()/10, args.bnds.z()/10});

  PMTGrid pmt_grid = GetPMTGrid(w_rat, CreateVBnds(args.bnds.Perp(), args.bnds.z()));


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
	// nTriggers = nTriggers > 1 ? 1 : nTriggers;

	for(auto iTrigger=0; iTrigger<nTriggers; iTrigger++){

	  // Get ID tag
	  const std::string tag = Form("Evt%dTrig%d", iEvt, iTrigger);
	  evt.iEvt = iEvt;
	  evt.iTrig = iTrigger;

	  // IF USE SPLITEVDAQ
	  // Get EV TrigTime
	  const auto TrigTime = w_rat.GetTriggerTime(iTrigger);

	  // Try to guess which particle is attached to this trigger
	  // Make sense only for IBD gen, when n capture will be >> in T that prompt event
	  // Note that also e+ is always first particle generated
	  auto nParticle = w_rat.GetNPrimaryParticle();
	  auto iParticle = nParticle > 1 ? (TrigTime > 1e3 ? 1 : 0) : 0;
	  // Skip if it is not prompt
	  // if(iParticle>0)
		// continue;

	  // Get True info to record
	  const auto PosTrue = w_rat.GetPosTrue(iParticle);
	  const auto DirTrue = w_rat.GetDirTrue(iParticle);

	  evt.MCPos = Vec(PosTrue);
	  evt.MCDir = Vec(DirTrue);
	  evt.MCT   = -TrigTime;
	  evt.ETrue = w_rat.GetETrue(iParticle);

	  const std::vector<double> xTrue = {PosTrue[0], PosTrue[1], PosTrue[2], -TrigTime};

	  // Get vector of hits
	  std::vector<Hit> vHits = w_rat.GetVHits(iTrigger);
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

	  // Get Seed
	  // Centroid
	  auto CentroidSeed = GetCentroidSeed(vHits, bnds, 4);
	  if(!bnds.IsIn(CentroidSeed))
		std::cerr << "NASTY Centroid seed" << std::endl;

	  // Record Centroid Seed
	  for(auto i=0; i<v_centroid_per_monitor.size(); i++)
		v_centroid_per_monitor[i].Fill(GetCentroidSeed(vHits, bnds, i), PosTrue);

	  // SEED Boundaries
	  Bnds localbnds = bnds;

	  std::vector<double> prout(2);
	  const double TSeed = TimePDF.GetTrigTime(CentroidSeed, prout);

	  std::vector<TVector3> vSeeds = GetVSeeds(vHits, hPDF_TRes, -TSeed, bnds, wPower, 5);

	  std::vector< std::vector<double> > vX;
	  vX.reserve(vSeeds.size());

	  for(auto &Seed: vSeeds){

		// Prep Recon
		ds.Reset();

		// Recon
		auto x = ReconPosTime(ds, localbnds, Seed, -TSeed);
		vX.emplace_back(x);

	  }

	  std::sort(vX.begin(), vX.end(), [](const std::vector<double>& v1, const std::vector<double>& v2){
		return v1[4] < v2[4];
	  });

	  // Recon
	  const TVector3 PosRec(vX.front()[0], vX.front()[1], vX.front()[2]);
	  recon_perf_monitor.Fill(PosRec, PosTrue);

	  evt.RecPos = Vec(vX.front());
	  evt.RecT   = vX.front()[3];
	  evt.Chi2   = vX.front()[4];
	  evt.NLOPT  = vX.front()[5];

	  tree.Fill();

	  // ###################################### //
	  // #### #### #### PLOTTING #### #### #### //
	  // ###################################### //

	  if(args.isDebug){

		//
		// EV DISPLAY
		//

		auto m_hits = w_rat.GetMHits(iTrigger);
		Save2ROOT( GetEventDisplay(w_rat.GetMHits(iTrigger), pmt_grid, tag, 4.0),
				   output );

		//
		// TRes Hist
		//

		// SEED
		// For debugging purposes
		const std::vector<double> xSeed = {CentroidSeed[0], CentroidSeed[1], CentroidSeed[2], -TSeed};
		Save2ROOT( GetTResHist(tag+"_Seed", vHits,
							   xTrue, xSeed,
							   wPower, hPDF_TRes), output);

		// REC
		Save2ROOT( GetTResHist(tag, vHits,
							   xTrue, vX.front(),
							   wPower, hPDF_TRes) , output);

		//
		// MAP NLL SPACE
		//

		// TSeed
		map_nll.Fill(vHits, -TSeed, hPDF_TRes, 1);
		Save2ROOT( GetMapPlots(map_nll.GetHGrid(), {PosTrue, PosRec, vSeeds.front()} ,
							   Form("cGridSeed_%s", tag.c_str())),
				   output );
		// HGrid GARBAGE COLLECTOR
		HGridGarbageCollector();
		// RESET
		map_nll.ResetGrid();

		// TTrue
		map_nll.Fill(vHits, -TrigTime, hPDF_TRes, 1);
		Save2ROOT( GetMapPlots(map_nll.GetHGrid(), {PosTrue, PosRec, vSeeds.front()} ,
							   Form("cGridTrue_%s", tag.c_str())),
				   output );
		// HGrid GARBAGE COLLECTOR
		HGridGarbageCollector();
		// RESET
		map_nll.ResetGrid();

		//
		// TCANVAS GARBAGE COLLECTOR
		//

		for (auto &&obj: *gROOT->GetListOfCanvases()) {
		  delete obj;
		}

	  }

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

  for(auto& fpm: v_centroid_per_monitor){
	auto cCentroid = fpm.GetPlot(true);
	cCentroid->Write();
	auto cDRho = fpm.GetDRhoPlot();
	cDRho->Write();
  }

  auto cFit = recon_perf_monitor.GetPlot(true);
  cFit->Write();
  auto cDRho = recon_perf_monitor.GetDRhoPlot();
  cDRho->Write();

  fOut.Close();

  return EXIT_SUCCESS;
}