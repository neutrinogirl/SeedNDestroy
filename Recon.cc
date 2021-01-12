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
#include "include/Recon.hh"

#include "Recon.hh"
#include "DrawNLLSpace.hh"

int main(int argc, char *argv[]){

  // ######################################## //
  // Create TApp
  TApplication theApp("App", &argc, argv);
  gROOT->SetBatch(true);


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
  auto hPDF = GetRootHisto<TH2D>(pdf.c_str(), "hCTVSTResPDF");
  TH1D *hPDF_TRes = hPDF->ProjectionX();


  // ######################################## //
  // Create wrapper object
  wRAT w_rat(args.filename);
  const unsigned long int nEvts = args.nEvts > 0 ? args.nEvts : w_rat.GetNEvts();
  const unsigned int wPower = args.wPower;


  // ######################################## //
  // Create structure holding data
  DataStruct1D ds = {hPDF_TRes, std::vector<Hit>(), wPower};
  TTree tree("T", "A flat tree");
  Event evt;
  SetTTree(tree, evt);


  // WM Boundaries
  Bnds WMBnds = {8.e3/*mm*/, 10./*ns*/};


  // ######################################## //
  // Loop and get vector of NHits
  ProgressBar progress_bar(nEvts, 70);
  signal(SIGINT, signal_handler);
  for(auto iEvt=0; iEvt<nEvts; iEvt++){

    // if ctrl+c exit loop
	if(gSignalStatus > 0)
	  break;

	// Record the tick
	++progress_bar;

	// Point to evt
	w_rat.SetEvt(iEvt);

	// Get number of trigger associated with an event
	// i.e, number of EV inside the rat DS
	auto nTriggers = w_rat.GetNTriggers();

	for(auto iTrigger=0; iTrigger<nTriggers; iTrigger++){

	  // IF USE SPLITEVDAQ
	  // Get EV TrigTime
	  const auto TrigTime = w_rat.GetTriggerTime(iTrigger);

	  // Try to guess which particle is attached to this trigger
	  // Make sense only for IBD gen, when n capture will be >> in T that prompt event
	  // Note that also e+ is always first particle generated
	  auto nParticle = w_rat.GetNPrimaryParticle();
	  auto iParticle = nParticle > 1 ? (TrigTime > 1e3 ? 1 : 0) : 0;

	  // Get True info to record
	  const auto PosTrue = w_rat.GetPosTrue(iParticle);
	  const auto DirTrue = w_rat.GetDirTrue(iParticle);
	  const auto TTrue = w_rat.GetTTrue(iParticle);

	  evt.MCPos = Vec(PosTrue);
	  evt.MCDir = Vec(DirTrue);
	  evt.MCT = TTrue;

	  const double TCor = TTrue - TrigTime;

	  // Get vector of hits
	  std::vector<Hit> vHits = w_rat.GetVHits(iTrigger);

	  // DO STUFF
	  auto PosTSeed = GetSeed(vHits, hPDF_TRes, 0., WMBnds.Pos, wPower);

	  ds.vHits.clear();
	  ds.vHits = vHits;

	  auto x = ReconPosTime(ds, WMBnds, PosTSeed.Pos, PosTSeed.T);

	  evt.RecPos = Vec(x);
	  evt.RecT = x[3]*1.e-2;

	  tree.Fill();

	  // ...

	}

	if(isVerbose)
	  progress_bar.display();

  }

  if(isVerbose)
	progress_bar.done();

  TFile fOut(output.c_str(), "RECREATE");
  tree.Write();
  fOut.Close();

  return EXIT_SUCCESS;
}