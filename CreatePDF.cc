//
// Created by zsoldos on 1/5/21.
//

#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

#include "CreatePDF.hh"
#include <Wrapper.hh>
#include <Hit.hh>
#include <ProgressBar.hpp>

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
  const std::string output = args.outname;


  // ######################################## //
  // Create wrapper object
  wRAT w_rat(args.filename);
  const unsigned long int nEvts = args.nEvts > 0 ? args.nEvts : w_rat.GetNEvts();
  const unsigned int wPower = args.wPower;


  // ######################################## //
  // #### #### #### HISTOGRAMS #### #### #### //
  // ######################################## //

  TH1D* hNHits = new TH1D("hNHits", "NHits per event ; NHits ; ",
						  1000, 0., 1000.);
  hNHits->SetDirectory(nullptr);

  TH2D* hTResVSCosT = new TH2D("hCTVSTResPDF", "T_{Res} VS Cos(#theta) ; T_{Res} [ns] ; Cos(#theta)",
									 600, -200., 400.,
									 24, -1., 1.);
  hTResVSCosT->SetDirectory(nullptr);


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

	for(auto iTrigger=0; iTrigger<nTriggers; iTrigger++){

	  // IF USE SPLITEVDAQ
	  // Get EV TrigTime
	  const auto TrigTime = w_rat.GetTriggerTime(iTrigger);

	  // Try to guess which particle is attached to this trigger
	  // Make sense only for IBD gen, when n capture will be >> in T that prompt event
	  // Note that also e+ is always first particle generated
	  auto nParticle = w_rat.GetNPrimaryParticle();
	  auto iParticle = nParticle > 1 ? (TrigTime > 1e3 ? 1 : 0) : 0;
	  // Skip if it is not prompt
	  if(iParticle>0)
		continue;

	  // Get True info to build PDFs
	  const auto PosTrue = w_rat.GetPosTrue(iParticle);
	  const auto DirTrue = w_rat.GetDirTrue(iParticle);
	  const auto TTrue = w_rat.GetTTrue(iParticle);

	  const double TCor = TTrue - TrigTime;

	  // Get vector of hits
	  std::vector<Hit> vHits = w_rat.GetVHits(iTrigger);

	  // DO STUFF

	  hNHits->Fill(GetNHits(vHits));
	  double AvgHits = 0;
	  unsigned int nPrompts = 0;
	  const double PromptCut = 4.0; //ns
	  auto IsPrompt = [&PromptCut](const Hit& hit){return hit.T < PromptCut;};

	  for(auto& hit: vHits){

		hTResVSCosT->Fill(hit.GetTRes(PosTrue, TTrue),
								hit.GetCosTheta(PosTrue, DirTrue),
								fweight(hit, wPower));
	  }

	  // ...

	}

	if(isVerbose)
	  progress_bar.display();

  }

  if(isVerbose)
	progress_bar.done();

  // #################################### //
  // #### #### #### FINISH #### #### #### //
  // #################################### //

  hTResVSCosT->Scale(1./static_cast<double>(nEvts));

  TFile fOut(output.c_str(), "RECREATE");
  hNHits->Write();
  hTResVSCosT->Write();
  hTResVSCosT->ProjectionX()->Write();
  hTResVSCosT->ProjectionY()->Write();
  fOut.Close();

  return EXIT_SUCCESS;
}
