//
// Created by zsoldos on 1/5/21.
//

#include <EventDisplay.hh>

#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TProfile.h>
#include <TVirtualFitter.h>

#include <Wrapper.hh>
#include <Hit.hh>
#include <ProgressBar.hpp>

#include "CreatePDF.hh"
#include "MathUtils.hh"

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


  // ######################################## //
  // Make PDFs for different weight
  std::vector<unsigned int> vPower = {0, 1, 2};


  // ######################################## //
  // #### #### #### HISTOGRAMS #### #### #### //
  // ######################################## //

  const zAxis axTRes(250, -50., 200.);
  const zAxis axCosT(24, -1., 1.);
  const zAxis axNHits(1000, 0., 2000.);

  TH1D* hNHits = new TH1D("hNHits", "NHits per event ; NHits ; ",
						  axNHits.nBins, axNHits.min, axNHits.max);

  TH1D* hN400 = new TH1D("hN400", "N_{400} per event ; N_{400} ; ",
						 axNHits.nBins, axNHits.min, axNHits.max);

  enum { kTHIT, kTOF };
  std::vector< std::vector<TH2D*> > vvHPDFs;
  vvHPDFs.reserve(vPower.size());

  for(const auto& wP : vPower){
	vvHPDFs.push_back(
		{
			new TH2D(Form("hCTVSTResPDF_THit_QW%d", wP), "T_{Res} VS Cos(#theta) ; T_{Res} [ns] ; Cos(#theta)",
					 axTRes.nBins, axTRes.min, axTRes.max,
					 axCosT.nBins, axCosT.min, axCosT.max),
			new TH2D(Form("hCTVSTResPDF_TTOF_QW%d", wP), "T_{Res} VS Cos(#theta) ; T_{Res} [ns] ; Cos(#theta)",
					 axTRes.nBins, axTRes.min, axTRes.max,
					 axCosT.nBins, axCosT.min, axCosT.max)
		}
	);
  }

  // ######################################## //
  // DET Boundaries
  bnds *b;
  if(args.isBox){
	b = new BoxBnds( { args.bnds[0] , args.bnds[1] , args.bnds[2] });
  } else{
	b = new CylBnds(args.bnds[0], args.bnds[1]);
  }
  if(args.isVerbose)
	b->Print();

  const double MaxDWall = b->GetMaxDWall();
  const double ScaleDWall = 1.5;

  auto hDWallVSTTime = new TH2D("hDWallVSTTime", "TRUE d_{Wall} vs T_{Trig} ; T_{Trig} [ns] ; d_{Wall} [mm]",
								20, b->vT.min, b->vT.max,
								20*ScaleDWall, 0., MaxDWall*ScaleDWall);


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

	// Assume only single particle are generated
	const int iParticle = 0;

	for(auto iTrigger=0; iTrigger<nTriggers; iTrigger++){

	  const auto TrigTime = w_rat.GetTriggerTime(iTrigger);

	  // Get True info to build PDFs
	  const auto PosTrue = w_rat.GetPosTrue(iParticle);
	  const auto DirTrue = w_rat.GetDirTrue(iParticle);
	  const auto TTrue = w_rat.GetTTrue(iParticle);

	  // Get vector of hits
	  std::vector<Hit> vHits = w_rat.GetVHits(iTrigger);
	  if(vHits.empty())
		continue;
	  std::sort(vHits.begin(), vHits.end());

	  //
	  // DO STUFF
	  //

	  // GetNHits per evt
	  hN400->Fill(GetNPrompts(vHits, 400));
	  hNHits->Fill(GetNPrompts(vHits));

	  for(auto iPower = 0; iPower<vPower.size(); iPower++){

		// Get PDF
		for(auto& hit: vHits){

		  const double QW = fweight(hit, vPower[iPower]);

		  vvHPDFs[iPower][kTHIT]->Fill(hit.GetTRes(PosTrue, TTrue),
									   hit.GetCosTheta(PosTrue, DirTrue),
									   QW);
		  vvHPDFs[iPower][kTOF]->Fill(hit.GetTRes(PosTrue, TTrue-TrigTime),
									  hit.GetCosTheta(PosTrue, DirTrue),
									  QW);

		}

	  }


	  double dWall = b->GetDWall(PosTrue);
	  hDWallVSTTime->Fill(TrigTime, dWall);

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


  TFile fOut(output.c_str(), "RECREATE");

  ScaleHist(hNHits, static_cast<double>(nEvts));
  hNHits->Write();

  ScaleHist(hN400, static_cast<double>(nEvts));
  hN400->Write();

  ScaleHist(hDWallVSTTime, static_cast<double>(nEvts));
  hDWallVSTTime->Write();

  for(auto& vHPDFs:vvHPDFs){

	ScaleHist(vHPDFs[kTHIT], static_cast<double>(nEvts));
	vHPDFs[kTHIT]->Write();
	vHPDFs[kTHIT]->ProjectionX()->Write();
	vHPDFs[kTHIT]->ProjectionY()->Write();

	ScaleHist(vHPDFs[kTOF], static_cast<double>(nEvts));
	vHPDFs[kTOF]->Write();
	vHPDFs[kTOF]->ProjectionX()->Write();
	vHPDFs[kTOF]->ProjectionY()->Write();

  }

  fOut.Close();

  theApp.Terminate();

  return EXIT_SUCCESS;
}
