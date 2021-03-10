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

#include <Wrapper.hh>
#include <Hit.hh>
#include <ProgressBar.hpp>

#include "CreatePDF.hh"
#include "MathUtils.hh"
#include "TriggerTimeMap.hh"
#include "Centroid.hh"

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

  const zAxis axTRes(150, -50., 100.);
  const zAxis axCosT(24, -1., 1.);
  const zAxis axNHits(1000, 0., 2000.);

  TH1D* hNHits = new TH1D("hNHits", "NHits per event ; NHits ; ",
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


  const std::vector<double> DetBnds = {args.bnds.x(), args.bnds.y(), args.bnds.z()};
  const std::vector<double> TBnds = {-10., args.bnds.Mag() / SOL};
  Bnds bnds = {DetBnds, TBnds};

  const int MaxRho = std::ceil(args.bnds.Perp()*1.e-3)*1.e3;
  const int MaxZ   = args.bnds.z();
  AxisGrid<int> agRho({0, MaxRho}, MaxRho/10);
  AxisGrid<int> agZ({0, MaxZ}, MaxZ/10);

  TrigTimePDF TTPDF(agRho.GetVCenters(), agZ.GetVCenters());

  auto hDWallVSTTime = new TH2D("hDWallVSTTime", "TRUE d_{Wall} vs T_{Trig} ; T_{Trig} [ns] ; d_{Wall} [mm]",
								20, 0., MaxZ / SOL,
								20, 0., MaxZ);

  auto hCentroidDWallVSTTime = new TH2D("hCentroidDWallVSTTime", "CENTROID d_{Wall} vs T_{Trig} ; T_{Trig} [ns] ; d_{Wall} [mm]",
										20, 0., MaxZ / SOL,
										20, 0., MaxZ);

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
	nTriggers = nTriggers > 1 ? 1 : nTriggers;

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

	  // Get vector of hits
	  std::vector<Hit> vHits = w_rat.GetVHits(iTrigger);
	  if(vHits.empty())
		continue;

	  std::sort(vHits.begin(), vHits.end());

	  //
	  // DO STUFF
	  //

	  // GetNHits per evt
	  hNHits->Fill(GetNPrompts(vHits, 400));

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


	  TTPDF.Fill(PosTrue, TrigTime);
	  hDWallVSTTime->Fill(TrigTime, GetDWall(PosTrue, 9000., 35000.));
	  hCentroidDWallVSTTime->Fill(TrigTime, GetDWall(GetCentroidSeed(vHits, bnds, 2), 9000., 35000.));

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

  ScaleHist(hDWallVSTTime, static_cast<double>(nEvts));
  hDWallVSTTime->Write();
  ScaleHist(hCentroidDWallVSTTime, static_cast<double>(nEvts));
  hCentroidDWallVSTTime->Write();

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

  for(auto& m : TTPDF.GetMT()){
	for(auto & mm : m.second){
	  mm.second->Write();
	}
  }

  TTPDF.Save();

  fOut.Close();

  return EXIT_SUCCESS;
}
