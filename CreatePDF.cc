//
// Created by zsoldos on 1/5/21.
//

#include <TFile.h>

#include <Wrapper.hh>
#include <Hit.hh>
#include <ProgressBar.hpp>

#include "CreatePDF.hh"
#include "MathUtils.hh"

int main(int argc, char *argv[]){

  // ######################################## //
  // Parse arguments
  // Simple struct containing filename and verbosity level
  CreatePDFArgs args;
  args(argc, argv);


  // ######################################## //
  // Create wrapper object
  wRAT w_rat(args.GetInput());
  const unsigned long nEvts = args.GetNEvts() > 0 ? args.GetNEvts() : w_rat.GetNEvts();


  // ######################################## //
  // Make PDFs for different weight
  std::vector<unsigned int> vPower = {0, 1, 2};


  // ######################################## //
  // #### #### #### HISTOGRAMS #### #### #### //
  // ######################################## //

  const zAxis axTRes(args.GetTResBins()[0], args.GetTResBins()[1], args.GetTResBins()[2]);
  const zAxis axCosT(12, -1., 1.);
  const zAxis axNHits(1000, 0., 1000.);

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
  bnds *b = new CylBnds(args.GetRadius(), args.GetHHeight());
  if(args.GetVerbose())
	b->Print();


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

	  // ...

	}

	if(args.GetVerbose())
	  progress_bar.display();

  }

  if(args.GetVerbose())
	progress_bar.done();

  // #################################### //
  // #### #### #### FINISH #### #### #### //
  // #################################### //


  TFile fOut(args.GetOutput(), "RECREATE");

  ScaleHist(hNHits, static_cast<double>(nEvts));
  hNHits->Write();

  ScaleHist(hN400, static_cast<double>(nEvts));
  hN400->Write();

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

  return EXIT_SUCCESS;
}
