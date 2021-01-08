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
  const std::string input = args.filename;
  const std::string output = args.outname;


  // ######################################## //
  // Create wrapper object
  wRAT w_rat(input);
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
  // Handle gen shift of ANNIE
  auto ANNIEShift = [](const TVector3& v){
    TVector3 vShifted;
    vShifted.SetX(v.X());
    vShifted.SetY(-1*(v.Z()-1724));
    vShifted.SetZ(v.Y()+133.3);
    return vShifted;
  };
  auto ANNIEDirShift = [](const TVector3& v){
    TVector3 vShifted;
    vShifted.SetX(v.X());
    vShifted.SetY(-1*v.Z());
    vShifted.SetZ(v.Y());
    return vShifted;
  };
  

  // ######################################## //
  // Loop and get vector of NHits
  ProgressBar progress_bar(nEvts, 70);
  for(auto iEvt=0; iEvt<nEvts; iEvt++){

    // Record the tick
    ++progress_bar;

    // Point to evt
    w_rat.SetEvt(iEvt);

    // Get True info to build PDFs
    const auto PosTrue = ANNIEShift(w_rat.GetPosTrue(0));
    const auto DirTrue = ANNIEDirShift(w_rat.GetDirTrue(0));
    const auto TTrue = w_rat.GetTTrue(0);

    // Get number of trigger associated with an event
    // i.e, number of EV inside the rat DS
    auto nTriggers = w_rat.GetNTriggers();

    for(auto iTrigger=0; iTrigger<nTriggers; iTrigger++){

      // Get vector of hits
      std::vector<Hit> vHits = w_rat.GetVHits(iTrigger);

      // DO STUFF

      hNHits->Fill(GetNHits(vHits));

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
