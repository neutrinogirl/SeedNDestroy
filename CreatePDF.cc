//
// Created by zsoldos on 1/5/21.
//

#include <EventDisplay.hh>

#include <iostream>
#include <numeric>

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
#include "TriggerTimeMap.hh"

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

  const int RoundedRho = static_cast<int>(std::round(b->GetTVector3().Perp()*1.e-3/2)*1.e3*2);
  const int RoundedZ   = static_cast<int>(std::round(b->GetTVector3().z()*1.e-3/2)*1.e3*2);
  AxisGrid<int> agRho({0, RoundedRho}, 1.e3);
  AxisGrid<int> agZ({0, RoundedZ}, 1.e3);

  TrigTimePDF TTPDF(agRho.GetVCenters(), agZ.GetVCenters());


  auto hDWallVSTTime = new TH2D("hDWallVSTTime", "TRUE d_{Wall} vs T_{Trig} ; T_{Trig} [ns] ; d_{Wall} [mm]",
				20, b->vT.min, b->vT.max,
				20, 0., MaxDWall);

  auto hRDWallVSTTime = new TH2D("hRDWallVSTTime", "TRUE d_{Wall} from R vs T_{Trig} ; T_{Trig} [ns] ; d_{Wall} [mm]",
				 20, b->vT.min, b->vT.max,
				 20, 0., MaxDWall);

  auto hZDWallVSTTime = new TH2D("hZDWallVSTTime", "TRUE d_{Wall} from Z vs T_{Trig} ; T_{Trig} [ns] ; d_{Wall} [mm]",
				 20, b->vT.min, b->vT.max,
				 20, 0., MaxDWall);

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

      double TrigTime = 0;

      // Get True info to build PDFs
      const auto PosTrue = ANNIEShift(w_rat.GetPosTrue(0));
      const auto DirTrue = ANNIEDirShift(w_rat.GetDirTrue(0));
      const auto TTrue = w_rat.GetTTrue(0);

      // Get vector of hits
      std::vector<Hit> vHits = w_rat.GetVHits(iTrigger);
      if(vHits.empty())
	continue;
      ReTriggerVHits(vHits, 2., TrigTime);
      // std::sort(vHits.begin(), vHits.end());

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

      std::size_t idx = 0;
      double dWall = GetDWall(PosTrue, 1516.6, 1973.8, idx);

      hDWallVSTTime->Fill(TrigTime, dWall);
      if(idx==0)
	hRDWallVSTTime->Fill(TrigTime, dWall);
      else if(idx==1)
	hZDWallVSTTime->Fill(TrigTime, dWall);


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

  //
  // ####
  //

  TF1 *fpol = new TF1("fpol", "pol1", -1, 1);
  fpol->SetLineWidth(2);
  fpol->SetParNames("B [ns]", "A [ns/mm]");

  auto hProf = hDWallVSTTime->ProfileY();
  hProf->Write();
  TFitResultPtr fr = hProf->Fit(fpol, "V");

  // Save fit in ttree
  typedef struct PolFitResults {
    double A, AErr, B, BErr, Chi2NDF;
  } PolFitResults;

  PolFitResults pfr;
  TTree T("pfr", "Linear fit results of d_Wall VS TTrig");
  T.Branch("A", &pfr.A, "A/D");
  T.Branch("AErr", &pfr.AErr, "AErr/D");
  T.Branch("B", &pfr.B, "B/D");
  T.Branch("BErr", &pfr.BErr, "BErr/D");
  T.Branch("Chi2NDF", &pfr.Chi2NDF, "Chi2NDF/D");
  pfr.Chi2NDF = fpol->GetChisquare() / fpol->GetNDF();
  pfr.B = fpol->GetParameter(0); pfr.BErr = fpol->GetParError(0);
  pfr.A = fpol->GetParameter(1); pfr.AErr = fpol->GetParError(1);
  T.Fill();
  T.Write();

  // /*Create a histogram to hold the confidence intervals*/
  // TH1D *hDWallVSTTime_Prof_Err = new TH1D("hDWallVSTTime_Prof_Err",
  // 						"TRUE d_{Wall} vs T_{Trig} ;  d_{Wall} [mm] ;",
  // 						hProf->GetNbinsX(), hProf->GetXaxis()->GetXmin(), hProf->GetXaxis()->GetXmax());
  // (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hDWallVSTTime_Prof_Err, 0.997);
  // hDWallVSTTime_Prof_Err->SetStats(kFALSE);
  // hDWallVSTTime_Prof_Err->SetFillColor(2);
  // hDWallVSTTime_Prof_Err->Write();


  hRDWallVSTTime->Write();
  hZDWallVSTTime->Write();

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

  theApp.Terminate();

  return EXIT_SUCCESS;
}
