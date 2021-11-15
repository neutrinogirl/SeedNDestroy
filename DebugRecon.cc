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

#include "Recon.hh"
#include "Centroid.hh"
#include "Multilateration.hh"

#include "DebugRecon.hh"
#include "FitPerfMonitor.hh"

#include "Output.hh"

int main(int argc, char *argv[]){

  // ######################################## //
  // Get ROOT Version
  const double ROOTVERSION = stod(std::string(gROOT->GetVersion()).substr(0, 4));
  std::cout << "USING ROOT V:" << ROOTVERSION << std::endl;


  // ######################################## //
  // Create TApp
  TApplication theApp("App", &argc, argv);
  gROOT->SetBatch(true);
  if(ROOTVERSION>6.)
    gStyle->SetPalette(kInvertedDarkBodyRadiator);


  // ######################################## //
  // Parse arguments
  // Simple struct containing filename and verbosity level
  Args args;
  ProcessArgs(theApp, args);
  const bool isVerbose = args.isVerbose;
  const std::string pdf = args.pdfname;
  const std::string output = args.outname;
  const unsigned int wPower = args.wPower;


  // ######################################## //
  // Load PDF
  auto hPDF = GetRootHisto<TH2D>(pdf.c_str(), Form("hCTVSTResPDF_TTOF_QW%d", wPower));
  TH1D *hPDF_TRes = hPDF->ProjectionX();
  if(args.isVerbose)
    std::cout << "PDF LOADED: " << hPDF_TRes->GetName() << std::endl;


  // ######################################## //
  // Create wrapper object
  wRAT w_rat(args.filename);
  const unsigned long int nEvts = args.nEvts > 0 ? args.nEvts : w_rat.GetNEvts();
  if(args.isVerbose)
    std::cout << "INPUT FILE LOADED: " << args.filename.front() << std::endl;


  // ######################################## //
  // Create structure holding data
  DataStruct1D ds = {hPDF_TRes, wPower};


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


  // ######################################## //
  // Create structure holding boundaries
  DetParams dp = {b, GetSOL()};
  const double MaxDWall = b->GetMaxDWall();
  if(args.isVerbose)
    std::cout << MaxDWall << std::endl;


  // ######################################## //
  // FitterDebbuger and co
  std::vector<FitPerfMonitor> v_centroid_per_monitor = {
    FitPerfMonitor("centroid_NoWeight",MaxDWall, MaxDWall/50.),
    FitPerfMonitor("centroid_Q",MaxDWall, MaxDWall/50.),
    FitPerfMonitor("centroid_Q2",MaxDWall, MaxDWall/50.),
    FitPerfMonitor("centroid_Q3",MaxDWall, MaxDWall/50.),
    FitPerfMonitor("centroid_Q4",MaxDWall, MaxDWall/50.)
  };
  FitPerfMonitor gridTDWall_seed_perf_monitor("gridTDWall_seed", MaxDWall, MaxDWall/50.);
  FitPerfMonitor seed_perf_monitor("seed", MaxDWall, MaxDWall/50.);
  FitPerfMonitor recon_perf_monitor("rec", MaxDWall, MaxDWall/50.);

  std::vector<double> DetBnds = b->GetVDetBnds();
  MapNLL map_nll(DetBnds, {DetBnds[0]/10, DetBnds[1]/10, DetBnds[2]/10});

  PMTGrid pmt_grid = GetPMTGrid(w_rat, b->GetTVector3());


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
    // nTriggers = nTriggers > 1 ? 1 : nTriggers;


    for(auto iTrig=0; iTrig<nTriggers; iTrig++){

      // Get ID tag
      const std::string tag = Form("Evt%dTrig%d", iEvt, iTrig);
      evt.iEvt = iEvt;
      evt.iTrig = iTrig;

      // DEBUG PRINTS
      if(args.isDDebug){
	std::cout << std::endl;
	std::cout << tag << std::endl;
	std::cout << std::endl;
      }

      // Get EV TrigTime
      const auto TrigTime = w_rat.GetTriggerTime(iTrig);
      // Get vector of hits
      std::vector<Hit> vHits = w_rat.GetVHits(iTrig);
      if(args.cvg > 0)
	SlimVHits(vHits, args.cvg);
      if(vHits.empty()){
	std::cout << "Trigger with no hits !?" << std::endl;
	continue;
      }

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
      evt.MCT   = TTrue-TrigTime;
      evt.ETrue = w_rat.GetETrue(iParticle);

      const std::vector<double> xTrue = {PosTrue[0], PosTrue[1], PosTrue[2], TTrue-TrigTime};

      // DEBUG PRINTS
      if(args.isDDebug){
	PrintDD("TRUTH", PosTrue, -(TTrue-TrigTime), b->GetDWall(PosTrue));
      }

      LoadVHits(evt, vHits);

      //
      // DO STUFF
      //

      // Prep Recon
      ds.vHits.clear();
      ds.vHits = vHits;

      // Record Centroid Seed
      for(auto i=0; i<v_centroid_per_monitor.size(); i++)
	v_centroid_per_monitor[i].Fill(GetCentroidSeed(vHits, *b, i), PosTrue);

      //
      // #### #### #### TIME SEEDING #### #### #### //
      //

      auto CentroidSeed = GetCentroidSeed(vHits, *b, 0);
      if(!b->IsInPos(CentroidSeed)){
	std::cerr << "NASTY Centroid seed" << std::endl;
	continue;
      }
      double DWallSeed = b->GetDWall(CentroidSeed);
      double TDWallSeed = DWallSeed/GetSOL();

      auto TDwallSeedBnds = GetTBndsLocal(TDWallSeed, {b->vT.min, b->vT.max});
      if(args.isDDebug) {
	std::cout << "[" << TDwallSeedBnds[0] << "," << TDwallSeedBnds[1] << "]: " << TDWallSeed << std::endl;
      }

      // DEBUG PRINTS
      if(args.isDDebug){
	PrintDD("CENTROID", CentroidSeed, TDWallSeed, DWallSeed,
		GetNLL(vHits, hPDF_TRes, CentroidSeed, -TDWallSeed, fweight, wPower, args.isUnbinned));
      }

      //
      // #### #### #### POS SEEDING #### #### #### //
      //

      const auto MaxSeeds = std::numeric_limits<unsigned int>::max();
      std::vector<PosT> vSeeds;
      // = GetVPosTSeeds(vHits, hPDF_TRes, *b, wPower, MaxSeeds);
      // if(!vSeeds.empty())
      // 	seed_perf_monitor.Fill(vSeeds.front().Pos, PosTrue);
      // else
      // 	std::cout << "Seeding algo empty" << std::endl;

      vSeeds.emplace_back(CentroidSeed, TDWallSeed);

      //
      // #### GRID SEARCH
      //

      if(args.useGridSearch){

	map_nll.Fill(vHits, -TDWallSeed, hPDF_TRes, wPower);
	vSeeds.emplace_back(map_nll.GetVMinNll(), TDWallSeed);
	gridTDWall_seed_perf_monitor.Fill(map_nll.GetVMinNll(), PosTrue);
	// RESET
	map_nll.ResetGrid();

	// DEBUG PRINTS
	if(args.isDDebug){
	  PrintDD("GRIDSEED", map_nll.GetVMinNll(), TDWallSeed, b->GetDWall(map_nll.GetVMinNll()),
		  GetNLL(vHits, hPDF_TRes, map_nll.GetVMinNll(), -TDWallSeed, fweight, wPower, args.isUnbinned));
	}

      }


      std::vector< std::vector<double> > vX;
      vX.reserve(vSeeds.size());

      for(auto &Seed: vSeeds){

	// DEBUG PRINTS
	if(args.isDDebug){
	  PrintDD("SEEDS", Seed.Pos, Seed.T, b->GetDWall(Seed.Pos),
		  GetNLL(vHits, hPDF_TRes, Seed.Pos, -Seed.T, fweight, wPower, args.isUnbinned));
	}

	// Prep Recon
	ds.Reset();

	// GetLocalBnds
	std::vector<double> vT = GetTBndsLocal(Seed.T, {b->vT.min, b->vT.max});
	bnds *localb;
	if(args.isBox)
	  localb = new BoxBnds(*b);
	else
	  localb = new CylBnds(*b);

	localb->vT.min = vT[0];
	localb->vT.max = vT[1];

	double NLLSeed = GetNLL(vHits, hPDF_TRes, Seed.Pos, -Seed.T, fweight, wPower, args.isUnbinned);

	// Recon
	// X = {XRec, YRec, ZRec, TRec, NLL, NLOPT::Results}
	auto x = ReconPosTime(ds, *localb, dp, Seed.Pos, -Seed.T, NLLSeed);

	// DEBUG PRINTS
	if(args.isDDebug){
	  PrintDD("REC", TVector3(x[0], x[1], x[2]), -x[3], b->GetDWall(TVector3(x[0], x[1], x[2])),
		  x[4]);
	  std::cout << "NLOPT RETURN: " << x[5] << std::endl;
	}

	delete localb;

	if(x[4] < NLLSeed)
	  vX.emplace_back(x);
	else
	  vX.emplace_back(std::vector<double>({Seed.Pos[0], Seed.Pos[1], Seed.Pos[2], -Seed.T, NLLSeed, 666.}));

      }

      std::sort(vX.begin(), vX.end(), [](const std::vector<double>& v1, const std::vector<double>& v2){
	  return v1[4] < v2[4];
	});

      if(!vX.empty()){

	// Recon
	const TVector3 PosRec(vX.front()[0], vX.front()[1], vX.front()[2]);
	recon_perf_monitor.Fill(PosRec, PosTrue);

	evt.RecPos = Vec(vX.front());
	evt.RecT   = vX.front()[3];
	evt.Chi2   = vX.front()[4];
	evt.NLOPT  = vX.front()[5];
	evt.dWall  = b->GetDWall(PosRec);

	tree.Fill();


	// ###################################### //
	// #### #### #### PLOTTING #### #### #### //
	// ###################################### //

	if(args.isDebug){

	  //
	  // EV DISPLAY
	  //

	  auto m_hits = w_rat.GetMHits(iTrig);
	  Save2ROOT( GetEventDisplay(w_rat.GetMHits(iTrig), pmt_grid, tag, 4.0),
		     output );

	  //
	  // TRes Hist
	  //

	  // Centroid
	  const std::vector<double> xCentroidSeed = {CentroidSeed[0], CentroidSeed[1], CentroidSeed[2], -TDWallSeed};
	  Save2ROOT( GetTResHist(tag+"_CentroidSeedTDWall", vHits,
				 xTrue, xCentroidSeed,
				 wPower, hPDF_TRes), output);

	  // Grid dWall
	  const std::vector<double> xGridSeedDWall = {vSeeds[vSeeds.size()-2].Pos[0], vSeeds[vSeeds.size()-2].Pos[1], vSeeds[vSeeds.size()-2].Pos[2], -TDWallSeed};
	  Save2ROOT( GetTResHist(tag+"_GridTDWallSeed", vHits,
				 xTrue, xGridSeedDWall,
				 wPower, hPDF_TRes), output);


	  // SEED
	  const std::vector<double> xSeed = {vSeeds.front().Pos[0], vSeeds.front().Pos[1], vSeeds.front().Pos[2], -vSeeds.front().T};
	  Save2ROOT( GetTResHist(tag+"_BestSeed", vHits,
				 xTrue, xSeed,
				 wPower, hPDF_TRes), output);

	  // REC
	  Save2ROOT( GetTResHist(tag, vHits,
				 xTrue, vX.front(),
				 wPower, hPDF_TRes) , output);

	  //
	  // MAP NLL SPACE
	  //

	  // RESET
	  map_nll.ResetGrid();

	  // TCentroid
	  map_nll.Fill(vHits, -TDWallSeed, hPDF_TRes, wPower);
	  Save2ROOT( GetMapPlots(map_nll.GetHGrid(), {PosTrue, PosRec, CentroidSeed} ,
				 Form("cMAP_CentroidSeedTDWall_%s", tag.c_str())),
		     output );
	  // HGrid GARBAGE COLLECTOR
	  HGridGarbageCollector();
	  // RESET
	  map_nll.ResetGrid();

	  // TCentroid
	  map_nll.Fill(vHits, -TDWallSeed, hPDF_TRes, wPower);
	  Save2ROOT( GetMapPlots(map_nll.GetHGrid(), {PosTrue, PosRec, TVector3(vSeeds[vSeeds.size()-2].Pos[0], vSeeds[vSeeds.size()-2].Pos[1], vSeeds[vSeeds.size()-2].Pos[2])} ,
				 Form("cMAP_GridTDWallSeed_%s", tag.c_str())),
		     output );
	  // HGrid GARBAGE COLLECTOR
	  HGridGarbageCollector();
	  // RESET
	  map_nll.ResetGrid();

	  // TSeed
	  map_nll.Fill(vHits, -vSeeds.front().T, hPDF_TRes, wPower);
	  Save2ROOT( GetMapPlots(map_nll.GetHGrid(), {PosTrue, PosRec, vSeeds.front().Pos} ,
				 Form("cMAP_BestSeed_%s", tag.c_str())),
		     output );
	  // HGrid GARBAGE COLLECTOR
	  HGridGarbageCollector();
	  // RESET
	  map_nll.ResetGrid();

	  // TTrue
	  map_nll.Fill(vHits, -TTrue, hPDF_TRes, wPower);
	  Save2ROOT( GetMapPlots(map_nll.GetHGrid(), {PosTrue, PosRec, vSeeds.front().Pos} ,
				 Form("cTrue_%s", tag.c_str())),
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

  auto cSeed = seed_perf_monitor.GetPlot(true);
  cSeed->Write();
  auto cDRhoSeed = seed_perf_monitor.GetDRhoPlot();
  cDRhoSeed->Write();

  if(args.useGridSearch){
    auto cGSeed = gridTDWall_seed_perf_monitor.GetPlot(true);
    cGSeed->Write();
    auto cDRhoGSeed = gridTDWall_seed_perf_monitor.GetDRhoPlot();
    cDRhoGSeed->Write();
  }

  auto cFit = recon_perf_monitor.GetPlot(true);
  cFit->Write();
  auto cDRho = recon_perf_monitor.GetDRhoPlot();
  cDRho->Write();

  fOut.Close();

  theApp.Terminate();

  return EXIT_SUCCESS;
}
