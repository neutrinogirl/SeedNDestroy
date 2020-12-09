//
// Created by zsoldos on 11/25/20.
//

#include "SeedNDestroy.hh"
#include "CreatePDF.hh"
#include "Multilateration.hh"
#include "Recon.hh"

#include "wRATter/include/Wrapper.hh"

#include <ProgressBar.hpp>

int main(int argc, char *argv){

  const char* filename;
  const char* pdfname;
  const char* cOutname = "fRecon.root";
  unsigned int nEvts = 0;
  int wPower = 1;
  bool isVerbose = false;

  // #################################################### //
  // #### #### #### OPEN FILE / READ TTREE #### #### #### //
  // #################################################### //

  wRAT w_rat(filename);
  nEvts = nEvts > 0 ? nEvts : w_rat.GetNEvts();
  ProgressBar progressBar(nEvts, 70);

  std::cout << "USE PDF: " << pdfname << std::endl;

  auto PDF = GetRootHisto<TH2D>(pdfname, "hCTVSTResPDF");
  auto PDF_TRes = PDF->ProjectionX();
  auto PDF_CT = PDF->ProjectionY();

  // Create structure holding data
  DataStruct1D ds = {PDF_TRes, std::vector<Hit>(), wPower};
  DataStructDir dsdir = {PDF_CT, std::vector<Hit>(), wPower, std::vector<double>()};

  // ################################## //
  // #### #### #### LOOP #### #### #### //
  // ################################## //

  ProgressBar progress_bar(nEvts, 70);

  // Original true Pos, Dir, T
  double mcx, mcy, mcz, mcT, mcdx, mcdy, mcdz; 
  // Recon Pos, T
  double recx, recy, recz, recT;
  
  TTree tree("T", "A zara tree");
  tree.Branch("mcx", &mcx, "mcx");
  tree.Branch("mcy", &mcy, "mcy");
  tree.Branch("mcz", &mcz, "mcz");
  tree.Branch("mcT", &mcT, "mcT");
  tree.Branch("mcdx", &mcdx, "mcdx");
  tree.Branch("mcdy", &mcdy, "mcdy");
  tree.Branch("mcdz", &mcdz, "mcdz");
  tree.Branch("recx", &recx, "recx");
  tree.Branch("recy", &recy, "recy");
  tree.Branch("recz", &recz, "recz");
  tree.Branch("recT", &recT, "recT");

  for(int iEvt=0; iEvt<nEvts; iEvt++) {
    ++progress_bar;
    w_rat.SetEvt(iEvt);

    auto nTriggers = w_rat.GetNTriggers();

    for (auto iTrigger = 0; iTrigger < nTriggers; iTrigger++) {

    // Get True Info
      const auto TrigTime = 0;
      
      auto nParticle = w_rat.GetNPrimaryParticle();
      auto iParticle = nParticle > 1 ? (TrigTime > 1e3 ? 1 : 0) : 0;

      if(iParticle>0)
	continue;

      if(iTrigger>0)
	continue;

      const auto PosTrue = w_rat.GetPosTrue(iParticle) - TVector3(0., -133.3, 1724.);
      const auto DirTrue = w_rat.GetDirTrue(iParticle);
      const auto TTrue = w_rat.GetTTrue(iParticle);

      mcx = PosTrue.x(); mcy = PosTrue.y(); mcz = PosTrue.z();
      mcT = TTrue;
      mcdx = DirTrue.x(); mcdy = DirTrue.y(); mcdz = DirTrue.z();

      auto vHits = w_rat.GetVHits(iTrigger);
      std::sort(vHits.begin(), vHits.end());

      // Get Seed
      auto PosTSeed = GetSeed(vHits, PDF_TRes, wPower);

      //
      // ####
      //

      ds.vHits.clear();
      ds.vHits = vHits;

      auto x = ReconPosTime(ds, PosTSeed.Pos, PosTSeed.T);
      std::cout << x[0] << "mm "
		<< x[1] << "mm "
		<< x[2] << "mm "
		<< x[3]*1.e-2 << "ns " << std::endl;

      recx = x[0]; recy = x[1]; recz = x[2];
      recT = x[3]*1.e-2;

      tree.Fill();

    }

    if(isVerbose)
      progress_bar.display();

  }

  if(isVerbose)
    progress_bar.done();

  TFile fOut(cOutname, "RECREATE");
  tree.Write();
  fOut.Close();

  return EXIT_SUCCESS;

}
