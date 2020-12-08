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
  const char* cOutname = "fRecon";
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

  for(int iEvt=0; iEvt<nEvts; iEvt++) {
    ++progress_bar;
    w_rat.SetEvt(iEvt);

    auto nTriggers = w_rat.GetNTriggers();

    for (auto iTrigger = 0; iTrigger < nTriggers; iTrigger++) {

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

    }

    if(isVerbose)
      progress_bar.display();

  }

  if(isVerbose)
    progress_bar.done();

  return EXIT_SUCCESS;

}
