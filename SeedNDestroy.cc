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

  const char* filename="test_Z.root";
  const char* pdfname="fout_Z.root";
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
//PDF_TRes->Rebin(10);

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
  tree.Branch("mcx", &mcx, "mcx/D");
  tree.Branch("mcy", &mcy, "mcy/D");
  tree.Branch("mcz", &mcz, "mcz/D");
  tree.Branch("mcT", &mcT, "mcT/D");
  tree.Branch("mcdx", &mcdx, "mcdx/D");
  tree.Branch("mcdy", &mcdy, "mcdy/D");
  tree.Branch("mcdz", &mcdz, "mcdz/D");
  tree.Branch("recx", &recx, "recx/D");
  tree.Branch("recy", &recy, "recy/D");
  tree.Branch("recz", &recz, "recz/D");
  tree.Branch("recT", &recT, "recT/D");

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

//      const auto PosTrue = w_rat.GetPosTrue(iParticle) - TVector3(0., -133.3, 1724.);
 //     const auto DirTrue = w_rat.GetDirTrue(iParticle);


      auto vHits = w_rat.GetVHits(iTrigger);
      //auto NHits = w_rat.GetNHits(iTrigger);

//std::cout<<"Hits"<<NHits<<"   "<<vHits.empty()<<std::endl;
if (vHits.empty()==1){
std::cout<<"Hits 0 "<<std::endl;
continue;
}

      TVector3  DefPosTrue = w_rat.GetPosTrue(iParticle);
      //const auto DefPosTrue = w_rat.GetPosTrue(iParticle);
	TVector3 NewPosTrue;
	NewPosTrue.SetX(DefPosTrue.X());
	NewPosTrue.SetY(-1*(DefPosTrue.Z()-1724));
	NewPosTrue.SetZ(DefPosTrue.Y()+133.3);

TVector3 DefDirTrue = w_rat.GetDirTrue(iParticle);

TVector3 NewDirTrue;
NewDirTrue.SetX(DefDirTrue.X());
NewDirTrue.SetY(-1*DefDirTrue.Z());
NewDirTrue.SetZ(DefDirTrue.Y());


      const auto TTrue = w_rat.GetTTrue(iParticle);

      mcx = NewPosTrue.X(); mcy = NewPosTrue.Y(); mcz = NewPosTrue.Z();
      mcT = TTrue;
      mcdx = NewDirTrue.X(); mcdy = NewDirTrue.Y(); mcdz = NewDirTrue.Z();

      std::sort(vHits.begin(), vHits.end());

      // Get Seed
      auto PosTSeed = TVector3(0,0,0);//GetCentroidSeed(vHits); //reverting to the previous seeding
      //auto PosTSeed = GetCentroidSeed(vHits); //reverting to the previous seeding
      //auto PosTSeed = GetSeed(vHits, PDF_TRes, wPower);  

      //
      // ####
      //

      ds.vHits.clear();
      ds.vHits = vHits;

      auto x = ReconPosTime(ds, PosTSeed, 0); //reverting to the previous seeding
      //auto x = ReconPosTime(ds, PosTSeed.Pos, PosTSeed.T);
      std::cout << x[0] << "mm "
		<< x[1] << "mm "
		<< x[2] << "mm "
		<< x[3] << "ns " << std::endl;

      recx = x[0]; recy = x[1]; recz = x[2];
      recT = x[3];

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
