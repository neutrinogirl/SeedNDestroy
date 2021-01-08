//
// Created by zsoldos on 1/7/21.
//

#ifndef SEEDNDESTROY_INCLUDE_DRAWNLLSPACE_HH_
#define SEEDNDESTROY_INCLUDE_DRAWNLLSPACE_HH_

#include <TStyle.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TH2D.h>

#include "MathUtils.hh"
#include "PathFit.hh"
#include "Centroid.hh"
#include "Multilateration.hh"

TCanvas *DrawNLLSpace(std::vector<Hit>& vHits, TH1D* hPDF,
					  const double& TTrue, const TVector3& PosTrue,
					  const std::string& tag,
					  const double& WidthBinsGrid = 500., const double HalfDetSize = 8000.){

  gStyle->SetOptStat(0);

  auto *c = new TCanvas(Form("c_%s", tag.c_str()), Form("c_%s", tag.c_str()),
						800, 800);
  auto *p1 = new TPad(Form("p_%s", tag.c_str()),
					  Form("p_%s", tag.c_str()), 0, 0, 1, 1);
  p1->SetFillStyle(0);
  p1->SetFillColor(0);
  p1->SetFrameFillStyle(0);

  const int nBinsGrid = std::ceil(2*HalfDetSize/WidthBinsGrid);

  TH2D *hNLLGrid = new TH2D(Form("hNLL_%s", tag.c_str()), "Grid NLL",
							nBinsGrid, -HalfDetSize, HalfDetSize,
							nBinsGrid, -HalfDetSize, HalfDetSize);

  // Get PMT hits position
  std::vector<TMarker*> vMHits(vHits.size());
  for(auto iHit=0; iHit<vHits.size(); iHit++){
	CylVec CylPosHit(vHits[iHit].PMTPos);
	vMHits[iHit] = new TMarker(CylPosHit.rho, CylPosHit.z, kFullCircle);
	vMHits[iHit]->SetMarkerSize(2);
	vMHits[iHit]->SetMarkerColor(kRed-4);
  }

  // Get True pos
  CylVec CylPosTrue(PosTrue);
  // Prepare marker
  auto TrueMarker = new TMarker(CylPosTrue.rho, CylPosTrue.z, kFullStar);
  TrueMarker->SetMarkerSize(2);
  TrueMarker->SetMarkerColor(kRed-4);

  // Integrate over theta
  const unsigned int nBinsTheta = 12;
  std::vector<double> vTheta(nBinsTheta);
  std::iota(vTheta.begin(), vTheta.end(), 0);
  std::transform(vTheta.begin(), vTheta.end(), vTheta.begin(),
				 [&nBinsTheta](const double& val){
				   return -1 + val / (double)(nBinsTheta);
				 }
  );

  for(int r=0; r<nBinsGrid; r++){

	// Get middle of the bin
	const double R = -HalfDetSize + (0.5 + r)*WidthBinsGrid;

	for(int z=0; z<nBinsGrid; z++){

	  // Get Middle of the bin
	  const double Z = -HalfDetSize + (0.5 + z)*WidthBinsGrid;

	  double AvNLL = 0.;

	  for(auto& t:vTheta){

		// Convert to cartesian
		TVector3 PosGuess(R*cos(t), R*sin(t), Z);
		AvNLL+=flatf(PosGuess, TTrue, vHits, hPDF, 1);

	  }

	  hNLLGrid->Fill(R, Z, AvNLL);

	}

  }

  auto vSeeds = GetVSeeds(vHits, TTrue, 0, hPDF, 1);
  auto nSeeds = vSeeds.size();
  nSeeds = nSeeds > 10 ? 10 : nSeeds;
  std::vector<TMarker*> vMSeed(nSeeds);
  for(auto iSeed=0; iSeed<nSeeds; iSeed++){
	CylVec CylPosSeed(vSeeds[iSeed]);
	vMSeed[iSeed] = new TMarker(CylPosSeed.rho, CylPosSeed.z, kOpenCross);
	vMSeed[iSeed]->SetMarkerSize(2);
	if(iSeed == 0)
	  vMSeed[iSeed]->SetMarkerSize(3);
	vMSeed[iSeed]->SetMarkerColor(kGreen+1);
  }

  auto CentroidSeed = GetCentroidSeed(vHits, 2);
  CentroidSeed.Print();
  CylVec CylCentroidSeed(CentroidSeed);
  CylCentroidSeed.Print();
  auto *mCentroid = new TMarker(CylCentroidSeed.rho, CylCentroidSeed.z, kOpenCross);
  mCentroid->SetMarkerSize(2);
  mCentroid->SetMarkerColor(kBlack);

  c->Draw();
  p1->Draw();
  hNLLGrid->Draw("COLZ");
  TrueMarker->Draw("");
  for(auto iHit=0; iHit<vHits.size(); iHit++){
	vMHits[iHit]->Draw("");
  }
  for(auto iSeed=0; iSeed<nSeeds; iSeed++){
	vMSeed[iSeed]->Draw("");
  }
  mCentroid->Draw("");
  return c;


}


#endif //SEEDNDESTROY_INCLUDE_DRAWNLLSPACE_HH_
