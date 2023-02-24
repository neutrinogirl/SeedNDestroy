//
// Created by Stephane Zsoldos on 7/6/22.
//

#include <csignal>

#include <TH2D.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMarker.h>

#include "Map.hh"
#include "Templates/TData.hh"

#include "Algo/VHits.hh"

#include <ROOT/Utils.hh>

MapAnalysis::MapAnalysis(const char *pdfname, const char *histname, const char* perpmthistname,
						 const double &R, const double &HH,
						 const int& nEvts,
						 const bool& applytrigger,
						 const int& nSBins, const int& nTBins,
						 const char *filename, const char *treename)
	: nMaxEvts(nEvts), isapplytrigger(applytrigger), nSpaceBins(nSBins), nTimeBins(nTBins){
  //
  hPDF = GetROOTObj<TH2D>(pdfname, histname)->ProjectionX("hPDF");
  std::cout << "Load PDF: " << hPDF->GetName() << std::endl;
  mPDF2D = GetROOTMObj<TH2D>(pdfname, perpmthistname, "TH2D");
  std::transform(
	  mPDF2D.begin(), mPDF2D.end(),
	  std::inserter(mPDF1D, mPDF1D.begin()),
	  [](const std::pair<int, TH2D*>& p){
		return std::make_pair(p.first, p.second->ProjectionX());
	  }
  );
  //
  Cyl = new Cylinder(R, HH);
  //
  OFile = new TFile(filename, "RECREATE");
  Tree = new TTree(treename, treename);
  RT.SetTree(Tree);
  //
  nMaxEvts = nMaxEvts < 0 ? std::numeric_limits<int>::max() : nMaxEvts;
  // Get cartesian coordinates of Cylinder from Radius R and half-heigh HH
  std::vector< std::vector<double> > vvEdges = {
	  linspace<double>(-Cyl->GetEdge().x(), Cyl->GetEdge().x(), nSpaceBins),
	  linspace<double>(-Cyl->GetEdge().y(), Cyl->GetEdge().y(), nSpaceBins),
	  linspace<double>(-Cyl->GetEdge().z(), Cyl->GetEdge().z(), nSpaceBins)
  };
  std::vector< std::vector<double> > vvCenter = vvEdges;
  std::for_each(
	  vvCenter.begin(),
	  vvCenter.end(),
	  [](std::vector<double>& v){
		std::transform(
			v.begin(),
			v.end() - 1,
			v.begin() + 1,
			v.begin(),
			[](double a, double b) {
			  return (a + b) / 2.0;
			}
		);
	  }
  );
  SG = new SpaceGrid(vvCenter);
  h3D = new TH3D("h3D", "h3D",
				 vvEdges[0].size() - 1, vvEdges[0].data(),
				 vvEdges[1].size() - 1, vvEdges[1].data(),
				 vvEdges[2].size() - 1, vvEdges[2].data());

  vTEdges = linspace<double>(-1.f, Cyl->GetTEdge(), nTimeBins);
  vTCenter = vTEdges;
  std::transform(
	  vTCenter.begin(),
	  vTCenter.end() - 1,
	  vTCenter.begin() + 1,
	  vTCenter.begin(),
	  [](double a, double b) {
		return (a + b) / 2.0;
	  }
  );

}

MapAnalysis::~MapAnalysis(){
  delete OFile;
  delete Cyl;
  delete hPDF;
  for(auto& p : mPDF2D)
	delete p.second;
  delete SG;
  // delete h3D;
}

void MapAnalysis::Do(void *Data) {

  // Get Data
  auto wData = static_cast<TData*>(Data);
  auto vHits = wData->GetVHits();
  //
  if(isapplytrigger){
	double T = GetFirstHitTime(vHits, 1.f);
	std::transform(
		vHits.begin(),
		vHits.end(),
		vHits.begin(),
		[T](const Hit& h){
		  return h-T;
		}
	);
  }
  //
  auto iEvt = wData->GetEventID();
  const char *tag = Form("Evt%d", iEvt);
  //
  if(iEvt > nMaxEvts)
	raise(SIGINT);

  // Init vSeeds with Centroid
  std::vector<PosT> vSeeds = {
	  GetCentroidBasedSeed(vHits, Cyl)
  };

  // Get LS seed
  vSeeds.emplace_back(GetLSBasedSeed(vHits, Cyl, vSeeds));

  TVector3 v3BF = ConvertTVector3Unit<double>(vSeeds.back().GetTVector3(), SpaceUnit::dm, SpaceUnit::mm);

  std::vector<TMarker*> vBFMarkers = {
	  new TMarker(v3BF.x(), v3BF.z(), kFullCross),
	  new TMarker(v3BF.y(), v3BF.z(), kFullCross),
	  new TMarker(v3BF.x(), v3BF.y(), kFullCross)
  };

  for(auto& m : vBFMarkers){
	m->SetMarkerColor(kGreen+1);
	m->SetMarkerSize(2);
  }

  TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 400);
  canvas->Divide(3,1);

  gSystem->Unlink(Form("map_%s.gif", tag));

  for(const auto& T: vTCenter){
	SG->ParallelWalk(*hPDF, vHits, T);

	for(int i=0; i<SG->GetNPts(); i++){
	  h3D->Fill(SG->vPts[i].x(), SG->vPts[i].y(), SG->vPts[i].z(), SG->vNLL[i]);
	}

	// Get the bin with the minimum content
	Int_t bin_min = h3D->GetMinimumBin();

	// Convert the bin number to x, y, z coordinates
	Int_t binx, biny, binz;
	h3D->GetBinXYZ(bin_min, binx, biny, binz);
	Double_t x_min = h3D->GetXaxis()->GetBinCenter(binx);
	Double_t y_min = h3D->GetYaxis()->GetBinCenter(biny);
	Double_t z_min = h3D->GetZaxis()->GetBinCenter(binz);

	// Print the coordinates and bin content
	// Double_t min_content = h3D->GetBinContent(bin_min);
	// std::cout << "Minimum bin content at (" << x_min << ", " << y_min << ", " << z_min << ") = " << min_content << " ";
	// std::cout << "LS: " << GetSum2Residuals(TVector3(x_min, y_min, z_min), T, vHits) << std::endl;

	std::vector<TMarker*> vMarkers = {
		new TMarker(x_min, z_min, kFullDiamond),
		new TMarker(y_min, z_min, kFullDiamond),
		new TMarker(x_min, y_min, kFullDiamond)
	};

	for(auto& m: vMarkers){
	  m->SetMarkerColor(kBlue);
	  m->SetMarkerSize(2);
	}

	auto hxz = h3D->Project3DProfile("xz");
	canvas->cd(1);
	hxz->Draw("colz");
	gStyle->SetOptStat(0);
	gPad->Update();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	vMarkers[0]->Draw("same");
	vBFMarkers[0]->Draw("same");
	auto hyz = h3D->Project3DProfile("yz");
	canvas->cd(2);
	hyz->Draw("colz");
	gStyle->SetOptStat(0);
	gPad->Update();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	vMarkers[1]->Draw("same");
	vBFMarkers[1]->Draw("same");
	auto hxy = h3D->Project3DProfile("xy");
	canvas->cd(3);
	hxy->Draw("colz");
	gStyle->SetOptStat(0);
	gPad->Update();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	vMarkers[2]->Draw("same");
	vBFMarkers[2]->Draw("same");
	canvas->Update();
	canvas->Print(Form("map_%s.gif+25", tag), "gif");

	OFile->cd();
	hxz->Write(Form("hxz_%s_%.1f", tag, T));
	hyz->Write(Form("hyz_%s_%.1f", tag, T));
	hxy->Write(Form("hxy_%s_%.1f", tag, T));

	h3D->Reset();

	delete hxz;
	delete hyz;
	delete hxy;

	for(auto& m: vMarkers) {
	  delete m;
	}

  }

  delete canvas;

  for(auto &m: vBFMarkers){
	delete m;
  }

}

void MapAnalysis::Export() const {
  OFile->cd();
  Tree->Write();
  OFile->Close();
}