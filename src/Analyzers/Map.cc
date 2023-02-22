//
// Created by Stephane Zsoldos on 7/6/22.
//

#include <csignal>

#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>

#include "Map.hh"
#include "Templates/TData.hh"

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
  auto iEvt = wData->GetEventID();
  auto iTrig = wData->GetTriggerID();
  const char *tag = Form("Evt%d_Trigger%d", iEvt, iTrig);
  //
  if(iEvt > nMaxEvts)
	raise(SIGINT);


  TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 400);
  canvas->Divide(3,1);

  gSystem->Unlink(Form("map_%s.gif", tag));

  for(const auto& T: vTCenter){
	SG->ParallelWalk(hPDF, vHits, T);

	for(int i=0; i<SG->GetNPts(); i++){
	  h3D->Fill(SG->vPts[i].x(), SG->vPts[i].y(), SG->vPts[i].z(), SG->vNLL[i]);
	}

	auto hxz = h3D->Project3D("xz");
	canvas->cd(1);
	hxz->Draw("colz");
	gStyle->SetOptStat(0);
	gPad->Update();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	auto hyz = h3D->Project3D("yz");
	canvas->cd(2);
	hyz->Draw("colz");
	gStyle->SetOptStat(0);
	gPad->Update();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	auto hxy = h3D->Project3D("xy");
	canvas->cd(3);
	hxy->Draw("colz");
	gStyle->SetOptStat(0);
	gPad->Update();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
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

  }

  delete canvas;

}

void MapAnalysis::Export() const {
  OFile->cd();
  Tree->Write();
  OFile->Close();
}