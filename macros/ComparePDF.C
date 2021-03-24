#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TColor.h>
#include <TMultiGraph.h>

#include <DrawingUtils.hh>

template <typename T>
T* GetRootHisto(const char* filename, const char* histname){
  auto f = TFile::Open(filename);
  // Check if key exist
  if(!f->GetListOfKeys()->Contains(histname))
	return nullptr;
  auto hist = dynamic_cast<T *>(f->Get(histname)->Clone());
  hist->SetDirectory(nullptr);
  f->Close();
  delete f;
  return hist;
}

void Compare(const std::vector<std::string>& vFilenames){

  const std::vector<std::string>& vHistnames = {
	  "hCTVSTResPDF_THit_px",
	  "hCTVSTResPDF_TTOF_px"
  };

  const std::vector<int> ColorPalette = {
	  TColor::GetColor("#4D9DE0"),
	  TColor::GetColor("#E15554"),
	  TColor::GetColor("#E1BC29"),
	  TColor::GetColor("#3BB273")
  };

  auto GetCol = [&ColorPalette](const unsigned int& i){
    return ColorPalette[i % ColorPalette.size()];
  };

  std::vector< std::vector< TH1D*> > vvHist( vFilenames.size() );

  std::vector<double> vMax(vHistnames.size(), std::numeric_limits<double>::min());
  std::vector<double> vMin(vHistnames.size(), std::numeric_limits<double>::max());

  unsigned int iFile = 0;
  unsigned int iHist = 0;
  for(auto& file: vFilenames){

    vvHist.emplace_back(std::vector<TH1D*>());

    for(auto& hist: vHistnames){

      // Get Histogram
      vvHist[iFile].emplace_back(GetRootHisto<TH1D>(file.c_str(), hist.c_str()));
      // Set Style
	  SetBasicTStyle(vvHist[iFile].back(), GetCol(iFile), 2, kSolid);
	  // Set Name (F you ROOT)
	  vvHist[iFile].back()->SetName(Form("%u_%s", iFile, hist.c_str()));
	  // Find min/max
	  vMax[iHist] = std::max(vMax[iHist], vvHist[iFile].back()->GetMaximum());
	  vMin[iHist] = std::min(vMin[iHist], vvHist[iFile].back()->GetMinimum());

	  iHist++;

    }

    iHist = 0;

    iFile++;

  }

  std::vector<TCanvas*> vC = {
  	PlotAHist(vvHist[0][0], "HIST"),
	PlotAHist(vvHist[0][1], "HIST")
  };

  vvHist[0][0]->GetYaxis()->SetRangeUser(vMin[0], vMax[0]);
  vvHist[0][1]->GetYaxis()->SetRangeUser(vMin[1], vMax[1]);

  for(iFile = 1; iFile<vvHist.size(); iFile++ ){
    vC[0]->cd();
    vvHist[iFile][0]->Draw("HISTSAME");
	vC[1]->cd();
	vvHist[iFile][1]->Draw("HISTSAME");
  }

}

void ComparePDF(){

  gStyle->SetOptStat(0);

  const std::vector<std::string> vPositronFilename = {
	  "PDFs/positrons_3.0MeV_90DegDirSmear.root",
	  "PDFs/positrons_3.0MeV_60DegDirSmear.root",
	  "PDFs/positrons_3.0MeV_45DegDirSmear.root",
	  "PDFs/positrons_3.0MeV_30DegDirSmear.root"
  };

  const std::vector<std::string> vGammaFilename = {
	  "PDFs/gamma_3.0MeV_90DegDirSmear.root",
	  "PDFs/gamma_3.0MeV_60DegDirSmear.root",
	  "PDFs/gamma_3.0MeV_45DegDirSmear.root",
	  "PDFs/gamma_3.0MeV_30DegDirSmear.root"
  };

  const std::vector<std::string> vElectronFilename = {
	  "PDFs/electrons_3.0MeV_90DegDirSmear.root",
	  "PDFs/electrons_3.0MeV_60DegDirSmear.root",
	  "PDFs/electrons_3.0MeV_45DegDirSmear.root",
	  "PDFs/electrons_3.0MeV_30DegDirSmear.root",
  };

  const std::vector<std::string> vLegName = {
	  "90^{#circ}",
	  "60^{#circ}",
	  "45^{#circ}",
	  "30^{#circ}"
  };

  const std::vector<int> ColorPalette = {
  	TColor::GetColor("#4D9DE0"),
	TColor::GetColor("#E15554"),
	TColor::GetColor("#E1BC29"),
	TColor::GetColor("#3BB273"),
  };

  const std::vector<int> StylePalette = {
	  kSolid,
	  kDotted,
	  kDashed
  };

  const std::size_t nCols   = ColorPalette.size();
  const std::size_t nStyles = StylePalette.size();
  auto GetIdx = [&nCols](unsigned int& i, const unsigned int& n){ return i % n;};

  const std::vector<std::vector<std::string>> vvFiles = {
  	// vPositronFilename,
  	// vGammaFilename,
  	vElectronFilename
  };

  auto ScaleHist = [](TH1D* h){h->Scale(1./h->GetMaximum());};
  unsigned int iFile = 0;
  unsigned int iStyle = 0;

  auto c = new TCanvas("c", "c", 800, 600);

  for(auto& vFile:vvFiles){
	for(auto& file:vFile){
	  auto h = GetRootHisto<TH1D>(file.c_str(), "hCTVSTResPDF_py");
	  h->SetName(Form("%s", vLegName[iFile].c_str()));
	  h->SetTitle(Form("%s ; cos(#theta) ; Counts (Normalized)", vLegName[iFile].c_str()));
	  ScaleHist(h);
	  SetBasicTStyle(h, ColorPalette[GetIdx(iFile, nCols)], 2, StylePalette[GetIdx(iStyle, nStyles)]);
	  h->Draw("HISTSAME");

	  // PlotAHist(h,"HIST");
	  iFile++;
	}
	iStyle++;
  }

  c->BuildLegend();

}