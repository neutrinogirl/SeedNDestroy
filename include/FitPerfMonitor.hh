//
// Created by zsoldos on 1/12/21.
//

#ifndef SEEDNDESTROY_INCLUDE_FITPERFMONITOR_HH_
#define SEEDNDESTROY_INCLUDE_FITPERFMONITOR_HH_

#include <TH1D.h>

#include <DrawingUtils.hh>

#include <FitResults.hh>

class FitPerfMonitor {

 private:
  const int nDim = 3;

  double dAxis = 8.e3;
  double wBin = 1.e2;
  int nBins = std::floor(2*dAxis/wBin);

  std::vector<int> ColorPalette = CreateBasicPalette();

  std::vector<TH1D*> hD;
  TH1D* hDRho;

  void CreateHistogram(const std::string& tag="hDAxis"){
	for(int iDim=0; iDim<3; iDim++){
	  hD.emplace_back(new TH1D(Form("%s%d", tag.c_str(), iDim), "", nBins, -dAxis, dAxis));
	  SetBasicTStyle(hD[hD.size()-1], ColorPalette[iDim], 2, kSolid);
	}
	hDRho = new TH1D(Form("hDRho_%s", tag.c_str()), "", nBins, 0., dAxis);
  }

 public:
  FitPerfMonitor(){
	CreateHistogram();
  }

  explicit FitPerfMonitor(const std::string& tag){
    CreateHistogram(tag);
  }

  FitPerfMonitor(const std::string& tag, const double& daxis, const double& wbin)
	  : dAxis(daxis), wBin(wbin) {
	CreateHistogram(tag);
  }

  virtual ~FitPerfMonitor() = default;

  void DeleteHist(){
	for(int iDim=0; iDim<3; iDim++) {
	  delete hD[iDim];
	}
	delete hDRho;
  }

  void Fill(const TVector3& PosGuess, const TVector3& PosTrue){
	auto diff = PosGuess - PosTrue;
	for(int iDim=0; iDim<3; iDim++) {
	  hD[iDim]->Fill(diff[iDim]);
	}
	hDRho->Fill(diff.Mag());
  }

  TCanvas *GetPlot(bool isFit = false){

    if(isFit){
      std::vector<GausFit::Res> vRes;
      for(auto& h:hD){
        GausFit::Res R = GausFit::GetFitRes(h);
		vRes.emplace_back(R);
		R.Sig.Print();
      }
    }

	std::cout << std::endl;

	auto c = PlotAHist(hD[kX], "HIST");
	hD[kY]->Draw("HISTSAME");
	hD[kZ]->Draw("HISTSAME");

	return c;

  }

  TCanvas *GetDRhoPlot(){
    auto c = PlotAHist(hDRho, "HIST");
    return c;
  }

  // ######################### //
  // #### GETTERS/SETTERS #### //
  // ######################### //

  double GetDAxis() const { return dAxis; }
  void SetDAxis(double d_axis) { dAxis = d_axis; }
  double GetWBin() const { return wBin; }
  void SetWBin(double w_bin) { wBin = w_bin; }
  int GetNBins() const { return nBins; }
  void SetNBins(int n_bins) { nBins = n_bins; }
  const std::vector<int> &GetColorPalette() const { return ColorPalette; }
  void SetColorPalette(const std::vector<int> &color_palette) { ColorPalette = color_palette; }

};

TCanvas* GetTResHist(const std::string& tag, const std::vector<Hit>& vHits,
					 const std::vector<double>& vTrue, const std::vector<double>& vRec,
					 const unsigned int wPower = 1,
					 TH1D* hPDF = nullptr){

  zAxis axTRes = hPDF ? zAxis(hPDF->GetXaxis()) : zAxis {600, -200, 400};

  std::vector<TH1D*> vHTRes = {
  	new TH1D(Form("hTrue%s", tag.c_str()), Form(" TRUE ; T_{Res} [ns] ; Hits"), axTRes.nBins, axTRes.min, axTRes.max),
	new TH1D(Form("hRec%s", tag.c_str()), " REC ; T_{Res} [ns] ; Hits", axTRes.nBins, axTRes.min, axTRes.max)
  };

  enum {
    HTRUE=0,
    HREC=1
  };

  SetBasicTStyle(vHTRes[HTRUE], kBlue-4);
  SetBasicTStyle(vHTRes[HREC], kRed-4);

  auto FillHist = [](TH1D*h, const Hit& hit, const std::vector<double> &v, const unsigned int& w){
	h->Fill(hit.GetTRes(TVector3(v[0], v[1], v[2]), v[3]), fweight(hit, w));
  };

  for(const auto& hit:vHits){
    FillHist(vHTRes[HTRUE], hit, vTrue, wPower);
	FillHist(vHTRes[HREC], hit, vRec, wPower);
  }

  if(hPDF)
	std::for_each(vHTRes.begin(), vHTRes.end(), [&hPDF](TH1D *h){h->Scale(hPDF->GetMaximum()/h->GetMaximum());});
  else
	std::for_each(vHTRes.begin(), vHTRes.end(), [](TH1D *h){h->Scale(1./h->GetMaximum());});

  auto c = PlotAHist(vHTRes[HTRUE], "HIST");
  vHTRes[HREC]->Draw("HISTSAME");
  if(hPDF)
    hPDF->Draw("HISTSAME");
  c->BuildLegend();

  return c;

}

TCanvas* GetTResHist(const std::vector<std::string>& vtag, const std::vector<Hit>& vHits,
					 const std::vector<std::vector<double>>& vv,
					 const unsigned int wPower = 1,
					 TH1D* hPDF = nullptr){

  std::vector<int> ColorPalette = CreateBasicPalette();

  zAxis axTRes = hPDF ? zAxis(hPDF->GetXaxis()) : zAxis {600, -200, 400};

  std::vector<TH1D*> vHTRes;
  const std::size_t nHits = vv.size();
  vHTRes.reserve(nHits);

  for(auto iHist=0; iHist<nHits; iHist++){
    std::string tag = iHist>vtag.size() ? Form("h%d", iHist) : vtag[iHist];
	vHTRes.emplace_back(GetTResHist(tag, vHits, vv[iHist], 1, axTRes.nBins, axTRes.min, axTRes.max));
	SetBasicTStyle(vHTRes.back(), ColorPalette[iHist % ColorPalette.size()]);
  }

  if(hPDF)
	std::for_each(vHTRes.begin(), vHTRes.end(), [&hPDF](TH1D *h){h->Scale(hPDF->GetMaximum()/h->GetMaximum());});
  else
	std::for_each(vHTRes.begin(), vHTRes.end(), [](TH1D *h){h->Scale(1./h->GetMaximum());});

  auto c = PlotAHist(vHTRes[0], "HIST");
  for(auto iHist=1; iHist<nHits; iHist++)
	vHTRes[iHist]->Draw("HISTSAME");
  if(hPDF)
	hPDF->Draw("HISTSAME");
  c->BuildLegend();

  return c;

}

#endif //SEEDNDESTROY_INCLUDE_FITPERFMONITOR_HH_
