//
// Created by zsoldos on 1/14/21.
//

#ifndef SEEDNDESTROY_INCLUDE_RECORDER_HH_
#define SEEDNDESTROY_INCLUDE_RECORDER_HH_

#include <numeric>
#include <vector>

#include <TROOT.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TMarker.h>

#include "NLL.hh"

typedef struct Recorder {

  unsigned long int iCall = 0;
  std::vector<TVector3> vPosGuess;
  std::vector<double> vTGuess;
  std::vector<double> vNLL;
  Recorder() = default;

  void Reset(){
	vPosGuess.clear();
	vTGuess.clear();
	vNLL.clear();
	iCall = 0;
  }

} Recorder;

typedef struct Grid1D {

  std::vector<double> vEdges;
  std::vector<double> vCenter;

  Grid1D(const double& bnds, const double& width) {
	const unsigned int nBins = std::ceil(2*bnds / width);
	vEdges.resize(nBins+1);
	vCenter.resize(nBins);
	std::iota(vEdges.begin(), vEdges.end(), 0);
	std::iota(vCenter.begin(), vCenter.end(), 0);

	std::transform(vEdges.begin(), vEdges.end(), vEdges.begin(),
				   [&](const double& val) {return -bnds + val*width;});
	std::transform(vCenter.begin(), vCenter.end(), vCenter.begin(),
				   [&](const double& val) {return -bnds + (val+0.5)*width;});
  }

} Grid1D;

static std::vector<TVector3> Get3DSmplVec(const double& width, const TVector3& orig = TVector3(0.,0.,0.)){

  const unsigned int nDim = 3;
  const std::size_t nPts = std::pow(2, nDim);
  std::vector<TVector3> vSmplVec; vSmplVec.reserve(nPts+1);
  vSmplVec.emplace_back(orig);

  auto Getm1 = [&width](const unsigned int &b){return !b ? -width/2 : width/2;};

  for(auto i=0; i<nPts;i++){
	std::bitset<3> word(i);
	const unsigned int x = word[2];
	const unsigned int y = word[1];
	const unsigned int z = word[0];
	vSmplVec.emplace_back(orig+TVector3(Getm1(x), Getm1(y), Getm1(z)));
  }

  return vSmplVec;

}

static std::vector<TVector3> Get3DSmplVec(const std::vector<double>& vWidth, const TVector3& orig = TVector3(0.,0.,0.)){

  const unsigned int nDim = 3;
  if(vWidth.size() != nDim)
    std::cerr << "static std::vector<TVector3> Get3DSmplVec ERROR specify vWidth for x y z" << std::endl;
  const std::size_t nPts = std::pow(2, nDim);
  std::vector<TVector3> vSmplVec; vSmplVec.reserve(nPts+1);
  vSmplVec.emplace_back(orig);

  auto Getm1 = [&vWidth](const unsigned int &b, const unsigned int& idx){return !b ? -vWidth[idx]/2 : vWidth[idx]/2;};

  for(auto i=0; i<nPts;i++){
	std::bitset<3> word(i);
	const unsigned int x = word[2];
	const unsigned int y = word[1];
	const unsigned int z = word[0];
	vSmplVec.emplace_back(orig+TVector3(Getm1(x, 0), Getm1(y, 1), Getm1(z, 2)));
  }

  return vSmplVec;

}


class MapNLL {
 private:

  enum {
	kX=0,
	kY,
	kZ
  };

  std::vector<Grid1D> vGrids;
  std::vector<TVector3> vSmplPts;

  TH3D *hGrid;
  TVector3 vMinNLL;
  double MinNLL;

 public:

  //
  // #### CONSTRUCTORS / DESTRUCTORS #### //
  //

  explicit MapNLL(std::vector<Grid1D> v_grids)
	  : vGrids(std::move(v_grids)){
    hGrid = CreateHGrid();
    vMinNLL = TVector3();
    MinNLL = std::numeric_limits<double>::max();
  }
  MapNLL(const double& bnds, const double& width){
	vGrids = std::vector<Grid1D>(3, Grid1D(bnds, width));
	vSmplPts = Get3DSmplVec(width);
	hGrid = CreateHGrid();
	vMinNLL = TVector3();
	MinNLL = std::numeric_limits<double>::max();
  }
  MapNLL(const std::vector<double>& bnds, const std::vector<double>& width){
	double b, w;
	BOOST_FOREACH(boost::tie(b, w), boost::combine(bnds, width)){
			vGrids.emplace_back(Grid1D(b, w));
		  };
	vSmplPts = Get3DSmplVec(width);
	hGrid = CreateHGrid();
	vMinNLL = TVector3();
	MinNLL = std::numeric_limits<double>::max();
  }

  virtual ~MapNLL(){
    delete hGrid;
  };

  void ResetGrid(){
    hGrid->Reset();
	vMinNLL = TVector3();
	MinNLL = std::numeric_limits<double>::max();
  }

  //
  // #### METHODS #### //
  //

  TH3D *CreateHGrid(const std::string& tag = "hGrid" ){
	return new TH3D(tag.c_str(), "NLL Phase Space ; x [mm] ; y[mm] ; z[mm]",
					vGrids[kX].vCenter.size(), &vGrids[kX].vEdges[0],
					vGrids[kY].vCenter.size(), &vGrids[kY].vEdges[0],
					vGrids[kZ].vCenter.size(), &vGrids[kZ].vEdges[0]);
  }

  void Fill(const std::vector<Hit>& vHits, const double& T,
			TH1D* hPDF, const unsigned int& wPower = 1){

	for(const auto& x:vGrids[kX].vCenter){
	  for(const auto& y:vGrids[kY].vCenter){
		for(const auto& z:vGrids[kZ].vCenter){
		  TVector3 v(x, y, z);
		  const unsigned int GlobalBin = hGrid->FindBin(x, y, z);

		  // double AvNLL = 0.;
		  // for(auto& vPt : vSmplPts){
		  //   vPt+=v;
		  //   AvNLL+= GetNLL(vHits, hPDF, v, T, fweight, wPower) / (double)(vSmplPts.size());
		  // }

		  double NLL = GetNLL(vHits, hPDF, v, T, fweight, wPower);

		  if(NLL<MinNLL){
		    MinNLL = NLL;
		    vMinNLL = v;
		  }

		  hGrid->SetBinContent(GlobalBin, NLL);

		}
	  }
	}

  }

  //
  // #### GETTERS / SETTERS #### //
  //

  const std::vector<Grid1D> &GetVGrids() const { return vGrids; }
  TH3D *GetHGrid() const { return hGrid; }
  const TVector3 &GetVMinNll() const { return vMinNLL; }
  double GetMinNll() const { return MinNLL; }

  // ...

};

static std::vector<TMarker*> GetVMarker(const TVector3& Pos,
										int Color = kRed - 4, int Style = kFullStar){

  auto SetMarker = [&Color](TMarker* M){
	M->SetMarkerSize(2);
	M->SetMarkerColor(Color);
  };

  auto XY = new TMarker(Pos[0], Pos[1], Style);
  SetMarker(XY);
  auto XZ = new TMarker(Pos[0], Pos[2], Style);
  SetMarker(XZ);
  auto YZ = new TMarker(Pos[1], Pos[2], Style);
  SetMarker(YZ);

  return {XY, XZ, YZ};

}

static std::string GetProfileAxis(const std::vector<std::string>& opt = {"x", "y"}){

  std::vector<std::string> vAxis = { "x", "y", "z" };

  const auto itIdx = std::find(vAxis.begin(), vAxis.end(), opt.front());
  vAxis.erase(itIdx);
  const auto itJdx = std::find(vAxis.begin(), vAxis.end(), opt.back());
  vAxis.erase(itJdx);

  return *vAxis.begin();

}

// Profile along axis
TH2D* ProfileH3(TH3D* h3,
				const std::vector<std::string>& opt = {"x", "y"}) {

  const std::string strAxeProfiled = GetProfileAxis(opt);

  std::map<std::string, int> mAxisIdx = {
	  std::make_pair("x", 0),
	  std::make_pair("y", 1),
	  std::make_pair("z", 2),
  };

  std::map<std::string, TAxis*> mAxis = {
	  std::make_pair("x", h3->GetXaxis()),
	  std::make_pair("y", h3->GetYaxis()),
	  std::make_pair("z", h3->GetZaxis()),
  };

  auto h2 = new TH2D(Form("%s_%s%s", h3->GetName(), opt.front().c_str(), opt.back().c_str()),
					 Form("Profile %s%s ; %s [mm] ; %s [mm]", opt.front().c_str(), opt.back().c_str(), opt.front().c_str(), opt.back().c_str()),
					 mAxis[ opt.front() ]->GetNbins(), mAxis[ opt.front() ]->GetXmin(), mAxis[ opt.front() ]->GetXmax(),
					 mAxis[ opt.back() ]->GetNbins(), mAxis[ opt.back() ]->GetXmin(), mAxis[ opt.back() ]->GetXmax());

  for (auto i0 = 1; i0 < mAxis[opt.front()]->GetNbins()+1; i0++) {
	for (auto i1 = 1; i1 < mAxis[opt.back()]->GetNbins()+1; i1++) {

	  double min = std::numeric_limits<double>::max();
	  auto SetMin = [&min](const double &val) { min = val < min ? val : min; };

	  for (auto i2 = 1; i2 < mAxis[strAxeProfiled]->GetNbins()+1; i2++) {

		switch (mAxisIdx[opt.front()]) {
		  case 0:
			switch (mAxisIdx[opt.back()]) {
			  case 1: SetMin(h3->GetBinContent(i0, i1, i2)); break;
			  case 2: SetMin(h3->GetBinContent(i0, i2, i1)); break;
			};
			break;
		  case 1:
			switch (mAxisIdx[opt.back()]) {
			  case 0: SetMin(h3->GetBinContent(i1, i0, i2)); break;
			  case 2: SetMin(h3->GetBinContent(i2, i0, i1)); break;
			};
			break;
		  case 2:
			switch (mAxisIdx[opt.back()]) {
			  case 0: SetMin(h3->GetBinContent(i1, i2, i0)); break;
			  case 1: SetMin(h3->GetBinContent(i2, i1, i0)); break;
			};
			break;
		} // END SWITCH

		// std::cout << min << std::endl;

		h2->SetBinContent(i0, i1, min);

	  } // END AX PROFILED

	} // END FOR i1

  } // END for i2

  return h2;

}

TCanvas *GetMapPlots(TH3D* hGrid, const std::vector<TVector3>& vPos,
					 const std::string& tag = "cGrid"){

  std::vector<std::vector<TMarker*>> vvM;
  vvM.reserve(vPos.size());

  std::size_t iM=0;
  std::vector<int> ColorPalette = {kBlue-4, kRed-4, kGreen+1};
  std::vector<int> StylePalette = {kFullCrossX, kOpenCrossX, kFullCross};

  for(auto& Pos: vPos){
	vvM.emplace_back(GetVMarker(Pos, ColorPalette[iM], StylePalette[iM]));
	iM++;
  }

  auto PlotMarker = [&vvM](const std::size_t idx){
	for(auto& vM: vvM)
	  vM[idx]->Draw("");
  };

  auto SetLegend = [](TH2D *hprof, const std::vector<std::string>& vs){
	hprof->GetXaxis()->SetName(Form("%s [mm]", vs[0].c_str()));
	hprof->GetYaxis()->SetName(Form("%s [mm]",  vs[1].c_str()));
  };

  auto c = new TCanvas(tag.c_str(), "", 1800, 600);
  c->Divide(3, 1);
  c->cd(1);
  // auto hXY = hGrid->Project3DProfile("xy");
  auto hXY = ProfileH3(hGrid, {"x", "y"});
  SetLegend(hXY, {"x", "z"});
  hXY->Draw("COLZ");
  PlotMarker(0);
  c->cd(2);
  // auto hXZ = hGrid->Project3DProfile("xz");
  auto hXZ = ProfileH3(hGrid, {"x", "z"});
  SetLegend(hXZ, {"x", "z"});
  hXZ->Draw("COLZ");
  PlotMarker(1);
  c->cd(3);
  // auto hYZ = hGrid->Project3DProfile("yz");
  auto hYZ = ProfileH3(hGrid, {"y", "z"});
  SetLegend(hYZ, {"y", "z"});
  hYZ->Draw("COLZ");
  PlotMarker(2);

  return c;

}

void HGridGarbageCollector(){
  for (auto &&obj: *gDirectory->GetList()) {
	if (!std::string(obj->GetName()).find("hGrid_")) {
	  // std::cout << obj->GetName() << std::endl;
	  delete obj;
	}
  }
}

#endif //SEEDNDESTROY_INCLUDE_RECORDER_HH_
