//
// Created by Stephane Zsoldos on 8/2/22.
//

#include <vector>
#include <map>

#include <SnD/Map.hh>
#include <SnD/GridWalker.hh>
#include <SnD/NLL.hh>

#include <TH3D.h>
#include <TH2D.h>
#include <TAxis.h>

//
static std::vector<double> GetSmplPts(double a=0.f, double b=1.f, int n=2){
  if(n<=2){
	return {a, b};
  } else {
	std::vector<double> x(n+1, 0.f);
	std::iota(x.begin(), x.end(), 0.f);
	std::transform(x.begin(), x.end(), x.begin(), [&](const double& val){return a + val*(b-a)/n;});
	return x;
  }
}

//
static std::string GetProfileAxis(const std::vector<std::string>& opt = {"x", "y"}){
  std::vector<std::string> vAxis = { "x", "y", "z" };

  const auto itIdx = std::find(vAxis.begin(), vAxis.end(), opt.front());
  vAxis.erase(itIdx);
  const auto itJdx = std::find(vAxis.begin(), vAxis.end(), opt.back());
  vAxis.erase(itJdx);

  return *vAxis.begin();
}
void SetLegend(TH2D *hprof, const std::vector<std::string>& vs){
  hprof->GetXaxis()->SetName(Form("%s [mm]", vs[0].c_str()));
  hprof->GetYaxis()->SetName(Form("%s [mm]", vs[1].c_str()));
};

// Profile along axis
static TH2D* ProfileH3(TH3D* h3,
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

//
std::vector< TCanvas *> GetMap(const std::vector<Hit> &vHits, TH1D *hPDF, Bnd *b){

  const int nT = 10, nAx = 10;

  std::vector<double> vT = GetSmplPts(0.f, b->GetTEdge(), nT);

  std::vector< std::vector<double> > vAxs = {
	  GetSmplPts(-b->GetEdge().x(), b->GetEdge().x(), nAx),
	  GetSmplPts(-b->GetEdge().y(), b->GetEdge().y(), nAx),
	  GetSmplPts(-b->GetEdge().z(), b->GetEdge().z(), nAx),
  };

  std::vector< TCanvas *> vCanvas;
  vCanvas.reserve(vT.size());

  for(const auto &t : vT){

	Grid G(vAxs);

	TH3D h(Form("hGrid_%f", t), "",
		   nAx, &vAxs[0][0],
		   nAx, &vAxs[1][0],
		   nAx, &vAxs[2][0]);
	h.SetDirectory(nullptr);

	std::vector<double> x(3, 0.f);

	while(G.Walk()){
	  G.GetGridPt(x);
	  double NLL = GetNLL(vHits, hPDF, TVector3(x[0], x[1], x[2]), t);
	  h.Fill(x[0], x[1], x[2], NLL);
	}

	vCanvas.push_back( new TCanvas(Form("c_%.1f",t), "", 1800, 600) );
	vCanvas.back()->Divide(3, 1);
	vCanvas.back()->cd(1);
	auto hXY = ProfileH3(&h, {"x", "y"});
	SetLegend(hXY, {"x", "y"});
	hXY->Draw("COLZ");
	vCanvas.back()->cd(2);
	auto hXZ = ProfileH3(&h, {"x", "z"});
	SetLegend(hXZ, {"x", "z"});
	hXZ->Draw("COLZ");
	vCanvas.back()->cd(3);
	auto hYZ = ProfileH3(&h, {"y", "z"});
	SetLegend(hYZ, {"y", "z"});
	hYZ->Draw("COLZ");

  }

  return vCanvas;

}