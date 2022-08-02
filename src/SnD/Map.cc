//
// Created by Stephane Zsoldos on 8/2/22.
//

#include <SnD/Map.hh>
#include <SnD/GridWalker.hh>
#include <SnD/NLL.hh>

#include <TH3D.h>

std::vector<double> GetSmplPts(double a=0.f, double b=1.f, int n=2){
  if(n<2){
	return {a, b};
  } else {
	std::vector<double> x(n+1, 0.f);
	std::iota(x.begin(), x.end(), 0.f);
	std::transform(x.begin(), x.end(), x.begin(), [&](const double& val){return a + val*(b-a)/n;});
	return x;
  }
}

void Map(const std::vector<Hit> &vHits, TH1D *hPDF, Bnd *b){

  std::vector<double> vT = GetSmplPts(0.f, b->GetTEdge(), 10);

  std::vector< std::vector<double> > vAxs = {
	  GetSmplPts(-b->GetEdge().x(), b->GetEdge().x(), 10),
	  GetSmplPts(-b->GetEdge().y(), b->GetEdge().y(), 10),
	  GetSmplPts(-b->GetEdge().z(), b->GetEdge().z(), 10),
  };

  for(const auto &t : vT){

	Grid G(vAxs);

	TH3D *h = new TH3D(Form("h_%f", t), "",
					   vAxs[0].size(), vAxs[0].data(),
					   vAxs[1].size(), vAxs[1].data(),
					   vAxs[2].size(), vAxs[2].data());

	std::vector<double> x(3, 0.f);
	while(G.Walk()){
	  G.GetGridPt(x);
	  double NLL = GetNLL(vHits, hPDF, TVector3(x[0], x[1], x[2]), t);
	  h->Fill(x[0], x[1], x[2], NLL);
	}

  }

}