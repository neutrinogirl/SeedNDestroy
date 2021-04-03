#include <iostream>
#include <vector>
#include <string>

#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include <MathUtils.hh>
#include <DrawingUtils.hh>
#include <TLine.h>

TGraph* GetInt(TH1D* h, const std::string& tag = ""){

  auto *gr = new TGraph();
  if(tag.empty())
	gr->SetName(Form("g_%s", h->GetName()));
  else
    gr->SetName(tag.c_str());

  const double NormInt = h->Integral();

  for(auto iBin=1; iBin<h->GetNbinsX()+1; iBin++){
	gr->SetPoint(gr->GetN(),
				 h->GetXaxis()->GetBinLowEdge(iBin),
				 h->Integral(0, iBin) / NormInt);
  }

  return gr;
}

void LLA(const std::string& histname = "nhits", const zAxis& Ax = {1000, 0., 10000.}){

  const std::string path = "outputs/electrons_5MeV_Fill_TrigThresh8_dWall_9000_35000_TTrigCut0_m50_p100/wGridSearch/";
  const std::string detname = "theia_90pct_9000_35000_wbls_3pct";

  const std::vector<std::string> vComponents {
	  "geoneutrino",
	  "reactorWorld",
	  "pmt_bi214*"
  };

  const std::vector<int> ColPalette {
    kBlue-4,
    kRed-4,
    kGreen+1
  };

  const std::string signame = vComponents[0];
  const std::vector< std::string > vSigNames = {signame};

  const int weight = 1;
  const std::string pdfweight = Form("_*_dWall_*_W%d_*.root", weight);

  std::map<std::string, TChain*> vCh;
  std::map<std::string, TH1D*> vHist;
  std::map<std::string, TGraph*> vGr;
  auto mg = new TMultiGraph("mg", "mg");

  int iCounter = 0;

  auto c =
	  new TCanvas(Form("c%s", histname.c_str()), Form("c%s", histname.c_str()), 800, 600);

  for(const auto& name:vComponents) {

    vHist[name] = new TH1D(Form("h_%s", name.c_str()),
						   Form("%s %s", name.c_str(), histname.c_str()),
						   Ax.nBins, Ax.min, Ax.max);
    SetBasicTStyle(vHist[name], ColPalette[iCounter]);

	std::string filenames = path + "/" + detname + "_" + name + pdfweight;
	vCh[name] = new TChain("T", "");
	vCh[name]->Add(filenames.c_str());

	std::cout << vCh[name]->GetNtrees() << std::endl;

	vCh[name]->Draw(Form("%s>>h_%s",histname.c_str(), name.c_str()), "itrig==0", "GOFF");
	vHist[name]->Scale(1./vCh[name]->GetNtrees());

	vGr[name] = GetInt(vHist[name]);
	SetBasicTStyle(vGr[name], ColPalette[iCounter]);

	mg->Add(vGr[name], "PL");

	iCounter++;
  }

  typedef std::pair<std::string, TH1D*> MyPairType;
  struct GetMin
  {
	bool operator()(const MyPairType& left, const MyPairType& right) const
	{
	  return left.second->GetMinimum() > right.second->GetMinimum();
	}
  };

  auto MaxHist = std::max_element(vHist.begin(), vHist.end(),GetMin());

  MaxHist->second->Draw("HIST");
  for(auto& pair: vHist){
    if(pair.first!=MaxHist->first)
      pair.second->Draw("HISTSAME");
  }

  c->BuildLegend();

  auto gDiscr = new TGraph();
  gDiscr->SetName(Form("gDiscr_%s", histname.c_str()));
  gDiscr->SetLineStyle(kDashed);
  gDiscr->SetLineWidth(2);
  const std::vector<double> vPts = Ax.GetStdVec();

  double maxDiscri = -std::numeric_limits<double>::max();
  double maxChi2Cut = -1;

  for(const auto& pt:vPts) {

	// BCKG
	double bckgfrac = 1.;
	for(const auto& name:vComponents){
	  for(const auto& sname:vSigNames){
		if(name != sname){
		  bckgfrac -= vGr[name]->Eval(pt);
		}
	  }
	}

	// SIG
	double sigfrac = 1-vGr[signame]->Eval(pt);

	// DEBUG
	// std::cout << bckgfrac << " " << sigfrac << std::endl;

	double disci = 0;
	const double discri = bckgfrac != 0 ? std::abs(sigfrac / bckgfrac) : bckgfrac;
	if(discri > maxDiscri){
	  maxDiscri = discri;
	  maxChi2Cut = pt;
	}

	gDiscr->SetPoint(gDiscr->GetN(), pt, discri);

  }

  // SetBasicTStyle(gDiscr, kBlack+1, 2, kDashed);

  // mg->Add(gDiscr, "C");


  auto lDiscri = new TLine(maxChi2Cut, 0., maxChi2Cut, MaxHist->second->GetMaximum());
  lDiscri->SetLineColor(kBlack);
  lDiscri->SetLineWidth(2);
  lDiscri->SetLineStyle(kDashed);

  lDiscri->Draw("SAME");


  c = new TCanvas(Form("cmg_%s", histname.c_str()), Form("cmg_%s", histname.c_str()), 800, 600);
  mg->Draw("A");
  c->BuildLegend();
  lDiscri->Draw("SAME");

  std::cout  << maxChi2Cut << " " << maxDiscri << std::endl;

}