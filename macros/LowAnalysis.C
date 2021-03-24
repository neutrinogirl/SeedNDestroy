#include <iostream>
#include <vector>
#include <string>

#include <TChain.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include <MathUtils.hh>
#include <DrawingUtils.hh>
#include <TLine.h>

TGraph *GetInt(const std::string& tag, TH2D* h){

  const double nEvts = h->GetEntries();

  auto g = new TGraph();
  g->SetName(Form("g_%s", tag.c_str()));
  g->SetTitle(Form("Fraction of events kept ; Chi^{2} cut ; Frac [%%] "));

  for(auto iBin=1; iBin<h->GetNbinsX()+1; iBin++){
	const double nevt = h->ProjectionX("hProj", 1, iBin+1)->GetEntries();
	g->SetPoint(g->GetN(), h->GetXaxis()->GetBinLowEdge(iBin+1), nevt / nEvts);
  }

  return g;

}

int main(){

  const std::string path = "outputs/electrons_5MeV_Fill_TrigThresh8_dWall_9000_35000_TTrigCut0_m50_p100";
  const std::string detname = "theia_90pct_9000_35000_wbls_3pct";

  const std::vector<std::string> vComponents {
	  "geoneutrino",
	  "pmt_*"
	  // "pmt_bi214",
	  // "pmt_tl208",
	  // "pmt_k40",
  };

  const std::string signame = vComponents[0];

  const int weight = 1;
  const std::string pdfweight = Form("_*_dWall_*_W%d_*.root", weight);

  std::map<std::string, TChain*> vCh;

  zAxis AxChi2 = {50, 0., 3.};
  zAxis AxRec = {50, -1.e4, 1.e4};

  const std::vector<std::string> vAxNames = {
	  "x",
	  "y",
	  "z"
  };

  std::map< std::string, std::map<std::string, TGraph*> > mG;

  for(const auto& axName: vAxNames){

	mG.insert(
		std::make_pair(axName, std::map<std::string, TGraph*>())
	);

    auto c =
		new TCanvas(Form("c%s", axName.c_str()), Form("c%s", axName.c_str()), 800, 600);

    auto mg = new TMultiGraph(Form("mg_%s", axName.c_str()),
							  Form("Fraction of events kept #vec{%s} ; Chi^{2} cut ; Frac [%%] ", axName.c_str()));

	std::vector<TH2D*> vHPerfs;
	vHPerfs.reserve(vComponents.size());
	for(const auto& name:vComponents){
	  vHPerfs.emplace_back(
		  new TH2D(Form("h_%s", name.c_str()),
				   Form("Rec - True %s VS Chi^{2} ; Chi^{2} ; ", name.c_str()),
				   AxChi2.nBins, AxChi2.min, AxChi2.max,
				   AxRec.nBins, AxRec.min, AxRec.max)

	  );
	  std::string filenames = path + "/" + detname + "_" + name + pdfweight;
	  vCh[name] = new TChain("T", "");
	  vCh[name]->Add(filenames.c_str());
	  vCh[name]->Draw(Form("mc%s-rec%s:chi2>>h_%s",axName.c_str(), axName.c_str(), name.c_str()), "", "GOFF");

	  auto g = GetInt(axName+name, vHPerfs.back());
	  if(name == signame)
		SetBasicTStyle(g, kBlue-4, 2, kSolid);
	  else
		SetBasicTStyle(g, kRed-4, 2, kSolid);

	  mg->Add(g, "PC");
	  mG[axName][name] = g;

	}

	auto gDiscr = new TGraph();
	gDiscr->SetName(Form("gDiscr_%s", axName.c_str()));
	const std::vector<double> vPts = AxChi2.GetStdVec();

	double maxDiscri = -std::numeric_limits<double>::max();
	double maxChi2Cut = -1;

	for(const auto& pt:vPts){
	  double bckgfrac = 0;
	  for(const auto& name:vComponents){
	    if(name != signame)
	      bckgfrac += mG[axName][name]->Eval(pt);
	  }
	  const double discri = mG[axName][signame]->Eval(pt) / bckgfrac;
	  if(discri > maxDiscri){
	    maxDiscri = discri;
	    maxChi2Cut = pt;
	  }
	  gDiscr->SetPoint(gDiscr->GetN(), pt, discri);
	}

	SetBasicTStyle(gDiscr, kGreen+1, 2, kDashed);

	mg->Add(gDiscr, "PC");
	mg->Draw("A");

	auto lDiscri = new TLine(maxChi2Cut, 0, maxChi2Cut, maxDiscri);
	lDiscri->SetLineColor(kBlack);
	lDiscri->SetLineWidth(2);
	lDiscri->SetLineStyle(kDashed);

	lDiscri->Draw("SAME");

	std::cout << axName << " " << maxChi2Cut << " " << maxDiscri << std::endl;

	c->Print(Form("%sW%d.png",c->GetName(), weight), "png");

  }

  return EXIT_SUCCESS;
}