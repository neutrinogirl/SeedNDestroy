//
// Created by Stephane Zsoldos on 7/10/22.
//

#include "ERecAnalysis.hh"
ERecAnalysis::ERecAnalysis() {

  vEdges = {0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0};

  for(const auto& E: vEdges){
	v2D.push_back(new TH2D(Form("h2D_E%.1f", E), Form("E=%.1f ; N_{hits} ; Q ", E), 50, 0., 2000., 50, 0., 2000.));
	v2DWall.push_back(new TH2D(Form("h2D_E%.1f_Wall", E), Form("E=%.1f ; N_{hits} ; Q ", E), 50, 0., 2000., 50, 0., 2000.));
  }

}
// 0 == underflow
// 1 == [0., 0.5[
// ...
// n == [n-1, n[
// size() == overflow
int ERecAnalysis::GetInterval(const double& E) {
  return std::upper_bound(vEdges.begin(), vEdges.end(), E) - vEdges.begin();
}

#include <SnD/RATData.hh>
void ERecAnalysis::Do(void *Data) {

  RATData *RData = static_cast<RATData*>(Data);

  const int EBin = GetInterval(RData->E);

  if(EBin > 0 && EBin <= vEdges.size()) {
	if(RData->Pos.Mag() < 3.e3)
	  v2D[EBin-1]->Fill(RData->NHits, RData->Q);
	else
	  v2DWall[EBin-1]->Fill(RData->NHits, RData->Q);
  }

}
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TF2.h>

// Calculate the error of a ratio of two double
static double GetErr(double X, double Y){
  return std::sqrt(X)/Y + X*std::sqrt(Y)/std::pow(Y, 2);
}

void ERecAnalysis::Export(const char *filename) {

  TGraphErrors gNHits;
  gNHits.SetName("gNHits");
  gNHits.SetLineWidth(2);
  gNHits.SetMarkerSize(2);
  gNHits.SetMarkerStyle(kPlus);
  gNHits.SetMarkerColor(kBlue-4);
  gNHits.SetLineColor(kBlue-4);
  TGraphErrors gQ;
  gQ.SetName("gQ");
  gQ.SetLineWidth(2);
  gQ.SetMarkerSize(2);
  gQ.SetMarkerStyle(kPlus);
  gQ.SetMarkerColor(kBlue-4);
  gQ.SetLineColor(kBlue-4);

  TFile *f = new TFile(filename, "RECREATE");
  for(auto i=0; i<vEdges.size(); i++){
	v2D[i]->Write();
	v2DWall[i]->Write();
	if(v2D[i]->GetEntries()>0){
	  TFitResultPtr r = v2D[i]->ProjectionX()->Fit("gaus", "SQ");
	  gNHits.SetPoint(gNHits.GetN(), vEdges[i], r->Parameter(2)/r->Parameter(1));
	  gNHits.SetPointError(gNHits.GetN()-1, 0., GetErr(r->Parameter(2), r->Parameter(1)));
	  r = v2D[i]->ProjectionY()->Fit("gaus", "SQ");
	  gQ.SetPoint(gQ.GetN(), vEdges[i], r->Parameter(2)/r->Parameter(1));
	  gQ.SetPointError(gQ.GetN()-1, 0., GetErr(r->Parameter(2), r->Parameter(1)));
	}
  }
  gNHits.Write();
  gQ.Write();
  f->Close();
}