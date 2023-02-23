//
// Created by Stephane Zsoldos on 11/16/22.
//

#include "MakePDF.hh"

#include "Templates/TData.hh"

void ShiftHistogram(TH2D* hist) {
  // find x position of maximum
  double x_max = hist->GetXaxis()->GetBinLowEdge(hist->ProjectionX()->GetMaximumBin());

  // get hist info
  int nbinsx = hist->GetXaxis()->GetNbins();
  int minx = hist->GetXaxis()->GetXmin();
  int maxx = hist->GetXaxis()->GetXmax();

  // shift the histogram
  hist->GetXaxis()->Set(nbinsx, minx-x_max, maxx-x_max);
}


MakePDF::MakePDF(const unsigned int& TResBins, const float& TResMin, const float& TResMax,
				 const bool &isshift,
				 const bool &applynorm,
				 const std::vector<float>& vPosShift)
	: isShift(isshift), isPosShifted(false), isApplyNorm(applynorm) {
  //
  const int NBinsNHits = 1000;
  const float MinNHits = 0.f;
  const float MaxNHits = 1000.f;
  const int NBinsCosT = 12;
  const float MinCosT = 0.f;
  const float MaxCosT = 1.f;
  //
  hNHits = new TH1D("hNHits", "NHits per event ; NHits ; ",
					NBinsNHits, MinNHits, MaxNHits);
  hNHits->SetDirectory(nullptr);
  //
  hN400 = new TH1D("hN400", "N_{400} per event ; N_{400} ; ",
				   NBinsNHits, MinNHits, MaxNHits);
  hN400->SetDirectory(nullptr);
  //
  std::vector<unsigned int> vPower = {0};
  vvHPDFs.reserve(vPower.size());
  for(const auto& wP : vPower){
	vvHPDFs.push_back(
		{
			new TH2D(Form("hCTVSTResPDF_TTOF_QW%d", wP), "T_{Res} VS Cos(#theta) ; T_{Res} [ns] ; Cos(#theta)",
					 TResBins, TResMin, TResMax,
					 NBinsCosT, MinCosT, MaxCosT)
		}
	);
  }
  // Check if vPosShift member are different from 0
  if(vPosShift[0] != 0 || vPosShift[1] != 0 || vPosShift[2] != 0){
	PosShift.SetXYZ(vPosShift[0], vPosShift[1], vPosShift[2]);
	isPosShifted = true;
  }
}

void MakePDF::Do(void *Data) {

  auto* wData = static_cast<TData*>(Data);

  auto vHits = wData->GetVHits();

  TVector3 Pos = wData->GetPosition();
  // Relative to ANNIE coordinate system
  if(isPosShifted){
	Pos.SetXYZ(wData->GetPosition().X(),
			   -1*(wData->GetPosition().Z()-PosShift.Z()),
			   wData->GetPosition().Y()-PosShift.Y());
  }
  TVector3 Dir = wData->GetDirection();
  double T     = 0.f;       // Timestamp of the event
  double TTrig = 0.f;       // Timestamp of the trigger
  double dT    = TTrig - T; // dT: Time between the vertex and the trigger
  // In real data, all hits are recorded with respect to the trigger time
  // Therefore, the time residuals must be corrected from the time between the event and the trigger.

  hN400->Fill(GetNPrompts(vHits, 400));
  hNHits->Fill(GetNPrompts(vHits, std::numeric_limits<double>::max()));

  std::vector<unsigned int> vPower = {0};
  enum { kTOF };

  for(auto iPower = 0; iPower<vPower.size(); iPower++){

	// Get PDF
	for(auto& hit: vHits){

	  const double QW = fWeight(hit, vPower[iPower]);

	  vvHPDFs[iPower][kTOF]->Fill(hit.GetTRes(Pos, dT),
								  hit.GetCosTheta(Pos, Dir),
								  QW);

	  if(!mPDFs[hit.ID]){
		mPDFs[hit.ID]=
			new TH2D(Form("hCTVSTResPDF_TTOF_QW%d_PMT%d", iPower, hit.ID),
					 "T_{Res} VS Cos(#theta) ; T_{Res} [ns] ; Cos(#theta)",
					 vvHPDFs[iPower][kTOF]->GetXaxis()->GetNbins(),
					 vvHPDFs[iPower][kTOF]->GetXaxis()->GetXmin(), vvHPDFs[iPower][kTOF]->GetXaxis()->GetXmax(),
					 vvHPDFs[iPower][kTOF]->GetYaxis()->GetNbins(),
					 vvHPDFs[iPower][kTOF]->GetYaxis()->GetXmin(), vvHPDFs[iPower][kTOF]->GetYaxis()->GetXmax());
	  } else {
		mPDFs[hit.ID]->Fill(hit.GetTRes(Pos, dT),
							hit.GetCosTheta(Pos, Dir),
							QW);
	  }

	}

  }

}

#include <TFile.h>
void MakePDF::Export(const char *filename) {

  TFile f(filename, "RECREATE");
  f.cd();
  hNHits->Write();
  hN400->Write();
  for(auto& vHPDF : vvHPDFs){
	for(auto& hPDF : vHPDF){
	  if(isApplyNorm)
		hPDF->Scale(1./hPDF->Integral());
	  if(isShift)
		ShiftHistogram(hPDF);
	  hPDF->Write();
	}
  }
  for(auto& m : mPDFs){
	if(isApplyNorm)
	  m.second->Scale(1./m.second->Integral());
	if(isShift)
	  ShiftHistogram(m.second);
	m.second->Write();
  }
  f.Close();
}
