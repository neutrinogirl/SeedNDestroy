//
// Created by Stephane Zsoldos on 11/16/22.
//

#include "MakePDF.hh"

#include "SnD/ZAxis.hh"

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
				 const bool &applytrigger,
				 const std::vector<float>& vPosShift)
	: isShift(isshift), isPosShifted(false), isApplyTrigger(applytrigger) {
  const zAxis axTRes(TResBins, TResMin, TResMax);
  const zAxis axCosT(12, -1., 1.);
  const zAxis axNHits(1000, 0., 1000.);
  hNHits = new TH1D("hNHits", "NHits per event ; NHits ; ",
					axNHits.nBins, axNHits.min, axNHits.max);
  hNHits->SetDirectory(nullptr);
  hN400 = new TH1D("hN400", "N_{400} per event ; N_{400} ; ",
				   axNHits.nBins, axNHits.min, axNHits.max);
  hN400->SetDirectory(nullptr);
  std::vector<unsigned int> vPower = {0};
  vvHPDFs.reserve(vPower.size());
  for(const auto& wP : vPower){
	vvHPDFs.push_back(
		{
			new TH2D(Form("hCTVSTResPDF_THit_QW%d", wP), "T_{Res} VS Cos(#theta) ; T_{Res} [ns] ; Cos(#theta)",
					 axTRes.nBins, axTRes.min, axTRes.max,
					 axCosT.nBins, axCosT.min, axCosT.max),
			new TH2D(Form("hCTVSTResPDF_TTOF_QW%d", wP), "T_{Res} VS Cos(#theta) ; T_{Res} [ns] ; Cos(#theta)",
					 axTRes.nBins, axTRes.min, axTRes.max,
					 axCosT.nBins, axCosT.min, axCosT.max)
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
  if(isApplyTrigger){
	std::sort(vHits.begin(), vHits.end(), [](const Hit& a, const Hit& b){return a.T < b.T;});
	auto T0 = vHits.begin()->T;
	std::transform(vHits.begin(), vHits.end(), vHits.begin(),
				   [T0](Hit& hit){hit.T -= T0; return hit;});
  }

  TVector3 Pos = wData->GetPosition();
  if(isPosShifted){
	Pos.SetXYZ(wData->GetPosition().X(),
			   -1*(wData->GetPosition().Z()-PosShift.Z()),
			   wData->GetPosition().Y()-PosShift.Y());
  }
  auto Dir = wData->GetDirection();
  auto T = wData->GetTime();
  auto TrigTime = 0;

  std::sort(vHits.begin(), vHits.end());

  hN400->Fill(GetNPrompts(vHits, 400));
  hNHits->Fill(GetNPrompts(vHits, std::numeric_limits<double>::max()));

  std::vector<unsigned int> vPower = {0};
  enum { kTHIT, kTOF };

  for(auto iPower = 0; iPower<vPower.size(); iPower++){

	// Get PDF
	for(auto& hit: vHits){

	  const double QW = fWeight(hit, vPower[iPower]);

	  vvHPDFs[iPower][kTHIT]->Fill(hit.GetTRes(Pos, T),
								   hit.GetCosTheta(Pos, Dir),
								   QW);
	  vvHPDFs[iPower][kTOF]->Fill(hit.GetTRes(Pos, T-TrigTime),
								  hit.GetCosTheta(Pos, Dir),
								  QW);

	  if(!mPDFs[hit.ID]){
		const zAxis axTRes(vvHPDFs[iPower][kTOF]->GetXaxis());
		const zAxis axCosT(vvHPDFs[iPower][kTOF]->GetYaxis());
		mPDFs[hit.ID]=
			new TH2D(Form("hCTVSTResPDF_TTOF_QW%d_PMT%d", iPower, hit.ID),
					 "T_{Res} VS Cos(#theta) ; T_{Res} [ns] ; Cos(#theta)",
					 axTRes.nBins, axTRes.min, axTRes.max,
					 axCosT.nBins, axCosT.min, axCosT.max)
			;
	  } else {
		mPDFs[hit.ID]->Fill(hit.GetTRes(Pos, T-TrigTime),
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
	  if(isShift)
		ShiftHistogram(hPDF);
	  hPDF->Write();
	}
  }
  for(auto& m : mPDFs){
	if(isShift)
	  ShiftHistogram(m.second);
	m.second->Write();
  }
  f.Close();
}
