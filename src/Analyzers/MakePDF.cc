//
// Created by Stephane Zsoldos on 11/16/22.
//

#include "MakePDF.hh"

#include "SnD/ZAxis.hh"

#include "Templates/TData.hh"

void ShiftHistogram(TH2D* hist) {
  // Find the bin with the maximum value in the x-axis
  int maxBin = hist->GetMaximumBin();
  double maxValue = hist->GetBinContent(maxBin);
  double maxX = hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(maxValue));

  // Shift the histogram so that the maximum bin x-value is 0
  hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin() - maxX, hist->GetXaxis()->GetXmax() - maxX);
}


MakePDF::MakePDF(const unsigned int& TResBins, const float& TResMin, const float& TResMax){
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
}

void MakePDF::Do(void *Data) {

  auto* wData = static_cast<TData*>(Data);

  auto vHits = wData->GetVHits();
  auto Pos = wData->GetPosition();
  auto Dir = wData->GetDirection();
  auto T = wData->GetTime();
  auto TrigTime = T;

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
			new TH2D(Form("hCTVSTResPDF_TTOF_QW%d_PMT%d", iPower, hit.ID), "T_{Res} VS Cos(#theta) ; T_{Res} [ns] ; Cos(#theta)",
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

  for(auto iPower = 0; iPower<vPower.size(); iPower++){
	ShiftHistogram(vvHPDFs[iPower][kTHIT]);
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
	  hPDF->Write();
	}
  }
  for(auto& m : mPDFs){
	m.second->Write();
  }
  f.Close();
}
