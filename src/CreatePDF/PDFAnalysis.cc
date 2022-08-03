//
// Created by Stephane Zsoldos on 7/6/22.
//

#include "PDFAnalysis.hh"

#include "SnD/ZAxis.hh"

Analysis::Analysis(const unsigned int& TResBins, const float& TResMin, const float& TResMax){
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

#include "SnD/RATData.hh"
void Analysis::Do(void *Data) {

  auto *RData = static_cast<RATData*>(Data);

  std::sort(RData->vHits.begin(), RData->vHits.end());

  hN400->Fill(GetNPrompts(RData->vHits, 400));
  hNHits->Fill(GetNPrompts(RData->vHits, std::numeric_limits<double>::max()));

  std::vector<unsigned int> vPower = {0};
  enum { kTHIT, kTOF };

  for(auto iPower = 0; iPower<vPower.size(); iPower++){

	// Get PDF
	for(auto& hit: RData->vHits){

	  const double QW = fWeight(hit, vPower[iPower]);

	  vvHPDFs[iPower][kTHIT]->Fill(hit.GetTRes(RData->Pos, RData->T),
								   hit.GetCosTheta(RData->Pos, RData->Dir),
								   QW);
	  vvHPDFs[iPower][kTOF]->Fill(hit.GetTRes(RData->Pos, RData->T-RData->TrigTime),
								  hit.GetCosTheta(RData->Pos, RData->Dir),
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
		mPDFs[hit.ID]->Fill(hit.GetTRes(RData->Pos, RData->T-RData->TrigTime),
							hit.GetCosTheta(RData->Pos, RData->Dir),
							QW);
	  }

	}

  }

}
#include <TFile.h>
void Analysis::Export(const char *filename) {
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
