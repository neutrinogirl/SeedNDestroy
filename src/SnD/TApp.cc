//
// Created by Stephane Zsoldos on 7/3/22.
//

#include "TApp.hh"
#include "ZAxis.hh"

#include <SnD/Hit.hh>

double GetNPrompts(const std::vector<Hit>& vHits, const double& T){
  double NPrompts = 0.;
  for(auto& hit: vHits){
	if(hit.T<T){
	  NPrompts++;
	}
  }
  return NPrompts;
};
double fWeight(const Hit& h, const int& P){
  return std::pow(h.Q, P);
}

class RATData {
 public:
  double TrigTime = 0.f;
  TVector3 Pos = TVector3(0.f, 0.f, 0.f);
  TVector3 Dir = TVector3(0.f, 0.f, 0.f);
  double T = 0.f;
  double E = 0.f;
  std::vector<Hit> vHits;
};

Analysis::Analysis(const TAppArgs &args) {
  const zAxis axTRes(args.GetTResBins()[0], args.GetTResBins()[1], args.GetTResBins()[2]);
  const zAxis axCosT(12, -1., 1.);
  const zAxis axNHits(1000, 0., 1000.);
  hNHits = new TH1D("hNHits", "NHits per event ; NHits ; ",
					axNHits.nBins, axNHits.min, axNHits.max);
  hNHits->SetDirectory(nullptr);
  hN400 = new TH1D("hN400", "N_{400} per event ; N_{400} ; ",
				   axNHits.nBins, axNHits.min, axNHits.max);
  hN400->SetDirectory(nullptr);
  std::vector<unsigned int> vPower = {0, 1, 2};
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
void Analysis::Do(void *Data) {

  auto *RData = reinterpret_cast<RATData*>(Data);

  std::sort(RData->vHits.begin(), RData->vHits.end());

  hN400->Fill(GetNPrompts(RData->vHits, 400));
  hNHits->Fill(GetNPrompts(RData->vHits, std::numeric_limits<double>::max()));

  std::vector<unsigned int> vPower = {0, 1, 2};
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

	}

  }

}

RATReader::RATReader(const TAppArgs &args) {
  w_rat.ReadFile(args.GetInput());
  w_rat.Set();
  progress_bar_.Set(w_rat.GetNEvts(), 70);
  verbose_ = args.GetVerbose();
}
bool RATReader::GetNextEvent() {
  w_rat.GetNextEvent();
}
bool RATReader::GetNextTrigger() {
  w_rat.GetNextTrigger();
}
void *RATReader::GetData() {

  auto *d = new RATData();
  w_rat.GetPrimaryParticleInfo(d->TrigTime, d->Pos, d->Dir, d->E, d->T);
  auto EV = w_rat.GetDS()->GetEV(w_rat.GetITrig());
  auto nPMTs = EV->GetPMTCount();

  for (auto iPMT = 0; iPMT < nPMTs; iPMT++) {

	auto PMT = EV->GetPMT(iPMT);
	auto ID = PMT->GetID();

	const auto PMTType = w_rat.GetRun()->GetPMTInfo()->GetType(ID);

	if (PMTType == 1) {
	  auto T = PMT->GetTime();
	  auto Pos = w_rat.GetRun()->GetPMTInfo()->GetPosition(ID);
	  auto QHit = PMT->GetCharge();
	  Hit hit(Pos, QHit, T);
	  d->vHits.emplace_back(hit);
	}

  }

  return d;
}

int main(int argc, char **argv) {

  // ######################################## //
  // Read arguments
  TAppArgs Args;
  Args.ProcessArgs(argc, argv);

  // ######################################## //
  // Create analysis class
  Analysis Ana(Args);

  // ######################################## //
  // Run analysis
  RATReader R(Args);
  R.Read(&Ana);


  return EXIT_SUCCESS;
}