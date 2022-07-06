//
// Created by Stephane Zsoldos on 7/5/22.
//

#include "RATReader.hh"

RATReader::RATReader(const char *filename, const bool &verbose) {
  w_rat.ReadFile(filename);
  w_rat.Set();
  progress_bar_.Set(w_rat.GetNEvts(), 70);
  verbose_ = verbose;
  d = new RATData();
}
bool RATReader::GetNextEvent() {
  w_rat.GetNextEvent();
}
bool RATReader::GetNextTrigger() {
  w_rat.GetNextTrigger();
}
void *RATReader::GetData() {

  d->Clear();
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
