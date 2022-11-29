//
// Created by Stephane Zsoldos on 11/16/22.
//

#include "Rat.hh"

RATReader::RATReader(const char *filename, const bool &verbose) {
  w_rat.ReadFile(filename);
  w_rat.Set();
  progress_bar_.Set(w_rat.GetNEvts(), 70);
  verbose_ = verbose;
}

bool RATReader::GetNextEvent() {
  w_rat.GetNextEvent();
}
bool RATReader::GetNextTrigger() {
  w_rat.GetNextTrigger();
}

RATData *RATReader::GetData() {

  data->GetVHits().clear();

  auto EV = w_rat.GetDS()->GetEV(w_rat.GetITrig());
  auto nPMTs = EV->GetPMTCount();

  for (auto iPMT = 0; iPMT < nPMTs; iPMT++) {

	auto PMT = EV->GetPMT(iPMT);
	auto ID = PMT->GetID();

	const auto PMTType = w_rat.GetRun()->GetPMTInfo()->GetType(ID);

	auto T = PMT->GetTime();
	auto Pos = w_rat.GetRun()->GetPMTInfo()->GetPosition(ID);
	auto QHit = PMT->GetCharge();
	Hit hit(Pos, QHit, T, ID);
	data->GetVHits().emplace_back(hit);

  }

  return data;
}
