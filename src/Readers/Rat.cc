//
// Created by Stephane Zsoldos on 11/16/22.
//

#include "Rat.hh"

void RATData::Update(wRAT *w_rat) {
  Pos = w_rat->GetPosTrue(0);
  Dir = w_rat->GetDirTrue(0);
  Energy = w_rat->GetETrue(0);
  Time = w_rat->GetTTrue(0);

  vHits.clear();
  auto EV = w_rat->GetDS()->GetEV(w_rat->GetITrig());
  auto nPMTs = EV->GetPMTCount();

  for (auto iPMT = 0; iPMT < nPMTs; iPMT++) {

	auto PMT = EV->GetPMT(iPMT);
	auto ID = PMT->GetID();

	const auto PMTType = w_rat->GetRun()->GetPMTInfo()->GetType(ID);

	auto T = PMT->GetTime();
	auto Pos = w_rat->GetRun()->GetPMTInfo()->GetPosition(ID);
	auto QHit = PMT->GetCharge();
	Hit hit(Pos, QHit, T, ID);
	vHits.emplace_back(hit);

  }

  EventID = w_rat->GetEvt();
  TriggerID = w_rat->GetITrig();

}

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

void *RATReader::GetData() {

  data->Update(&w_rat);
  return data;
}
