//
// Created by Stephane Zsoldos on 11/16/22.
//

#include "NTuple.hh"

#include <range/v3/all.hpp>

Flat::Flat(std::vector< TTreeReader* > &vTreeReaders) : many(
	{
		new TTreeReaderValue<double>(*vTreeReaders[kTree], "mcx"),
		new TTreeReaderValue<double>(*vTreeReaders[kTree], "mcy"),
		new TTreeReaderValue<double>(*vTreeReaders[kTree], "mcz"),
		new TTreeReaderValue<double>(*vTreeReaders[kTree], "mcu"),
		new TTreeReaderValue<double>(*vTreeReaders[kTree], "mcv"),
		new TTreeReaderValue<double>(*vTreeReaders[kTree], "mcw"),
		new TTreeReaderValue<double>(*vTreeReaders[kTree], "mcke"),
		new TTreeReaderValue<int>(*vTreeReaders[kTree], "evid"),
		new TTreeReaderValue<int>(*vTreeReaders[kTree], "subev"),
		new TTreeReaderValue<std::vector<int>>(*vTreeReaders[kTree], "hitPMTID"),
		new TTreeReaderValue<std::vector<double>>(*vTreeReaders[kTree], "hitPMTTime"),
		new TTreeReaderValue<std::vector<double>>(*vTreeReaders[kTree], "hitPMTCharge"),
		new TTreeReaderValue<double>(*vTreeReaders[kTree], "triggerTime")
	}
){
  TTreeReaderValue<std::vector<int>> pmtId(*vTreeReaders[kMeta], "pmtId");
  TTreeReaderValue<std::vector<double>> pmtX(*vTreeReaders[kMeta], "pmtX");
  TTreeReaderValue<std::vector<double>> pmtY(*vTreeReaders[kMeta], "pmtY");
  TTreeReaderValue<std::vector<double>> pmtZ(*vTreeReaders[kMeta], "pmtZ");
  vTreeReaders[kMeta]->Next();
  auto combinedRange = ranges::views::zip(*pmtId, *pmtX, *pmtY, *pmtZ);
  for(auto tup : combinedRange) {
	int id;
	double x, y, z;
	std::tie(id, x, y, z) = tup;
	mPMTPos[id] = Vector3(x, y, z, SpaceUnit::mm);
  }
}


Vector3 Flat::GetPosition() {
  auto x = std::any_cast<TTreeReaderValue<double>*>(many[0])->Get();
  auto y = std::any_cast<TTreeReaderValue<double>*>(many[1])->Get();
  auto z = std::any_cast<TTreeReaderValue<double>*>(many[2])->Get();
  return {*x, *y, *z, SpaceUnit::mm};
}

Vector3 Flat::GetDirection() {
  auto x = std::any_cast<TTreeReaderValue<double>*>(many[3])->Get();
  auto y = std::any_cast<TTreeReaderValue<double>*>(many[4])->Get();
  auto z = std::any_cast<TTreeReaderValue<double>*>(many[5])->Get();
  return {*x, *y, *z, SpaceUnit::mm};
}

double Flat::GetEnergy() {
  return *std::any_cast<TTreeReaderValue<double>*>(many[6])->Get();
}

double Flat::GetTime() {
  return *std::any_cast<TTreeReaderValue<double>*>(many[12])->Get();
}

int Flat::GetEventID() {
  return *std::any_cast<TTreeReaderValue<int>*>(many[7])->Get();
}

int Flat::GetTriggerID() {
  return *std::any_cast<TTreeReaderValue<int>*>(many[8])->Get();
}

std::vector<Hit> Flat::GetVHits() {
  std::vector<Hit> vHits;
  auto pmtId = std::any_cast<TTreeReaderValue<std::vector<int>>*>(many[9]);
  auto pmtTime = std::any_cast<TTreeReaderValue<std::vector<double>>*>(many[10]);
  auto pmtCharge = std::any_cast<TTreeReaderValue<std::vector<double>>*>(many[11]);
  auto combinedRange = ranges::views::zip(**pmtId, **pmtTime, **pmtCharge);
  for (auto tup : combinedRange) {
	int ID;
	double T, Q;
	std::tie(ID, T, Q) = tup;
	vHits.emplace_back(mPMTPos[ID], Q, T, ID);
  }
  return vHits;
}

// #### #### #### #### #### #### #### #### #### #### #### #### //
FlatReader::FlatReader(const char *filename,
					   const char *treename, const char *metaname,
					   const bool &verbose) {
  //
  f = new TFile(filename);
  vTreeReaders = {
	  new TTreeReader(treename, f),
	  new TTreeReader(metaname, f)
  };
  //
  fFlat = new Flat(vTreeReaders);
  //
  progress_bar_.Set(vTreeReaders[kTree]->GetEntries(), 70);
  verbose_ = verbose;
}

FlatReader::~FlatReader() {
  delete f;
  delete vTreeReaders[kTree];
  delete vTreeReaders[kMeta];
  delete fFlat;
}

bool FlatReader::GetNextEvent() {
  return vTreeReaders[kTree]->Next();
}

void *FlatReader::GetData() {
  return fFlat;
}