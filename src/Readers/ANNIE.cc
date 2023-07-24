//
// Created by Stephane Zsoldos on 7/24/23.
//

#include "ANNIE.hh"

#include <range/v3/all.hpp>

ANNIE::ANNIE(std::vector<TTreeReader *> &vTreeReaders) : many(
	{
	  new TTreeReaderValue<std::vector<double>>(*vTreeReaders[0], "hitX"),
	  new TTreeReaderValue<std::vector<double>>(*vTreeReaders[0], "hitY"),
	  new TTreeReaderValue<std::vector<double>>(*vTreeReaders[0], "hitZ"),
	  new TTreeReaderValue<std::vector<double>>(*vTreeReaders[0], "hitT"),
	  new TTreeReaderValue<std::vector<double>>(*vTreeReaders[0], "hitQ"),
	  new TTreeReaderValue<std::vector<int>>(*vTreeReaders[0], "hitDetID"),
	  new TTreeReaderValue<unsigned long long>(*vTreeReaders[0], "eventTimeTank"),
	  new TTreeReaderValue<int>(*vTreeReaders[0], "eventNumber")
	}
){}

double ANNIE::GetTime() {
  return *std::any_cast<TTreeReaderValue<unsigned long long>*>(many[6])->Get();
}

int ANNIE::GetEventID() {
  return *std::any_cast<TTreeReaderValue<int>*>(many[7])->Get();
}

std::vector<Hit> ANNIE::GetVHits() {
  std::vector<Hit> vHits;
  auto hitX = std::any_cast<TTreeReaderValue<std::vector<double>>*>(many[0]);
  auto hitY = std::any_cast<TTreeReaderValue<std::vector<double>>*>(many[1]);
  auto hitZ = std::any_cast<TTreeReaderValue<std::vector<double>>*>(many[2]);
  auto hitT = std::any_cast<TTreeReaderValue<std::vector<double>>*>(many[3]);
  auto hitQ = std::any_cast<TTreeReaderValue<std::vector<double>>*>(many[4]);
  auto hitDetID = std::any_cast<TTreeReaderValue<std::vector<int>>*>(many[5]);
  auto combinedRange = ranges::views::zip(**hitX, **hitY, **hitZ, **hitT, **hitQ, **hitDetID);
  for (auto tup : combinedRange) {
	int ID;
	double X, Y, Z, T, Q;
	std::tie(X, Y, Z, T, Q, ID) = tup;
	vHits.emplace_back(Vector3(X, Y, Z, SpaceUnit::mm), Q, T * 1.e-3, ID);
  }
  return vHits;
}

// #### #### #### #### #### #### #### #### #### #### #### #### //
ANNIEReader::ANNIEReader(const char *filename,
						 const char *treename,
						 const bool &verbose) {
  //
  f = new TFile(filename);
  vTreeReaders = {
	  new TTreeReader(treename, f),
  };
  //
  fFlat = new ANNIE(vTreeReaders);
  //
  progress_bar_.Set(vTreeReaders[0]->GetEntries(), 70);
  verbose_ = verbose;
}

ANNIEReader::~ANNIEReader() {
  delete f;
  delete vTreeReaders[0];
  delete fFlat;
}

bool ANNIEReader::GetNextEvent() {
  return vTreeReaders[0]->Next();
}

void *ANNIEReader::GetData() {
  return fFlat;
}