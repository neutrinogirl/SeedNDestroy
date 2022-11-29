//
// Created by Stephane Zsoldos on 11/16/22.
//

#include "NTuple.hh"

#include <boost/range/combine.hpp>

Flat::Flat(std::vector< TTreeReader* > &vTreeReaders) {
  many = {
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
	  new TTreeReaderValue<std::vector<double>>(*vTreeReaders[kTree], "hitPMTCharge")
  };
  TTreeReaderValue<std::vector<int>> pmtId(*vTreeReaders[kMeta], "pmtId");
  TTreeReaderValue<std::vector<double>> pmtX(*vTreeReaders[kMeta], "pmtX");
  TTreeReaderValue<std::vector<double>> pmtY(*vTreeReaders[kMeta], "pmtY");
  TTreeReaderValue<std::vector<double>> pmtZ(*vTreeReaders[kMeta], "pmtZ");
  vTreeReaders[kMeta]->Next();
  for(auto tup : boost::combine(*pmtId, *pmtX, *pmtY, *pmtZ)) {
	int id;
	double x, y, z;
	boost::tie(id, x, y, z) = tup;
	mPMTPos[id] = TVector3(x, y, z);
  }

}

TVector3 Flat::GetPosition() {
  auto x = boost::any_cast<TTreeReaderValue<double>*>(many[0])->Get();
  auto y = boost::any_cast<TTreeReaderValue<double>*>(many[1])->Get();
  auto z = boost::any_cast<TTreeReaderValue<double>*>(many[2])->Get();
  return TVector3(*x,
				  *y,
				  *z);
}

TVector3 Flat::GetDirection() {
  auto x = boost::any_cast<TTreeReaderValue<double>*>(many[3])->Get();
  auto y = boost::any_cast<TTreeReaderValue<double>*>(many[4])->Get();
  auto z = boost::any_cast<TTreeReaderValue<double>*>(many[5])->Get();
  return TVector3(*x,
				  *y,
				  *z);
}

double Flat::GetEnergy() {
  return *boost::any_cast<TTreeReaderValue<double>*>(many[6])->Get();
}

double Flat::GetTime() {
  return 0;
}

int Flat::GetEventID() {
  return *boost::any_cast<TTreeReaderValue<int>*>(many[7])->Get();
}

int Flat::GetTriggerID() {
  return *boost::any_cast<TTreeReaderValue<int>*>(many[8])->Get();
}

std::vector<Hit> Flat::GetVHits() {
  std::vector<Hit> vHits;
  auto pmtId = boost::any_cast<TTreeReaderValue<std::vector<int>>*>(many[9]);
  auto pmtTime = boost::any_cast<TTreeReaderValue<std::vector<double>>*>(many[10]);
  auto pmtCharge = boost::any_cast<TTreeReaderValue<std::vector<double>>*>(many[11]);
  for (auto tup : boost::combine(**pmtId, **pmtTime, **pmtCharge)) {
	int ID;
	double T, Q;
	boost::tie(ID, T, Q) = tup;
	vHits.emplace_back(mPMTPos[ID], T, Q, ID);
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
  //
  iTrig = -1;
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

bool FlatReader::GetNextTrigger() {
  if (++iTrig < 1) {
	return true;
  } else {
	iTrig = -1;
	return false;
  }
}

void *FlatReader::GetData() {
  return fFlat;
}