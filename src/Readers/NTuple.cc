//
// Created by Stephane Zsoldos on 11/16/22.
//

#include "NTuple.hh"

// #include <boost/range/combine.hpp>
// MetaNTuple::MetaNTuple(TTreeReader *Reader) {
//   TTreeReaderValue<std::vector<int>> pmtId(*Reader, "pmtId");
//   TTreeReaderValue<std::vector<double>> pmtX(*Reader, "pmtX");
//   TTreeReaderValue<std::vector<double>> pmtY(*Reader, "pmtY");
//   TTreeReaderValue<std::vector<double>> pmtZ(*Reader, "pmtZ");
//   Reader->Next();
//   for(auto tup : boost::combine(*pmtId, *pmtX, *pmtY, *pmtZ)) {
// 	int id;
// 	double x, y, z;
// 	boost::tie(id, x, y, z) = tup;
// 	mPMTPos[id] = TVector3(x, y, z);
//   }
// }
//
// TVector3 MetaNTuple::GetPMTPosition(const int& PMTId) {
//   return mPMTPos[PMTId];
// }
//
// MetaNTuple::~MetaNTuple() {
//   mPMTPos.clear();
// }
//
// // #### #### #### #### #### #### #### #### #### #### #### #### //
// NTuple::NTuple(TTreeReader *Reader, TTreeReader *mReader) {
//   SetReader(Reader);
//   meta = new MetaNTuple(mReader);
// }
//
// NTuple::~NTuple() {
//   delete hitPMTID, delete hitPMTTime, delete hitPMTCharge;
// }
//
// void NTuple::SetReader(TTreeReader *Reader) {
//   hitPMTID = new TTreeReaderValue<std::vector<int>>(*Reader, "hitPMTID");
//   hitPMTTime = new TTreeReaderValue<std::vector<double>>(*Reader, "hitPMTTime");
//   hitPMTCharge = new TTreeReaderValue<std::vector<double>>(*Reader, "hitPMTCharge");
// }
//
// std::vector<Hit> NTuple::GetVHits() {
//   std::vector<Hit> vHits;
//   for (auto tup :
// 	  boost::combine(**hitPMTID, **hitPMTTime, **hitPMTCharge)) {
// 	int ID;
// 	double T, Q;
// 	boost::tie(ID, T, Q) = tup;
// 	vHits.emplace_back(meta->GetPMTPosition(ID), T, Q, ID);
//   }
//   return vHits;
// }

Flat::Flat(TTreeReader *Reader) {
  many = {
	  TTreeReaderValue<double>(*Reader, "mcx"),
	  TTreeReaderValue<double>(*Reader, "mcy"),
	  TTreeReaderValue<double>(*Reader, "mcz"),
	  TTreeReaderValue<double>(*Reader, "mcu"),
	  TTreeReaderValue<double>(*Reader, "mcv"),
	  TTreeReaderValue<double>(*Reader, "mcw"),
	  TTreeReaderValue<double>(*Reader, "mcke"),
	  TTreeReaderValue<int>(*Reader, "evid"),
	  TTreeReaderValue<int>(*Reader, "subev")
  };
}
TVector3 Flat::GetPosition() {
  return {*boost::any_cast<TTreeReaderValue<double>>(many[0]),
		  *boost::any_cast<TTreeReaderValue<double>>(many[1]),
		  *boost::any_cast<TTreeReaderValue<double>>(many[2])};
}

TVector3 Flat::GetDirection() {
  return {*boost::any_cast<TTreeReaderValue<double>>(many[3]),
		  *boost::any_cast<TTreeReaderValue<double>>(many[4]),
		  *boost::any_cast<TTreeReaderValue<double>>(many[5])};
}

double Flat::GetEnergy() {
  return *boost::any_cast<TTreeReaderValue<double>>(many[6]);
}

int Flat::GetEventID() {
  return *boost::any_cast<TTreeReaderValue<int>>(many[7]);
}

int Flat::GetSubEventID() {
  return *boost::any_cast<TTreeReaderValue<int>>(many[8]);
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
  fFlat = new Flat(vTreeReaders[kTree]);
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