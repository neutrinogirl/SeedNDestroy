//
// Created by Stephane Zsoldos on 11/16/22.
//

#include "NTuple.hh"

#include <boost/range/combine.hpp>

MetaNTuple::MetaNTuple(TTreeReader *Reader) {
  TTreeReaderValue<std::vector<int>> pmtId(*Reader, "pmtId");
  TTreeReaderValue<std::vector<double>> pmtX(*Reader, "pmtX");
  TTreeReaderValue<std::vector<double>> pmtY(*Reader, "pmtY");
  TTreeReaderValue<std::vector<double>> pmtZ(*Reader, "pmtZ");
  Reader->Next();
  for(auto tup : boost::combine(*pmtId, *pmtX, *pmtY, *pmtZ)) {
	int id;
	double x, y, z;
	boost::tie(id, x, y, z) = tup;
	mPMTPos[id] = TVector3(x, y, z);
  }
}

TVector3 MetaNTuple::GetPMTPosition(const int& PMTId) {
  return mPMTPos[PMTId];
}

MetaNTuple::~MetaNTuple() {
  mPMTPos.clear();
}

// #### #### #### #### #### #### #### #### #### #### #### #### //
NTuple::NTuple(TTreeReader *Reader, TTreeReader *mReader) {
  SetReader(Reader);
  meta = new MetaNTuple(mReader);
}

NTuple::~NTuple() {
  delete hitPMTID, delete hitPMTTime, delete hitPMTCharge;
}

void NTuple::SetReader(TTreeReader *Reader) {
  hitPMTID = new TTreeReaderValue<std::vector<int>>(*Reader, "hitPMTID");
  hitPMTTime = new TTreeReaderValue<std::vector<double>>(*Reader, "hitPMTTime");
  hitPMTCharge = new TTreeReaderValue<std::vector<double>>(*Reader, "hitPMTCharge");
}

std::vector<Hit> NTuple::GetVHits() {
  std::vector<Hit> vHits;
  for (auto tup :
	  boost::combine(**hitPMTID, **hitPMTTime, **hitPMTCharge)) {
	int ID;
	double T, Q;
	boost::tie(ID, T, Q) = tup;
	vHits.emplace_back(meta->GetPMTPosition(ID), T, Q, ID);
  }
  return vHits;
}

// #### #### #### #### #### #### #### #### #### #### #### #### //
NTupleReader::NTupleReader(const char *filename,
						   const char *treename, const char *metaname,
						   const bool &verbose) {
  f = new TFile(filename);
  t = new TTreeReader(treename, f);
  m = new TTreeReader(metaname, f);
  data = new NTuple(t, m);
  iTrig = -1;
  progress_bar_.Set(t->GetEntries(), 70);
  verbose_ = verbose;
}

NTupleReader::~NTupleReader() {
  delete f;
  delete t;
  delete m;
  delete data;
}

bool NTupleReader::GetNextEvent() {
  return t->Next();
}

bool NTupleReader::GetNextTrigger() {
  if (++iTrig < 1) {
	return true;
  } else {
	iTrig = -1;
	return false;
  }
}

TData *NTupleReader::GetData() {
  return data;
}
