//
// Created by Stephane Zsoldos on 7/6/22.
//

#include "TAppAnalysis.hh"

// #### #### #### #### #### #### #### #### #### #### #### #### //
MetaNTuple::MetaNTuple(TTreeReader *Reader) {
  SetReader(Reader);
  Reader->Next();
}

MetaNTuple::~MetaNTuple() {
  delete pmtId, delete pmtX, delete pmtY, delete pmtZ;
}

void MetaNTuple::SetReader(TTreeReader *Reader) {
  // MetaNTuple::~MetaNTuple();
  pmtId = new TTreeReaderValue<std::vector<int>>(*Reader, "pmtId");
  pmtX = new TTreeReaderValue<std::vector<double>>(*Reader, "pmtX");
  pmtY = new TTreeReaderValue<std::vector<double>>(*Reader, "pmtY");
  pmtZ = new TTreeReaderValue<std::vector<double>>(*Reader, "pmtZ");
}

TVector3 MetaNTuple::GetPMTPosition(int PMTId) {
  auto ID = (*pmtId->Get())[PMTId];
  return {(*pmtX->Get())[ID], (*pmtY->Get())[ID], (*pmtZ->Get())[ID]};
}

// #### #### #### #### #### #### #### #### #### #### #### #### //
NTuple::NTuple(TTreeReader *Reader) {
  SetReader(Reader);
}

NTuple::~NTuple() {
  delete hitPMTID, delete hitPMTTime, delete hitPMTCharge;
}

void NTuple::SetReader(TTreeReader *Reader) {
  // NTuple::~NTuple();
  hitPMTID = new TTreeReaderValue<std::vector<int>>(*Reader, "hitPMTID");
  hitPMTTime = new TTreeReaderValue<std::vector<double>>(*Reader, "hitPMTTime");
  hitPMTCharge = new TTreeReaderValue<std::vector<double>>(*Reader, "hitPMTCharge");
}

#include <boost/range/combine.hpp>
std::vector<Hit> NTuple::GetVHits() {
  std::vector<Hit> vHits;
  for (auto tup :
	  boost::combine(**hitPMTID, **hitPMTTime, **hitPMTCharge)) {
	int ID;
	double T, Q;
	boost::tie(ID, T, Q) = tup;
	vHits.emplace_back(Meta->GetPMTPosition(ID), T, Q, ID);
  }
  return vHits;
}

// #### #### #### #### #### #### #### #### #### #### #### #### //
NTupleReader::NTupleReader(const char *filename,
						   const char *treename, const char *metaname) {
  f = new TFile(filename);
  t = new TTreeReader(treename, f);
  data = new NTuple(t);
  m = new TTreeReader(metaname, f);
  meta = new MetaNTuple(m);
  iTrig = 1;
}

NTupleReader::~NTupleReader() {
  delete f;
  delete t;
  delete data;
  delete m;
  delete meta;
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

// #### #### #### #### #### #### #### #### #### #### #### #### //
void TAppAnalysis::Do(TData *Data) {

  for(const auto& hit : Data->GetVHits()) {
	hit.PMTPos.Print();
  }

}

#include <TFile.h>
void TAppAnalysis::Export(const char *filename) {
  TFile f(filename, "RECREATE");
  f.cd();
  f.Close();
}
