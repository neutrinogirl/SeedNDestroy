//
// Created by Stephane Zsoldos on 7/5/22.
//

#include "SnD/NTupleReader.hh"

NTupleReader::NTupleReader(const char *filename, const char *treename,
						   const bool &verbose) {

  f = new TFile(filename, "READ");
  t = new TTreeReader(treename, f);
  iTrig = -1;
  d = new NTupleData(t);
  progress_bar_.Set(t->GetEntries(), 70);
  verbose_ = verbose;
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

void *NTupleReader::GetData() {
  return d;
}