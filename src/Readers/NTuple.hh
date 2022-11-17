//
// Created by Stephane Zsoldos on 11/16/22.
//

#ifndef SND_SRC_READERS_NTUPLE_HH_
#define SND_SRC_READERS_NTUPLE_HH_

#include <Templates/TReader.hh>

#include <map>

#include <TFile.h>
#include <TTreeReader.h>

class MetaNTuple {
 private:
  std::map<int, TVector3> mPMTPos;
 public:
  explicit MetaNTuple(TTreeReader *Reader);
  ~MetaNTuple();
  TVector3 GetPMTPosition(const int& PMTID);
};

class NTuple : public TData {
 private:
  MetaNTuple *meta;
  TTreeReaderValue<std::vector<int>> *hitPMTID;
  TTreeReaderValue<std::vector<double>> *hitPMTTime;
  TTreeReaderValue<std::vector<double>> *hitPMTCharge;
 public:
  NTuple() = default;
  ~NTuple();
  NTuple(TTreeReader *Reader, TTreeReader *mReader);
  void SetReader(TTreeReader *Reader);
  std::vector<Hit> GetVHits() override;
};

class NTupleReader : public TReader {
 private:
  TFile *f;
  TTreeReader *t;
  TTreeReader *m;
  NTuple *data;
  int iTrig;
  ProgressBar progress_bar_;
  bool verbose_;
 public:
  explicit NTupleReader(const char *filename,
						const char *treename="output", const char *metaname="meta",
						const bool &verbose=false);
  ~NTupleReader();
  bool GetNextEvent() override;
  bool GetNextTrigger() override;
  TData *GetData() override;
  ProgressBar *GetProgressBar() override { return &progress_bar_; }
  bool GetVerbosity() override { return verbose_; }
};

#endif //SND_SRC_READERS_NTUPLE_HH_
