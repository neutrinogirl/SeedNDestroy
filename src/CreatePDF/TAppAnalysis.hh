//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_SRC_SND_TAPPANALYSIS_HH_
#define SND_SRC_SND_TAPPANALYSIS_HH_

// #### #### #### #### #### #### #### #### #### #### #### #### //
// This goes into "SnD/TData.hh"
#include <SnD/Hit.hh>

class TData{
 protected:
 public:
  virtual std::vector<Hit> GetVHits() = 0;
};

// #### #### #### #### #### #### #### #### #### #### #### #### //
// This goes into "SnD/TAnalysis.hh"
class TAAnalysis {
 protected:
 public:
  virtual void Do(TData *Data) = 0;
};

// #### #### #### #### #### #### #### #### #### #### #### #### //
// This goes into "SnD/TReader.hh"
#include "ProgressBar/ProgressBar.hpp"

#include <csignal>
#include <cstdio>

class TRReader {
 protected:
  virtual bool GetNextEvent() = 0;
  virtual bool GetNextTrigger() = 0;
  virtual TData *GetData() = 0;
  virtual ProgressBar *GetProgressBar() = 0;
  virtual bool GetVerbosity() = 0;
 public:
  static volatile sig_atomic_t gSignalStatus;
  static void MyHandler(int sig){
	TRReader::gSignalStatus = sig;
  };
  virtual void Read(TAAnalysis *Ana){
	signal(SIGINT, MyHandler);
	while(this->GetNextEvent()){
	  if(TRReader::gSignalStatus == SIGINT)
		break;
	  ++(*GetProgressBar());
	  while(this->GetNextTrigger()){
		Ana->Do(this->GetData());
	  }
	  if(GetVerbosity())
		GetProgressBar()->display();
	}
	GetProgressBar()->done();
  }
};

// #### #### #### #### #### #### #### #### #### #### #### #### //
//   #### #### #### #### DERIVED CLASSES #### #### #### ####   //
// #### #### #### #### #### #### #### #### #### #### #### #### //
#include <TTreeReader.h>
class MetaNTuple {
 private:
  std::map<int, TVector3> mPMTPos;
 public:
  MetaNTuple(TTreeReader *Reader);
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

#include <TFile.h>
class NTupleReader : public TRReader {
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

class TAppAnalysis : public TAAnalysis {
 private:
 public:
  TAppAnalysis() = default;
  void Do(TData* Data) override;
  void Export(const char *filename);
};

#endif //SND_SRC_SND_PDFANALYSIS_HH_
