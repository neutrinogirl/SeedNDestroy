//
// Created by Stephane Zsoldos on 11/16/22.
//

#ifndef SND_SRC_READERS_NTUPLE_HH_
#define SND_SRC_READERS_NTUPLE_HH_

//
#include <Templates/TReader.hh>
//
#include <map>
//
#include <TFile.h>
#include <TTreeReader.h>
//
#include <boost/any.hpp>
//
#include "Templates/TData.hh"

//
enum ETreeReaders{
  kTree,
  kMeta
};

//
class Flat : public TData {
 private:
  //
  std::map<int, TVector3> mPMTPos;
  //
  std::vector<boost::any> many;
 public:
  explicit Flat(std::vector< TTreeReader* > &vTreeReaders);
  TVector3 GetPosition() override;
  TVector3 GetDirection() override;
  double GetEnergy() override;
  double GetTime() override;
  std::vector<Hit> GetVHits() override;
  int GetEventID() override;
  int GetTriggerID() override;
};

//
class FlatReader : public TReader {
 private:
  //
  Flat *fFlat;
  //
  TFile *f;
  std::vector< TTreeReader* > vTreeReaders;
  //
  int iTrig;
  //
  ProgressBar progress_bar_;
  bool verbose_;
 public:
  explicit FlatReader(const char *filename,
						const char *treename="output", const char *metaname="meta",
						const bool &verbose=false);
  ~FlatReader();
  bool GetNextEvent() override;
  bool GetNextTrigger() override;
  void *GetData() override;
  ProgressBar *GetProgressBar() override { return &progress_bar_; }
  bool GetVerbosity() override { return verbose_; }
};

#endif //SND_SRC_READERS_NTUPLE_HH_
