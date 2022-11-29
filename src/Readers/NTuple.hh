//
// Created by Stephane Zsoldos on 11/16/22.
//

#ifndef SND_SRC_READERS_NTUPLE_HH_
#define SND_SRC_READERS_NTUPLE_HH_

#include <Templates/TReader.hh>

#include <map>

#include <TFile.h>
#include <TTreeReader.h>
#include <TVector3.h>

#include <boost/any.hpp>

class Flat {
 private:
  std::vector<boost::any> many;
 public:
  explicit Flat(TTreeReader *Reader);
  TVector3 GetPosition();
  TVector3 GetDirection();
  double GetEnergy();
  int GetEventID();
  int GetSubEventID();
};

class FlatReader : public TReader {
 private:
  //
  enum ETreeReaders{
	kTree,
	kMeta
  };
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
