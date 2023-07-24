//
// Created by Stephane Zsoldos on 11/16/22.
//

#ifndef SND_SRC_READERS_ANNIE_HH_
#define SND_SRC_READERS_ANNIE_HH_

//
#include <Templates/TReader.hh>
//
#include <map>
#include <any>
//
#include <TFile.h>
#include <TTreeReader.h>
//
//
#include "Templates/TData.hh"

//
class ANNIE : public TData {
 private:
  //
  std::map<int, Vector3> mPMTPos;
  //
  std::vector<std::any> many;
 public:
  explicit ANNIE(std::vector< TTreeReader* > &vTreeReaders);
  Vector3 GetPosition() override {return Vector3();}
  Vector3 GetDirection() override {return Vector3();}
  double GetEnergy() override {return 0.0;}
  double GetTime() override;
  std::vector<Hit> GetVHits() override;
  int GetEventID() override;
  int GetTriggerID() override {return 0;}
};

//
class ANNIEReader : public TReader {
 private:
  //
  ANNIE *fFlat;
  //
  TFile *f;
  std::vector< TTreeReader* > vTreeReaders;
  //
  ProgressBar progress_bar_;
  bool verbose_;
 public:
  explicit ANNIEReader(const char *filename,
					   const char *treename="phaseIITankClusterTree",
					   const bool &verbose=false);
  ~ANNIEReader();
  bool GetNextEvent() override;
  void *GetData() override;
  ProgressBar *GetProgressBar() override { return &progress_bar_; }
  bool GetVerbosity() override { return verbose_; }
};

#endif //SND_SRC_READERS_ANNIE_HH_
