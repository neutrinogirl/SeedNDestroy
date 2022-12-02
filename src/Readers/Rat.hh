//
// Created by Stephane Zsoldos on 11/16/22.
//

#ifndef SND_SRC_READERS_RAT_HH_
#define SND_SRC_READERS_RAT_HH_

#include "wRATter/Wrapper.hh"

#include "Templates/TData.hh"

class RATData : public TData {
 private:
  TVector3 Pos;
  TVector3 Dir;
  double Energy;
  double Time;
  std::vector<Hit> vHits;
  int EventID;
  int TriggerID;
 public:
  RATData() = default;
  void Update(wRAT *w_rat);
  TVector3 GetPosition() override {return Pos;};
  TVector3 GetDirection() override {return Dir;};
  double GetEnergy() override {return Energy;};
  double GetTime() override {return Time;};
  std::vector<Hit> GetVHits() override {return vHits;};
  int GetEventID() override {return EventID;};
  int GetTriggerID() override {return TriggerID;};
};


#include <Templates/TReader.hh>

class RATReader : public TReader{
 private:
  wRAT w_rat;
  RATData *data;
  ProgressBar progress_bar_;
  bool verbose_;
 protected:
  bool GetNextEvent() override;
  bool GetNextTrigger() override;
  void* GetData() override;
  ProgressBar *GetProgressBar() override { return &progress_bar_; }
  bool GetVerbosity() override { return verbose_; }
 public:
  explicit RATReader(const char *filename, const bool &verbose = false);
  ~RATReader() { delete data; };

};

#endif //SND_SRC_READERS_RAT_HH_
