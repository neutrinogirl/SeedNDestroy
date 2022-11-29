//
// Created by Stephane Zsoldos on 11/16/22.
//

#ifndef SND_SRC_READERS_RAT_HH_
#define SND_SRC_READERS_RAT_HH_

#include "wRATter/Wrapper.hh"

class RATData : public TData {
 private:
  std::vector<Hit> vHits;
 public:
  RATData() = default;
  explicit RATData(const std::vector<Hit> &v_hits)
	  : vHits(v_hits) {}
  std::vector<Hit> GetVHits() override { return vHits; };
};


#include "Templates/TReader.hh"

class RATReader : public TReader{
 private:
  wRAT w_rat;
  RATData *data;
  ProgressBar progress_bar_;
  bool verbose_;
 protected:
  bool GetNextEvent() override;
  bool GetNextTrigger() override;
  RATData* GetData() override;
  ProgressBar *GetProgressBar() override { return &progress_bar_; }
  bool GetVerbosity() override { return verbose_; }
 public:
  explicit RATReader(const char *filename, const bool &verbose = false);
  ~RATReader() { delete d; };

};

#endif //SND_SRC_READERS_RAT_HH_
