//
// Created by Stephane Zsoldos on 7/5/22.
//

#ifndef SND_SRC_SND_RATANALYSIS_HH_
#define SND_SRC_SND_RATANALYSIS_HH_

#include <wRATter/Wrapper.hh>
#include "TReader.hh"
#include "RATData.hh"
class RATReader : public TReader {
 private:
  wRAT w_rat;
  RATData *d;
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
  ~RATReader() { delete d; };
};

#endif //SND_SRC_SND_RATANALYSIS_HH_
