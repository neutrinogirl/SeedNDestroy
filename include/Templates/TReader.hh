//
// Created by Stephane Zsoldos on 11/16/22.
//

#ifndef SND_INCLUDE_TEMPLATES_TREADER_HH_
#define SND_INCLUDE_TEMPLATES_TREADER_HH_

#include "TAnalysis.hh"

#include "ProgressBar/ProgressBar.hpp"

#include <csignal>
#include <cstdio>

class TReader {
 protected:
  virtual bool GetNextEvent() = 0;
  virtual void *GetData() = 0;
  virtual ProgressBar *GetProgressBar() = 0;
  virtual bool GetVerbosity() = 0;
 public:
  static volatile sig_atomic_t gSignalStatus;
  static void MyHandler(int sig);
  virtual void Read(TAnalysis *Ana);
};

#endif //SND_INCLUDE_TEMPLATES_TREADER_HH_
