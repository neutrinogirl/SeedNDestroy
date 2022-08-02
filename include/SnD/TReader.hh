//
// Created by Stephane Zsoldos on 7/5/22.
//

#ifndef SND_SRC_SND_TREADER_HH_
#define SND_SRC_SND_TREADER_HH_

#include "TAnalysis.hh"

#include "ProgressBar/ProgressBar.hpp"

#include <csignal>
#include <cstdio>

class TReader {
 protected:
  virtual bool GetNextEvent() = 0;
  virtual bool GetNextTrigger() = 0;
  virtual void *GetData() = 0;
  virtual ProgressBar *GetProgressBar() = 0;
  virtual bool GetVerbosity() = 0;
 public:
  static volatile sig_atomic_t gSignalStatus;
  static void MyHandler(int sig){
	TReader::gSignalStatus = sig;
  };
  virtual void Read(TAnalysis *Ana){
	signal(SIGINT, MyHandler);
	while(this->GetNextEvent()){
	  if(TReader::gSignalStatus == SIGINT)
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


#endif //SND_SRC_SND_TREADER_HH_
