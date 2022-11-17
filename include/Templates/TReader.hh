//
// Created by Stephane Zsoldos on 11/16/22.
//

#ifndef SND_INCLUDE_TEMPLATES_TREADER_HH_
#define SND_INCLUDE_TEMPLATES_TREADER_HH_

#include "TData.hh"
#include "TAnalysis.hh"

#include "ProgressBar/ProgressBar.hpp"

#include <csignal>
#include <cstdio>

class TReader {
 protected:
  virtual bool GetNextEvent() = 0;
  virtual bool GetNextTrigger() = 0;
  virtual TData *GetData() = 0;
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

#endif //SND_INCLUDE_TEMPLATES_TREADER_HH_
