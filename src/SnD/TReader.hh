//
// Created by Stephane Zsoldos on 7/5/22.
//

#ifndef SND_SRC_SND_TREADER_HH_
#define SND_SRC_SND_TREADER_HH_

#include "TAnalysis.hh"

#include <ProgressBar/ProgressBar.hpp>

class TReader {
 protected:
  virtual bool GetNextEvent() = 0;
  virtual bool GetNextTrigger() = 0;
  virtual void *GetData() = 0;
  virtual ProgressBar *GetProgressBar() = 0;
  virtual bool GetVerbosity() = 0;
 public:
  virtual void Read(TAnalysis *Ana){
	while(this->GetNextEvent()){
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
