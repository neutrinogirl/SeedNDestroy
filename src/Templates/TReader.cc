//
// Created by Stephane Zsoldos on 2/22/23.
//

#include "Templates/TReader.hh"

volatile sig_atomic_t TReader::gSignalStatus = 0;

void TReader::MyHandler(int sig){
  TReader::gSignalStatus = sig;
};
void TReader::Read(TAnalysis *Ana){
  signal(SIGINT, MyHandler);
  while(this->GetNextEvent()){
	if(TReader::gSignalStatus == SIGINT)
	  break;
	++(*GetProgressBar());
	Ana->Do(this->GetData());
	if(GetVerbosity())
	  GetProgressBar()->display();
  }
  GetProgressBar()->done();
}