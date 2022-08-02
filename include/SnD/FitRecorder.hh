//
// Created by Stephane Zsoldos on 8/2/22.
//

#ifndef SND_INCLUDE_SND_FITRECORDER_HH_
#define SND_INCLUDE_SND_FITRECORDER_HH_

typedef struct StepRecorder {
  double X=0.f, Y=0.f, Z=0.f, T=0.f, NLL=0.f;
  int i=0;
  void Reset(){
	X=0.f; Y=0.f; Z=0.f; T=0.f; NLL=0.f;
	i=0;
  }
} StepRecorder;

#endif //SND_INCLUDE_SND_FITRECORDER_HH_
