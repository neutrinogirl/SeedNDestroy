//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_SRC_SND_RATDATA_HH_
#define SND_SRC_SND_RATDATA_HH_

#include <TVector3.h>

#include <SnD/Hit.hh>

class RATData {
 public:
  double TrigTime = 0.f;
  TVector3 Pos = TVector3(0.f, 0.f, 0.f);
  TVector3 Dir = TVector3(0.f, 0.f, 0.f);
  double T = 0.f;
  double E = 0.f;
  std::vector<Hit> vHits;
  void Clear(){
	TrigTime = 0.f;
	Pos = TVector3(0.f, 0.f, 0.f);
	Dir = TVector3(0.f, 0.f, 0.f);
	T = 0.f;
	E = 0.f;
	vHits.clear();
  }
};

#endif //SND_SRC_SND_RATDATA_HH_
