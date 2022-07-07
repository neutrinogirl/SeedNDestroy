//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_SRC_SND_RATDATA_HH_
#define SND_SRC_SND_RATDATA_HH_

#include <TTree.h>

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
  void SetTTree(TTree *Tree){
	Tree->Branch("TrigTime", &this->TrigTime, "TrigTime/D");
	Tree->Branch("Pos", &this->Pos, "Pos[3]/D");
	Tree->Branch("Dir", &this->Dir, "Dir[3]/D");
	Tree->Branch("T", &this->T, "T/D");
	Tree->Branch("E", &this->E, "E/D");
  }
};

#endif //SND_SRC_SND_RATDATA_HH_
