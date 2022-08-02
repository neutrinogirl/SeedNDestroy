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
  double Q = 0.f;
  double NHits = 0.f;
  void Clear(){
	TrigTime = 0.f;
	Pos = TVector3(0.f, 0.f, 0.f);
	Dir = TVector3(0.f, 0.f, 0.f);
	T = 0.f;
	E = 0.f;
	Q = 0.f;
	NHits = 0.f;
	vHits.clear();
  }
  void SetTTree(TTree *Tree){
	Tree->Branch("TrigTime", &this->TrigTime, "TrigTime/D");
	Tree->Branch("PosTrue", &this->Pos, "PosTrue[3]/D");
	Tree->Branch("DirTrue", &this->Dir, "DirTrue[3]/D");
	Tree->Branch("DirX", &this->Dir[0], "DirX/D");
	Tree->Branch("DirY", &this->Dir[1], "DirY/D");
	Tree->Branch("DirZ", &this->Dir[2], "DirZ/D");
	Tree->Branch("TTrue", &this->T, "T/D");
	Tree->Branch("ETrue", &this->E, "E/D");
	Tree->Branch("Q", &this->Q, "Q/D");
	Tree->Branch("NHits", &this->NHits, "NHits/D");
  }
};

#endif //SND_SRC_SND_RATDATA_HH_
