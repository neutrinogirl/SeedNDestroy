//
// Created by Stephane Zsoldos on 7/7/22.
//

#ifndef SND_INCLUDE_SND_POST_HH_
#define SND_INCLUDE_SND_POST_HH_

#include <iostream>

#include <TVector3.h>
#include <TTree.h>

class PosT{
 public:
  TVector3 Pos = TVector3(0.f, 0.f, 0.f);
  double T = 0.f;
  PosT() = default;
  PosT(const TVector3& Pos, const double& T) : Pos(Pos), T(T) {}
  void Clear(){
	Pos = TVector3(0.f, 0.f, 0.f);
	T = 0.f;
  }
  void SetTree(TTree *Tree){
	Tree->Branch("Pos", &this->Pos, "Pos[3]/D");
	Tree->Branch("T", &this->T, "T/D");
  }
  void Print(){
	std::cout << "Pos: " << Pos.X() << " " << Pos.Y() << " " << Pos.Z() << std::endl;
	std::cout << "T: " << T << std::endl;
  }
};

bool operator==(const PosT& s1, const PosT& s2);

#endif //SND_INCLUDE_SND_POST_HH_
