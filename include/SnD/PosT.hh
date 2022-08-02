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
	Tree->Branch("PosX", &this->Pos[0], "PosX/D");
	Tree->Branch("PosY", &this->Pos[1], "PosY/D");
	Tree->Branch("PosZ", &this->Pos[2], "PosZ/D");
	Tree->Branch("T", &this->T, "T/D");
  }
  void Print() const{
	std::cout << "Pos: " << Pos.X() << " " << Pos.Y() << " " << Pos.Z() << std::endl;
	std::cout << "T: " << T << std::endl;
  }
};

bool operator==(const PosT& s1, const PosT& s2);

#endif //SND_INCLUDE_SND_POST_HH_
