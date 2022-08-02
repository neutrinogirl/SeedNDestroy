//
// Created by Stephane Zsoldos on 7/7/22.
//

#ifndef SND_INCLUDE_SND_POST_HH_
#define SND_INCLUDE_SND_POST_HH_

#include <iostream>

#include <TVector3.h>
#include <TTree.h>

class PosT{
 protected:
 public:
  TVector3 Pos;
  double T=0.f;
  //
  PosT() : Pos(0.f, 0.f, 0.f), T(0.f) {}
  //
  PosT(const TVector3& P, const double& t){
	Pos = P;
	T = t;
  }
  explicit PosT(const std::vector<double> &x)
	  : Pos(x[0], x[1], x[2]), T(x[3]) {}
  PosT(const double &x1, const double &x2, const double &x3, const double &x4)
	  : Pos(x1, x2, x3), T(x4) {}
  //
  void Clear(){
	Pos = TVector3(0.f, 0.f, 0.f);
	T = 0.f;
  }
  virtual void SetTree(TTree *Tree){
	Tree->Branch("X", &this->Pos[0], "X/D");
	Tree->Branch("Y", &this->Pos[1], "Y/D");
	Tree->Branch("Z", &this->Pos[2], "Z/D");
	Tree->Branch("T", &this->T, "T/D");
  }
  virtual void Print() const{
	std::cout << "Pos: " << Pos.X() << " " << Pos.Y() << " " << Pos.Z() << std::endl;
	std::cout << "T: " << T << std::endl;
  }
  std::vector<double> GetStdVec() const{
	return {Pos.X(), Pos.Y(), Pos.Z(), T};
  }
};

bool operator==(const PosT& s1, const PosT& s2);

class RecT : public PosT {
 public:
  double NLL=0.f;
  //
  RecT() : PosT(), NLL(0.f) {}
  RecT(const PosT &P, const double& minf)
	  : PosT(P), NLL(minf) {}
  RecT(const TVector3& Pos, const double& T, const double& minf)
	  : PosT(Pos, T), NLL(minf) {}
  explicit RecT(const std::vector<double> &x) : PosT(x), NLL(x[4]) {}
  RecT(const double& x1, const double& y1, const double& z1,
	   const double& t1,
	   const double& minf)
	  : PosT(x1, y1, z1, t1), NLL(minf) {}
  //
  void SetTree(TTree *Tree) override{
	PosT::SetTree(Tree);
	Tree->Branch("NLL", &this->NLL, "NLL/D");
  }
};

#endif //SND_INCLUDE_SND_POST_HH_
