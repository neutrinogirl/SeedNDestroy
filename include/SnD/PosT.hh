//
// Created by Stephane Zsoldos on 7/7/22.
//

#ifndef SND_INCLUDE_SND_POST_HH_
#define SND_INCLUDE_SND_POST_HH_

#include <iostream>

#include <TVector3.h>
#include <TTree.h>

#include "SnD/ZVector.hh"

class PosT{
 protected:
 public:
  double X=0.f, Y=0.f, Z=0.f, T=0.f;
  //
  PosT() = default;
  PosT(const PosT &rhs) = default;
  PosT(const double &x, const double &y, const double &z, const double &t)
	  : X(x), Y(y), Z(z), T(t) {}
  PosT(const TVector3 &v, const double &t)
	  : X(v.x()), Y(v.y()), Z(v.z()), T(t) {}
  PosT(const Vector3<double> &v, const double& t)
	  : X(v.GetX()), Y(v.GetY()), Z(v.GetZ()), T(t) {}
  //
  //
  void Clear(){
	X = 0.f;
	Y = 0.f;
	Z = 0.f;
	T = 0.f;
  }
  virtual void SetTree(TTree *Tree){
	Tree->Branch("X", &this->X, "X/D");
	Tree->Branch("Y", &this->Y, "Y/D");
	Tree->Branch("Z", &this->Z, "Z/D");
	Tree->Branch("T", &this->T, "T/D");
  }
  virtual void Print() const{
	std::cout << "Pos: " << X << " " << Y << " " << Z << std::endl;
	std::cout << "T: " << T << std::endl;
  }
  std::vector<double> GetStdVec() const{
	return {X, Y, Z, T};
  }
  TVector3 GetTVector3() const{
	return TVector3(X, Y, Z);
  }
};

bool operator==(const PosT& s1, const PosT& s2);

class RecT : public PosT {
 public:
  double NLL=0.f;
  RecT() = default;
  RecT(const RecT &rhs) = default;
  RecT(const double &x, const double &y, const double &z,
	   const double &t,
	   const double &nll)
	  : PosT(x, y, z, t), NLL(nll) {}
  //
  //
  void SetTree(TTree *Tree) override{
	PosT::SetTree(Tree);
	Tree->Branch("NLL", &this->NLL, "NLL/D");
  }
  //
  PosT GetPosT() const{
	return PosT(X, Y, Z, T);
  }
  void Print() const override{
	PosT::Print();
	std::cout << "NLL: " << NLL << std::endl;
  }
};

#endif //SND_INCLUDE_SND_POST_HH_
