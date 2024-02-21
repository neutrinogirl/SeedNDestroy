//
// Created by Stephane Zsoldos on 6/25/23.
//

#ifndef SND_INCLUDE_SND_COORD_HH_
#define SND_INCLUDE_SND_COORD_HH_

#include "Vector3.hh"

#include <TTree.h>

class Coord : public Vector3 {
 protected:
  double T{};
 public:
  Coord() = default;
  Coord(double x, double y, double z, SpaceUnit unit, double t)
	  : Vector3(x, y, z, unit), T(t) {}
  virtual void SetTree(TTree *Tree){
	Tree->Branch("X", &GetXRef(), "X/D");
	Tree->Branch("Y", &GetYRef(), "Y/D");
	Tree->Branch("Z", &GetZRef(), "Z/D");
	Tree->Branch("T", &this->T, "T/D");
  }
  friend std::ostream& operator<<(std::ostream& os, const Coord& coord) {
	os << static_cast<const Vector3&>(coord);
	os << " T: " << coord.T;
	return os;
  }
  // Get T
  [[nodiscard]] double GetT() const { return T; }
};

class RecCoord : public Coord {
 protected:
  double NLL{};
  double mcx{}, mcy{}, mcz{};
 public:
  RecCoord() = default;
  RecCoord(double x, double y, double z, SpaceUnit unit, double t, double nll)
	  : Coord(x, y, z, unit, t), NLL(nll) {}
  void SetTree(TTree* Tree) override {
	Coord::SetTree(Tree);  // Call the base class implementation
	Tree->Branch("NLL", &this->NLL, "NLL/D");
	Tree->Branch("mcx", &this->mcx, "mcx/D");
	Tree->Branch("mcy", &this->mcy, "mcy/D");
	Tree->Branch("mcz", &this->mcz, "mcz/D");
  }
  void SetMC(const Vector3& v){
	mcx = v.GetX();
	mcy = v.GetY();
	mcz = v.GetZ();
  }
  friend std::ostream& operator<<(std::ostream& os, const RecCoord& reccoord) {
	os << static_cast<const Coord&>(reccoord);
	os << " NLL: " << reccoord.NLL;
	return os;
  }
  // Get NLL
  [[nodiscard]] double GetNLL() const { return NLL; }
};

#endif //SND_INCLUDE_SND_COORD_HH_
