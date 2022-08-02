//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_SRC_SND_GEOM_HH_
#define SND_SRC_SND_GEOM_HH_

#include <TVector3.h>
#include "SnD/Utils.hh"

class Bnd {
 public:
  virtual double GetDWall(const TVector3& pos) = 0;
  virtual double GetTWall(const TVector3& pos) = 0;
  virtual bool IsInside(const TVector3& pos) = 0;
  virtual TVector3 GetEdge() = 0;
  virtual double GetTEdge() = 0;
};

class Cylinder : public Bnd {
 private:
  double R, HH;
  double T;
 public:
  Cylinder(double r, double hh) : R(r), HH(hh) { T = Cylinder::GetTWall(TVector3(0, 0, 0)); }
  double GetDWall(const TVector3& pos) override {
	std::min(R - pos.Perp(), HH - std::abs(pos.z()));
  };
  double GetTWall(const TVector3& pos) override {
	return GetDWall(pos) / Csts::GetSoL();
  }
  bool IsInside(const TVector3& pos) override {
	return pos.Perp() < R && std::abs(pos.z()) < HH;
  }
  TVector3 GetEdge() override {
	return {R, R, HH};
  }
  double GetTEdge() override {
	return T;
  }
};

#endif //SND_SRC_SND_GEOM_HH_
