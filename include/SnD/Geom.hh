//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_INCLUDE_SND_GEOM_HH_
#define SND_INCLUDE_SND_GEOM_HH_

#include "SnD/Vector3.hh"
#include "SnD/Utils.hh"

struct Edges {
  virtual double GetDWall(const Vector3& pos) = 0;
  virtual double GetTWall(const Vector3& pos) = 0;
  virtual bool   IsInside(const Vector3& pos) = 0;
};

struct CylEdges : public Edges {
  double radius;
  double halfheight;
  double T;
  SpaceUnit unit;
  CylEdges(double r, double hh, SpaceUnit u)
	  : radius(r), halfheight(hh), unit(u) {
	T = GetTWall(Vector3(0, 0, 0, unit));
  };
  [[nodiscard]] double GetDWall(const Vector3& pos) override {
	return std::min(radius - pos.Get(unit).GetPerp(),
					halfheight - std::abs(pos.Get(unit).GetZ()));
  };
  [[nodiscard]] double GetTWall(const Vector3& pos) override {
	return GetDWall(pos) / Csts::GetSoL();
  }
  [[nodiscard]] bool IsInside(const Vector3& pos) override {
	return GetDWall(pos) > 0;
  }
};

#endif //SND_INCLUDE_SND_GEOM_HH_
