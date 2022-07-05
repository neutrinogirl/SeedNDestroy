//
// Created by Stephane Zsoldos on 7/2/22.
//

#ifndef SND_INCLUDE_SND_VHITS_HH_
#define SND_INCLUDE_SND_VHITS_HH_

#include "Hit.hh"

class VHits {
 private:
  std::vector<Hit> vHit;
 public:
  VHits() = default;
  VHits(const std::vector<Hit> &vHit) : vHit(vHit) {}
  ~VHits() = default;
  void AddHit(const Hit &hit) { vHit.emplace_back(hit); }
  void AddHits(const std::vector<Hit> &vHit) {
	for (auto &hit: vHit) {
	  this->vHit.emplace_back(hit);
	}
  }
  const std::vector<Hit> &GetHits() const { return vHit; }
  std::vector<Hit> &GetHits() { return vHit; }
  void Clear() { vHit.clear(); }
  void Print() const {
	for (auto &hit: vHit) {
	  hit.Print();
	}
  }
  // Calculate the centroid of a vector of hits
  TVector3 GetCentroid() const {
	TVector3 centroid(0., 0., 0.);
	for (auto &hit: vHit) {
	  centroid += hit.PMTPos;
	}
	centroid.SetMag(centroid.Mag()/vHit.size());
  };
};
#endif //SND_INCLUDE_SND_VHITS_HH_
