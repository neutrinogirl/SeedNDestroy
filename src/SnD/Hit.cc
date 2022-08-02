//
// Created by Stephane Zsoldos on 7/6/22.
//

#include <SnD/Hit.hh>

// Generate comparison between two hits based on T
bool operator<(const Hit& h1, const Hit& h2){
  return h1.T < h2.T;
}

double GetNPrompts(const std::vector<Hit>& vHits, const double& T){
  double NPrompts = 0.;
  for(auto& hit: vHits){
	if(hit.T<T){
	  NPrompts++;
	}
  }
  return NPrompts;
};

double fWeight(const Hit& h, const int& P){
  return std::pow(h.Q, P);
}

// Calculate centroid of a vector of hits
TVector3 GetCentroid(const std::vector<Hit>& vHits){
  TVector3 centroid(0., 0., 0.);
  const auto NHits = static_cast<double>(vHits.size());
  for(const auto& hit: vHits){
	centroid += hit.PMTPos * (1./NHits);
  }
  return centroid;
}