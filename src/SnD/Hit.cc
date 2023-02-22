//
// Created by Stephane Zsoldos on 7/6/22.
//

#include <SnD/Hit.hh>

// Generate comparison between two hits based on T
bool operator<(const Hit& h1, const Hit& h2){
  return h1.T < h2.T;
}
bool operator==(const Hit& h1, const Hit& h2){
  return h1.T == h2.T;
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

double GetFirstHitTime(const std::vector<Hit>& vHits, const double& threshold){
  double T = 0.;
  for(const auto& hit: vHits){
	if(hit.Q>threshold){
	  T = hit.T;
	  break;
	}
  }
  return T;
}

double GetFirstHitTime(const std::vector<Hit>& vHits){
  double T = std::numeric_limits<double>::max();
  for(const auto& hit: vHits){
	if(hit.T<T){
	  T = hit.T;
	}
  }
  return T;
}

double GetWindowHitTime(const std::vector<Hit>& vHits, const double& threshold, const int& windowsize){
  double T = 0.;
  for(int i=0; i<vHits.size()-windowsize; i++){
	double Q = 0.;
	for(int j=0; j<windowsize; j++) {
	  Q += vHits[i + j].Q;
	}
	if(Q>threshold){
	  T = vHits[i].T;
	  break;
	}
  }
  return T;
}

double GetMaxHitTime(const std::vector<Hit>& vHits){
  double T = 0.;
  double Q = 0.;
  for(const auto& hit: vHits){
	if(hit.Q>Q){
	  T = hit.T;
	  Q = hit.Q;
	}
  }
  return T;
}