//
// Created by zsoldos on 8/12/20.
//

#ifndef _CENTROID_HH_
#define _CENTROID_HH_

#include <numeric>

#include <Hit.hh>
#include "MathUtils.hh"

TVector3 GetCentroidSeed(const std::vector<Hit>& vHits, const bnds& b,
						 const unsigned int& wPower = 0){

  TVector3 seed = std::accumulate(vHits.begin(), vHits.end(),
								  TVector3(),
								  [](const TVector3 &lhs, const Hit &rhs){
									return lhs+rhs.PMTPos;
								  }
  );
  double norm = std::accumulate(vHits.begin(), vHits.end(),
								0.,
								[&wPower](const double &lhs, const Hit &rhs){
								  return lhs+fweight(rhs, wPower);
								}
  );

  seed *= 1 / norm;

  return b.IsInPos(seed) ? seed : b.GetTVector3();

}

#endif //_CENTROID_HH_
