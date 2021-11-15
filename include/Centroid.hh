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

  TVector3 Seed(0.,0.,0.);

  double NormQ = 0.;
  for(const auto& hit:vHits)
    NormQ += fweight(hit, wPower);
  double Norm = NormQ;

  // ##################################### //
  // ### #### ### FIND SEED #### #### #### //
  // ##################################### //

  for(const auto& hit:vHits)
    Seed+=hit.PMTPos*fweight(hit, wPower);
  Seed *= 1 / Norm;

  return b.IsInPos(Seed) ? Seed : b.GetTVector3();

}

#endif //_CENTROID_HH_
