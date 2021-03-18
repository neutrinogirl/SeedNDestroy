//
// Created by zsoldos on 8/12/20.
//

#ifndef _CENTROID_HH_
#define _CENTROID_HH_

#include <numeric>

#include <Hit.hh>
#include "MathUtils.hh"

TVector3 GetCentroidSeed(const std::vector<Hit>& vHits, const bnds& b,
						 const unsigned int& wPower = 2){

  TVector3 Seed(0.,0.,0.);
  double xMean=0; double yMean=0; double zMean=0;

  double NormQ = 0.;
  for(const auto& hit:vHits)
    NormQ += fweight(hit, wPower);
  double Norm = NormQ;

  // ##################################### //
  // ### #### ### FIND SEED #### #### #### //
  // ##################################### //

  for(const auto& hit:vHits){

	xMean+=hit.PMTPos.x()*fweight(hit, wPower)/Norm;
	yMean+=hit.PMTPos.y()*fweight(hit, wPower)/Norm;
	zMean+=hit.PMTPos.z()*fweight(hit, wPower)/Norm;

  }

  Seed = TVector3(xMean, yMean, zMean);

  return b.IsInPos(Seed) ? Seed : b.GetTVector3();

}

#endif //_CENTROID_HH_
