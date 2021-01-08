//
// Created by zsoldos on 8/12/20.
//

#ifndef _CENTROID_HH_
#define _CENTROID_HH_

#include "../wRATter/include/Hit.hh"
#include "MathUtils.hh"

static double GetTotQ(const std::vector<Hit>& vHits){
  double Q=0.;
  for(const auto& hit:vHits)
	Q+=hit.Q;
  return Q;
}

static double GetTotQ2(const std::vector<Hit>& vHits){
  double Q2=0.;
  for(const auto& hit:vHits)
	Q2+=hit.Q*hit.Q;
  return Q2;
}

TVector3 GetCentroidSeed(std::vector<Hit>& vHits, const double& bnds,
								const int& weightPower = 2){

  TVector3 Seed(0.,0.,0.);
  double xMean=0; double yMean=0; double zMean=0;

  std::sort(vHits.begin(), vHits.end());
  auto NHits = vHits.size();

  auto wQ = [&weightPower](const Hit& hit){
	return weightPower > 0 ? std::pow(hit.Q, weightPower) : 1;
	// return weightPower > 0 ? 1 - exp(-std::pow(hit.Q, weightPower)) : 1;
  };

  auto QNorm = GetTotQ(vHits);
  auto Q2Norm = GetTotQ2(vHits);

  Hit hCut;
  hCut.T = 0;
  const int TPrompt = 1000/*ns*/;
  const auto T0Hit = std::lower_bound(vHits.begin(), vHits.end(), hCut);
  const double T0 = T0Hit->T; // Get first hit > 0
  // std::cout << T0 << std::endl;
  const double TCut = T0+TPrompt; // upper bound cut above scintillation
  const double alpha = 1e1; // arbitral weight

  auto isAbovePreTrigger = [&T0](const Hit& h){return h.T > T0;};
  auto isPrompt = [&isAbovePreTrigger, &TCut](const Hit& h){return isAbovePreTrigger(h) && h.T < TCut;};

  // auto weight = [](const double& T){return T;};
  // auto weight = [&T0](const double& T){return T0/T;};
  // auto weight = [&T0, &alpha](const double& T){return exp(T0/(alpha*T));};
  // auto weight = [&T0, &TCut](const double& T){return T0 * (1 - (T-T0)/(TCut-T0));};

  // ##################################### //
  // ### #### ### FIND SEED #### #### #### //
  // ##################################### //

  for(const auto& hit:vHits){

	// if(!isAbovePreTrigger(hit)){

	// hit.PMTPos.Print();
	// std::cout << hit.T << " " << hit.Q << std::endl;

	xMean+=hit.PMTPos.x()*wQ(hit)/static_cast<double>(NHits);
	yMean+=hit.PMTPos.y()*wQ(hit)/static_cast<double>(NHits);
	zMean+=hit.PMTPos.z()*wQ(hit)/static_cast<double>(NHits);

	// }


  }

  // std::cout << std::endl;

  Seed = TVector3(xMean, yMean, zMean);

  // ######################################### //
  // RESCALE seed vector if outside boundaries //
  // ######################################### //

  auto isFV = [&bnds](const TVector3& Pos){
	return abs(Pos.x()) < bnds && abs(Pos.y()) < bnds && abs(Pos.z()) < bnds;
  };

  if(isFV(Seed)){
	return Seed;
  } else {
	Seed.SetMag(bnds);
	return Seed;
  }
}

#endif //_CENTROID_HH_
