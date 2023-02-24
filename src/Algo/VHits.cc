//
// Created by Stephane Zsoldos on 2/23/23.
//

#include "Algo/VHits.hh"
#include "Algo/WOpt.hh"
#include "SnD/ZVector.hh"

PosT GetCentroidBasedSeed(const std::vector<Hit>& vHits, Bnd *b){
  // Get centroid seed
  TVector3 Centroid = GetCentroid(vHits);
  // Get time seed
  double TSeed = b->GetTWall(Centroid);
  //
  return PosT(ConvertTVector3Unit<double>(Centroid, SpaceUnit::mm, SpaceUnit::dm), TSeed);

}

boost::optional<PosT> GetMLATSeed(const std::vector<Hit>& vHits, Bnd *b){
  try{

	TVector3 MLAT = GetMLAT(vHits);

	if(b->IsInside(MLAT)){
	  // Get time seed
	  double TSeed = b->GetTWall(MLAT);
	  return PosT(ConvertTVector3Unit<double>(MLAT, SpaceUnit::mm, SpaceUnit::dm), TSeed);
	}

  } catch (std::exception& e){

  }

}

PosT GetLSBasedSeed(const std::vector<Hit>& vHits, Bnd *b, std::vector<PosT> vSeeds){
  RecT RT = Recon((void*)&vHits, b, vSeeds, GetAlgo(2), fLS, {SetBounds, SetPars, SetInequalityConstraint});
  return RT.GetPosT();
}
