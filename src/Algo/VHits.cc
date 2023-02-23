//
// Created by Stephane Zsoldos on 2/23/23.
//

#include "Algo/VHits.hh"
#include "Algo/WOpt.hh"

PosT GetCentroidBasedSeed(const std::vector<Hit>& vHits, Bnd *b){
  // Get centroid seed
  TVector3 Centroid = GetCentroid(vHits);
  Vector3<double> v3Centroid(Centroid.x(), Centroid.y(), Centroid.z(), SpaceUnit::mm);
  Vector3<double> vSeed = v3Centroid.ConvertTo(SpaceUnit::dm);

  // Get time seed
  double TSeed = b->GetTWall(Centroid);

  return PosT(vSeed, TSeed);

}

boost::optional<PosT> GetMLATSeed(const std::vector<Hit>& vHits, Bnd *b){
  try{

	TVector3 MLAT = GetMLAT(vHits);

	if(b->IsInside(MLAT)){
	  Vector3<double> v3MLAT(MLAT.x(), MLAT.y(), MLAT.z(), SpaceUnit::mm);
	  Vector3<double> vSeed = v3MLAT.ConvertTo(SpaceUnit::dm);

	  // Get time seed
	  double TSeed = b->GetTWall(MLAT);
	  return PosT(vSeed, TSeed);
	}

  } catch (std::exception& e){

  }

}

boost::optional<PosT> GetLSBasedSeed(const std::vector<Hit>& vHits, Bnd *b, std::vector<PosT> vSeeds){
  try{
	RecT RT = Recon((void*)&vHits, b, vSeeds, GetAlgo(2), fLS, {SetBounds, SetPars, SetInequalityConstraint});
	return RT.GetPosT();
  } catch (std::exception& e){

  }
}
