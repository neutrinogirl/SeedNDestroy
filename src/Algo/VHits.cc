//
// Created by Stephane Zsoldos on 2/23/23.
//

#include "Algo/VHits.hh"
#include "Algo/WOpt.hh"

Coord GetCentroidBasedSeed(const std::vector<Hit>& vHits, Bnd *b){
  // Get centroid seed
  TVector3 Centroid = GetCentroid(vHits);
  // Get time seed
  double TSeed = b->GetTWall(Centroid);
  //
  return Coord(ConvertTVector3Unit<double>(Centroid, SpaceUnit::mm, SpaceUnit::dm), TSeed);

}

boost::optional<Coord> GetMLATSeed(const std::vector<Hit>& vHits, Bnd *b){
  try{

	TVector3 MLAT = GetMLAT(vHits);

	if(b->IsInside(MLAT)){
	  // Get time seed
	  double TSeed = b->GetTWall(MLAT);
	  return Coord(ConvertTVector3Unit<double>(MLAT, SpaceUnit::mm, SpaceUnit::dm), TSeed);
	}

  } catch (std::exception &e){

	// std::cerr << "MLAT failed: " << e.what() << std::endl;

  }

  return boost::none;

}

Coord GetLSBasedSeed(const std::vector<Hit>& vHits, Bnd *b, std::vector<Coord> vSeeds){
  RecT RT = Recon((void*)&vHits, b, vSeeds, nlopt::LN_SBPLX, fLS, {SetBounds, SetPars, SetInequalityConstraint});
  return RT.GetPosT();
}

std::vector<PosT> GetSeeds(std::vector<Hit> vHits, Bnd *b){
  // Init vSeeds with Centroid
  std::vector<PosT> vSeeds = {
	  GetCentroidBasedSeed(vHits, b)
  };

  // Sort all hits from closest to farthest from Centroid position
  auto Centroid = GetCentroid(vHits);
  SortHitsFromPos(vHits, Centroid);

  // Init vSeeds with MLAT
  auto MLAT = GetMLATSeed(vHits, b);
  if(MLAT.is_initialized())
	vSeeds.emplace_back(*MLAT);

  // Get subsets of seeds
  auto mHits = GetSubsets(vHits, Centroid, 1e1);
  for (auto it = mHits.begin(); it != mHits.end(); ++it){
	vSeeds.emplace_back(GetCentroidBasedSeed(it->second, b));
	SortHitsFromPos(it->second, GetCentroid(it->second));
	auto MLAT = GetMLATSeed(it->second, b);
	if(MLAT.is_initialized()){
	  vSeeds.emplace_back(*MLAT);
	}
  }

  return vSeeds;
}